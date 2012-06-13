#include "common_neig.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MAXNEAR   25
#define MAXCOMMON 20
#define MAXLINE   512

enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};
enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ----------------------------------------------------------------------
 * Reference: Comp. Phys. Comm. 177:518, (2007).
 * ---------------------------------------------------------------------- */

ComputeCNAAtom::ComputeCNAAtom(const int job, const int ntm, int **list, int **pair, double **pos, double *box, FILE *fp)
{
  x = pos;
  natom   = ntm;
  nearest = list;
  bonded  = pair;
  L[0] = box[0]; L[1] = box[1]; L[2] = box[2];
  for (int i=0; i<3; i++) hL[i] = 0.5*L[i];

  memory = new Memory();
  pattern  = memory->create(pattern, natom+1, "pattern");

  // to the real job
  if (job == 1) compute_cna();
  else          compute_cnp();

  output(fp);

return;
}

/* ---------------------------------------------------------------------- */

ComputeCNAAtom::~ComputeCNAAtom()
{
  memory->destroy(pattern);
  nearest = NULL;
  bonded = NULL;
  x = NULL;

  delete memory;
}

/* ----------------------------------------------------------------------
 * Common Neighbor Analysis
 * ---------------------------------------------------------------------- */
void ComputeCNAAtom::compute_cna()
{
  int i,j,k,jj,kk,m,n,inum,jnum,inear,jnear;
  int ncommon,nbonds,maxbonds,minbonds;
  int nfcc,nhcp,nbcc4,nbcc6,nico,cj,ck,cl,cm;
  int MaxNear = 15, MaxComm = 20;

  int **cna, *common, *bonds;

  cna = memory->create(cna,MaxNear+1,4,"cna");
  bonds = memory->create(bonds, MaxComm, "bonds");
  common = memory->create(common, MaxComm, "common");

  // compute CNA for each atom in group
  // only performed if # of nearest neighbors = 12 or 14 (fcc,hcp)

  for (i = 1; i <= natom; i++) {
    if (nearest[0][i] != 12 && nearest[0][i] != 14) {
      pattern[i] = OTHER;
      continue;
    }

    // loop over near neighbors of I to build cna data structure
    // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
    // cna[k][NBONDS] = # of bonds between those common neighbors
    // cna[k][MAXBOND] = max # of bonds of any common neighbor
    // cna[k][MINBOND] = min # of bonds of any common neighbor
    /* not necessary, since only considers 12 and 14
    if (nearest[0][i] > MaxNear){
      MaxNear = nearest[0][i];
      cna = memory->grow(cna,MaxNear+1,4,"cna");
    } */

    for (m = 1; m <= nearest[0][i]; m++) {
      j = nearest[m][i];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

	   ncommon = 0;
	   for (inear = 1; inear <= nearest[0][i]; inear++)
	   for (jnear = 1; jnear <= nearest[0][j]; jnear++){
	     if (nearest[inear][i] == nearest[jnear][j]) {
          if (ncommon >= MaxComm){
            MaxComm = ncommon + 5;
            bonds = memory->grow(bonds,MaxComm,"bonds");
            common = memory->grow(common,MaxComm,"common");
          }
	       common[ncommon++] = nearest[inear][i];
        }
	   }

      cna[m][NCOMMON] = ncommon;

      // calculate total # of bonds between common neighbor atoms
      // also max and min # of common atoms any common atom is bonded to
      // bond = pair of atoms within cutoff

      for (n = 0; n < ncommon; n++) bonds[n] = 0;

      nbonds = 0;
      for (jj = 0; jj < ncommon; jj++) {
	     j = common[jj];
	     for (kk = jj+1; kk < ncommon; kk++) {
	       k = common[kk];
	       if (bonded[j][k]) {
	         nbonds++;
	         bonds[jj]++;
	         bonds[kk]++;
	       }
	     }
      }
      cna[m][NBOND] = nbonds;

      maxbonds = 0;
      minbonds = MaxComm;
      for (n = 0; n < ncommon; n++) {
	     maxbonds = MAX(bonds[n],maxbonds);
	     minbonds = MIN(bonds[n],minbonds);
      }      
      cna[m][MAXBOND] = maxbonds;
      cna[m][MINBOND] = minbonds;
    }

    // detect CNA pattern of the atom

    nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;
    pattern[i] = OTHER;

    if (nearest[0][i] == 12) {
      for (inear = 0; inear < 12; inear++) {
        cj = cna[inear+1][NCOMMON];
        ck = cna[inear+1][NBOND];
        cl = cna[inear+1][MAXBOND];
        cm = cna[inear+1][MINBOND];
        if (cj == 4 && ck == 2 && cl == 1 && cm == 1) nfcc++;
        else if (cj == 4 && ck == 2 && cl == 2 && cm == 0) nhcp++;
        else if (cj == 5 && ck == 5 && cl == 2 && cm == 2) nico++;
      }
      if (nfcc == 12) pattern[i] = FCC;
      else if (nfcc == 6 && nhcp == 6) pattern[i] = HCP;
      else if (nico == 12) pattern[i] = ICOS;
      
    } else if (nearest[0][i] == 14) {
      for (inear = 0; inear < 14; inear++) {
        cj = cna[inear+1][NCOMMON];
        ck = cna[inear+1][NBOND];
        cl = cna[inear+1][MAXBOND];
        cm = cna[inear+1][MINBOND];
        if (cj == 4 && ck == 4 && cl == 2 && cm == 2) nbcc4++;
        else if (cj == 6 && ck == 6 && cl == 2 && cm == 2) nbcc6++;
      }
      if (nbcc4 == 6 && nbcc6 == 8) pattern[i] = BCC;
    }
  }

  memory->destroy(cna);
  memory->destroy(bonds);
  memory->destroy(common);
}

/* ----------------------------------------------------------------------
 * Common Neighborhood Parameter
 * ---------------------------------------------------------------------- */
void ComputeCNAAtom::compute_cnp()
{
  for (int i = 1; i <= natom; i++) {
    pattern[i] = 0.;
    for (int m = 1; m <= nearest[0][i]; m++) {
      int j = nearest[m][i];

      // common = list of neighbors common to atom I and atom J
      double Rij[3];
      Rij[0] = Rij[1] = Rij[2] = 0.;

	   for (int inear = 1; inear <= nearest[0][i]; inear++)
	   for (int jnear = 1; jnear <= nearest[0][j]; jnear++){
	     if (nearest[inear][i] == nearest[jnear][j]) {
          int k = nearest[inear][i];
          for (int idim=0; idim<3; idim++){
            double xik = x[k][idim] - x[i][idim];
            double xjk = x[k][idim] - x[j][idim];
            while (xik > hL[idim]) xik -= L[idim];
            while (xik <-hL[idim]) xik += L[idim];
            while (xjk > hL[idim]) xjk -= L[idim];
            while (xjk <-hL[idim]) xjk += L[idim];
            Rij[idim] += xik + xjk;
          }
        }
	   }
      pattern[i] += Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2];
    }
    if (nearest[0][i] > 0) pattern[i] /= double(nearest[0][i]);
  }

return;
}

/* ----------------------------------------------------------------------
 * write Common Neighbor Parameter result
 * ---------------------------------------------------------------------- */
void ComputeCNAAtom::output(FILE *fp)
{
  for (int i=1; i<=natom; i++) fprintf(fp,"%d %lg %lg %lg %lg\n", i, x[i][0], x[i][1], x[i][2], pattern[i]);
return;
}
