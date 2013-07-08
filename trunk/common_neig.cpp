#include "common_neig.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MAXNEAR   25
#define MAXCOMMON 20
#define MAXLINE   512
#define ZERO 1.e-10

enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};
enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ----------------------------------------------------------------------
 * Reference: Comp. Phys. Comm. 177:518, (2007).
 * ----------------------------------------------------------------------
 * job  (in)  : 1 for cna; else for cnp
 * ntm  (in)  : # of atoms
 * list (in)  : nearest neighbor list
 * pos  (in)  : atomic positions in cartesian
 * box  (in)  : box dimensions
 * fp   (in)  : FILE to output the result
 * ---------------------------------------------------------------------- */
ComputeCNAAtom::ComputeCNAAtom(const int job, const int ntm, int **list, double **pos, double *box, FILE *fp)
{
  x = pos;
  natom   = ntm;
  nearest = list;
  L[0] = box[0]; L[1] = box[1]; L[2] = box[2];
  xy = box[3]; xz = box[4]; yz = box[5];
  for (int i=0; i<3; i++) hL[i] = 0.5*L[i];

  non_ortho_box = 0;
  if (xy*xy+xz*xz+yz*yz > ZERO) non_ortho_box = 1;

  memory = new Memory();
  pattern  = memory->create(pattern, natom+1, "pattern");

  // to the real job
  if (job == 1) compute_cna();
  else if (job == 2) compute_cnp();
  else centro_atom(-job);

  output(fp);

return;
}

/* ----------------------------------------------------------------------
 * Deconstuctor, free memory
 * ---------------------------------------------------------------------- */
ComputeCNAAtom::~ComputeCNAAtom()
{
  memory->destroy(pattern);
  nearest = NULL;
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
	       if (bonded(j,k)) {
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
      double Rij[3], xik[3], xjk[3];
      Rij[0] = Rij[1] = Rij[2] = 0.;

	   for (int inear = 1; inear <= nearest[0][i]; inear++)
	   for (int jnear = 1; jnear <= nearest[0][j]; jnear++){
	     if (nearest[inear][i] == nearest[jnear][j]) {
          int k = nearest[inear][i];
          double xik[3], xjk[3];
          
          for (int idim=0; idim<3; idim++){
            xik[idim] = x[k][idim] - x[i][idim];
            xjk[idim] = x[k][idim] - x[j][idim];
          }
          // apply pbc
          apply_pbc(xik[0], xik[1], xik[2]);
          apply_pbc(xjk[0], xjk[1], xjk[2]);
            
          for (int idim=0; idim<3; idim++) Rij[idim] += xik[idim] + xjk[idim];
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

/* ----------------------------------------------------------------------
 * Private method, to check if two atoms are bonded to each other
 * ---------------------------------------------------------------------- */
int ComputeCNAAtom::bonded(int id, int jd)
{
  int ni = nearest[0][id];
  for (int jj=1; jj <= ni; jj++){
    if (nearest[jj][id] == jd) return 1;
  }

return 0;
}

/* ----------------------------------------------------------------------
 * Private method, Copied from LAMMPS compute_centro_atom
 * ---------------------------------------------------------------------- */
#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

void ComputeCNAAtom::select(int k, int n, double *arr)
{
  int i,ir,j,l,mid;
  double a,tmp;

  arr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir])
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir])
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      ISWAP(iarr[mid],iarr[l+1])
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir])
        ISWAP(iarr[l+1],iarr[ir])
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1])
        ISWAP(iarr[l],iarr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j])
        ISWAP(iarr[i],iarr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
 * Central symmetry parameter
 * ---------------------------------------------------------------------- */
void ComputeCNAAtom::centro_atom(const int nnn)
{

  int nhalf = nnn/2;
  int npairs = nnn * (nnn-1) / 2;
  double *pairs = memory->create(pairs, npairs, "pairs");

  double *distsq;
  int *neighb;
  int maxneigh = nnn+nnn;
  distsq = memory->create(distsq,maxneigh,"centro/atom:distsq");
  neighb = memory->create(neighb,maxneigh,"centro/atom:neighb");

  // compute centro-symmetry parameter for each atom in group, use full neighbor list

  for (int i = 1; i <= natom; i++) {
    pattern[i] = 0.;
    int jnum = nearest[0][i];

    // insure distsq and nearest arrays are long enough
    if (jnum > maxneigh) {
      maxneigh = jnum;
      distsq = memory->grow(distsq,maxneigh,"centro/atom:distsq");
      neighb = memory->grow(neighb,maxneigh,"centro/atom:neighb");
    }

    // loop over list of all neighbors within force cutoff
    // distsq[] = distance sq to each
    // nearest[] = atom indices of neighbors

    int n = 0;
    for (int jj = 1; jj <= jnum; jj++) {
      int j = nearest[jj][i];

      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      apply_pbc(delx, dely, delz);

      double rsq = delx*delx + dely*dely + delz*delz;

      distsq[n] = rsq;
      neighb[n++] = j;
    }

    // if not nnn neighbors, centro = 0.0

    if (n < nnn) {
      pattern[i] = 0.;
      continue;
    }

    // store nnn nearest neighs in 1st nnn locations of distsq and nearest

    select2(nnn,n,distsq,neighb);

    // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
    // pairs = squared length of each R
    n = 0;
    for (int jj = 0; jj < nnn; jj++) {
      int j = neighb[jj];
      for (int kk = jj+1; kk < nnn; kk++) {
        int k = neighb[kk];
        double xij[3], xik[3];
        for (int idim=0; idim<3; idim++){
          xij[idim] = x[j][idim] - x[i][idim];
          xik[idim] = x[k][idim] - x[i][idim];
        }
        double delx = xij[0] + xik[0];
        double dely = xij[1] + xik[1];
        double delz = xij[2] + xik[2];
        apply_pbc(delx, dely, delz);

        pairs[n++] = delx*delx + dely*dely + delz*delz;
      }
    }

    // store nhalf smallest pair distances in 1st nhalf locations of pairs

    select(nhalf,npairs,pairs);

    // centrosymmetry = sum of nhalf smallest squared values

    double value = 0.;
    for (int jj = 0; jj < nhalf; jj++) value += pairs[jj];
    pattern[i] = value;
  }

  memory->destroy(pairs);
  memory->destroy(distsq);
  memory->destroy(neighb);
return;
}

/* ----------------------------------------------------------------------
 * To apply PBC
 * ---------------------------------------------------------------------- */
void ComputeCNAAtom::apply_pbc(double &xtmp, double &ytmp, double &ztmp)
{
  if (non_ortho_box){
    while (ztmp > hL[2]){
      xtmp -= xz;
      ytmp -= yz;
      ztmp -= L[2];
    }
    while (ztmp <-hL[2]){
      xtmp += xz;
      ytmp += yz;
      ztmp += L[2];
    }
  
    while (ytmp > hL[1]){
      xtmp -= xy;
      ytmp -= L[1];
    }
    while (ytmp <-hL[1]){
      xtmp += xy;
      ytmp += L[1];
    }
  
    while (xtmp > hL[0]) xtmp -= L[0];
    while (xtmp <-hL[0]) xtmp += L[0];
  
  } else {

    while (xtmp > hL[0]) xtmp -= L[0];
    while (ytmp > hL[1]) ytmp -= L[1];
    while (ztmp > hL[2]) ztmp -= L[2];
    while (xtmp <-hL[0]) xtmp += L[0];
    while (ytmp <-hL[1]) ytmp += L[1];
    while (ztmp <-hL[2]) ztmp += L[2];
  
  }

return;
}
/*------------------------------------------------------------------------------ */
