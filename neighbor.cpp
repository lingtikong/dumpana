#include "atom.h"
#include "time.h"
#include "random.h"
#include "global.h"

/*------------------------------------------------------------------------------
 * Method to calculate the neighbor list based on distances
 *------------------------------------------------------------------------------
 * Parameters:
 *  cutoffs  : square of cutoff distances for each pair of atom types
 *----------------------------------------------------------------------------*/
void DumpAtom::ComputeNeiList(double **rcut_sq)
{
  if (neilist) memory->destroy(neilist);
  memory->create(neilist, MaxNei+1, natom+1, "neilist");

  for (int i = 0; i <= natom; ++i) neilist[0][i] = 0;

  for (int i = 1; i <= natom; ++i){
      int ip = attyp[i];
      for (int j = i+1; j <= natom; ++j){
          int jp = attyp[j];

          double xij = atpos[j][0] - atpos[i][0];
          double yij = atpos[j][1] - atpos[i][1];
          double zij = atpos[j][2] - atpos[i][2];

          ApplyPBC(xij, yij, zij);

          double r2 = get_dist2(xij, yij, zij);

          if (r2 > rcut_sq[ip][jp]) continue;

          int ni = ++neilist[0][i];
          int nj = ++neilist[0][j];
          int nmax = MAX(ni, nj);
          if (nmax > MaxNei){
             MaxNei = nmax + 2;
             memory->grow(neilist, MaxNei+1, natom+1, "neilist");
          }

          neilist[ni][i] = j;
          neilist[nj][j] = i;
      }
  }
return;
}
