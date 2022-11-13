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
void DumpAtom::ComputeNeiList(double ***rcut_sq)
{
  if (neilist) memory->destroy(neilist);
  memory->create(neilist, MaxNei+1, natom+1, "neilist");

  for (int id = 0; id <= natom; ++id) neilist[0][id] = 0;

  for (int id = 1; id <= natom; ++id){
      int ip = attyp[id];
      for (int jd = id+1; jd <= natom; ++jd){
          int jp = attyp[jd];

          double xij = atpos[jd][0] - atpos[id][0];
          double yij = atpos[jd][1] - atpos[id][1];
          double zij = atpos[jd][2] - atpos[id][2];

          ApplyPBC(xij, yij, zij);

          double r2 = get_dist2(xij, yij, zij);

          if (r2 < rcut_sq[ip][jp][0] || r2 > rcut_sq[ip][jp][1]) continue;

          int ni = ++neilist[0][id];
          int nj = ++neilist[0][jd];
          int nmax = MAX(ni, nj);
          if (nmax > MaxNei){
             MaxNei = nmax + 2;
             memory->grow(neilist, MaxNei+1, natom+1, "neilist");
          }

          neilist[ni][id] = jd;
          neilist[nj][jd] = id;
      }
  }
return;
}
