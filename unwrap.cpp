#include "driver.h"
#include "math.h"

/*------------------------------------------------------------------------------
 * Method to unwrap atoms so that the direct distance between bonded atoms are
 * at the expected relative positions.
 * 
 * It should work well for molecule system, but not for bulk system where all atoms
 * form a network.
 *----------------------------------------------------------------------------*/
void Driver::unwrap()
{
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("   To unwrap the atomic positions for the selected frames.");

  one = all[istr];
  int ntype = one->ntype;
  choose_neighbor_method(0);

  // now to do the real method
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // ntype of different frames must be the same
    if (one->ntype != ntype) continue;

    // get neighbor list
    if (neighbor_method == 1) one->ComputeVoro(mins);
    else one->ComputeNeiList(r2cuts);

    // set local variables
    int *unwraped = new int[one->natom+1];
    int *parents = new int[one->natom+1];
    for (int i = 1; i <= one->natom; ++i) unwraped[i] = 0;

    for (int id = 1; id <= one->natom; ++id){
        if (unwraped[id] == 0) unwraped[id] = 1;
        parents[0] = id;
        unwrap_neighbors(id, unwraped, parents, 0);
    }

    if (min_mem) one->FreeVoro();
    delete unwraped;
  }
  
  printf("\nSelected frames are now unwrapped, you can use option 11 to output them.\n");
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Recursive method to unwrap the neighbors of current central atom.
 *----------------------------------------------------------------------------*/
void Driver::unwrap_neighbors(int id, int *status, int *parents, int level)
{
   for (int jj = 1; jj <= one->neilist[0][id]; ++jj){
       int jd = one->neilist[jj][id];
       if (status[jd] == 0){
          double xij = one->atpos[jd][0] - one->atpos[id][0];
          double yij = one->atpos[jd][1] - one->atpos[id][1];
          double zij = one->atpos[jd][2] - one->atpos[id][2];
          one->ApplyPBC(xij, yij, zij);

          one->atpos[jd][0] = one->atpos[id][0] + xij;
          one->atpos[jd][1] = one->atpos[id][1] + yij;
          one->atpos[jd][2] = one->atpos[id][2] + zij;

          status[jd] = 1;
       }
       int done = 0;
       for (int il = 0; il < level; il ++) done += (jd == parents[il]);
       if (done == 0){
          parents[level+1] = jd;
          unwrap_neighbors(jd, status, parents, level+1);
       }
   }

   return;
}
