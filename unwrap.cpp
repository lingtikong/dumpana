#include "driver.h"
#include "math.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to unwrap atoms so that the direct distance between bonded atoms are
 * at the expected relative positions.
 * 
 * It should work well for molecule system, but not for bulk system where all atoms
 * form a network.
 *----------------------------------------------------------------------------*/
void Driver::unwrap()
{
  printf("\n"); for (int i = 0; i < 20; ++i) printf("====");
  printf("\nTo unwrap the atomic positions for the selected frames, so that bonded atoms\n");
  printf("are neighboring to each other, instead of crossing the PBC boundary. It might\n");
  printf("only work for molecular systems.\n");

  one = all[istr];
  int ntype = one->ntype;
  choose_neighbor_method(0);

  Timer * timer = new Timer();
  // now to do the real method
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // ntype of different frames must be the same
    if (one->ntype != ntype) continue;

    // Must be in fraction coordinate
    one->car2dir();

    // get neighbor list
    if (neighbor_method == 1) one->ComputeVoro(mins);
    else one->ComputeNeiList(r2cuts);

    // set local variables
    int *visited = new int[one->natom+1];
    for (int id = 0; id <= one->natom; ++id) visited[id] = 0;

    for (int id = 1; id <= one->natom; ++id)
        if (visited[id] == 0){
           for (int idim = 0; idim < 3; ++idim){
              while (one->atpos[id][idim] >= 0.5) one->atpos[id][idim] -= 1.;
              while (one->atpos[id][idim] < -0.5) one->atpos[id][idim] += 1.;
           }

           unwrap_neighbors(id, visited);
        }

    one->dir2car();
    if (min_mem) one->FreeVoro();

    delete visited;
  }
  
  printf("\nSelected frames are now unwrapped, now you can output them.\n");
  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Recursive method to unwrap the neighbors of current central atom.
 *----------------------------------------------------------------------------*/
void Driver::unwrap_neighbors(int id, int *visited)
{
   visited[id] = 1;

   for (int jj = 1; jj <= one->neilist[0][id]; ++jj){
       int jd = one->neilist[jj][id];
       if (visited[jd] == 0){
          double xij = one->atpos[jd][0] - one->atpos[id][0];
          double yij = one->atpos[jd][1] - one->atpos[id][1];
          double zij = one->atpos[jd][2] - one->atpos[id][2];
          one->ApplyPBC(xij, yij, zij);

          one->atpos[jd][0] = one->atpos[id][0] + xij;
          one->atpos[jd][1] = one->atpos[id][1] + yij;
          one->atpos[jd][2] = one->atpos[id][2] + zij;

          unwrap_neighbors(jd, visited);
       }
   }

   return;
}
