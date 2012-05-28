#include "driver.h"
#include "voro++.hh"
#include "math.h"

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to evaluate the chemical short range order based on the
 * voronoi neighbors for selected frames.
 * The chemical short range order is defined as (Warren-Cowley):
 *  eta = 1. - n_j2i/(N_i*c_j)
 * Where n_j2i is the number of type j neighbors of type i, N_i is the total #
 * of neighbors for type i, c_j is the concentration of type j
 *----------------------------------------------------------------------------*/
void Driver::csro()
{
  char str[MAXLINE];
  int job = 2;
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. CSRO based on direct Voronoi info;\n");
  printf("  2. CSRO based on refined Voronoi info, skip tiny surfaces;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 2){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  double surf_min = 1.e-4;
  int refine = 0;

  if (job == 2){
    printf("Please input your criterion for tiny surfaces [%g]: ", surf_min);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) surf_min = atof(ptr);
    printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", surf_min);

    if (surf_min > 0.) refine = 1;
  }

  one = all[istr];
  int ntype = one->ntype;
  bigint *NumType, **NumNei;
  NumType= memory->create(NumType,ntype+1,"csro:NumType");
  NumNei = memory->create(NumNei,ntype+1,ntype+1,"csro:NumNei");
  for (int i=0; i<= ntype; i++){
    NumType[i] = 0;
    for (int j=0; j<= ntype; j++) NumNei[i][j] = 0;
  }

  // now to do the real job
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // not possible to evaluate voro info for triclinic box
    if (one->triclinic) continue;

    // ntype of different frames must be the same
    if (one->ntype != ntype) continue;

    // set local variables
    double xlo = one->xlo, xhi = one->xhi, lx = one->lx;
    double ylo = one->ylo, yhi = one->yhi, ly = one->ly;
    double zlo = one->zlo, zhi = one->zhi, lz = one->lz;
    double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;

    int n = one->natom;
    int *attyp = one->attyp;
    double **atpos = one->atpos;
    for (int i=1; i<=ntype; i++) NumType[i] += one->numtype[i];

    // need cartesian coordinates
    one->dir2car();

    // compute optimal size for container, and then contrust it
    double l = pow(double(n)/(5.6*lx*ly*lz), 1./3.);
    int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);
    voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

    // put atoms into the container
    for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor cell;
    voro::c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(cell,cl)){
      int id = cl.pid();
      int ip = one->attyp[id];

      std::vector<int> neigh;    // neigh list
      cell.neighbors(neigh);
      int nf = neigh.size();

      if (refine) { // based on refined voro info
        double fcut = surf_min * cell.surface_area();
        std::vector<double> fs;    // face areas
        cell.face_areas(fs);
        for (int ii=0; ii<nf; ii++){
          if (fs[ii] < fcut) continue;
          int jd = neigh[ii];
          int jp = one->attyp[jd];
          NumNei[ip][jp]++;
        }

      } else { // based on direct voro info
        for (int ii=0; ii<nf; ii++){
          int jd = neigh[ii];
          int jp = one->attyp[jd];
          NumNei[ip][jp]++;
        }
      }

    } while (cl.inc());
  }
  
  bigint ntotal = 0;
  double concentration[ntype+1];
  for (int i=1; i<= ntype; i++) ntotal += NumType[i];
  for (int i=1; i<= ntype; i++) concentration[i] = double(NumType[i])/double(ntotal);

  printf("\n"); for (int i=0; i<20; i++) printf("____");
  printf("\nType         :"); for (int ip=1; ip<= ntype; ip++) printf("    %2d    ",ip);
  printf("\n"); for (int i=0; i<20; i++) printf("----");
  printf("\nConcentration:"); for (int ip=1; ip<= ntype; ip++) printf("%10.6f", concentration[ip]);
  printf("\n"); for (int i=0; i<20; i++) printf("____");
  printf("\nCSRO         :"); for (int jp=1; jp<= ntype; jp++) printf("    %2d    ",jp);
  printf("\n"); for (int i=0; i<20; i++) printf("----");
  for (int ip=1; ip<= ntype; ip++){
    printf("\n      %2d     :", ip);
    ntotal = 0;
    for (int jp=1; jp<= ntype; jp++) ntotal += NumNei[ip][jp];
    for (int jp=1; jp<= ntype; jp++) printf("%10.6f", 1.-double(NumNei[ip][jp])/double(ntotal)/concentration[jp]);
  }
  printf("\n"); for (int i=0; i<20; i++) printf("===="); printf("\n");

  memory->destroy(NumType);
  memory->destroy(NumNei);
return;
}
