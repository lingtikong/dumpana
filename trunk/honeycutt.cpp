#include "driver.h"
#include "voro++.hh"
#include "math.h"
#include <vector>

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to compute the Honeycutt-Andersen index based on the  voronoi neighbors
 * for selected frames. Only bonded pairs are analysed.
 *
 * The first index 1 means the selected pair are bonded, the second gives the
 * number of common neighbors between these two, and the third gives the number
 * of bonds formed between the common neighbors; the fourth are generally 1, while
 * for 144?, if all common neighbors are linked together, it is set to 2.
 *----------------------------------------------------------------------------*/
void Driver::honeycutt_andersen()
{
  char str[MAXLINE];
  int job = 2;
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. HA index based on direct Voronoi info;\n");
  printf("  2. HA index based on refined Voronoi info, skip tiny surfaces;\n");
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

  printf("Please input the file name to output the HA bond index info [ha.dat]: ");
  fgets(str, MAXLINE, stdin);
  char *fname = strtok(str, " \n\t\r\f");
  if (fname == NULL){
    strcpy(str,"ha.dat\n");
    fname = strtok(str, " \n\t\r\f");
  }

  FILE *fp = fopen(fname, "w");
  fprintf(fp, "# id  jd index\n"); fflush(fp);

  one = all[istr];
  int neimax = 24;
  int **bonded, **neilist;
  int natom = one->natom;
  bonded  = memory->create(bonded,natom+1,natom+1,"bonded");
  neilist = memory->create(neilist,neimax+1,natom+1,"neilist");
  for (int i=0; i<= natom; i++)
  for (int j=0; j<= natom; j++) bonded[i][j] = 0;

  // now to loop over all asked images
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // not possible to evaluate voro info for triclinic box
    if (one->triclinic) continue;

    // set local variables
    double xlo = one->xlo, xhi = one->xhi, lx = one->lx;
    double ylo = one->ylo, yhi = one->yhi, ly = one->ly;
    double zlo = one->zlo, zhi = one->zhi, lz = one->lz;
    double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;

    int n = one->natom;
    int *attyp = one->attyp;
    double **atpos = one->atpos;

    // reallocated bonded and neilist if different
    if (n != natom){
      natom = n;
      memory->destroy(bonded);
      memory->destroy(neilist);

      bonded  = memory->create(bonded,natom+1,natom+1,"bonded");
      neilist = memory->create(neilist,neimax+1,natom+1,"neilist");
    }

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

      if (nf > neimax){
        neimax = nf + 2;
        neilist = memory->grow(neilist,neimax+1,natom+1,"neilist");
      }

      if (refine) { // based on refined voro info
        double fcut = surf_min * cell.surface_area();
        std::vector<double> fs;    // face areas
        cell.face_areas(fs);
        int nnei = 0;
        for (int ii=0; ii<nf; ii++){
          if (fs[ii] < fcut) continue;
          int jd = neigh[ii];
          neilist[++nnei][id] = jd;
        }
        neilist[0][id] = nnei;

      } else { // based on direct voro info
        neilist[0][id] = nf;
        for (int ii=0; ii<nf; ii++) neilist[ii+1][id] = neigh[ii];
      }

    } while (cl.inc());

    // to find out if any two atoms are bonded or not
    for (int i=0; i<= natom; i++)
    for (int j=0; j<= natom; j++) bonded[i][j] = 0;

    for (int id=1; id<= natom; id++){
      int nni = neilist[0][id];
      for (int ii=1; ii<= nni; ii++){
        int jd = neilist[ii][id];
        bonded[id][jd] = 1;
      }
    }
    fprintf(fp,"# frame number: %d\n", img);
    // now to analyse the Honeycutt-Andersen bond type info
    std::vector<int> comms;
    for (int id=1; id<= natom; id++){
      int nni = neilist[0][id];
      for (int ii=1; ii<= nni; ii++){
        int jd = neilist[ii][id];
        if (jd >= id) continue;
        int nnj = neilist[0][jd];

        comms.clear();
        for (int jj=1; jj<=nnj; jj++)
        for (int kk=1; kk<=nni; kk++) if (neilist[jj][jd] == neilist[kk][id]) comms.push_back(neilist[jj][jd]);

        int ncomm = comms.size();
        int nbond = 0;
        for (int mm=0; mm<ncomm; mm++)
        for (int nn=mm+1; nn<ncomm; nn++) nbond += bonded[comms[mm]][comms[nn]];

        int nconf = 1;
        // needs to distinct same ncomm-nbond for 144
        if (ncomm == 4 && nbond == 4){
          int ned[4]; ned[0] = ned[1] = ned[2] = ned[3] = 0;
          for (int mm=0; mm<ncomm; mm++)
          for (int nn=mm+1; nn<ncomm; nn++){
            int md = comms[mm], nd = comms[nn];
            ned[mm] += bonded[md][nd];
            ned[nn] += bonded[md][nd];
          }
          for (int mm=0; mm<ncomm; mm++) if (ned[mm] > 2) nconf = 2;
        }

        fprintf(fp,"%d %d 1%d%d%d\n", id, jd, ncomm, nbond, nconf);
      }
    }
  }
  printf("\n"); for (int i=0; i<20; i++) printf("===="); printf("\n");
  fclose(fp);

  memory->destroy(bonded);
  memory->destroy(neilist);

return;
}
