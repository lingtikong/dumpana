#include "driver.h"
#include "voro++.hh"
#include "common_neig.h"
#include "math.h"
#include <vector>

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Driver to execute the Common neighbor analysis or common neighbor parameter
 * computations;
 * Ref: Comp. Phys. Comm. 177:518, (2007).
 *----------------------------------------------------------------------------*/
void Driver::Compute_CNACNP()
{
  char str[MAXLINE];
  int job = 1;
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. CNA based on direct Voronoi info;\n");
  printf("  2. CNA based on refined Voronoi info, skip tiny surfaces;\n");
  printf("  3. CNP based on direct Voronoi info;\n");
  printf("  4. CNP based on refined Voronoi info, skip tiny surfaces;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 4){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  double surf_min = 1.e-4;
  int refine = 1-job%2;
  job = (job+1)/2; // 1, cna; 2, cnp

  if (refine){
    printf("Please input your criterion for tiny surfaces [%g]: ", surf_min);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) surf_min = atof(ptr);
    printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", surf_min);

    if (surf_min > 0.) refine = 1;
  }
  
  char *fname = new char[8];
  if (job == 1) strcpy(fname, "cna.dat");
  else strcpy(fname, "cnp.dat");

  printf("Please input the file name to output the CNA/CNP info[%s]: ", fname);
  fgets(str, MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr){
    delete []fname;
    fname = new char[strlen(ptr)+1];
    strcpy(fname, ptr);
  }

  FILE *fp = fopen(fname, "w");
  if (job == 1)
    fprintf(fp,"# CNA: 1, FCC; 2, HCP; 3, BCC; 4, ICOS; 4, OTHER; 5, UNKNOWN.\n# id x y z cna\n");
  else fprintf(fp, "# id x y z cnp\n");
  fflush(fp);

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

    // set local variables
    double xlo = one->xlo, xhi = one->xhi;
    double ylo = one->ylo, yhi = one->yhi;
    double zlo = one->zlo, zhi = one->zhi;
    double lx  = one->lx,  ly  = one->ly,  lz = one->lz;
    double xy  = one->xy,  xz  = one->xz,  yz = one->yz;
    double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;

    int n = one->natom;
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

    if (one->triclinic){ // non-orthogonal box

      voro::container_periodic con(lx,xy,ly,xz,yz,lz,nx,ny,nz,8);
  
      // put atoms into the container
      for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
  
      // loop over all particles and compute their voronoi cell
      voro::voronoicell_neighbor cell;
      voro::c_loop_all_periodic cl(con);
      if (cl.start()) do if (con.compute_cell(cell,cl)){
        int id = cl.pid();
  
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

    } else { // orthogonal box

      voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);
  
      // put atoms into the container
      for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
  
      // loop over all particles and compute their voronoi cell
      voro::voronoicell_neighbor cell;
      voro::c_loop_all cl(con);
      if (cl.start()) do if (con.compute_cell(cell,cl)){
        int id = cl.pid();
  
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
    }

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
    // now to compute the CNA/CNP info
    ComputeCNAAtom *cna = new ComputeCNAAtom(job, natom, neilist, bonded, atpos, one->box, fp);
    delete cna;

    printf("Frame %d done, CNA/CNP info written to: %s\n", img+1, fname);
  }
  printf("\n"); for (int i=0; i<20; i++) printf("===="); printf("\n");
  fclose(fp);

  memory->destroy(bonded);
  memory->destroy(neilist);
  delete []fname;

return;
}
