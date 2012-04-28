#include "driver.h"
#include "voro++.hh"
#include "math.h"

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to compute the voro index info for the selected frames
 * Resultant voro index info will be written to files: voro%{tstep}
 *----------------------------------------------------------------------------*/
void Driver::voro()
{
  // ask for frames
  setrange();

  // ask for threshold
  char str[MAXLINE];
  double fmin = 0.001;
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nPlease input a threshold for the evaluation of voro neighbors, if the\n");
  printf("surface area corresponding to an atom takes less ratio of the whole surface\n");
  printf("than the defined threshold, it will not be seen as a neighbor [%lg]: ", fmin);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) fmin = fabs(atof(ptr));
  printf("\nA threshold of %g will be used.\n\n", fmin);

  // now to do the real job
  for (img = istr; img <= iend; img += inc){
    one = all[img];

    double xlo = one->xlo, xhi = one->xhi, lx = one->lx;
    double ylo = one->ylo, yhi = one->yhi, ly = one->ly;
    double zlo = one->zlo, zhi = one->zhi, lz = one->lz;
    double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;
    int n = one->natom;
    // need cartesian coordinates
    one->dir2car();

    // compute optimal size for container, and then contrust it
    double l = pow(double(n)/(5.6*lx*ly*lz), 1./3.);
    int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);
    voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

    // put atoms into the container
    con.clear();
    for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);

    // open file for output
    sprintf(str,"voro_%d.dat ", one->tstep);
    ptr = strtok(str," \n\t\r\f");
    FILE *fp = fopen(ptr, "w");
    fprintf(fp,"#Box info: %lg %lg %lg %lg %lg %lg\n", xlo, xhi, ylo, yhi, zlo, zhi);
    fprintf(fp,"# 1 2    3 4 5 6   7         8    9    10\n");
    fprintf(fp,"#id type x y z vol voroindex f5%% NNei NeiList surfaceareas\n");
    // loop over all atoms and compute their voronoi cell
    voro::voronoicell_neighbor c, c2;
    voro::c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(c,cl)){
      int id;
      double x, y, z;
      std::vector<int> neigh;
      std::vector<double> fs;
       
      cl.pos(x,y,z);
      id = cl.pid();
      c.neighbors(neigh);
      c.face_areas(fs);

      // initialize a second cell
      c2.init(-lx,lx,-ly,ly,-lz,lz);

      // add condition on surface
      double fcut = fmin * c.surface_area();
      for (int i=0; i<(signed int) fs.size(); i++){
        if (fs[i] > fcut){
          int j = neigh[i];

          double xij = one->atpos[j][0]-x;
          while (xij > hx) xij -= lx;
          while (xij <-hx) xij += lx;
          double yij = one->atpos[j][1]-y;
          while (yij > hy) yij -= ly;
          while (yij <-hy) yij += ly;
          double zij = one->atpos[j][2]-z;
          while (zij > hz) zij -= lz;
          while (zij <-hz) zij += lz;

          c2.nplane(xij,yij,zij,j);
        }
      }

      // output voro index info
      std::vector<int> ff;
      c2.face_areas(fs);
      c2.neighbors(neigh);
      c2.face_freq_table(ff);
      int index[7];
      int nn = ff.size();
      for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
      for (int i=MIN(6,nn)+1; i<= 6; i++) index[i] = 0;
      nn = fs.size();
      double wf = double(index[5])/double(nn)*100.;
      fprintf(fp,"%d %d %lg %lg %lg %lg %d,%d,%d,%d %g %d", id, one->attyp[id], x, y, z, c2.volume(), index[3], index[4], index[5], index[6], wf, nn);
      for (int i=0; i<nn; i++) fprintf(fp," %d", neigh[i]);
      for (int i=0; i<nn; i++) fprintf(fp," %lg", fs[i]);
      fprintf(fp,"\n");

    } while (cl.inc());

    fclose(fp);
    printf("Frame %d done, voro info written to: %s\n", img+1, ptr);
  }
  for (int i=0; i<20; i++) printf("====");  printf("\n");
return;
}
/*----------------------------------------------------------------------------*/
