#include "driver.h"
#include "voro++.hh"
#include "math.h"

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to do voronoi diagram analysis for the selected frames
 * Resultant voro index info will be written to files: voro%{tstep}
 *----------------------------------------------------------------------------*/
void Driver::voro()
{
  char str[MAXLINE];
  int job = 4;
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. Voronoi index info;\n");
  printf("  2. Refined Voronoi info, skip tiny surfaces;\n");
  printf("  3. Refined Voronoi info, skip ultra short edges;\n");
  printf("  4. Refined Voronoi info, skip tiny surfaces and short edges;\n");
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

  double surf_min = 1.e-4, edge_min = 1.e-4;
  int nminnei = 0;

  if (job == 2 || job == 4){
    printf("Please input your criterion for tiny surfaces [%g]: ", surf_min);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) surf_min = atof(ptr);
    printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", surf_min);

    printf("Sometimes it might be desirable to keep a minimum # of neighbors when refining\n");
    printf("the Voronoi index, for example, keep at least 14 for a bcc lattice, 12 for hcp\n");
    printf("or fcc. If you prefer to do so, input a positive number now [%d]: ", nminnei);
    if (count_words(fgets(str,MAXLINE, stdin)) > 0){
      nminnei = atoi(strtok(str, " \n\t\r\f"));
      if (nminnei < 1) nminnei = 0;
      printf("\nA minimum number of %d neighobrs will be kept no matter how tiny the surface is.\n", nminnei);
    }
  }

  if (job == 3 || job == 4){
    printf("Please input your criterion for ultra short   [%g]: ", edge_min);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) edge_min = atof(ptr);
    printf("Edges whose length takes less ratio than %lg will be skiped!\n\n", edge_min);
  }

  printf("Please input the prefix for output files [voro]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"voro"); ptr = strtok(str, " \n\t\r\f");}
  char *prefix = new char[strlen(ptr)+1];
  strcpy(prefix, ptr);
  
  FILE *fpsurf, *fpedge;
  if (flag_out & OutSurf) fpsurf = fopen("surf_ratio.dat", "w");
  if (flag_out & OutEdge) fpedge = fopen("edge_ratio.dat", "w");

  // now to do the real job
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

    int *attyp = one->attyp;
    double **atpos = one->atpos;

    // need cartesian coordinates
    one->dir2car();

    // open file for output
    sprintf(str,"%s_%d.dat ", prefix, one->tstep);
    ptr = strtok(str," \n\t\r\f");
    FILE *fp = fopen(ptr, "w");
    fprintf(fp,"#Box info: %lg %lg %lg %lg %lg %lg %d %lg %lg %lg\n", xlo, xhi, ylo, yhi, zlo, zhi, n, xy, xz, yz);
    fprintf(fp,"# 1  2    3  4  5  6   7         8    9    10\n");
    fprintf(fp,"# id type x  y  z  vol voroindex f5%%  NNei NeiList surfaceareas\n");

    // compute optimal size for container, and then contrust it
    double l = pow(double(n)/(5.6*lx*ly*lz), 1./3.);
    int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);

    if (one->triclinic){

      voro::container_periodic con(lx,xy,ly,xz,yz,lz,nx,ny,nz,8);
      // put atoms into the container
      for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
  
      // loop over all particles and compute their voronoi cell
      voro::voronoicell_neighbor c1, c2, *cell;
      voro::c_loop_all_periodic cl(con);
      if (cl.start()) do if (con.compute_cell(c1,cl)){
        int id;
        double x, y, z, vol;
        std::vector<int> ff;       // face_freq
        std::vector<int> neigh;    // neigh list
        std::vector<double> fs;    // face areas
        int index[7];
        for (int i=0; i<7; i++) index[i] = 0;
         
        cl.pos(x,y,z);
        id = cl.pid();
        c1.neighbors(neigh);
        c1.face_areas(fs);
        cell = &c1;
  
        // refine the voronoi cell if asked by removing tiny surfaces
        if ((job == 2 || job == 4) && surf_min > 0.){
          c2.init(-lx,lx,-ly,ly,-lz,lz);
    
          int nf = fs.size();
          // sort neighbors by area if asked to keep a minimum # of neighbors
          if (nminnei > 0){
            for (int i=0; i<nf; i++)
            for (int j=i+1; j<nf; j++){
              if (fs[j] > fs[i]){
                double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
                int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
              }
            }
          }
  
          // add condition on surface
          double fcut = surf_min * cell->surface_area();
          for (int i=0; i<nf; i++){
            if (flag_out & OutSurf) fprintf(fpsurf, "%lg\n", fs[i]/cell->surface_area());
            if (i < nminnei || fs[i] > fcut){
              int j = neigh[i];
    
              // apply pbc
              double xij = one->atpos[j][0]-x;
              double yij = one->atpos[j][1]-y;
              double zij = one->atpos[j][2]-z;

              while (zij > hz){
                xij -= xz;
                yij -= yz;
                zij -= lz;
              }
              while (zij <-hz){
                xij += xz;
                yij += yz;
                zij += lz;
              }

              while (yij > hy){
                xij -= xy;
                yij -= ly;
              }
              while (yij <-hy){
                xij += xy;
                yij += ly;
              }
  
              while (xij > hx) xij -= lx;
              while (xij <-hx) xij += lx;
    
              c2.nplane(xij,yij,zij,j);
            }
          }
          c2.face_areas(fs);
          c2.neighbors(neigh);
          cell = &c2;
        }
  
        vol = cell->volume();
        cell->face_freq_table(ff);
        int nn = ff.size()-1;
        for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
      
        // refine the voronoi cell if asked by skipping ultra short edges
        if ((job == 3 || job == 4) && edge_min > 0.){
          std::vector<double> vpos;
          std::vector<int>    vlst;
          double lcut2 = cell->total_edge_distance()*edge_min;
          lcut2 = lcut2*lcut2;
  
          cell->vertices(vpos);
          cell->face_vertices(vlst);
  
          int nf = fs.size();
          int ford[nf];
          int k = 0, iface = 0;
          while (k < vlst.size()){
            int ned = vlst[k++];
            int nuc = 0;
            for (int ii=0; ii<ned; ii++){
              int jj = (ii+1)%ned;
              int v1 = vlst[k+ii], v2 = vlst[k+jj];
              double dx = vpos[v1*3]   - vpos[v2*3];
              double dy = vpos[v1*3+1] - vpos[v2*3+1];
              double dz = vpos[v1*3+2] - vpos[v2*3+2];
              double r2 = dx*dx+dy*dy+dz*dz;
              if ((flag_out & OutEdge) && (v1 > v2)) fprintf(fpedge, "%lg\n", sqrt(r2)/cell->total_edge_distance());
              if (r2 <= lcut2) nuc++;
            }
            ford[iface++] = ned - nuc;
            k += ned;
          }
  
          for (int i=3; i<7; i++) index[i] = 0;
          for (int i=0; i<nf; i++){
            if (ford[i] < 7) index[ford[i]] += 1;
          }
        }
  
        // output voro index info
        int nf = fs.size();
        double wf = double(index[5])/double(nf)*100.;
        fprintf(fp,"%d %d %lg %lg %lg %lg %d,%d,%d,%d %g %d", id, one->attyp[id],
        x, y, z, vol, index[3], index[4], index[5], index[6], wf, nf);
        for (int i=0; i<nf; i++) fprintf(fp," %d", neigh[i]);
        for (int i=0; i<nf; i++) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
  
      } while (cl.inc());

    } else {  // orthogonal box

      voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);
  
      // put atoms into the container
      for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
  
      // loop over all particles and compute their voronoi cell
      voro::voronoicell_neighbor c1, c2, *cell;
      voro::c_loop_all cl(con);
      if (cl.start()) do if (con.compute_cell(c1,cl)){
        int id;
        double x, y, z, vol;
        std::vector<int> ff;       // face_freq
        std::vector<int> neigh;    // neigh list
        std::vector<double> fs;    // face areas
        int index[7];
        for (int i=0; i<7; i++) index[i] = 0;
         
        cl.pos(x,y,z);
        id = cl.pid();
        c1.neighbors(neigh);
        c1.face_areas(fs);
        cell = &c1;
  
        // refine the voronoi cell if asked by removing tiny surfaces
        if ((job == 2 || job == 4) && surf_min > 0.){
          c2.init(-lx,lx,-ly,ly,-lz,lz);
    
          int nf = fs.size();
          // sort neighbors by area if asked to keep a minimum # of neighbors
          if (nminnei > 0){
            for (int i=0; i<nf; i++)
            for (int j=i+1; j<nf; j++){
              if (fs[j] > fs[i]){
                double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
                int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
              }
            }
          }
  
          // add condition on surface
          double fcut = surf_min * cell->surface_area();
          for (int i=0; i<nf; i++){
            if (flag_out & OutSurf) fprintf(fpsurf, "%lg\n", fs[i]/cell->surface_area());
            if (i < nminnei || fs[i] > fcut){
              int j = neigh[i];
    
              // apply pbc
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
          c2.face_areas(fs);
          c2.neighbors(neigh);
          cell = &c2;
        }
  
        vol = cell->volume();
        cell->face_freq_table(ff);
        int nn = ff.size()-1;
        for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
      
        // refine the voronoi cell if asked by skipping ultra short edges
        if ((job == 3 || job == 4) && edge_min > 0.){
          std::vector<double> vpos;
          std::vector<int>    vlst;
          double lcut2 = cell->total_edge_distance()*edge_min;
          lcut2 = lcut2*lcut2;
  
          cell->vertices(vpos);
          cell->face_vertices(vlst);
  
          int nf = fs.size();
          int ford[nf];
          int k = 0, iface = 0;
          while (k < vlst.size()){
            int ned = vlst[k++];
            int nuc = 0;
            for (int ii=0; ii<ned; ii++){
              int jj = (ii+1)%ned;
              int v1 = vlst[k+ii], v2 = vlst[k+jj];
              double dx = vpos[v1*3]   - vpos[v2*3];
              double dy = vpos[v1*3+1] - vpos[v2*3+1];
              double dz = vpos[v1*3+2] - vpos[v2*3+2];
              double r2 = dx*dx+dy*dy+dz*dz;
              if ((flag_out & OutEdge) && (v1 > v2)) fprintf(fpedge, "%lg\n", sqrt(r2)/cell->total_edge_distance());
              if (r2 <= lcut2) nuc++;
            }
            ford[iface++] = ned - nuc;
            k += ned;
          }
  
          for (int i=3; i<7; i++) index[i] = 0;
          for (int i=0; i<nf; i++){
            if (ford[i] < 7) index[ford[i]] += 1;
          }
        }
  
        // output voro index info
        int nf = fs.size();
        double wf = double(index[5])/double(nf)*100.;
        fprintf(fp,"%d %d %lg %lg %lg %lg %d,%d,%d,%d %g %d", id, one->attyp[id],
        x, y, z, vol, index[3], index[4], index[5], index[6], wf, nf);
        for (int i=0; i<nf; i++) fprintf(fp," %d", neigh[i]);
        for (int i=0; i<nf; i++) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
  
      } while (cl.inc());
    }

    fclose(fp);
    printf("Frame %d done, voro info written to: %s\n", img+1, ptr);
  }

  if (flag_out & OutSurf) fclose(fpsurf);
  if (flag_out & OutEdge) fclose(fpedge);

  delete []prefix;
  for (int i=0; i<20; i++) printf("===="); printf("\n");
return;
}
