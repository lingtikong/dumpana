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
void Driver::voro(const int job)
{
  // ask for frames
  setrange();

  char str[MAXLINE];
  printf("\n"); for (int i=0; i<20; i++) printf("====");

  if (job == 1){  // compute refined voronoi index info
    // ask for threshold
    double fmin = 0.001;
    printf("\nPlease input a threshold for the evaluation of voro neighbors, if the\n");
    printf("surface area corresponding to an atom takes less ratio of the whole surface\n");
    printf("than the defined threshold, it will not be seen as a neighbor [%lg]: ", fmin);
    fgets(str,MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) fmin = atof(ptr);
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
      //con.clear();
      for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
  
      // open file for output
      sprintf(str,"voro_%d.dat ", one->tstep);
      ptr = strtok(str," \n\t\r\f");
      FILE *fp = fopen(ptr, "w");
      fprintf(fp,"#Box info: %lg %lg %lg %lg %lg %lg %d\n", xlo, xhi, ylo, yhi, zlo, zhi, n);
      fprintf(fp,"# 1  2    3  4  5  6   7         8    9    10\n");
      fprintf(fp,"# id type x  y  z  vol voroindex f5%%  NNei NeiList surfaceareas\n");
      //FILE *fp2 = fopen("c2.gnu", "w");
      // loop over all atoms and compute their voronoi cell
      voro::voronoicell_neighbor c, c2;
      voro::c_loop_all cl(con);
      if (cl.start()) do if (con.compute_cell(c,cl)){
        int id;
        double x, y, z, vol;
        std::vector<int> ff;       // face_freq
        std::vector<int> neigh;    // neigh list
        std::vector<double> fs;    // face areas
        int index[7];
        for (int i=0; i<7; i++) index[i] = 0;
         
        cl.pos(x,y,z);
        id = cl.pid();
        c.neighbors(neigh);
        c.face_areas(fs);
  
        if (fmin > 0.){ // if fmin <= 0., no need to refine
          // initialize a second cell
          c2.init(-lx,lx,-ly,ly,-lz,lz);
    
          // add condition on surface
          double fcut = fmin * c.surface_area();
          int nf = fs.size();
          for (int i=0; i<nf; i++){
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
          c2.face_areas(fs);
          c2.neighbors(neigh);
          c2.face_freq_table(ff);
          vol = c2.volume();

          std::vector<double> vcord;
          std::vector<int>    vlist;
          double lcut = c2.total_edge_distance()*1.e-5;
          c2.vertices(vcord);
          c2.face_vertices(vlist);
          int nv = vcord.size()/3;
          int count[nv][nv];
          for (int i=0; i<nv; i++)
          for (int j=0; j<nv; j++) count[i][j] = 0;

          for (int i=0; i<nv; i++)
          for (int j=i+1; j<nv; j++){
            int ip = i*3;
            int jp = j*3;
            double dx = vcord[jp]-vcord[ip];
            double dy = vcord[jp+1]-vcord[ip+1];
            double dz = vcord[jp+2]-vcord[ip+2];
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            if (r < lcut) count[i][j] = count[j][i] = 1;
          }
          nf = fs.size();
          int ford[nf];
          for (int i=0; i<nf; i++) ford[i] = 0;
          int k = 0, iface = 0;
          while (k < vlist.size()){
            int ned = vlist[k++];
            int nuc = 0;
            for (int ii=0; ii<ned; ii++)
            for (int jj=ii+1; jj<ned; jj++){
              int b1 = vlist[k+ii], b2 = vlist[k+jj];
              nuc += count[b1][b2];
            }
            ford[iface++] = ned - nuc;
            k += ned;
          }
          for (int i=0; i<nf; i++){
            if (ford[i] < 7) index[ford[i]] +=1;
          }

        } else {
          c.face_freq_table(ff);
          vol = c.volume();
        }
          
  
        // output voro index info
        int nn = ff.size()-1;
        //for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
        int nf = fs.size();
        double wf = double(index[5])/double(nf)*100.;
        fprintf(fp,"%d %d %lg %lg %lg %lg %d,%d,%d,%d %g %d", id, one->attyp[id],
        x, y, z, vol, index[3], index[4], index[5], index[6], wf, nf);
        for (int i=0; i<nf; i++) fprintf(fp," %d", neigh[i]);
        for (int i=0; i<nf; i++) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
  
      } while (cl.inc());
  
      fclose(fp);
      printf("Frame %d done, voro info written to: %s\n", img+1, ptr);
    }

  } else if (job == 3){ // compute voro surface areas

    printf("\nPlease input the file name to output the voronoi surface areas[voro_areas.dat]: ");
    if (count_words(fgets(str,MAXLINE, stdin)) < 1) strcpy(str,"voro_areas.dat");

    char *ptr = strtok(str, " \n\t\r\f");
    FILE *fp = fopen(ptr,"w");

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
      for (int i=1; i<= n; i++) con.put(i,one->atpos[i][0],one->atpos[i][1],one->atpos[i][2]);
  
      // loop over all atoms and compute their voronoi cell
      voro::voronoicell_neighbor c;
      voro::c_loop_all cl(con);
      if (cl.start()) do if (con.compute_cell(c,cl)){
        std::vector<double> fs;
        c.face_areas(fs);
  
        for (int i=0; i<(signed int) fs.size(); i++) fprintf(fp,"%lg\n", fs[i]);
      } while (cl.inc());
    }
    fclose(fp);
    printf("Job done, the results are written to: %s\n", ptr);

  } else if (job == 2){ // use qhull to compute the voro index
    // ask for threshold
    double fmin = 0.001;
    printf("\nPlease input a threshold for the evaluation of voro neighbors, if the\n");
    printf("surface area corresponding to an atom takes less ratio of the whole surface\n");
    printf("than the defined threshold, it will not be seen as a neighbor [%lg]: ", fmin);
    fgets(str,MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) fmin = atof(ptr);
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
      //con.clear();
      for (int i=1; i<= n; i++) con.put(i,one->atpos[i][0],one->atpos[i][1],one->atpos[i][2]);
  
      // open file for output
      sprintf(str,"voro_%d.dat ", one->tstep);
      ptr = strtok(str," \n\t\r\f");
      FILE *fp = fopen(ptr, "w");
      fprintf(fp,"#Box info: %lg %lg %lg %lg %lg %lg %d\n", xlo, xhi, ylo, yhi, zlo, zhi, n);
      fprintf(fp,"# 1  2    3  4  5  6   7         8    9    10\n");
      fprintf(fp,"# id type x  y  z  vol voroindex f5%%  NNei NeiList surfaceareas\n");

      // loop over all atoms and compute their voronoi cell
      voro::voronoicell_neighbor c;
      voro::c_loop_all cl(con);
      if (cl.start()) do if (con.compute_cell(c,cl)){
        int id;
        double x, y, z, vol;
        std::vector<int> ff;       // face_freq
        std::vector<int> neigh;    // neigh list
        std::vector<double> fs;    // face areas
         
        cl.pos(x,y,z);
        id = cl.pid();
        c.neighbors(neigh);
        c.face_areas(fs);
  
        int index[7];
        for (int i=3; i<7; i++) index[i] = 0;

        if (fmin > 0.){ // if fmin <= 0., no need to refine
          std::vector<double> dx, dy, dz;
          dx.clear(); dy.clear(); dz.clear();
          dx.push_back(0.); dy.push_back(0.); dz.push_back(0.);
          // add condition on surface
          double fcut = fmin * c.surface_area();
          int nf = fs.size();
          for (int i=0; i<nf; i++){
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
    
              dx.push_back(xij); dy.push_back(yij); dz.push_back(zij);
            }
          }
          FILE *fpt = fopen(".qconvex_input_tmp","w");
          fprintf(fpt,"3 %d\n", (int)dx.size());
          for (int i=0; i<dx.size(); i++) fprintf(fpt,"%lg %lg %lg\n", dx[i], dy[i],dz[i]);
          fclose(fpt);
          //system("qvoronoi p QzV0 < .qconvex_input_tmp > .qconvex_output_intermediate");
          //fgets(str,MAXLINE,fpt); int ns = atoi(strtok(str, " \n\t\r\f"));

          system("qvoronoi p QzV0 < .qconvex_input_tmp |qconvex Fv FS > .qconvex_output_tmp");
          fpt = fopen(".qconvex_output_tmp", "r");
          fgets(str,MAXLINE,fpt); int ns = atoi(strtok(str, " \n\t\r\f"));
          for (int i=0; i<ns; i++){
            fgets(str,MAXLINE,fpt);
            int ned = atoi(strtok(str, " \n\t\r\f"));
            if (ned < 7) index[ned]++;
          }
          fgets(str,MAXLINE,fpt); fgets(str,MAXLINE,fpt);
          strtok(str, " \n\t\r\f");
          strtok(NULL," \n\t\r\f");
          vol = atof(strtok(NULL," \n\t\r\f"));
          fclose(fpt);
        } else {
          c.face_freq_table(ff);
          vol = c.volume();
          int nn = ff.size()-1;
          for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
          for (int i=MIN(6,nn)+1; i<= 6; i++) index[i] = 0;
        }
  
        // output voro index info
        int nf = fs.size();
        double wf = double(index[5])/double(nf)*100.;
        fprintf(fp,"%d %d %lg %lg %lg %lg %d,%d,%d,%d %g %d",id,one->attyp[id],
        x, y, z, vol, index[3], index[4], index[5], index[6], wf, nf);
        for (int i=0; i<nf; i++) fprintf(fp," %d", neigh[i]);
        for (int i=0; i<nf; i++) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
  
      } while (cl.inc());
  
      fclose(fp);
      printf("Frame %d done, voro info written to: %s\n", img+1, ptr);
    }
  }

  for (int i=0; i<20; i++) printf("====");  printf("\n");
return;
}
/*----------------------------------------------------------------------------*/
