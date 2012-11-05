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
  printf("  1. Overall CSRO based on direct Voronoi info;\n");
  printf("  2. Overall CSRO based on refined Voronoi info, skip tiny surfaces;\n");
  printf("  3. Per atom CSRO based on refined Voronoi info;");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 3){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  // refinement info
  double surf_min = 1.e-4;
  int refine = 0;
  if (job >= 2){
    printf("Please input your criterion for tiny surfaces [%g]: ", surf_min);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) surf_min = atof(ptr);
    printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", surf_min);

    if (surf_min > 0.) refine = 1;
  }

  one = all[istr];
  int ntype = one->ntype;

  // output file name for per atom CSRO
  FILE *fp; fp = NULL;
  if (job == 3){
    printf("Please input the output file name [csro.dat]: ");
    if (count_words(fgets(str,MAXLINE, stdin)) < 1) strcpy(str,"csro.dat");
    ptr = strtok(str, " \n\t\r\f");
    fp = fopen(ptr, "w");
  }

  // work spaces
  bigint *NumType, **NumNei, *nnei;
  nnei   = memory->create(nnei, ntype+1, "csro:nnei"); // per atom info
  NumType= memory->create(NumType,ntype+1,"csro:NumType");
  NumNei = memory->create(NumNei,ntype+1,ntype+1,"csro:NumNei");
  for (int i=0; i<= ntype; i++){
    NumType[i] = 0;
    for (int j=0; j<= ntype; j++) NumNei[i][j] = 0;
  }
  // chemical composition for each frame
  double *cc;
  cc = new double [ntype+1];

  // now to do the real job
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // ntype of different frames must be the same
    if (one->ntype != ntype) continue;

    if (job == 3){
      fprintf(fp, "# Per atom CSRO info for frame %d; istep = %d\n# id type", img+1, one->tstep);
      for (int jp=1; jp<=ntype; jp++) fprintf(fp," csro-%d", jp);
      fprintf(fp,"\n");
    }

    // set local variables
    double xlo = one->xlo, xhi = one->xhi;
    double ylo = one->ylo, yhi = one->yhi;
    double zlo = one->zlo, zhi = one->zhi;
    double lx  = one->lx,  ly  = one->ly,  lz = one->lz;
    double xy  = one->xy,  xz  = one->xz,  yz = one->yz;

    int n = one->natom;
    int *attyp = one->attyp;
    double **atpos = one->atpos;
    for (int i=1; i<=ntype; i++) NumType[i] += one->numtype[i];
    for (int i=1; i<=ntype; i++) cc[i] = double(one->numtype[i])/double(one->natom);

    // need cartesian coordinates
    one->dir2car();

    // compute optimal size for container, and then contrust it
    double l = pow(double(n)/(5.6*lx*ly*lz), 1./3.);
    int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);

    if (one->triclinic){

      voro::container_periodic con(lx,xy,ly,xz,yz,lz,nx,ny,nz,8);
  
      // put atoms into the container
      for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
  
      // loop over all particles and compute their voronoi cell
      voro::voronoicell_neighbor cell;
      voro::c_loop_all_periodic cl(con);
      if (cl.start()) do if (con.compute_cell(cell,cl)){
        int id = cl.pid();
        int ip = one->attyp[id];
  
        std::vector<int> neigh;    // neigh list
        cell.neighbors(neigh);
        int nf = neigh.size();
  
        // zero per atom info
        for (int jp=1; jp<= ntype; jp++) nnei[jp] = 0;

        if (refine) { // based on refined voro info
          double fcut = surf_min * cell.surface_area();
          std::vector<double> fs;    // face areas
          cell.face_areas(fs);
          for (int ii=0; ii<nf; ii++){
            if (fs[ii] < fcut) continue;
            int jd = neigh[ii];
            int jp = one->attyp[jd];
            nnei[jp]++;
          }
  
        } else { // based on direct voro info
          for (int ii=0; ii<nf; ii++){
            int jd = neigh[ii];
            int jp = one->attyp[jd];
            nnei[jp]++;
          }
        }
  
        // accumulate overall info
        for (int jp=1; jp<= ntype; jp++) NumNei[ip][jp] += nnei[jp];

        // output per atom info if needed
        if (job == 3){
          fprintf(fp,"%d %d", id, ip);
          int nn = 0;
          for (int jp=1; jp<= ntype; jp++) nn += nnei[jp];
          for (int jp=1; jp<= ntype; jp++) fprintf(fp," %g", 1.-double(nnei[jp])/(double(nn)*cc[jp]));
          fprintf(fp,"\n");
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
        int ip = one->attyp[id];
  
        std::vector<int> neigh;    // neigh list
        cell.neighbors(neigh);
        int nf = neigh.size();
  
        // zero per atom info
        for (int jp=1; jp<= ntype; jp++) nnei[jp] = 0;

        if (refine) { // based on refined voro info
          double fcut = surf_min * cell.surface_area();
          std::vector<double> fs;    // face areas
          cell.face_areas(fs);
          for (int ii=0; ii<nf; ii++){
            if (fs[ii] < fcut) continue;
            int jd = neigh[ii];
            int jp = one->attyp[jd];
            nnei[jp]++;
          }
  
        } else { // based on direct voro info
          for (int ii=0; ii<nf; ii++){
            int jd = neigh[ii];
            int jp = one->attyp[jd];
            nnei[jp]++;
          }
        }
  
        // accumulate overall info
        for (int jp=1; jp<= ntype; jp++) NumNei[ip][jp] += nnei[jp];

        // output per atom info if needed
        if (job == 3){
          fprintf(fp,"%d %d", id, ip);
          int ntot = 0;
          for (int jp=1; jp<= ntype; jp++) ntot += nnei[jp]; ntot = MAX(1,ntot);
          for (int jp=1; jp<= ntype; jp++) fprintf(fp," %g", 1.-double(nnei[jp])/(double(ntot)*cc[jp]));
          fprintf(fp,"\n");
        }

      } while (cl.inc());
    }
  }
  
  if (job == 3) fclose(fp);

  bigint ntotal = 0;
  double concentration[ntype+1];
  for (int i=1; i<= ntype; i++) ntotal += NumType[i];
  for (int i=1; i<= ntype; i++) concentration[i] = double(NumType[i])/double(ntotal);

  printf("\n"); for (int i=0; i<20; i++) printf("____"); printf("\nType         :");
  if (type2atnum == NULL){
    for (int ip=1; ip<= ntype; ip++) printf("    %2d    ",ip);
    printf("\n"); for (int i=0; i<20; i++) printf("----");
    printf("\nConcentration:"); for (int ip=1; ip<= ntype; ip++) printf("%10.6f", concentration[ip]);
    printf("\n"); for (int i=0; i<20; i++) printf("____");
    printf("\nCSRO         :"); for (int jp=1; jp<= ntype; jp++) printf("    %2d    ",jp);
  } else {
    char ename[3];
    for (int ip=1; ip<= ntype; ip++){
      element->Num2Name(type2atnum[ip],ename);
      printf("    %2s    ", ename);
    }
    printf("\n"); for (int i=0; i<20; i++) printf("----");
    printf("\nConcentration:"); for (int ip=1; ip<= ntype; ip++) printf("%10.6f", concentration[ip]);
    printf("\n"); for (int i=0; i<20; i++) printf("____"); printf("\nCSRO         :");
    for (int jp=1; jp<= ntype; jp++){
      element->Num2Name(type2atnum[jp],ename);
      printf("    %2s    ", ename);
    }
  }

  printf("\n"); for (int i=0; i<20; i++) printf("----");
  for (int ip=1; ip<= ntype; ip++){
    if (type2atnum == NULL) printf("\n      %2d     :", ip);
    else {
      char ename[3]; element->Num2Name(type2atnum[ip], ename);
      printf("\n      %2s     :", ename);
    }
    ntotal = 0;
    for (int jp=1; jp<= ntype; jp++) ntotal += NumNei[ip][jp];
    for (int jp=1; jp<= ntype; jp++) printf("%10.6f", 1.-double(NumNei[ip][jp])/double(ntotal)/concentration[jp]);
  }
  printf("\n"); for (int i=0; i<20; i++) printf("===="); printf("\n");

  memory->destroy(nnei);
  memory->destroy(NumType);
  memory->destroy(NumNei);
return;
}
