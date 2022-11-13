#include "driver.h"
#include <complex>
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the static structure factor of the system
 *----------------------------------------------------------------------------*/
void Driver::strfac()
{
  char str[MAXLINE], *ptr;
  double dq[3], qmax, rdq;
  int nq[3], nbin = 101;

  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("    Static  Structure  Factor   ");
  for (int i = 0; i < 6; ++i) printf("====");

  // ask for the max value of k's
  qmax = 10.;
  printf("\nPlease input the upper bound of the q-vectors [%g]: ", qmax);
  input->read_stdin(str);
  if (count_words(str) >= 1){
    ptr = strtok(str, " \n\t\r\f");
    qmax = atof(ptr);
  }
  
  // ask for the # of q-points along each direction
  nq[0] = nq[1] = nq[2] = 31;

  printf("\nThe computation for 3D is rather expensive, be patient if you use a large q-mesh.\n");
  printf("Please input the # of q-points along each direction [%d %d %d]: ", nq[0], nq[1], nq[2]);
  input->read_stdin(str);
  if (count_words(str) >= 3){
    ptr = strtok(str, " \n\t\r\f");
    for (int i = 0; i < 2; ++i){
      nq[i] = abs(atoi(ptr));
      ptr = strtok(NULL, " \n\t\r\f");
    }
    nq[2] = abs(atoi(ptr));
  }

  // ask for the # of bins for S(q)
  printf("Please input the # of bins for S(q) [%d]: ", nbin);
  input->read_stdin(str);
  if (count_words(str) > 0){
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) nbin = atoi(ptr);
    nbin = MAX(2,nbin);
  }
  rdq = double(nbin-1)/qmax;
  
  for (int i = 0; i < 3; ++i) dq[i] = qmax/MAX(1.,double(nq[i]-1));

  // selection commands for atoms
  one = all[istr];
  char selcmd[MAXLINE];
  while (1){
    printf("\nPlease input the atom selection command, `h` for help [all]: ");
    input->read_stdin(str);
    if (count_words(str) > 0){
      strcpy(selcmd, str);
      ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(selcmd,"all\n");
  
    // check the selection command on the first frame
    one->selection(selcmd); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      input->read_stdin(str);
      if (count_words(str) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  // S(kx,ky,kz) accumulator
  double ***sqall;
  memory->create(sqall, 2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "sqall");
  for (int ix = 0; ix < 2*nq[0]+1; ++ix)
  for (int iy = 0; iy < 2*nq[1]+1; ++iy)
  for (int iz = 0; iz <   nq[2]+1; ++iz) sqall[ix][iy][iz] = 0.;

  // prepare for loop
  int nused = 0;
  const std::complex<double> I0 = std::complex<double>(0,-1.);

  // timer
  Timer * timer = new Timer();
  printf("\nComputing, takes some time, please be patient...\n");

  double q[3];
  // loop over all frames
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // select atoms as source
    one->selection(selcmd);
    if (min_mem) one->FreeVoro();
    if (one->nsel < 1) continue;

    printf("  Now to process frame %d... ", img+1); fflush(stdout);

    // cartesian coordinate needed
    one->dir2car();

    double inv_n = 1./double(one->nsel);
    // loop over q-point
#pragma omp parallel for default(shared)  private(q)
    for (int qx = -nq[0]; qx <= nq[0]; ++qx)
    for (int qy = -nq[1]; qy <= nq[1]; ++qy)
    for (int qz =      0; qz <= nq[2]; ++qz){
      int ix = qx + nq[0];
      int iy = qy + nq[1];
      int iz = qz;
      q[0] = qx*dq[0];
      q[1] = qy*dq[1];
      q[2] = qz*dq[2];
 
      std::complex<double> qrnow = std::complex<double>(0.,0.);
      // loop over all atoms
      for (int id = 1; id <= one->natom; ++id){
        if (one->atsel[id] == 0) continue;

        double qr = q[0]*one->atpos[id][0] + q[1]*one->atpos[id][1] + q[2]*one->atpos[id][2];
        qrnow += exp(I0*qr);
      }
      sqall[ix][iy][iz] += real(conj(qrnow)*qrnow) * inv_n;
    }

    ++nused;
    printf("Done! Time used: %g seconds.\n", timer->sincelast());
  }
  timer->stop();
  printf("Total CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  if (nused < 1) return;

  double inv_n = 1./double(nused);
  for (int ix = 0; ix < 2*nq[0]+1; ++ix)
  for (int iy = 0; iy < 2*nq[1]+1; ++iy)
  for (int iz = 0; iz <   nq[2]+1; ++iz) sqall[ix][iy][iz] *= inv_n;

  // output the result, and compute S(q)
  double **sq;
  memory->create(sq,  nbin, 2, "sq");
  for (int i = 0; i < nbin; ++i) sq[i][0] = sq[i][1]  = 0.;

  printf("\nIf you want to output S(kx,ky,kz), input a file name, enter to skip: ");
  input->read_stdin(str);
  if (count_words(str) > 0){
    ptr = strtok(str, " \n\t\r\f");
    FILE *fp = fopen(ptr,"w");
    fprintf(fp,"# kx ky kz  S(k)\n");

    for (int qx = -nq[0]; qx <= nq[0]; ++qx)
    for (int qy = -nq[1]; qy <= nq[1]; ++qy){
      for (int qz = -nq[2]; qz < 0; ++qz){
        int ix = -qx + nq[0];
        int iy = -qy + nq[1];
        int iz = -qz;
        q[0] = qx*dq[0];
        q[1] = qy*dq[1];
        q[2] = qz*dq[2];

        fprintf(fp,"%lg %lg %lg %lg\n", q[0], q[1], q[2], sqall[ix][iy][iz]);
      }
      for (int qz = 0; qz <= nq[2]; ++qz){
        int ix = qx + nq[0];
        int iy = qy + nq[1];
        int iz = qz;
        q[0] = qx*dq[0];
        q[1] = qy*dq[1];
        q[2] = qz*dq[2];

        fprintf(fp,"%lg %lg %lg %lg\n", q[0], q[1], q[2], sqall[ix][iy][iz]);
      }
    }
    fclose(fp);
  }

  // now to compute s(q)
  for (int qx = -nq[0]; qx <= nq[0]; ++qx)
  for (int qy = -nq[1]; qy <= nq[1]; ++qy){
    for (int qz = -nq[2]; qz < 0; ++qz){
      int ix = -qx + nq[0];
      int iy = -qy + nq[1];
      int iz = -qz;
      q[0] = qx*dq[0];
      q[1] = qy*dq[1];
      q[2] = qz*dq[2];
  
      int ibin = int(sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])*rdq + 0.5);
      if (ibin < nbin){
        sq[ibin][0] += sqall[ix][iy][iz];
        sq[ibin][1] += 1.;
      }
    }
    for (int qz = 0; qz <= nq[2]; ++qz){
      int ix = qx + nq[0];
      int iy = qy + nq[1];
      int iz = qz;
      q[0] = qx*dq[0];
      q[1] = qy*dq[1];
      q[2] = qz*dq[2];
  
      int ibin = int(sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])*rdq + 0.5);
      if (ibin < nbin){
        sq[ibin][0] += sqall[ix][iy][iz];
        sq[ibin][1] += 1.;
      }
    }
  }

  printf("Please input the file name to output S(q) [sq.dat]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL){
     strcpy(str, "sq.dat");
     ptr = strtok(str, " \n\t\r\f");
  }
  char *fname = new char [strlen(ptr)+1];
  strcpy(fname, ptr);
  ConfirmOverwrite(fname);
  FILE *fp = fopen(fname,"w");
  fprintf(fp,"# q S(q)\n");

  dq[0] = 1./rdq;
  q[0]  = dq[0];
  for (int i = 1; i < nbin; ++i){
    fprintf(fp,"%lg %lg\n", q[0], sq[i][0]/MAX(1.,sq[i][1]));
    q[0] += dq[0];
  }
  fclose(fp);
  
  memory->destroy(sqall);
  memory->destroy(sq);
  printf("\n%d images were used in the evaluation of S(k), which is written to %s\n", nused, fname);
  delete []fname;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
