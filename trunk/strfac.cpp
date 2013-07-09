#include "driver.h"
#include <complex>
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the static structure factor of the system
 *
 * Currently it works for orthogonal box only.
 *----------------------------------------------------------------------------*/
void Driver::strfac()
{
  char str[MAXLINE], *ptr;
  double kmax[3], dk[3], qmax, rdq;
  int nk[3], nbin = 101;

  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("    Static  Structure  Factor   ");
  for (int i=0; i<6; i++) printf("====");

  // ask for the max value of k's
  kmax[0] = kmax[1] = kmax[2] = 10.;
  printf("\nPlease input the upper bound of the k-vectors [%g %g %g]: ", kmax[0], kmax[1], kmax[2]);
  if (count_words(fgets(str,MAXLINE, stdin)) >= 3){
    ptr = strtok(str, " \n\t\r\f");
    for (int i=0; i<2; i++){
      kmax[i] = atof(ptr);
      ptr = strtok(NULL, " \n\t\r\f");
    }
    kmax[2] = atof(ptr);
  }
  
  // ask for the # of bins of k along each direction
  nk[0] = nk[1] = nk[2] = 20;
  for (int i=0; i<3; i++) if (abs(kmax[i]) < ZERO) nk[i] = 1;

  printf("\nThe computation for 3D is rather expensive, be patient if you use a large k-mesh.\n");
  printf("Please input the # of k-points along each direction [%d %d %d]: ", nk[0], nk[1], nk[2]);
  if (count_words(fgets(str,MAXLINE, stdin)) >= 3){
    ptr = strtok(str, " \n\t\r\f");
    for (int i=0; i<2; i++){
      nk[i] = atoi(ptr);
      ptr = strtok(NULL, " \n\t\r\f");
    }
    nk[2] = atoi(ptr);
  }
  for (int i=0; i<3; i++){
    if (abs(kmax[i]) < ZERO) nk[i] = 1;
    nk[i] = MAX(1,nk[i]);
  }

  // ask for the # of bins for S(q)
  printf("Please input the # of bins for S(k) [%d]: ", nbin);
  if (count_words(fgets(str,MAXLINE, stdin)) > 0){
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) nbin = atoi(ptr);
    nbin = MAX(2,nbin);
  }
  qmax = MAX(kmax[0], MAX(kmax[1], kmax[2]));
  rdq = double(nbin-1)/qmax;
  
  for (int i=0; i<3; i++) dk[i] = kmax[i]/MAX(1.,double(nk[i]-1));

  // selection commands for atoms
  one = all[istr];
  char selcmd[MAXLINE];
  while (1){
    printf("\nPlease input the atom selection command, `h` for help [all]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(selcmd, str);
      ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(selcmd,"all\n");
  
    // check the selection command on the first frame
    one->selection(selcmd); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  // S(kx,ky,kz) accumulator
  double ***skall;
  memory->create(skall, nk[0], nk[1], nk[2], "skall");
  for (int i= 0; i< nk[0]; i++)
  for (int j= 0; j< nk[1]; j++)
  for (int k= 0; k< nk[2]; k++) skall[i][j][k] = 0.;

  // space for exp(-i*K_a*r_a)
  complex<double> **kxrx, **kyry, **kzrz;
  memory->create(kxrx, nk[0], one->nsel+1, "kxrx");
  memory->create(kyry, nk[1], one->nsel+1, "kyry");
  memory->create(kzrz, nk[2], one->nsel+1, "kzrz");
  int nprev = one->nsel;

  // S(k) is only valid for k > pi/L_min; here we use LMax instead for practical reason
  double LMax = 0.;

  // prepare for loop
  int nused = 0, nnorm = 0;
  const complex<double> I0 = complex<double>(0,1.);

  // timer
  Timer * timer = new Timer();
  printf("\nStructure factor calculation takes time, please be patient.\n");

  // loop over all frames
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // does not work for non-orthogonal box
    if (one->triclinic) continue;

    // select atoms as source
    one->selection(selcmd);
    if (one->nsel < 1) continue;

    printf("  Now to process frame %d... ", img+1); fflush(stdout);

    // factional coordinate needed
    one->car2dir();
    LMax += MAX(one->box[0], MAX(one->box[1], one->box[2]));

    if (one->nsel != nprev){
      memory->destroy(kxrx);
      memory->destroy(kyry);
      memory->destroy(kzrz);

      memory->create(kxrx, nk[0], one->nsel, "kxrx");
      memory->create(kyry, nk[1], one->nsel, "kyry");
      memory->create(kzrz, nk[2], one->nsel, "kzrz");
      nprev = one->nsel;
    }

    // loops over k
    complex<double> dq[3];
    for (int idim=0; idim<3; idim++) dq[idim] = dk[idim]*one->box[idim]*I0;

    int inext = 0;
    for (int ii=1; ii<= one->natom; ii++){
      if (one->atsel[ii] == 0) continue;

      for (int i = 0; i < nk[0]; i++){
        complex<double> kx = double(i) * dq[0];

        kxrx[i][inext] = exp(kx*one->atpos[ii][0]);
      }

      for (int j = 0; j < nk[1]; j++){
        complex<double> ky = double(j) * dq[1];

        kyry[j][inext] = exp(ky*one->atpos[ii][1]);
      }

      for (int k = 0; k < nk[2]; k++){
        complex<double> kz = double(k) * dq[2];

        kzrz[k][inext] = exp(kz*one->atpos[ii][2]);
      }

      inext++;
    }

    for (int i = 0; i < nk[0]; i++)
    for (int j = 0; j < nk[1]; j++)
    for (int k = 0; k < nk[2]; k++){
      complex<double> skone = complex<double>(0.,0.);
      for (int ii = 0; ii < one->nsel; ii++) skone += kxrx[i][ii]*kyry[j][ii]*kzrz[k][ii];

      skall[i][j][k] += real(skone*conj(skone));
    }

    nused++; nnorm += one->nsel;
    printf("Done! Time used: %g seconds.\n", timer->sincelast());
  }
  memory->destroy(kxrx);
  memory->destroy(kyry);
  memory->destroy(kzrz);
  timer->stop();
  printf("Total CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  if (nused < 1) return;

  double fac = 1./double(nnorm);
  // output the result, and compute S(k)
  int *hit; double *Sk;
  memory->create(Sk,  nbin, "Sk");
  memory->create(hit, nbin, "hit");
  for (int i=0; i< nbin; i++) Sk[i]  = 0.;
  for (int i=0; i< nbin; i++) hit[i] = 0;

  printf("\nPlease input the file name to output S(kx,ky,kz) [skv.dat]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "skv.dat");
  ptr = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(ptr,"w");
  fprintf(fp,"# kx ky kz  S(k)\n");
  double q[3];
  for (int i = 0; i < nk[0]; i++){ q[0] = double(i)*dk[0];
  for (int j = 0; j < nk[1]; j++){ q[1] = double(j)*dk[1];
  for (int k = 0; k < nk[2]; k++){ q[2] = double(k)*dk[2];

    double sknow = skall[i][j][k]*fac;
    fprintf(fp,"%g %g %g %lg\n", q[0], q[1], q[2], sknow);

    int ibin = int(sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])*rdq + 0.5);
    if (ibin < nbin){
      Sk[ibin]  += sknow;
      hit[ibin] ++;
    }
  } fprintf(fp,"\n"); }}
  fclose(fp);

  printf("Please input the file name to output S(k) [sk.dat]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "sk.dat");
  ptr = strtok(str, " \n\t\r\f");
  fp = fopen(ptr,"w");
  fprintf(fp,"# k S(k)\n");

  int nmin = int(double(nused)*4.*atan(1.)/LMax*rdq + 0.5);
  double dq = 1./rdq, qr = double(nmin)*dq;
  for (int i=nmin; i<nbin; i++){
    if (hit[i] > 0) fprintf(fp,"%lg %lg\n", qr, Sk[i]/double(hit[i]));
    qr += dq;
  }
  fclose(fp);
  
  memory->destroy(skall);
  memory->destroy(hit);
  memory->destroy(Sk);
  printf("\n%d images were used in the evaluation of S(k), which is written to %s\n", nused, ptr);

  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
