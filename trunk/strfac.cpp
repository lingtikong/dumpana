#include "driver.h"
#include "math.h"
#include <complex>

#define MAXLINE 512
#define ZERO 1.e-8
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to compute the static structure factor of the system
 *
 * Currently it works for orthogonal box only.
 *----------------------------------------------------------------------------*/
void Driver::strfac()
{
  char str[MAXLINE], header[MAXLINE], *ptr;
  double kmax[3], dk[3], qmax, rdq;
  int nk[3], nbin = 101;

  kmax[0] = kmax[1] = kmax[2] = 15.;
  printf("Please input the upper bound of the k-vectors [15 15 15]: ");
  if (count_words(fgets(str,MAXLINE, stdin)) >= 3){
    ptr = strtok(str, " \n\t\r\f");
    for (int i=0; i<2; i++){
      kmax[i] = atof(ptr);
      ptr = strtok(NULL, " \n\t\r\f");
    }
    kmax[2] = atof(ptr);
  }
  
  nk[0] = nk[1] = nk[2] = 15;
  for (int i=0; i<3; i++) if (abs(kmax[i]) < ZERO) nk[i] = 1;

  printf("\nThe computation for 3D is rather expensive, be patient if you use a large k-mesh.\n");
  printf("Please input the # of k-points along each direction [%d %d %d]: ", nk[0], nk[1], nk[2]);
  if (count_words(fgets(str,MAXLINE, stdin)) >= 3){
    ptr = strtok(str, " \n\t\r\f");
    for (int i=0; i<2; i++){
      nk[i] = atoi(ptr)-1;
      ptr = strtok(NULL, " \n\t\r\f");
    }
    nk[2] = atoi(ptr)-1;
  }
  for (int i=0; i<3; i++){
    if (abs(kmax[i]) < ZERO) nk[i] = 0;
    nk[i] = MAX(0,nk[i]);
  }

  printf("Please input the # of bins for S(k) [%d]: ", nbin);
  if (count_words(fgets(str,MAXLINE, stdin)) > 0){
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) nbin = atoi(ptr);
    nbin = MAX(2,nbin);
  }
  qmax = MAX(kmax[0], MAX(kmax[1], kmax[2]));
  rdq = double(nbin-1)/qmax;
  
  for (int i=0; i<3; i++) dk[i] = kmax[i]/MAX(1.,double(nk[i]));

  int *hit;
  double ***skall, *Sk;
  skall = memory->create(skall, 2*nk[0]+1, 2*nk[1]+1, 2*nk[2]+1, "skall");
  Sk  = memory->create(Sk,  nbin, "Sk");
  hit = memory->create(hit, nbin, "hit");
  for (int i= -nk[0]; i<= nk[0]; i++)
  for (int j= -nk[1]; j<= nk[1]; j++)
  for (int k= -nk[2]; k<= nk[2]; k++) skall[i+nk[0]][j+nk[1]][k+nk[2]] = 0.;

  for (int i=0; i< nbin; i++) Sk[i]  = 0.;
  for (int i=0; i< nbin; i++) hit[i] = 0;
  
  one = all[0];
  char selcmd[MAXLINE];
  
  // selection commands for atoms
  // atoms as source
  while (1){
    printf("\nPlease input the selection command for atoms, `h` for help [all]:");
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

  complex<double> **kxrx, **kyry, **kzrz;
  kxrx = memory->create(kxrx, one->nsel, 2*nk[0]+1, "kxrx");
  kyry = memory->create(kyry, one->nsel, 2*nk[1]+1, "kyry");
  kzrz = memory->create(kzrz, one->nsel, 2*nk[2]+1, "kzrz");
  int nprev = one->nsel;

  int nused = 0, nnorm = 0;
  const complex<double> I0 = complex<double>(0,1.);

  char flag[4];
  flag[0] = '-'; flag[1] = '\\'; flag[2] = '|'; flag[3] = '/';
  printf("\nComputing, it takes time ... "); fflush(stdout);

  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // select atoms as source
    one->selection(selcmd);
    if (one->nsel < 1) continue;

    printf("\b%c", flag[nused%4]); fflush(stdout);

    // cartesian coordinate needed
    one->car2dir();

    if (one->nsel != nprev){
      kxrx = memory->grow(kxrx, one->nsel, 2*nk[0]+1, "kxrx");
      kyry = memory->grow(kyry, one->nsel, 2*nk[1]+1, "kyry");
      kzrz = memory->grow(kzrz, one->nsel, 2*nk[2]+1, "kzrz");
      nprev = one->nsel;
    }

    // loops over k
    int inext = 0;
    complex<double> dq[3];
    for (int idim=0; idim<3; idim++) dq[idim] = dk[idim]*one->box[idim]*I0;

    for (int ii=1; ii<= one->natom; ii++){
      if (one->atsel[ii] == 0) continue;

      for (int ix = -nk[0]; ix <= nk[0]; ix++){
        complex<double> kx = double(ix) * dq[0];
        int x = ix + nk[0];

        kxrx[inext][x] = exp(-kx*one->atpos[ii][0]);
      }

      for (int iy = -nk[1]; iy <= nk[1]; iy++){
        complex<double> ky = double(iy) * dq[1];
        int y = iy + nk[1];

        kyry[inext][y] = exp(-ky*one->atpos[ii][1]);
      }

      for (int iz = -nk[2]; iz <= nk[2]; iz++){
        complex<double> kz = double(iz) * dq[2];
        int z = iz + nk[2];

        kzrz[inext][z] = exp(-kz*one->atpos[ii][2]);
      }

      inext++;
    }

    for (int ix = -nk[0]; ix <= nk[0]; ix++){ int x = ix + nk[0];
    for (int iy = -nk[1]; iy <= nk[1]; iy++){ int y = iy + nk[1];
    for (int iz = -nk[2]; iz <= nk[2]; iz++){ int z = iz + nk[2];
      complex<double> skone = complex<double>(0.,0.);
      for (int i=0; i< one->nsel; i++) skone += kxrx[i][x]*kyry[i][y]*kzrz[i][z];
      skall[x][y][z] += real(skone * conj(skone));
    }}}

    nused++; nnorm += one->nsel;
  }
  memory->destroy(kxrx);
  memory->destroy(kyry);
  memory->destroy(kzrz);

  double fac = 1./double(nnorm);
  // output the result
  printf("\bDone!\n\nPlease input the file to output S(kx,ky,kz) [Skv.dat]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "Skv.dat");
  ptr = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(ptr,"w");
  fprintf(fp,"# kx ky kz  S(k)\n");
  for (int ix = -nk[0]; ix <= nk[0]; ix++){
  for (int iy = -nk[1]; iy <= nk[1]; iy++){
  for (int iz = -nk[2]; iz <= nk[2]; iz++){
    int x = ix + nk[0];
    int y = iy + nk[1];
    int z = iz + nk[2];
    double q[3];
    q[0] = ix*dk[0]; q[1] = iy*dk[1]; q[2] = iz*dk[2];
    double sknow = skall[x][y][z]*fac;

    fprintf(fp,"%g %g %g %lg\n", q[0], q[1], q[2], sknow);

    int ibin = int(sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])*rdq + 0.5);
    if (ibin < nbin){ Sk[ibin] += sknow; hit[ibin]++;}
  } fprintf(fp,"\n"); }}
  fclose(fp);

  printf("Please input the file to output S(k) [Sk.dat]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "Sk.dat");
  ptr = strtok(str, " \n\t\r\f");
  fp = fopen(ptr,"w");
  fprintf(fp,"# k S(k)\n");
  double dq = 1./rdq, q = dq;
  for (int i=1; i<nbin; i++){
    if (hit[i] > 0) fprintf(fp,"%lg %lg\n", q, Sk[i]/double(hit[i]));
    else fprintf(fp,"%lg %lg\n", q, 0.);
    q += dq;
  }
  fclose(fp);
  
  memory->destroy(skall);
  memory->destroy(hit);
  memory->destroy(Sk);
  printf("\n%d images were used in the evaluation of S(k), which is written to %s\n", nused, ptr);

return;
}
