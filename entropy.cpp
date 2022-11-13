#include "driver.h"
#include "math.h"

/*------------------------------------------------------------------------------
 * Method to compute the configurational entropy of mixing.
 *
 * Ref:
 *    JY Qin, Acta Phys.-Chim. Sin. 28(7):1586-1592, 2012.
 *----------------------------------------------------------------------------*/
void Driver::compute_smix()
{
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 4; ++i) printf("====");
  printf(" Configurational entropy of mixing calculation  ");
  for (int i = 0; i < 4; ++i) printf("===="); printf("\n");

  // voro refinement info
  set_cutoffs(0);

  // output info
  printf("\nPlease input the file to output the entropy info [smix.dat]: ");
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "smix.dat");
  ptr = strtok(str, " \n\t\r\f");
  char *fname = new char [strlen(ptr)+1];
  strcpy(fname, ptr);
  ConfirmOverwrite(fname);
  FILE *fp = fopen(fname,"w");
  fprintf(fp,"# Configurational entropy of mixing, in unit of kB/mole.\n");
  fprintf(fp,"# Ref: JY Qin, Acta Phys.-Chim. Sin. 28(7):1586-1592, 2012.\n# MD-step smix\n");

  double smix_tot = 0.;
  int nused = 0;
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    if (one->flag_smix == 0){
      one->ComputeVoro(mins);

      int **nei_count;
      memory->create(nei_count, one->ntype+1, one->ntype+2, "nei_count");
      for (int ip = 0; ip <= one->ntype;   ++ip)
      for (int jp = 0; jp <= one->ntype+1; ++jp) nei_count[ip][jp] = 0;

#pragma omp parallel for default(shared)
      for (int id = 1; id <= one->natom; ++id){
        int ip = one->attyp[id];
        int ni = one->neilist[0][id];
        for (int jj = 1; jj <= ni; ++jj){
          int jd = one->neilist[jj][id];
          int jp = one->attyp[jd];

#pragma omp atomic
          ++nei_count[ip][jp];
        }
#pragma omp atomic
        ++nei_count[ip][0];
#pragma omp atomic
        nei_count[ip][one->ntype+1] += ni;
      }

      double *ni, *xi, **fij;
      ni = new double[one->ntype+1];
      xi = new double[one->ntype+1];
      memory->create(fij, one->ntype+1, one->ntype+1, "fij");
      for (int ip = 1; ip <= one->ntype; ++ip){
        ni[ip] = double(nei_count[ip][one->ntype+1])/double(MAX(1,nei_count[ip][0]));
        xi[ip] = double(nei_count[ip][0])/double(one->natom);

        for (int jp = 1; jp <= one->ntype; ++jp) fij[ip][jp] = double(nei_count[ip][jp])/double(MAX(1,nei_count[ip][one->ntype+1]));
      }

      one->smix = 0.;
      double nxi = 0.;
      for (int jp = 1; jp <= one->ntype; ++jp){
        nxi += ni[jp] * xi[jp];

        for (int kp = 1; kp <= one->ntype; ++kp){
          double nxf_lk = 0.;
          for (int lp = 1; lp <= one->ntype; ++lp) nxf_lk += ni[lp] * xi[lp] * fij[lp][kp];
          double nxf_jk = ni[jp] * xi[jp] * fij[jp][kp];
          one->smix += nxf_jk * log(nxf_lk / nxf_jk);
        }
      }
      one->smix /= nxi;

      delete []ni;
      delete []xi;
      memory->destroy(fij);
      memory->destroy(nei_count);
    }

    smix_tot += one->smix;
    fprintf(fp, "%d %lg\n", one->tstep, one->smix);
    if (min_mem) one->FreeVoro();

    ++nused;
  }
  smix_tot /= double(MAX(1, nused));
  fprintf(fp, "# Average: %lg kB/mole\n", smix_tot);
  fclose(fp);
  
  printf("\n%d images were used in computing smix (ave: %g kB), written to %s\n", nused, smix_tot, fname);
  delete []fname;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------*/
