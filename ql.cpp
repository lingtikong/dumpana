#include "spherical.h"
#include "driver.h"

/*------------------------------------------------------------------------------
 * Driver to compute the local order parameber based on the Spherical Haromics.
 * Ref:
 * 1)  Phys. Rev. B 28, 784 (1983).
 * 2)  J Chem Phys 122:214722 (2005).
 *----------------------------------------------------------------------------*/
void Driver::compute_sh()
{
  char str[MAXLINE], *ptr;
  printf("\n"); for (int i = 0; i < 4; ++i) printf("====");
  printf("   Spherical Harmonics  Local Order Parameter   ");
  for (int i = 0; i < 4; ++i) printf("====");

  printf("\nReference: Phys. Rev. B 28:784, 1983.\n");
  printf("             J. Chem. Phys. 122:214722, 2005.\n");
  for (int i = 0; i < 20; ++i) printf("----");

  int L = 6;
  // ask for l
  printf("\nPlease input the value of l [6]: ");
  if (count_words(fgets(str, MAXLINE, stdin)) > 0){
    ptr = strtok(str, " \n\t\r\f");
    L = atoi(ptr);

    if (L < 0) L = 6;
  }
  printf("q%dq%d will be evaluated for your selected frames.\n", L, L);

  char prefix[MAXLINE];
  sprintf(prefix, "q%dq%d", L, L);
  // prefix to output files
  printf("\nPlease input the prefix to output the results [%s]: ", prefix);
  if (count_words(fgets(str, MAXLINE, stdin)) > 0){
    ptr = strtok(str, " \n\t\r\f");

    strcpy(prefix, ptr);
  }

  // identify local environment based on a threshold value
  int flag_env = 0;
  printf("\nWould you like to discriminate the local environment based on the measured\n");
  printf("  1.  q%d;\n  2.  w%d;\n  3.  q%dq%d;\n  0.  None.\nYour choice [0]: ", L, L, L, L);
  if (count_words(fgets(str, MAXLINE, stdin)) > 0){
    ptr = strtok(str, " \n\t\r\f");

    flag_env = atoi(ptr);
  }
  double thr = 0.;
  if (flag_env >= 1 && flag_env <= 3){
    printf("Please input the threshold value for that: ");
    if (count_words(fgets(str, MAXLINE, stdin)) > 0){
      ptr = strtok(str, " \n\t\r\f");

      thr = atof(ptr);
    }
    if (thr <= 0.) flag_env = 0;
  } else flag_env = 0;

  // thresholds for surface and edges
  set_cutoffs(0);
  one = all[istr];

  const double facql = sqrt(16.*atan(1.)/double(L+L+1));
  double **qlm, **qw; qlm = qw = NULL;
  SphericalHarmonics *sh = new SphericalHarmonics();

  // now to loop over all asked images
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // Compute the Voronoi info
    one->ComputeVoro(mins);
    one->dir2car();

    memory->grow(qw, 3, one->natom+1, "qw");
    memory->grow(qlm, one->natom+1, L+L+1, "qlm");
    // Loop over all atoms
#pragma omp parallel for default(shared)
    for (int id = 1; id <= one->natom; ++id){
      double norm2 = 0.;
      
      for (int m = -L; m <= L; ++m){
        qlm[id][m+L] = 0.;

        int ni = one->neilist[0][id];
        for (int jj = 1; jj <= ni; ++jj){
          int jd = one->neilist[jj][id];
          double rij[3];
          for (int idim = 0; idim < 3; ++idim) rij[idim] = one->atpos[jd][idim] - one->atpos[id][idim];

          one->ApplyPBC(rij[0], rij[1], rij[2]);

          qlm[id][m+L] += sh->Y(L, m, rij);
        }
        if (ni > 0) qlm[id][m+L] /= double(ni);
        norm2 += qlm[id][m+L] * qlm[id][m+L];
      }
      double norm = sqrt(norm2);

      qw[0][id] = norm*facql;
      if (norm < ZERO) norm = 1.;
      for (int i = 0; i <= L+L; ++i) qlm[id][i] /= norm;

      qw[1][id] = 0.;
      for (int m1 = -L; m1 <= L; ++m1)
      for (int m2 = -L; m2 <= L; ++m2)
      for (int m3 = -L; m3 <= L; ++m3){
        if ((m1+m2+m3) != 0) continue;

        double w3j = sh->w3j(L, m1, m2, m3);
        qw[1][id] += w3j * qlm[id][m1+L] * qlm[id][m2+L] * qlm[id][m3+L];
      }
    }

    // Now to evaluate qlql
    for (int id = 1; id <= one->natom; ++id){
      qw[2][id] = 0.;

      int ni = one->neilist[0][id];
      for (int jj = 1; jj <= ni; ++jj){
        int jd = one->neilist[jj][id];

        for (int i = 0; i <= L+L; ++i) qw[2][id] += qlm[id][i]*qlm[jd][i];
      }
      if (ni > 0) qw[2][id] /= double(ni);
    }

    // identify local environment
    if (flag_env){
      if (one->prop) memory->destroy(one->prop);
      one->prop = qw[flag_env-1];

      one->identify_env(thr);
    }
    // Now to output the results
    sprintf(str, "%s_%d.dat", prefix, one->iframe);
    ptr = strtok(str, " \n\t\r\f");
    FILE *fp = fopen(ptr, "w");
    fprintf(fp, "#Voronoi refinement info: surf_min = %g, edge_min = %g, nei_min = %d\n", mins[0], mins[2], int(mins[1]));
    fprintf(fp, "#%s local order parameter for frame %d of %s\n", prefix, one->iframe, one->fname);
    if (flag_env){
      fprintf(fp, "# 1  2    3 4 5 6   7   8   9\n");
      fprintf(fp, "# id type x y z q%d w%d q%dq%d env\n", L, L, L, L);
      for (int id = 1; id <= one->natom; ++id){
        fprintf(fp, "%d %d %lg %lg %lg %lg %lg %lg %d\n", id, one->attyp[id], one->atpos[id][0], one->atpos[id][1],
        one->atpos[id][2], qw[0][id], qw[1][id], qw[2][id], one->env[id]);
      }
      one->prop = NULL;

    } else {
      fprintf(fp, "# 1  2    3 4 5 6   7   8\n");
      fprintf(fp, "# id type x y z q%d w%d q%dq%d\n", L, L, L, L);
      for (int id = 1; id <= one->natom; ++id){
        fprintf(fp, "%d %d %lg %lg %lg %lg %lg %lg\n", id, one->attyp[id], one->atpos[id][0], one->atpos[id][1],
        one->atpos[id][2], qw[0][id], qw[1][id], qw[2][id]);
      }
    }
    fclose(fp);

    if (min_mem) one->FreeVoro();
    printf(" Frame %d done, the results are written to: %s\n", img+1, ptr);
  }

  if (qw)  memory->destroy(qw);
  if (qlm) memory->destroy(qlm);
  delete sh;
  printf("\n"); for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------ */
