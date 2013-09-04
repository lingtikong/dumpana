#include "driver.h"
#include "common_neig.h"

/*------------------------------------------------------------------------------
 * Driver to execute the Common neighbor analysis or common neighbor parameter
 * computations;
 * Ref: Comp. Phys. Comm. 177:518, (2007).
 *----------------------------------------------------------------------------*/
void Driver::Compute_CNACNP()
{
  char str[MAXLINE], jobstr[MAXLINE];
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("   Common  Neighber  Analysis   ");
  for (int i = 0; i < 6; ++i) printf("====");

  int job = 1;
  printf("\nReference: Comp. Phys. Comm. 177:518, 2007.\n");
  for (int i = 0; i < 20; ++i) printf("----");
  printf("\nNow please select your desired job:\n");
  printf("  1. CNA based on Voronoi info;\n");
  printf("  2. CNP based on Voronoi info;\n");
  printf("  3. Centro-symmetry parameter;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 3){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  // thresholds for surface and edges
  set_cutoffs(0);

  // Show relevant info if Weighted Voronoi is used
  one = all[istr]; ShowRadius4Voro();

  char *fname = new char[8];
  if (job == 1){
    strcpy(fname, "cna.dat");
    strcpy(jobstr, "CNA");
  } else if (job == 2){
    strcpy(fname, "cnp.dat");
    strcpy(jobstr, "CNP");
  } else {
    strcpy(fname, "csp.dat");
    strcpy(jobstr, "CSP");
  }

  printf("\nPlease input the file name to output the %s info [%s]: ", jobstr, fname);
  fgets(str, MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr){
    delete []fname;
    fname = new char[strlen(ptr)+1];
    strcpy(fname, ptr);
  }

  FILE *fp = fopen(fname, "w");
  fprintf(fp,"#Voronoi refinement info: surf_min = %g, edge_min = %g, nei_min = %d\n", mins[0], mins[1], int(mins[2]));
  if (job == 1) fprintf(fp,"# CNA: 1, FCC; 2, HCP; 3, BCC; 4, ICOS; 4, OTHER; 5, UNKNOWN.\n# id x y z cna\n");
  else fprintf(fp, "# id x y z %s\n", jobstr);
  fflush(fp);

  if (job == 3){ // job will carry the # of nearest neighbor info now
    while ( 1 ){
      printf("\nPlease input the # of nearest neighbors for your reference lattice [8]: ");
      fgets(str, MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr){
        job = -atoi(ptr);
        if (job%2 == 1 || job > 0){
          printf("The # of nearest neighbors must be an even number!\n");
          continue;
        }
      } else job = -8;

      break;
    }
  }

  // now to loop over all asked images
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // Compute the Voronoi info
    one->ComputeVoro(mins, weighted);
    one->dir2car();

    fprintf(fp,"# frame number: %d\n", img);
    // now to compute the CNA/CNP info
    ComputeCNAAtom *cna = new ComputeCNAAtom(job, one->natom, one->neilist, one->atpos, one->box, fp);
    delete cna;

    if (min_mem) one->FreeVoro();
    printf(" Frame %d done, %s info written to: %s\n", img+1, jobstr, fname);
  }

  fclose(fp);
  delete []fname;

  printf("\n"); for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------ */
