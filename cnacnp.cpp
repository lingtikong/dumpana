#include "driver.h"
#include "common_neig.h"

/*------------------------------------------------------------------------------
 * Driver to execute the Common neighbor analysis or common neighbor parameter
 * computations;
 * Ref: Comp. Phys. Comm. 177:518, (2007).
 *----------------------------------------------------------------------------*/
void Driver::Compute_CNACNP()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("   Common  Neighber  Analysis   ");
  for (int i=0; i<6; i++) printf("====");

  int job = 1;
  printf("\nReference: Comp. Phys. Comm. 177:518, 2007.\n");
  for (int i=0; i<20; i++) printf("----");
  printf("\nNow please select your desired job:\n");
  printf("  1. CNA based on Voronoi info;\n");
  printf("  2. CNP based on Voronoi info;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 2){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  // thresholds for surface and edges
  double mins[3];
  mins[0] = 1.e-2; mins[1] =-1.; mins[2] = 0.;

  // ask for Voronoi refinement parameters
  printf("Please input your criterion for tiny surfaces, 0 to keep all [%g]: ", mins[0]);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[0] = atof(ptr);
  printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", mins[0]);

  char *fname = new char[8];
  if (job == 1) strcpy(fname, "cna.dat");
  else strcpy(fname, "cnp.dat");

  printf("Please input the file name to output the CNA/CNP info [%s]: ", fname);
  fgets(str, MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr){
    delete []fname;
    fname = new char[strlen(ptr)+1];
    strcpy(fname, ptr);
  }

  FILE *fp = fopen(fname, "w");
  fprintf(fp,"#Voronoi refinement info: surf_min%% = %g, edge_min%% = %g, nei_min = %d\n", mins[0], mins[1], int(mins[2]));
  if (job == 1) fprintf(fp,"# CNA: 1, FCC; 2, HCP; 3, BCC; 4, ICOS; 4, OTHER; 5, UNKNOWN.\n# id x y z cna\n");
  else fprintf(fp, "# id x y z cnp\n");
  fflush(fp);

  // now to loop over all asked images
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // Compute the Voronoi info
    one->ComputeVoro(mins);
    one->dir2car();

    fprintf(fp,"# frame number: %d\n", img);
    // now to compute the CNA/CNP info
    ComputeCNAAtom *cna = new ComputeCNAAtom(job, one->natom, one->neilist, one->atpos, one->box, fp);
    delete cna;

    printf(" Frame %d done, CNA/CNP info written to: %s\n", img+1, fname);
  }

  fclose(fp);
  delete []fname;

  printf("\n"); for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------ */
