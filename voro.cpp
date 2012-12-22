#include "driver.h"

/*------------------------------------------------------------------------------
 * Method to do voronoi diagram analysis for the selected frames
 * Resultant voro index info will be written to files: voro%{tstep}
 *----------------------------------------------------------------------------*/
void Driver::voro()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("   Voronoi  Diagram  Analysis   ");
  for (int i=0; i<6; i++) printf("====");

  // thresholds for surface and edges
  double mins[3];
  mins[0] = 1.e-2; mins[1] = 2.e-3; mins[2] = 0.;

  printf("\nPlease input your criterion for tiny surfaces, 0 to keep all [%g]: ", mins[0]);
  fgets(str,MAXLINE, stdin);
  char * ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[0] = atof(ptr);
  printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", mins[0]);

  printf("Sometimes it might be desirable to keep a minimum # of neighbors when refining\n");
  printf("the Voronoi index, for example, keep at least 14 for a bcc lattice, 12 for hcp\n");
  printf("or fcc. If you prefer to do so, input a positive number now [%d]: ", int(mins[2]));
  if (count_words(fgets(str,MAXLINE, stdin)) > 0){
    mins[2] = atof(strtok(str, " \n\t\r\f"));
    if (mins[2] < 1.) mins[2] = 0.;
    else printf("\nA minimum number of %d neighobrs will be kept no matter how tiny the surface is.\n", int(mins[2]));
  }

  printf("Please input your criterion for ultra short edges, 0 to keep all [%g]: ", mins[1]);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[1] = atof(ptr);
  printf("Edges whose lengths take less ratio than %lg will be skipped!\n", mins[1]);

  // show relevant info if weighted Voronoi is used
  one = all[istr]; ShowRadius4Voro();

  printf("\nPlease input the prefix for output files [voro]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"voro"); ptr = strtok(str, " \n\t\r\f");}
  char *prefix = new char[strlen(ptr)+1];
  strcpy(prefix, ptr);
  
  FILE *fpsurf, *fpedge; fpsurf = fpedge = NULL;
  if (flag_out & OutSurf) fpsurf = fopen("surf_ratio.dat", "w");
  if (flag_out & OutEdge) fpedge = fopen("edge_ratio.dat", "w");

  // now to do the real job
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // open file for output
    sprintf(str,"%s_%d.dat ", prefix, one->tstep);
    ptr = strtok(str," \n\t\r\f");
    FILE *fp = fopen(ptr, "w");
    fprintf(fp,"#Box info: %lg %lg %lg %lg %lg %lg %d %lg %lg %lg\n",
    one->xlo, one->xhi, one->ylo, one->yhi, one->zlo, one->zhi, one->natom, one->xy, one->xz, one->yz);
    fprintf(fp,"#Voronoi refinement info: surf_min%% = %g, edge_min%% = %g, nei_min = %d\n", mins[0], mins[1], int(mins[2]));
    fprintf(fp,"# 1  2    3  4  5  6   7         8    9    10\n");
    fprintf(fp,"# id type x  y  z  vol voroindex f5%%  NNei NeiList surfaceareas\n");

    // ask for Voronoi info
    one->ComputeVoro(mins, fp, fpsurf, fpedge, weighted);

    fclose(fp);
    printf("  Frame %d done, voro info written to: %s\n", img+1, ptr);
  }

  if (fpsurf) fclose(fpsurf);
  if (fpedge) fclose(fpedge);

  delete []prefix;
  for (int i=0; i<20; i++) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------*/
