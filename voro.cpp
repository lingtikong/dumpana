#include "driver.h"

/*------------------------------------------------------------------------------
 * Method to do voronoi diagram analysis for the selected frames
 * Resultant voro index info will be written to files: voro%{tstep}
 *----------------------------------------------------------------------------*/
void Driver::voro()
{
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("   Voronoi  Diagram  Analysis   ");
  for (int i = 0; i < 6; ++i) printf("====");

  // thresholds for surface and edges
  set_cutoffs(1);
  one = all[istr];

  printf("\nPlease input the prefix for output files [voro]: ");
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
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
    fprintf(fp,"#Voronoi refinement info: surf_min = %g, edge_min = %g, nei_min = %d\n", mins[0], mins[2], int(mins[1]));
    fprintf(fp,"# 1  2    3  4  5  6   7         8  9    10\n");
    fprintf(fp,"# id type x  y  z  vol voroindex n5 NNei NeiList surfaceareas\n");

    // ask for Voronoi info
    one->ComputeVoro(mins, fp, fpsurf, fpedge);

    fclose(fp);
    if (min_mem) one->FreeVoro();
    printf("  Frame %d done, voro info written to: %s\n", img+1, ptr);
  }

  if (fpsurf) fclose(fpsurf);
  if (fpedge) fclose(fpedge);

  delete []prefix;
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------*/
