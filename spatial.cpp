#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the spatial distribution of atoms
 * Currently it works for orthogonal box only.
 *----------------------------------------------------------------------------*/
void Driver::spatial()
{
  char str[MAXLINE], *ptr;
  double ***hits;
  int nbin[3];
  double lo[3], hi[3], ds[3], inv_ds[3];

  nbin[0] = nbin[1] = nbin[2] = 1;
  lo[0] = lo[1] = lo[2] = 0.;
  hi[0] = hi[1] = hi[2] = 1.;

  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("  Spatial distribution of atoms ");
  for (int i=0; i<6; i++) printf("====");

  // ask for the # of bins along each direction
  printf("\nThe spatial distribution of atoms will be evaluated. You can define the\n");
  printf("way to discretize the space by input your command as `dir lo hi n` where\n");
  printf("`dir` is `x`, `y`, or `z`; `lo` and `hi` are the lower and upper bounds\n");
  printf("along the corresponding direction in fractional coordinate ([0-1]); `n`\n");
  printf("is the number of bins between `lo` and `hi`. A sequential input of this\n");
  printf("command for x, y, and z will discretize the space into 3D mesh. Any\n");
  printf("combination of `x`, `y`, and/or `z` is allowed. If not specified, the\n");
  printf("whole box along that direction is seen as a block.\n");
  printf("\nNow, please input your command: ");
  if (count_words(fgets(str,MAXLINE, stdin)) < 4){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }
  // analyse the input
  ptr = strtok(str, " \n\t\r\f");
  while (ptr){
    int idir = 0;
    if (strcmp(ptr,"x") == 0 || strcmp(ptr, "X") == 0){
      idir = 0;
    } else if (strcmp(ptr,"y") == 0 || strcmp(ptr, "Y") == 0){
      idir = 1;
    } else if (strcmp(ptr,"z") == 0 || strcmp(ptr, "Z") == 0){
      idir = 2;
    } else {
      printf("\nWrong input: %s\n", ptr);
      for (int i=0; i<20; i++) printf("===="); printf("\n");
      return;
    }
    ptr = strtok(NULL, " \n\t\r\f"); if (ptr == NULL) break;
    lo[idir] = atof(ptr);
    lo[idir] = MAX(0., MIN(1.,lo[idir]));

    ptr = strtok(NULL, " \n\t\r\f"); if (ptr == NULL) break;
    hi[idir] = atof(ptr);
    hi[idir] = MAX(0., MIN(1.,hi[idir]));

    if (lo[idir] >= hi[idir]){
      printf("\nWrong input: lo[%c] >= hi[%c]!\n", idir+'x', idir+'x');
      for (int i=0; i<20; i++) printf("===="); printf("\n");
      return;
    }

    ptr = strtok(NULL, " \n\t\r\f"); if (ptr == NULL) break;
    nbin[idir] = atoi(ptr);
    if (nbin[idir] < 1) nbin[idir] = 1;

    ptr = strtok(NULL, " \n\t\r\f");
  }
  int nb = nbin[0]*nbin[1]*nbin[2];
  if (nb <= 1){
    printf("\nThe whole simulation box is set as one bin, no computation needed!\n");
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }

  // print out the dicretization info
  printf("\nSpace discretization info based on your input:\n");
  for (int i=0; i<20; i++) printf("----");
  printf("\ndirection    lower-bound    upper-bound   #-of-bins");
  printf("\n"); for (int i=0; i<20; i++) printf("----");
  for (int idim=0; idim<3; idim++){
    printf("\n    %c   %14g %14g %8d", idim+'x', lo[idim], hi[idim], nbin[idim]);
  }
  printf("\n"); for (int i=0; i<20; i++) printf("----"); printf("\n");

  // assign the step size
  for (int i=0; i<3; i++){
    ds[i] = (hi[i]-lo[i])/double(nbin[i]);
    inv_ds[i] = 1./ds[i];
  }

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

  // allocate the memory for accumulation
  hits = memory->create(hits, nbin[0], nbin[1], nbin[2], "hits");
  for (int i= 0; i< nbin[0]; i++)
  for (int j= 0; j< nbin[1]; j++)
  for (int k= 0; k< nbin[2]; k++) hits[i][j][k] = 0.;

  int nused = 0;

  // loop over all frames
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // does not work for non-orthogonal box
    if (one->triclinic) continue;

    // select atoms as source
    one->selection(selcmd);
    if (one->nsel < 1) continue;

    // factional coordinate needed
    one->car2dir();

    for (int ii=1; ii<= one->natom; ii++){
      if (one->atsel[ii] == 0) continue;
      int idx[3];
      int outside = 0;
      for (int idim=0; idim<3; idim++){
        idx[idim] = (one->atpos[ii][idim] - lo[idim])*inv_ds[idim];
        if (idx[idim] < 0 || idx[idim] >= nbin[idim]) outside = 1;
      }

      if (outside == 0) hits[idx[0]][idx[1]][idx[2]] += 1.;
    }

    nused++;
  }
  if (nused < 1) return;

  // normalize the data
  double fac = 1./double(nused);
  for (int ii = 0; ii < nbin[0]; ii++)
  for (int jj = 0; jj < nbin[0]; jj++)
  for (int kk = 0; kk < nbin[0]; kk++) hits[ii][jj][kk] *= fac;

  // output the result
  int fastz = 1;
  printf("\nPlease input the fastest direction to output the data (x/z)[z]: ");
  fgets(str,MAXLINE, stdin); ptr = strtok(str, " \n\t\r\f");
  if (ptr != NULL && (strcmp(ptr,"x") == 0 || strcmp(ptr, "X") == 0)) fastz = 0;

  // output a blank line after each cycle of the fatest direction
  int blank = 0;
  printf("\nWould you like to write a blank line after each cycle of the fatest direction (y/n)[n]: ");
  fgets(str,MAXLINE, stdin); ptr = strtok(str, " \n\t\r\f");
  if (ptr != NULL && (strcmp(ptr,"y") == 0 || strcmp(ptr, "Y") == 0)) blank = 1;

  printf("\nPlease input the file name to output the result [spatial.dat]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "spatial.dat");
  ptr = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(ptr,"w");
  fprintf(fp,"# Spatial distribution for atoms: %s\n", selcmd);
  fprintf(fp,"# Similar to fix-ave-spatial format; mesh size: %d x %d x %d.\n", nbin[0], nbin[1], nbin[2]);
  fprintf(fp,"# index sx sy sz #atoms\n");
  fprintf(fp,"# 0 %d\n", nbin[0]*nbin[1]*nbin[2]);
  int index = 0;
  if (fastz){
    for (int ii=0; ii< nbin[0]; ii++)
    for (int jj=0; jj< nbin[1]; jj++){
      for (int kk=0; kk< nbin[2]; kk++){
        fprintf(fp, "%d %lg %lg %lg %lg\n", index++, (double(ii)+0.5)*ds[0]+lo[0],
        (double(jj)+0.5)*ds[1]+lo[1], (double(kk)+0.5)*ds[2]+lo[2],hits[ii][jj][kk]);
      }
      if (blank) fprintf(fp,"\n");
    }
  } else {
    for (int kk=0; kk< nbin[2]; kk++)
    for (int jj=0; jj< nbin[1]; jj++){
      for (int ii=0; ii< nbin[0]; ii++){
        fprintf(fp, "%d %lg %lg %lg %lg\n", index++, (double(ii)+0.5)*ds[0]+lo[0],
        (double(jj)+0.5)*ds[1]+lo[1], (double(kk)+0.5)*ds[2]+lo[2], hits[ii][jj][kk]);
      }
      if (blank) fprintf(fp,"\n");
    }
  }
  fclose(fp);

  memory->destroy(hits);
  printf("\n%d images were used in the evaluation and the result is written to %s\n", nused, ptr);

  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
