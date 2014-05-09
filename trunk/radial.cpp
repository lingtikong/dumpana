#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the radial distribution of atoms
 *----------------------------------------------------------------------------*/
void Driver::radial()
{
  char str[MAXLINE], *ptr, *prefix;
  double cpos[3], dr, inv_dr;
  int ic = 0;

  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("  Radial distribution of atoms  ");
  for (int i = 0; i < 6; ++i) printf("====");

  one = all[istr];

  // ask for the # of bins along each direction
  printf("The subroutine will help to analyze the radial distribution of atoms,\n");
  printf("the related info will be written for each frame, and you can further\n");
  printf("analyze the results later on.\n");
  printf("You need to define the location of the center by `type [options] step`,\n");
  printf("where `type` can be `id`, `x`, or `X`. `step` is the bin size along the\n");
  printf("the radial direction. If `type` == `id`, the only option will be the\n");
  printf("ID of the atom that will be seen as the center; if `type` == `x` or `r`,\n");
  printf("then three numbers specifing the fractional coordinate of the center\n");
  printf("point is expected. While if `type` == `X` or `R`, three numbers giving\n");
  printf("the cartesian coordinate of the center should be followed.\n");
  printf("\nNow, please input your command: ");
  if (count_words(fgets(str,MAXLINE, stdin)) < 2){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  // analyse the input
  ptr = strtok(str, " \n\t\r\f");
  if (strcmp(ptr,"id") == 0){
    ptr = strtok(NULL, " \n\t\r\f");
    if (ptr == NULL){
      printf("\nInsufficient input!\n");
      for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
      return;
    }
    ic = atoi(ptr);
    if (ic < 1){
      printf("\nWrong input! The center ID must be a positive number!\n");
      for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
      return;

    } else if (ic > one->natom){
      printf("\nWrong input! The center ID does not exist!\n");
      for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
      return;
    }

  } else if (strcmp(ptr,"x") == 0 || strcmp(ptr, "r") == 0){
    for (int idim = 0; idim < 3; ++idim){
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL){
        printf("\nInsufficient input!\n");
        for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
        return;
      }
      cpos[idim] = atof(ptr);
      if (cpos[idim] < -1. || cpos[idim] > 1){
        printf("\nFractional coordinate expected!\n");
        for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
        return;
      }
    }
    ic = -1;

  } else if (strcmp(ptr,"X") == 0 || strcmp(ptr, "R") == 0){
    for (int idim = 0; idim < 3; ++idim){
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL){
        printf("\nInsufficient input!\n");
        for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
        return;
      }
      cpos[idim] = atof(ptr);
    }
    ic = 0;

  } else {
    printf("\nWrong input: %s\n", ptr);
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }

  ptr = strtok(NULL, " \n\t\r\f");
  if (ptr == NULL){
    printf("\nInsufficient input!\n");
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  dr = atof(ptr);
  inv_dr = 1./dr;
  
  // output filename prefix
  printf("Please input the prefix for output file names [radial]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str, "radial");
  ptr = strtok(str," \n\t\r\f");
  prefix = new char[strlen(ptr)+1];
  strcpy(prefix, ptr);

  // loop over all frames
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    one->dir2car();

    double c[3];
    // define center position
    if (ic > 0){
      for (int id = 1; id <= one->natom; ++id){
        if (ic == id){
          for (int idim = 0; idim < 3; ++idim) c[idim] = one->atpos[id][idim];
          break;
        }
      }

    } else if (ic == 0){
      for (int idim = 0; idim < 3; ++idim) c[idim] = cpos[idim];

    } else {
      c[0] = cpos[0]*one->lx + cpos[1]*one->xy + cpos[2]*one->xz + one->xlo;
      c[1] = cpos[1]*one->ly + cpos[2]*one->yz + one->ylo;
      c[2] = cpos[2]*one->lz + one->zlo;
    }

    // open file, and write the header
    sprintf(str, "%s_%d", prefix, one->iframe);
    FILE *fp = fopen(str, "w");
    fprintf(fp, "# Frame %d from file %s\n", one->iframe, one->fname);
    fprintf(fp, "# 1  2 3 4  5    6\n");
    fprintf(fp, "# id x y z ibin d2c\n");

    for (int id = 1; id <= one->natom; ++id){
      double xij = one->atpos[id][0] - c[0];
      double yij = one->atpos[id][1] - c[1];
      double zij = one->atpos[id][2] - c[2];

      one->ApplyPBC(xij, yij, zij);

      double r = sqrt(xij*xij + yij*yij + zij*zij);
      int ibin = int(r*inv_dr + 0.5);
      fprintf(fp, "%d %lg %lg %lg %d %lg\n", id, one->atpos[id][0], one->atpos[id][1], one->atpos[id][2], ibin, r);
    }
    fclose(fp);

    printf("Frame %d done, and the results are written to %s\n", one->iframe, str);
  }
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
