#include "driver.h"
#include "rmsd.h"

/*------------------------------------------------------------------------------
 * Method to compare the RMSD between frames
 *----------------------------------------------------------------------------*/
void Driver::compare_rmsd()
{
  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("  RMSD between different frames ");
  for (int i=0; i<6; i++) printf("====");

  int iref = 0;
  char str[MAXLINE];
  printf("\nTotal number of frames read: %d\n", nframe);
  printf("Please input your desired frame that will be set as a reference [1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) iref = atoi(strtok(str, " \n\t\r\f"))-1;
  printf("Frame %d will be set as a reference.\n", iref+1);

  printf("\nPlease input the frame, or frame range that will be compared to frame %d: ", iref+1);
  int nw = count_words(fgets(str,MAXLINE,stdin));
  if (nw < 1){
    istr = 0; iend = 0; inc = 1;

  } else if (nw == 1){
    istr = iend = atoi(strtok(str, " \n\t\r\f"))-1;
    inc = 1;

  } else if (nw == 2){
    istr = atoi(strtok(str, " \n\t\r\f"))-1;
    iend = atoi(strtok(NULL," \n\t\r\f"))-1;
    inc = 1;

  } else {
    istr = atoi(strtok(str, " \n\t\r\f"))-1;
    iend = atoi(strtok(NULL," \n\t\r\f"))-1;
    inc  = atoi(strtok(NULL," \n\t\r\f"));
  }

  istr = MAX(0,istr);
  iend = MIN(iend,nframe-1);
  inc = MAX(1,inc);

  printf("Frames from No. %d to No. %d with increment of %d will be compared to %d.\n", istr+1, iend+1, inc, iref);

  one = all[iref];
  int nref = one->natom;
  double pos_ref[nref][3], pos_one[nref][3];
  double mcom[3], v2ref[3], U[3][3];
  one->dir2car();

  for (int i = 0; i < nref; ++i)
  for (int j = 0; j < 3; ++j) pos_ref[i][j] = one->atpos[i][j];
  
  RMSD *rmsd = new RMSD();
  
  // compute rmsd one by one
  printf("\n"); for (int i = 0; i < 80; ++i) printf("-");
  printf("\nRMSD with respect to frame %d read from: %s\n", iref+1, one->fname);
  printf("Frames with different # of atom to %d will be skipped.\n", iref+1);
  for (int i = 0; i < 80; ++i) printf("-");  printf("\n");
  printf("  Frame         RMSD (Ang)    Source\n");
  for (int i = 0; i < 80; ++i) printf("-");  printf("\n");

  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];
    if (one->natom != nref) continue;
    one->dir2car();

    for (int i = 0; i < nref; ++i)
    for (int j = 0; j < 3; ++j) pos_one[i][j] = one->atpos[i][j];
    
    // to compute the rmsd
    double res;
    rmsd->calculate_rotation_rmsd(pos_ref, pos_one, nref, mcom, v2ref, U, &res);

    printf("  %5d  %18.8f   %s\n", img+1, res, one->fname);
  }

  delete rmsd;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
