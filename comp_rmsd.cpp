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
  iref = MAX(0, MIN(iref, nframe-1));
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

  FILE *fp1 = NULL, *fp2 = NULL;
  printf("Frames from No. %d to No. %d with increment of %d will be compared to %d.\n", istr+1, iend+1, inc, iref);
  printf("\nIf you want to output the results to a file, input the filename now: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str, " \n\t\r\f");
    fp1 = fopen(ptr, "w");
    if (fp1) printf("The related info will be written to file: %s\n", ptr);
  }
  printf("If you want to output per-atom displacement, input the filename now: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str, " \n\t\r\f");
    fp2 = fopen(ptr, "w");
    if (fp2) printf("The related info will be written to file: %s\n", ptr);
  }
  double disp_tol = 0.;
  if (fp2){
    printf("If you want to count  the displaced atoms, input a threshold  now: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) disp_tol = atof(strtok(str, " \n\t\r\f"));
    if (disp_tol > 0.) printf("Atoms with a displacment > %g will be counted as moved.\n", disp_tol);
    else printf("All atoms will be counted as moved.\n");
    disp_tol *= disp_tol;
  }

  one = all[iref];
  int nref = one->natom;
  double pos_ref[nref][3], pos_one[nref][3];
  one->dir2car();

  for (int i = 0; i < nref; ++i){
    for (int j = 0; j < 3; ++j) pos_ref[i][j] = one->atpos[i+1][j] - one->atpos[1][j];
    one->ApplyPBC(pos_ref[i][0], pos_ref[i][1], pos_ref[i][2]);
  }
  
  RMSD *rmsd = new RMSD();
  
  // compute rmsd one by one
  printf("\n"); for (int i = 0; i < 80; ++i) printf("-");
  printf("\nRMSD with respect to frame %d read from: %s\n", iref+1, one->fname);
  printf("Frames with different # of atom to %d will be skipped.\n", iref+1);
  for (int i = 0; i < 80; ++i) printf("-");  printf("\n");
  printf("   Frame   timestep        RMSD (Ang)   Source\n");
  for (int i = 0; i < 80; ++i) printf("-");  printf("\n");
  if (fp1){
    fprintf(fp1,"\n#"); for (int i = 0; i < 80; ++i) fprintf(fp1,"-");
    fprintf(fp1,"\n#RMSD with respect to frame %d read from: %s\n#", iref+1, one->fname);
    for (int i = 0; i < 80; ++i) fprintf(fp1,"-");  fprintf(fp1,"\n");
    fprintf(fp1,"#  Frame   timestep        RMSD (Ang)   Source\n#");
    for (int i = 0; i < 80; ++i) fprintf(fp1,"-");  fprintf(fp1,"\n");
  }

  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];
    if (one->natom != nref) continue;
    one->dir2car();

    // make sure the configurations are comparable to each other
    for (int i = 0; i < nref; ++i){
      for (int j = 0; j < 3; ++j) pos_one[i][j] = one->atpos[i+1][j] - one->atpos[1][j];
      one->ApplyPBC(pos_one[i][0], pos_one[i][1], pos_one[i][2]);
      double dx = pos_one[i][0] - pos_ref[i][0];
      double dy = pos_one[i][1] - pos_ref[i][1];
      double dz = pos_one[i][2] - pos_ref[i][2];
      one->ApplyPBC(dx, dy, dz);

      pos_one[i][0] = pos_ref[i][0] + dx;
      pos_one[i][1] = pos_ref[i][1] + dy;
      pos_one[i][2] = pos_ref[i][2] + dz;
    }
    
    // to compute the rmsd
    double mcom[3], v2ref[3], U[3][3], res;
    rmsd->calculate_rotation_rmsd(pos_ref, pos_one, nref, mcom, v2ref, U, &res);

    printf("  %5d %10d %18.8f   %s\n", img+1, one->tstep, res, one->fname);
    if (fp1) fprintf(fp1,"  %5d %10d %18.8f   %s\n", img+1, one->tstep, res, one->fname);

    // output the per-atom info for each frame, similar to dump atomic style
    if (fp2){
      double rotated[3], disp[nref], dr_all = 0.;
      int nmoved = 0;
      for (int i = 0; i < nref; ++i){
        rotated[0] = rotated[1] = rotated[2] = 0.;
        for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) rotated[j] += U[j][k]*pos_one[i][k];
        disp[i] = 0.;
        for (int j = 0; j < 3; ++j){
          pos_one[i][j] = rotated[j] - pos_ref[i][j];
          disp[i] += pos_one[i][j] * pos_one[i][j];
        }
        dr_all += disp[i];
        if (disp[i] >= disp_tol) ++nmoved;
      }

      fprintf(fp2, "# natom timestep source frame-id\n");
      fprintf(fp2, "# %d %d %s %d\n", one->natom, one->tstep, one->fname, img+1);
      fprintf(fp2, "# Total rmsd: %lg, rmsd-per-atom: %g; total # of atoms moved: %d\n", sqrt(dr_all), sqrt(dr_all/double(nref)), nmoved);
      fprintf(fp2, "# 1  2   3    4  5   6 \n");
      fprintf(fp2, "# id ip  dx  dy  dz  dr\n");
      for (int i = 0; i < nref; ++i){
        fprintf(fp2, "%d %d %lg %lg %lg %lg\n", i+1, one->attyp[i+1], pos_one[i][0], pos_one[i][1], pos_one[i][2], sqrt(disp[i]));
      }
    }
  }
  if (fp1) fclose(fp1);
  if (fp2) fclose(fp2);

  delete rmsd;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
