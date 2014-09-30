#include "driver.h"
#include "math.h"
#include <vector>

/*------------------------------------------------------------------------------
 * Method to compute the Mean square displacement of selected atoms.
 *
 * Note that in this case the dump file read must contain the image info.
 *----------------------------------------------------------------------------*/
void Driver::compute_msd()
{
  char str[MAXLINE], header[MAXLINE];
  printf("\n"); for (int i = 0; i < 5; ++i) printf("====");
  printf("  Mean-squared displacement calculation ");
  for (int i = 0; i < 5; ++i) printf("===="); printf("\n");

  // selection commands; the selections will be intersected among frames.
  char selcmd[MAXLINE];
  one = all[istr];
  int natom = one->natom;
  int *insel = memory->create(insel, natom+1, "insel");
  for (int id = 1; id <= natom; ++id) insel[id] = 1;

  int nused = 0;
  // selection commands for atoms
  while (1){
    printf("\nPlease input the selection command for atoms, `h` for help [all]: ");
  
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(selcmd, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }

    } else strcpy(selcmd,"all\n");
    
    if (strcmp(selcmd, "all") == 0) break;

    // determine if selection is based on first frame or all frames
    int flag_sel = 2;
    printf("Should the selection be applied to:\n");
    printf("  1. the first frame only;\n");
    printf("  2. all frames and get the common selection.\n");
    printf("Your choice [2]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr = strtok(str," \n\t\r\f");
      flag_sel = MIN(2, MAX(atoi(ptr),1));
    }

    // check the selection command
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      if (one->natom != natom){
        printf("\nError: number of atoms from frame %d differs from previous ones!\n", img+1);
        for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
        return;
      }
      if (one->image == NULL){
        printf("\nError: image info not available for frame %d! MSD calculation is not possible.\n", img+1);
        for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
        return;
      }
      if (flag_sel == 2 || img == istr){
        one->selection(selcmd);
        for (int id = 1; id <= natom; ++id) insel[id] &= one->atsel[id];
      }

      // cartesian coordinate needed for MSD calculations
      one->dir2car();

      ++nused;
    }
    int nsel = 0;
    for (int id = 1; id <= natom; ++id) nsel += insel[id];
    if (nsel < 1){
      printf("It seems that no atom is selected, do you wish to make another selection? (y/n)[n]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")== 0 || strcmp(ptr,"Y")==0) continue;
      }
      for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
      return;
    }

    break;
  }
  
  double **msd = memory->create(msd, nused, 6, "msd");
  for (int img = 0; img < nused; ++img)
  for (int idim = 0; idim < 6; ++idim) msd[img][idim] = 0.;

  for (int ir = istr; ir <= iend; ir += inc){
    one = all[ir];
    int istep = 0;
    for (int in = ir+inc; in <= iend; in += inc){
      DumpAtom *now = all[in];
      msd[++istep][0] = double(now->tstep - one->tstep);

      if (one->triclinic == 0) {
        for (int id = 1; id <= natom; ++id){
          if (insel[id] == 0) continue;

          int xbox = one->image[id][0];
          int ybox = one->image[id][1];
          int zbox = one->image[id][2];
          double x0 = one->atpos[id][0] + xbox*one->lx;
          double y0 = one->atpos[id][1] + ybox*one->ly;
          double z0 = one->atpos[id][2] + zbox*one->lz;

          xbox = now->image[id][0];
          ybox = now->image[id][1];
          zbox = now->image[id][2];
          double dx = now->atpos[id][0] - x0 + xbox*now->lx;
          double dy = now->atpos[id][1] - y0 + ybox*now->ly;
          double dz = now->atpos[id][2] - z0 + zbox*now->lz;

          msd[istep][1] += dx*dx;
          msd[istep][2] += dy*dy;
          msd[istep][3] += dz*dz;
          msd[istep][4] += dx*dx + dy*dy + dz*dz;
          msd[istep][5] += 1.;
        }

      } else {
        for (int id = 1; id <= natom; ++id){
          if (insel[id] == 0) continue;

          int xbox = one->image[id][0];
          int ybox = one->image[id][1];
          int zbox = one->image[id][2];
          double x0 = one->atpos[id][0] + xbox * one->lx + ybox * one->xy + zbox * one->xz;
          double y0 = one->atpos[id][1] + ybox * one->ly + zbox * one->yz;
          double z0 = one->atpos[id][2] + zbox * one->lz;

          xbox = now->image[id][0];
          ybox = now->image[id][1];
          zbox = now->image[id][2];
          double dx = now->atpos[id][0] - x0 + xbox * now->lx + ybox * now->xy + zbox * now->xz;
          double dy = now->atpos[id][1] - y0 + ybox * now->ly + zbox * now->yz;
          double dz = now->atpos[id][2] - z0 + zbox * now->lz;

          msd[istep][1] += dx*dx;
          msd[istep][2] += dy*dy;
          msd[istep][3] += dz*dz;
          msd[istep][4] += dx*dx + dy*dy + dz*dz;
          msd[istep][5] += 1.;
        }
      }
    }
  }

  // average
  for (int it = 1; it < nused; ++it){
    if (msd[it][5] > 0.)
      for (int idim = 1; idim <= 4; ++idim) msd[it][idim] /= msd[it][5];
  }

  // output the result
  printf("Please input the file to output the MSD info [msd.dat]: ");
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "msd.dat");
  ptr = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(ptr,"w");
  fprintf(fp,"# MSD for: %s", selcmd);
  fprintf(fp,"# d-step x y z r\n");
  for (int it = 0; it < nused; ++it){
    fprintf(fp,"%d %lg %lg %lg %lg\n", int(msd[it][0]), msd[it][1], msd[it][2], msd[it][3], msd[it][4]);
  }
  fclose(fp);
  
  memory->destroy(msd);
  printf("\n%d images were used in the evaluation of MSD, which is written to %s\n", nused, ptr);

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------*/
