#include "driver.h"
#include "math.h"
#include <set>

/*------------------------------------------------------------------------------
 * Method to compute the bond length/angles; Voronoi neighbors are seen as bonded
 *----------------------------------------------------------------------------*/
void Driver::bonds()
{
  char str[MAXLINE]; int job = 1;
  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("    Bond length / bond angle    ");
  for (int i=0; i<6; i++) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. bond lengths between selected atoms;\n");
  printf("  2. bond angles  between selected atoms;\n");
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

  // voro refinement info
  double mins[3];
  mins[0] = 1.e-2; mins[1] = 1.e-4; mins[2] = 0.;

  printf("\nRefined Voronoi tesselation will be needed to procceed.\n");
  printf("Now please input your criterion for tiny surfaces, 0 to keep all [%g]: ", mins[0]);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[0] = atof(ptr);
  printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", mins[0]);

  printf("Sometimes it might be desirable to keep a minimum # of neighbors when refining\n");
  printf("the Voronoi index, for example, keep at least 14 for a bcc lattice, 12 for hcp\n");
  printf("or fcc. If you prefer to do so, input a positive number now [%d]: ", int(mins[2]));
  if (count_words(fgets(str,MAXLINE, stdin)) > 0){
    double dum = atof(strtok(str, " \n\t\r\f"));
    if (dum > 0.) mins[2] = dum;
    printf("\nA minimum number of %d neighobrs will be kept no matter how tiny the surface is.\n", int(mins[2]));
  }

  printf("\nPlease input your criterion for ultra short   [%g]: ", mins[1]);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[1] = atof(ptr);
  printf("Edges whose length takes less ratio than %lg will be skipped!\n\n", mins[1]);

  char isel[MAXLINE], jsel[MAXLINE], ksel[MAXLINE];
  // selection commands for atoms
  one = all[istr]; WeightVoro();
  one->ComputeVoro(mins, weighted);

  if (job == 1) printf("\nA pair of atoms defines a bond, first define one end of the bond.\n");
  else printf("\nThree atoms define an angle, now define the central atom.\n");
  while (1){
    printf("Please input the atom selection command, `h` for help [all]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(isel, str);
      ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(isel,"all\n");
  
    // check the selection command on the first frame
    one->selection(isel); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  if (job == 1) printf("\nNow to define the other end of the bond.\n");
  else printf("\nNow to define the first end of the bond angle.\n");
  while (1){
    printf("Please input the atom selection command, `h` for help [all]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(jsel, str);
      ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(jsel,"all\n");
  
    // check the selection command on the first frame
    one->selection(jsel); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  if (job == 2){
    printf("\nNow to define the second end of the bond angle.\n");
    while (1){
      printf("Please input the atom selection command, `h` for help [all]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        strcpy(ksel, str);
        ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
      } else strcpy(ksel,"all\n");
    
      // check the selection command on the first frame
      one->selection(ksel); one->SelInfo();
      if (one->nsel < 1){
        printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
        if (count_words(fgets(str,MAXLINE,stdin)) > 0){
          char *ptr = strtok(str," \n\t\r\f");
          if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
        }
      }
      break;
    }
  }

  // output file name
  char *fname;
  printf("\nNow please input the output file name ");
  if (job == 1){
    printf("[bondlen.dat]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str, "bondlen.dat");
  } else {
    printf("[bondang.dat]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str, "bondang.dat");
  }
  ptr = strtok(str," \n\t\r\f");
  fname = new char[strlen(ptr)+1];
  strcpy(fname, ptr);

  // write header of output
  FILE *fp = fopen(fname, "w");
  if (job == 1){
    fprintf(fp,"# selection command for the 1st end of bonds: %s", isel);
    fprintf(fp,"# selection command for the 2nd end of bonds: %s", jsel);
    fprintf(fp,"# frame index id jd ip jp bondlength\n");
  } else {
    fprintf(fp,"# selection command for the center  of angles: %s", isel);
    fprintf(fp,"# selection command for the 1st end of angles: %s", jsel);
    fprintf(fp,"# selection command for the 2nd end of angles: %s", ksel);
    fprintf(fp,"# frame index jd id kd theta cos(theta)\n");
  }

  // working space for selection info
  const double rad2ang = 45./atan(1.);
  int natom = one->natom;
  int *is, *js, *ks;
  is = memory->create(is, natom+1, "is");
  js = memory->create(js, natom+1, "is");
  ks = NULL;

  int nused = 0;
  // loop over all frames
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // Compute the Voronoi neighbors
    one->ComputeVoro(mins, weighted);

    // make the first selection
    one->selection(isel);
    if (one->nsel < 1) continue;
    if (one->natom > natom){
      natom = one->natom;
      is = memory->grow(is, natom+1, "is:grow");
      js = memory->grow(js, natom+1, "is:grow");
    }
    for (int id=1; id<= one->natom; id++) is[id] = one->atsel[id];

    // make the second selection
    one->selection(jsel);
    if (one->nsel < 1) continue;
    for (int id=1; id<= one->natom; id++) js[id] = one->atsel[id];

    // make the third selection if needed
    if (job == 2){
      one->selection(ksel);
      if (one->nsel < 1) continue;
      ks = one->atsel;
    }

    // use fractional coordinate
    one->car2dir();
    double xij[3], xik[3], rij, rik;
    double lx = one->lx, ly = one->ly, lz = one->lz;
    double xy = one->xy, xz = one->xz, yz = one->yz;

    set<int> counted;

    // now to do the real job
    if (job == 1){ // bond length

      int ibond = 0;
      set<int> counted; counted.clear();

      for (int id=1; id<= one->natom; id++){
        if (is[id] == 0) continue;
        int ip = one->attyp[id];
  
        int ni = one->neilist[0][id];
        for (int jj=1; jj<= ni; jj++){
          int jd = one->neilist[jj][id];
          int index = MIN(id,jd)*(one->natom+1)+MAX(id,jd);
          if (js[jd] == 0 || counted.count(index)) continue;
          counted.insert(index);
          int jp = one->attyp[jd];
  
          for (int idim=0; idim<3; idim++){
            xij[idim] = one->atpos[jd][idim] - one->atpos[id][idim];
            while (xij[idim] > 0.5) xij[idim] -= 1.;
            while (xij[idim] <-0.5) xij[idim] += 1.;
          }
          xij[0] = xij[0] * lx + xij[1] * xy + xij[2] * xz;
          xij[1] = xij[1] * ly + xij[2] * yz;
          xij[2] = xij[2] * lz;
          rij = sqrt(xij[0]*xij[0] + xij[1]*xij[1] + xij[2]*xij[2]);
  
          fprintf(fp,"%d %d %d %d %d %d %lg\n", img+1, ++ibond, id, jd, ip, jp, rij);
        }
      }
      counted.clear();

    } else {       // bond angle

      int iangle = 0;

      for (int id=1; id<= one->natom; id++){
        if (is[id] == 0) continue;

        int ni = one->neilist[0][id];
        for (int jj=1; jj< ni; jj++){
          int jd = one->neilist[jj][id];
          if (js[jd] == 0 && ks[jd] == 0) continue;
  
          for (int idim=0; idim<3; idim++){
            xij[idim] = one->atpos[jd][idim] - one->atpos[id][idim];
            while (xij[idim] > 0.5) xij[idim] -= 1.;
            while (xij[idim] <-0.5) xij[idim] += 1.;
          }
          xij[0] = xij[0] * lx + xij[1] * xy + xij[2] * xz;
          xij[1] = xij[1] * ly + xij[2] * yz;
          xij[2] = xij[2] * lz;
          rij = xij[0]*xij[0] + xij[1]*xij[1] + xij[2]*xij[2];
  
          for (int kk=jj+1; kk<= ni; kk++){
            int kd = one->neilist[kk][id];
            if ( (js[jd] & ks[kd]) || (js[kd] & ks[jd]) ){
              for (int idim=0; idim<3; idim++){
                xik[idim] = one->atpos[kd][idim] - one->atpos[id][idim];
                while (xik[idim] > 0.5) xik[idim] -= 1.;
                while (xik[idim] <-0.5) xik[idim] += 1.;
              }
              xik[0] = xik[0] * lx + xik[1] * xy + xik[2] * xz;
              xik[1] = xik[1] * ly + xik[2] * yz;
              xik[2] = xik[2] * lz;
              rik = xik[0]*xik[0] + xik[1]*xik[1] + xik[2]*xik[2];
  
              double cos = (xij[0]*xik[0] + xij[1]*xik[1] + xij[2]*xik[2]) / sqrt(rij*rik);
              double ang = acos(cos)*rad2ang;
  
              fprintf(fp,"%d %d %d %d %d %g %lg\n", img+1, ++iangle, jd, id, kd, ang, cos);
            }
          }
        }
      }
    }

    nused++;
  }
  fclose(fp);
  memory->destroy(is); memory->destroy(js); ks = NULL;
  printf("\n%d images were used in counting bonds, the results is written to %s\n", nused, fname);
  delete []fname;

  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
