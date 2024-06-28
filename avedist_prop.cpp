#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the average distance for pairs of same property
 *----------------------------------------------------------------------------*/
void Driver::ave_dist_same_property()
{
  char str[MAXLINE], header[MAXLINE];
  double rmin = 0., rmax = 0.5*MIN(MIN(all[0]->lx,all[0]->ly), all[0]->lz);

  // menu
  int job = 1;
  printf("\n"); for (int i = 0; i < 5; ++i) printf("===");
  printf(" Average distance for atom pairs of same property ");
  for (int i = 0; i < 5; ++i) printf("===");
  printf("\nPlease select your desired job:\n");
  for (int i = 0; i < 20; ++i) printf("----"); printf("\n");
  printf("  1. Average distance for atom pairs of same property within selection;\n");
  printf("  2. Average distance for atom pairs of same property between selections;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  input->read_stdin(str);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 2){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  // selection commands
  char srcsel[MAXLINE], dessel[MAXLINE];
  int *insrc; insrc = NULL;
  
  one = all[istr];
  if (one->prop_label.size() < 1){
    printf("\nNo property read in the dump file, the calculation cannot be proceeded.\n");
    return;
  }

  // selection commands for atoms
  while (1){
    if (job == 1) printf("\nPlease input the selection command for atoms, `h` for help [all]: ");
    else printf("\nPlease input the selection command for source atoms, `h` for help [all]: ");

    input->read_stdin(str);
    if (count_words(str) > 0){
      strcpy(srcsel, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(srcsel,"all\n");
  
    // check the selection command on the first frame
    one->selection(srcsel); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      input->read_stdin(str);
      if (count_words(str) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  if (job == 1) sprintf(header,"# Average distance for atom pairs selected by: %s", srcsel);
  else {
    // atoms as neighbors
    while (1){
      printf("\nPlease input the selection command for neighbors, `h` for help [all]: ");
      input->read_stdin(str);
      if (count_words(str) > 0){
        strcpy(dessel, str);
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
      } else strcpy(dessel,"all\n");

      // check the selection command on the first frame
      one->selection(dessel); one->SelInfo();
      if (one->nsel < 1){
        printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
        input->read_stdin(str);
        if (count_words(str) > 0){
          char *ptr = strtok(str," \n\t\r\f");
          if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
        }
      }
      break;
    }
    memory->create(insrc, one->natom+1, "insrc");
    sprintf(header,"# Average distance for pairs between two selections:\n# source atoms: %s# neighbors: %s", srcsel, dessel);
  }
  
  // property 
  int pid = 0;
  printf("\nThe available properties are: ");
  for (int i = 0; i < one->prop_label.size(); ++i) printf("%d) %s; ", i+1, one->prop_label[i].c_str());
  printf("\nPlease input the id of the property to be used [%d]: ", pid+1);
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) pid = atoi(ptr) - 1;
  pid = MIN(MAX(pid, 0), one->prop_label.size());
  printf("The property chosen is: %s\n", one->prop_label[pid].c_str());

  // output average distance or for each frame
  printf("\nWould you like to output:\n");
  printf("  1. the average distance for all frames used;\n");
  printf("  2. the average distance for each frame;\n");
  printf("Your choice [1]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  int flag_each = 0;
  if (ptr) flag_each = atoi(ptr)-1;

  // output file name or prefix
  printf("\nPlease input the file to output the result [avedist.dat]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "avedist.dat");

  ptr = strtok(str, " \n\t\r\f");
  char *fname = new char [strlen(ptr)+1];
  strcpy(fname, ptr);

  FILE *fp = fopen(fname,"w");
  fprintf(fp,"# Average distance for atom pairs with same %s\n", one->prop_label[pid].c_str());

  // timer
  Timer * timer = new Timer();

  // working space
  double total_dist = 0.;
  int npairs = 0;

  int nused = 0;
  // now to compute the average distance
  if (job == 1){  // average distance for atom pairs within selected atoms
    fprintf(fp, "# within selection: %s", srcsel);
    if (flag_each == 1){
       fprintf(fp, "# img  pairs average-distance\n");
    } else {
       fprintf(fp, "# frames-used pairs average-distance\n");
    }

    for (int img = istr; img <= iend; img += inc){ // loop over frames
      one = all[img];
      one->selection(srcsel);
      if (one->nsel < 1) continue;

      // need fractional coordinates
      one->car2dir();
  
      // set local variables
#pragma omp parallel for default(shared)
      for (int i = 1; i < one->natom; ++i){
        if (one->atsel[i] == 0) continue;

        for (int j = i+1; j <= one->natom; ++j){
          if (one->atsel[j] == 0) continue;
          if (fabs(one->atprop[i][pid] - one->atprop[j][pid]) > ZERO) continue;

          double dx[3];
          for (int idim = 0; idim < 3; ++idim){
            dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
            while (dx[idim] > 0.5) dx[idim] -= 1.;
            while (dx[idim] <-0.5) dx[idim] += 1.;
          }
          dx[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
          dx[1] = dx[1]*one->ly + dx[2]*one->yz;
          dx[2] = dx[2]*one->lz;
          double r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

#pragma omp atomic
          total_dist += r;
#pragma omp atomic
          ++npairs;
        }
      }
      ++nused;

      if (min_mem) one->FreeVoro();
      if (flag_each == 1){
         fprintf(fp, "%d %d %lg\n", img, npairs, total_dist/MAX(1, npairs));
         total_dist = 0.;
         npairs = 0;
      }
    } // end loop over frames
    if (flag_each != 1) fprintf(fp, "%d %d %lg\n", nused, npairs, total_dist/MAX(1, npairs));

  } else {        // Average distance of atom pairs of same property between two selections

    fprintf(fp, "# between selections: %s# and %s", srcsel, dessel);
    if (flag_each == 1){
       fprintf(fp, "# img  pairs average-distance\n");
    } else {
       fprintf(fp, "# frames-used pairs average-distance\n");
    }
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
  
      // select atoms as source
      one->selection(srcsel);
      if (one->nsel < 1) continue;
      memory->grow(insrc, one->natom+1, "insrc");
      for (int ii = 1; ii <= one->natom; ++ii) insrc[ii] = one->atsel[ii];
      int nsrc = one->nsel;

      // select atoms as neighbors
      one->selection(dessel);
      if (one->nsel < 1) continue;

      // need fractional coordinates
      one->car2dir();

      // set local variables
#pragma omp parallel for default(shared)
      for (int i = 1; i <= one->natom; ++i){
        if (insrc[i] == 0) continue;

        for (int j = 1; j <= one->natom; ++j){
          if (one->atsel[j] == 0 || i==j) continue;
          if (fabs(one->atprop[i][pid] - one->atprop[j][pid]) > ZERO) continue;

          double dx[3];
          for (int idim = 0; idim < 3; ++idim){
            dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
            while (dx[idim] > 0.5) dx[idim] -= 1.;
            while (dx[idim] <-0.5) dx[idim] += 1.;
          }
          dx[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
          dx[1] = dx[1]*one->ly + dx[2]*one->yz;
          dx[2] = dx[2]*one->lz;
          double r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

#pragma omp atomic
          total_dist += r;
#pragma omp atomic
          ++npairs;
        }
      }

      ++nused;
      if (min_mem) one->FreeVoro();
      if (flag_each == 1){
         fprintf(fp, "%d %d %lg\n", img, npairs, total_dist/MAX(1, npairs));
         total_dist = 0.;
         npairs = 0;
      }
    }
    if (flag_each != 1) fprintf(fp, "%d %d %lg\n", nused, npairs, total_dist/MAX(1, npairs));
  }

  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  memory->destroy(insrc);
  if (flag_each == 0) printf("\n%d images were used in the evaluation, which is written to %s\n", nused, fname);
  delete []fname;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------ */
