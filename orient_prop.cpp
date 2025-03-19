#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the orientation for pairs of same property
 *----------------------------------------------------------------------------*/
void Driver::orient_same_property()
{ 
  const double rad2deg = 45. / atan(1.);
  char str[MAXLINE], header[MAXLINE];
  double rmin = 0., rmax = 0.5*MIN(MIN(all[0]->lx,all[0]->ly), all[0]->lz);

  // menu
  int job = 1;
  printf("\n"); printf("==");
  printf(" Orientation (cosine wrt. reference vector) for atom pairs of same property ==");
  printf("\nPlease select your desired job:\n");
  for (int i = 0; i < 20; ++i) printf("----"); printf("\n");
  printf("  1. Orientation for atom pairs of same property within selection;\n");
  printf("  2. Orientation for atom pairs of same property between selections;\n");
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

  // define reference orientation
  double refOrient[3];
  int nW = 0;
  while (nW < 3){
    printf("\nPlease input the reference vector (cartesian) to measure the orientation: ");
    input->read_stdin(str);
    ptr = strtok(str, " \n\t\r\f");
    while (ptr && nW < 3){
      refOrient[nW++] = atof(ptr);
      ptr = strtok(NULL, " \n\t\r\f");
    } 
  }
  double invNorm = 1./sqrt(refOrient[0]*refOrient[0] + refOrient[1]*refOrient[1] + refOrient[2]*refOrient[2]);
  for (int i = 0; i < 3; ++i) refOrient[i] *= invNorm;

  int flagHalfPi = 0;
  printf("\nDo you prefer to have the absolute values of the cosines measured? (y/n)[n]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr){
     if (strcmp(ptr, "y") == 0 || strcmp(ptr,"Y") == 0) flagHalfPi = 1;
  }
  
  
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

  if (job == 1) sprintf(header,"# Cosine wrt [%f %f %f] for atom pairs selected by: %s", refOrient[0], refOrient[1], refOrient[2], srcsel);
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
    sprintf(header,"# Cosine wrt [%f %f %f] for atom pairs between two selections:\n# source atoms: %s# neighbors: %s", refOrient[0], refOrient[1], refOrient[2], srcsel, dessel);
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

  // output file name or prefix
  printf("\nPlease input the file to output the result [orient.dat]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "orient.dat");

  ptr = strtok(str, " \n\t\r\f");
  char *fname = new char [strlen(ptr)+1];
  strcpy(fname, ptr);

  FILE *fp = fopen(fname,"w");
  fprintf(fp,"# Orientation for atom pairs with same %s\n", one->prop_label[pid].c_str());
  fprintf(fp,header);

  // timer
  Timer * timer = new Timer();

  // working space
  double com[3];
  int npairs = 0;

  int nused = 0;
  // now to compute the cosine with respect to the reference vector
  if (job == 1){  // cosine for atom pairs within selected atoms

    for (int img = istr; img <= iend; img += inc){ // loop over frames
      one = all[img];
      one->selection(srcsel);
      if (one->nsel < 1) continue;

      fprintf(fp, "# Frame ID: %d\n", img);
      fprintf(fp, "# atom-1 atom-2 com_x com_y com_z cosine angle\n");

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
          double rij = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

          com[0] = one->atpos[i][0];
          com[1] = one->atpos[i][1];
          com[2] = one->atpos[i][2];
          com[0] = com[0]*one->lx + com[1]*one->xy + com[2]*one->xz + dx[0];
          com[1] = com[1]*one->ly + com[2]*one->yz + dx[1];
          com[2] = com[2]*one->lz + dx[2];

          double cosine = (dx[0]*refOrient[0] + dx[1]*refOrient[1] + dx[2]*refOrient[2])/rij;
          if (flagHalfPi) cosine = fabs(cosine);
          double angle = acos(cosine);
          fprintf(fp, "%d %d %f %f %f %g %f\n", i, j, com[0], com[1], com[2], cosine, angle*rad2deg);

#pragma omp atomic
          ++npairs;
        }
      }
      ++nused;

      if (min_mem) one->FreeVoro();
    } // end loop over frames

  } else {        // cosine for atom pairs of same property between two selections

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

      fprintf(fp, "# Frame ID: %d\n", img);
      fprintf(fp, "# atom-1 atom-2 com_x com_y com_z cosine angle\n");
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
          double rij = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

          com[0] = one->atpos[i][0];
          com[1] = one->atpos[i][1];
          com[2] = one->atpos[i][2];
          com[0] = com[0]*one->lx + com[1]*one->xy + com[2]*one->xz + dx[0];
          com[1] = com[1]*one->ly + com[2]*one->yz + dx[1];
          com[2] = com[2]*one->lz + dx[2];

          double cosine = (dx[0]*refOrient[0] + dx[1]*refOrient[1] + dx[2]*refOrient[2])/rij;
          if (flagHalfPi) cosine = fabs(cosine);
          double angle = acos(cosine);

          fprintf(fp, "%d %d %f %f %f %g %f\n", i, j, com[0], com[1], com[2], cosine, angle*rad2deg);

#pragma omp atomic
          ++npairs;
        }
      }

      ++nused;
      if (min_mem) one->FreeVoro();
    }
  }

  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  memory->destroy(insrc);
  printf("\n%d images with %d pairs were used and the results are written to %s\n", nused, npairs, fname);
  fclose(fp);
  delete []fname;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------ */
