#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the pair correlation function of the system
 *----------------------------------------------------------------------------*/
void Driver::paircorr()
{
  char str[MAXLINE], header[MAXLINE];
  double rmin = 0., rmax = 0.5*MIN(MIN(all[0]->lx,all[0]->ly), all[0]->lz);
  int nbin = 201;

  // menu
  int job = 1;
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("  Pair  Correlation  Functions  ");
  for (int i = 0; i < 6; ++i) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i = 0; i < 20; ++i) printf("----"); printf("\n");
  printf("  1. g(r) of selected atoms;\n");
  printf("  2. g(r) between two selections;\n");
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
  int *insrc;
  
  one = all[istr];
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

  if (job == 1) sprintf(header,"# g(r) for atoms selected by: %s", srcsel);
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
    sprintf(header,"# g(r) between two selections:\n# source atoms: %s# neighbors: %s", srcsel, dessel);
  }
  
  // bounds of g(r)
  printf("\nPlease input the lower bound of g(r) [%g]: ", rmin);
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) rmin = atof(ptr);
  rmin = MAX(rmin, 0.);

  printf("Please input the upper bound of g(r) [%g]: ", rmax);
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) rmax = atof(ptr);
  rmax = MAX(rmax,rmin+1.);

  printf("Please input the # of points of g(r) [%d]: ", nbin);
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) nbin = atoi(ptr);
  nbin = MAX(2,nbin);
  
  // output average g(r) or for each frame
  printf("\nWould you like to output:\n");
  printf("  1. the average g(r) for all frames used;\n");
  printf("  2. a separate  g(r) for each frame;\n");
  printf("Your choice [1]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  int flag_each = 0;
  if (ptr) flag_each = atoi(ptr)-1;

  // output file name or prefix
  if (flag_each == 0){
    printf("\nPlease input the file to output g(r) [gr.dat]: ");
    input->read_stdin(str);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr == NULL) strcpy(str, "gr.dat");

  } else {
    printf("\nPlease input the prefix of the output g(r) files [gr]: ");
    input->read_stdin(str);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr == NULL) strcpy(str, "gr");
  }
  ptr = strtok(str, " \n\t\r\f");
  char *fprefix = new char [strlen(ptr)+1];
  strcpy(fprefix, ptr);

  // derived variables
  double delr = (rmax-rmin)/double(nbin-1);
  double rdr = 1./delr;
  
  // working space
  double **gr;
  memory->create(gr, nbin, 2, "gr");
#pragma omp parallel for default(shared)
  for (int i = 0; i < nbin; ++i) gr[i][0] = gr[i][1] = 0.;
  
  const double tpi = 8.*atan(1.);
  int nused = 0;

  // timer
  Timer * timer = new Timer();

  // now to compute g(r)
  if (job == 1){  // g(r) for selected atoms
    for (int img = istr; img <= iend; img += inc){ // loop over frames
      one = all[img];
      one->selection(srcsel);

      double dg = one->vol/(tpi*delr*one->nsel*one->nsel);
      double dh = 2./double(one->nsel);
  
      // need fractional coordinates
      one->car2dir();
  
      // set local variables
#pragma omp parallel for default(shared)
      for (int i = 1; i < one->natom; ++i){
        if (one->atsel[i] == 0) continue;

        for (int j = i+1; j <= one->natom; ++j){
          if (one->atsel[j] == 0) continue;
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
          int ibin = (r-rmin)*rdr;

          if (ibin >= 0 && ibin < nbin){
#pragma omp atomic
            gr[ibin][0] += dg;
#pragma omp atomic
            gr[ibin][1] += dh;
          }
        }
      }
      ++nused;

      if (min_mem) one->FreeVoro();
      if (flag_each == 1){
         sprintf(str, "%s_%d.dat", fprefix, img+1);
         char *fname = new char [strlen(str)+1];
         strcpy(fname, ptr);
         write_gr(rmin, delr, nbin, gr, nused, fname, header);
         nused = 0;
         printf("  PDF of frame %d evaluated and written to %s\n", img+1, fname);
#pragma omp parallel for default(shared)
         for (int i = 0; i < nbin; ++i) gr[i][0] = gr[i][1] = 0.;

         delete []fname;
      }
    } // end loop over frames

  } else {        // g(r) between two selections

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

      double dg = one->vol/(2.*tpi*delr*nsrc*one->nsel);
      double dh = 1./double(nsrc);

      // need fractional coordinates
      one->car2dir();

      // set local variables
#pragma omp parallel for default(shared)
      for (int i = 1; i <= one->natom; ++i){
        if (insrc[i] == 0) continue;

        for (int j = 1; j <= one->natom; ++j){
          if (one->atsel[j] == 0 || i==j) continue;
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
          int ibin = (r-rmin)*rdr;

          if (ibin >= 0 && ibin < nbin){
#pragma omp atomic
            gr[ibin][0] += dg;
#pragma omp atomic
            gr[ibin][1] += dh;
          }
        }
      }

      ++nused;
      if (min_mem) one->FreeVoro();
      if (flag_each == 1){
         sprintf(str, "%s_%d.dat", fprefix, img+1);
         char *fname = new char [strlen(str)+1];
         strcpy(fname, ptr);
         write_gr(rmin, delr, nbin, gr, nused, fname, header);
         nused = 0;
         printf("  PDF of frame %d evaluated and written to %s\n", img+1, fname);
#pragma omp parallel for default(shared)
         for (int i = 0; i < nbin; ++i) gr[i][0] = gr[i][1] = 0.;

         delete []fname;
      }
    }
  }

  if (flag_each == 0) write_gr(rmin, delr, nbin, gr, nused, fprefix, header);

  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  memory->destroy(gr);
  if (flag_each == 0) printf("\n%d images were used in the evaluation of g(r), which is written to %s\n", nused, fprefix);
  delete []fprefix;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 *
 *------------------------------------------------------------------------------*/
void Driver::write_gr(double rmin, double delr, int nbin, double **gr, int nused, char *fname, char *header)
{
  // normalize the g(r)
  double r = rmin - 0.5*delr;
  for (int i = 0; i < nbin; ++i){
    r += delr;
    gr[i][0] /= r*r*MAX(1,nused);
    gr[i][1] /= MAX(1,nused);
  }
  
  // output the result
  ConfirmOverwrite(fname);
  FILE *fp = fopen(fname,"w");
  fprintf(fp,"%s", header);
  fprintf(fp,"# r  g(r) int\n");
  r = rmin - 0.5*delr;
  double nsum = 0.;
  for (int i = 0; i < nbin; ++i){
    r += delr; nsum += gr[i][1];
    fprintf(fp,"%lg %lg %g\n", r, gr[i][0], nsum);
  }
  fclose(fp);

return;
}
/*------------------------------------------------------------------------------*/
