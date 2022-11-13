#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the property pair correlation function of the system
 *----------------------------------------------------------------------------*/
void Driver::property_pc()
{
  char str[MAXLINE], header[MAXLINE];
  double rmin = 0., rmax = 0.5*MIN(MIN(all[0]->lx,all[0]->ly), all[0]->lz);
  int nbin = 201;

  // menu
  int job = 1;
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("  Property Pair  Correlation  Functions  ");
  for (int i = 0; i < 6; ++i) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i = 0; i < 20; ++i) printf("----"); printf("\n");
  printf("  1. C(r) of selected atoms;\n");
  printf("  2. C(r) between two selections;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
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

    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(srcsel, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(srcsel,"all\n");
  
    // check the selection command on the first frame
    one->selection(srcsel); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  if (job == 1) sprintf(header,"# C(r) of properties for atoms selected by: %s", srcsel);
  else {
    // atoms as neighbors
    while (1){
      printf("\nPlease input the selection command for neighbors, `h` for help [all]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        strcpy(dessel, str);
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
      } else strcpy(dessel,"all\n");

      // check the selection command on the first frame
      one->selection(dessel); one->SelInfo();
      if (one->nsel < 1){
        printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
        if (count_words(fgets(str,MAXLINE,stdin)) > 0){
          char *ptr = strtok(str," \n\t\r\f");
          if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
        }
      }
      break;
    }
    memory->create(insrc, one->natom+1, "insrc");
    sprintf(header,"# C(r) between two selections:\n# source atoms: %s# neighbors: %s", srcsel, dessel);
  }
  
  // property 
  int pid = 0;
  printf("\nThe available properties are: ");
  for (int i = 0; i < one->prop_label.size(); ++i) printf("%d) %s; ", i+1, one->prop_label[i].c_str());
  printf("\nPlease input the id of the property to be used [%d]: ", pid+1);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) pid = atoi(ptr) - 1;
  pid = MIN(MAX(pid, 0), one->prop_label.size());
  printf("The property chosen is: %s\n", one->prop_label[pid].c_str());

  // bounds of g(r)
  printf("\nPlease input the lower bound of C(r) [%g]: ", rmin);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) rmin = atof(ptr);
  rmin = MAX(rmin, 0.);

  printf("Please input the upper bound of C(r) [%g]: ", rmax);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) rmax = atof(ptr);
  rmax = MAX(rmax,rmin+1.);

  printf("Please input the # of points of C(r) [%d]: ", nbin);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) nbin = atoi(ptr);
  nbin = MAX(2,nbin);
  
  // output average C(r) or for each frame
  printf("\nWould you like to output:\n");
  printf("  1. the average C(r) for all frames used;\n");
  printf("  2. a separate  C(r) for each frame;\n");
  printf("Your choice [1]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  int flag_each = 0;
  if (ptr) flag_each = atoi(ptr)-1;

  // output file name or prefix
  if (flag_each == 0){
    printf("\nPlease input the file to output C(r) [pc.dat]: ");
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr == NULL) strcpy(str, "pc.dat");

  } else {
    printf("\nPlease input the prefix of the output C(r) files [pc]: ");
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr == NULL) strcpy(str, "pc");
  }
  ptr = strtok(str, " \n\t\r\f");
  char *fprefix = new char [strlen(ptr)+1];
  strcpy(fprefix, ptr);

  // derived variables
  double delr = (rmax-rmin)/double(nbin-1);
  double rdr = 1./delr;
  
  // working space
  double **gr, *one_prop, *ave_prop, *std_prop;
  memory->create(gr, nbin, 2, "gr");
  one_prop = ave_prop = std_prop = NULL;

#pragma omp parallel for default(shared)
  for (int i = 0; i < nbin; ++i) gr[i][0] = gr[i][1] = 0.;
  
  const double tpi = 8.*atan(1.);
  int nused = 0;

  // timer
  Timer * timer = new Timer();

  // now to compute C(r)
  if (job == 1){  // C(r) for selected atoms
    for (int img = istr; img <= iend; img += inc){ // loop over frames
      one = all[img];
      one->selection(srcsel);
      if (one->nsel < 1) continue;

      // prepare the property
      memory->grow(one_prop, one->natom+1, "one_prop");
      memory->grow(ave_prop, one->ntype+1, "ave_prop");
      memory->grow(std_prop, one->ntype+1, "std_prop");
      for (int ip = 0; ip <= one->ntype; ++ip) ave_prop[ip] = std_prop[ip] = 0.;
      for (int i = 1; i < one->natom+1; ++i){
        int ip = one->attyp[i];
        ave_prop[ip] += one->atprop[i][pid];
      }
      for (int ip = 1; ip <= one->ntype; ++ip) ave_prop[ip] /= MAX(1., one->numtype[ip]);
      for (int i = 1; i < one->natom+1; ++i){
        int ip = one->attyp[i];
        one_prop[i] = one->atprop[i][pid] - ave_prop[ip];
        std_prop[ip] += one_prop[i] * one_prop[i];
      }
      for (int ip = 1; ip <= one->ntype; ++ip) std_prop[ip] = sqrt(std_prop[ip]/MAX(1.,one->numtype[ip]));

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
            gr[ibin][0] += dg * one_prop[i] * one_prop[j] / (std_prop[one->attyp[i]] * std_prop[one->attyp[j]]);
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
         strcpy(fname, str);
         write_prop_pc(rmin, delr, nbin, gr, nused, fname, header);
         nused = 0;
         printf("  PDF of frame %d evaluated and written to %s\n", img+1, fname);
         delete []fname;
#pragma omp parallel for default(shared)
         for (int i = 0; i < nbin; ++i) gr[i][0] = gr[i][1] = 0.;
      }
    } // end loop over frames

  } else {        // C(r) between two selections

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

      // prepare the property
      memory->grow(one_prop, one->natom+1, "one_prop");
      memory->grow(ave_prop, one->ntype+1, "ave_prop");
      memory->grow(std_prop, one->ntype+1, "std_prop");
      for (int ip = 0; ip <= one->ntype; ++ip) ave_prop[ip] = std_prop[ip] = 0.;
      for (int i = 1; i <= one->natom; ++i){
        int ip = one->attyp[i];
        ave_prop[ip] += one->atprop[i][pid];
      }
      for (int ip = 1; ip <= one->ntype; ++ip) ave_prop[ip] /= MAX(1.,one->numtype[ip]);
      for (int i = 1; i < one->natom+1; ++i){
        int ip = one->attyp[i];
        one_prop[i] = one->atprop[i][pid] - ave_prop[ip];
        std_prop[ip] += one_prop[i] * one_prop[i];
      }
      for (int ip = 1; ip <= one->ntype; ++ip) std_prop[ip] = sqrt(std_prop[ip]/MAX(1.,one->numtype[ip]));

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
            gr[ibin][0] += dg * one_prop[i] * one_prop[j] / (std_prop[one->attyp[i]] * std_prop[one->attyp[j]]);
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
         strcpy(fname, str);
         write_prop_pc(rmin, delr, nbin, gr, nused, fname, header);
         nused = 0;
         printf("  PDF of frame %d evaluated and written to %s\n", img+1, fname);
         delete []fname;
#pragma omp parallel for default(shared)
         for (int i = 0; i < nbin; ++i) gr[i][0] = gr[i][1] = 0.;
      }
    }
  }

  if (flag_each == 0) write_prop_pc(rmin, delr, nbin, gr, nused, fprefix, header);

  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  memory->destroy(gr);
  memory->destroy(insrc);
  memory->destroy(one_prop);
  memory->destroy(ave_prop);
  memory->destroy(std_prop);

  if (flag_each == 0) printf("\n%d images were used in the evaluation of C(r), which is written to %s\n", nused, fprefix);
  delete []fprefix;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 *
 *------------------------------------------------------------------------------*/
void Driver::write_prop_pc(double rmin, double delr, int nbin, double **gr, int nused, char *fname, char *header)
{
  // normalize the C(r)
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
  fprintf(fp,"# r  C(r) int_gr\n");
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
