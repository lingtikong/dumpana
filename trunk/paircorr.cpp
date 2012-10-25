#include "driver.h"
#include "math.h"

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to compute the pair correlation function of the system
 *----------------------------------------------------------------------------*/
void Driver::paircorr()
{
  char str[MAXLINE], header[MAXLINE];
  double rmin = 0., rmax = 0.5*MIN(MIN(all[0]->lx,all[0]->ly), all[0]->lz);
  int nbin = 201;

  int job = 1;
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  while (1){
    printf("\nPlease select your desired job:\n");
    for (int i=0; i<20; i++) printf("----"); printf("\n");
    printf("  1. g(r) of the system;\n");
    printf("  2. g(r) of an atomic type;\n");
    printf("  3. g(r) based on an atomic type;\n");
    printf("  4. g(r) between two atomic types;\n");
    printf("  5. g(r) between two selections;\n");
    printf("  0. Return;\nYour choice [%d]: ", job);
    fgets(str,MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) job = atoi(ptr);
    printf("Your selection : %d\n", job);
    if (job < 1 || job > 5){
      for (int i=0; i<20; i++) printf("===="); printf("\n");
      return;
    }
    printf("\n");
  
    printf("Please input the lower bound of g(r) [%g]: ", rmin);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) rmin = atof(ptr);
    rmin = MAX(rmin, 0.);
  
    printf("Please input the upper bound of g(r) [%g]: ", rmax);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) rmax = atof(ptr);
    rmax = MAX(rmax,rmin+1.);
  
    printf("Please input the # of points of g(r) [%d]: ", nbin);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) nbin = atoi(ptr);
    nbin = MAX(2,nbin);
  
    double delr = (rmax-rmin)/double(nbin-1);
    double rdr = 1./delr;
  
    double *gr;
    gr = memory->create(gr,nbin,"gr");
    for (int i=0; i<nbin; i++) gr[i] = 0.;
  
    const double tpi = 8.*atan(1.);
    int nused = 0;
  
    // now to do the real job
    if (job == 1){
  
      sprintf(header,"# g(r) for all, from %g to %g with %d points.\n", rmin, rmax, nbin);
      for (int img = istr; img <= iend; img += inc){
        one = all[img];
    
        // need fractional coordinates
        one->car2dir();
  
        nused++;
        const double dg = one->vol/(tpi*delr*one->natom*one->natom);
  
        // compute G(r)
        for (int i=1; i<= one->natom; i++)
        for (int j=i+1; j<= one->natom; j++){
          double dx[3], dr[3];
          for (int idim=0; idim<3; idim++){
            dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
            while (dx[idim] > 0.5) dx[idim] -= 1.;
            while (dx[idim] <-0.5) dx[idim] += 1.;
          }
          dr[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
          dr[1] = dx[1]*one->ly + dx[2]*one->yz;
          dr[2] = dx[2]*one->lz;
          double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
          int ibin = (r-rmin)*rdr;
          if (ibin >= 0 && ibin<nbin) gr[ibin] += dg;
        }
      }
  
    } else if (job == 2){
  
      // ask for the atomic type to get g(r)
      int srctype = 1;
      printf("\nThere are %d atomic types in the first image,\n", all[0]->ntype);
      printf("please input the atomic type # to get g(r) [%d]: ", srctype);
      fgets(str,MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) srctype = atoi(ptr);
      srctype = MIN(MAX(1,srctype), all[0]->ntype);
  
      if (type2atnum == NULL)
        sprintf(header,"# partial g(r) for type %d from %g to %g with %d points.\n", srctype, rmin, rmax, nbin);
      else {
        char ename[3];
        element->Num2Name(type2atnum[srctype], ename);
        sprintf(header,"# %s-%s partial g(r) from %g to %g with %d points.\n", ename, ename, rmin, rmax, nbin);
      }
  
      for (int img = istr; img <= iend; img += inc){
        one = all[img];
    
        // need fractional coordinates
        one->car2dir();
  
        nused++;
        const double dg = one->vol/(tpi*delr*one->numtype[srctype]*one->numtype[srctype]);
  
        // set local variables
        for (int i=1; i<= one->natom; i++){
          if (one->attyp[i] != srctype) continue;
  
          for (int j=i+1; j<= one->natom; j++){
            if (one->attyp[j] != srctype||j == i) continue;
            double dx[3], dr[3];
            for (int idim=0; idim<3; idim++){
              dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
              while (dx[idim] > 0.5) dx[idim] -= 1.;
              while (dx[idim] <-0.5) dx[idim] += 1.;
            }
            dr[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
            dr[1] = dx[1]*one->ly + dx[2]*one->yz;
            dr[2] = dx[2]*one->lz;
            double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
            int ibin = (r-rmin)*rdr;
            if (ibin >= 0 && ibin<nbin) gr[ibin] += dg;
          }
        }
      }
  
    } else if (job == 3){
  
      // ask for the atomic type to get g(r)
      int srctype = 1;
      printf("\nThere are %d atomic types in the first image,\n", all[0]->ntype);
      printf("please input the source type to get g(r) [%d]: ", srctype);
      fgets(str,MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) srctype = atoi(ptr);
      srctype = MIN(MAX(1,srctype), all[0]->ntype);
  
      if (type2atnum == NULL)
        sprintf(header,"# partial g(r) centered on type %d from %g to %g with %d points.\n", srctype, rmin, rmax, nbin);
      else {
        char ename[3];
        element->Num2Name(type2atnum[srctype], ename);
        sprintf(header,"# %s-all partial g(r) from %g to %g with %d points.\n", ename, rmin, rmax, nbin);
      }
  
      for (int img = istr; img <= iend; img += inc){
        one = all[img];
    
        // need fractional coordinates
        one->car2dir();
  
        nused++;
        const double dg = one->vol/(2.*tpi*delr*one->numtype[srctype]*one->natom);
  
        // set local variables
        for (int i=1; i<= one->natom; i++){
          if (one->attyp[i] != srctype) continue;
  
          for (int j=1; j<= one->natom; j++){
            if (j == i) continue;
            double dx[3], dr[3];
            for (int idim=0; idim<3; idim++){
              dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
              while (dx[idim] > 0.5) dx[idim] -= 1.;
              while (dx[idim] <-0.5) dx[idim] += 1.;
            }
            dr[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
            dr[1] = dx[1]*one->ly + dx[2]*one->yz;
            dr[2] = dx[2]*one->lz;
            double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
            int ibin = (r-rmin)*rdr;
            if (ibin >= 0 && ibin<nbin) gr[ibin] += dg;
          }
        }
      }
  
    } else if (job == 4){
  
      // ask for the atomic type to get g(r)
      int srctype = 1, destype = 1;
      printf("\nThere are %d atomic types in the first image,\n", all[0]->ntype);
      printf("please input the source type to get g(r) [%d]: ", srctype);
      fgets(str,MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) srctype = atoi(ptr);
      srctype = MIN(MAX(1,srctype), all[0]->ntype);
  
      printf("please input the atomic type of the neighbors [%d]: ", destype);
      fgets(str,MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) destype = atoi(ptr);
      destype = MIN(MAX(1,destype), all[0]->ntype);
  
      if (type2atnum == NULL)
        sprintf(header,"# partial g(r) for pair %d-%d from %g to %g with %d points.\n", srctype, destype, rmin, rmax, nbin);
      else {
        char sname[3], dname[3];
        element->Num2Name(type2atnum[srctype], sname);
        element->Num2Name(type2atnum[destype], dname);
        sprintf(header,"# %s-%s partial g(r) from %g to %g with %d points.\n", sname, dname, rmin, rmax, nbin);
      }
  
      for (int img = istr; img <= iend; img += inc){
        one = all[img];
    
        // need fractional coordinates
        one->car2dir();
  
        nused++;
        const double dg = one->vol/(2.*tpi*delr*one->numtype[srctype]*one->numtype[destype]);
  
        // set local variables
        for (int i=1; i<= one->natom; i++){
          if (one->attyp[i] != srctype) continue;
  
          for (int j=1; j<= one->natom; j++){
            if (one->attyp[j] != destype || i==j) continue;
            double dx[3], dr[3];
            for (int idim=0; idim<3; idim++){
              dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
              while (dx[idim] > 0.5) dx[idim] -= 1.;
              while (dx[idim] <-0.5) dx[idim] += 1.;
            }
            dr[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
            dr[1] = dx[1]*one->ly + dx[2]*one->yz;
            dr[2] = dx[2]*one->lz;
            double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
            int ibin = (r-rmin)*rdr;
            if (ibin >= 0 && ibin<nbin) gr[ibin] += dg;
          }
        }
      }
    } else if (job == 5){ // first selection as source, second selection as neighbors
  
      one = all[0];
  
      // local common variables
      char srcsel[MAXLINE], dessel[MAXLINE];
  
      // selection commands for atoms
      // atoms as source
      while (1){
        printf("\nPlease input the selection command for source atoms, `h` for help [all]:");
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
      // atoms as neighbors
      while (1){
        printf("\nPlease input the selection command for neighbors, `h` for help [all]:");
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
      int *insrc = memory->create(insrc, one->natom+1, "insrc");
  
      sprintf(header,"# g(r) between two selections:\n# source atoms: %s# neighbors: %s", srcsel, dessel);
  
      for (int img = istr; img <= iend; img += inc){
        one = all[img];
    
        // select atoms as source
        one->selection(srcsel);
        if (one->nsel < 1) continue;
        insrc = memory->grow(insrc, one->natom+1, "insrc");
        for (int ii=1; ii<= one->natom; ii++) insrc[ii] = one->atsel[ii];
        int nsrc = one->nsel;
  
        // select atoms as neighbors
        one->selection(dessel);
        if (one->nsel < 1) continue;
  
        nused++;
        const double dg = one->vol/(2.*tpi*delr*nsrc*one->nsel);
  
        // need fractional coordinates
        one->car2dir();
  
        // set local variables
        for (int i=1; i<= one->natom; i++){
          if (insrc[i] == 0) continue;
  
          for (int j=1; j<= one->natom; j++){
            if (one->atsel[j] == 0 || i==j) continue;
            double dx[3], dr[3];
            for (int idim=0; idim<3; idim++){
              dx[idim] = one->atpos[j][idim] - one->atpos[i][idim];
              while (dx[idim] > 0.5) dx[idim] -= 1.;
              while (dx[idim] <-0.5) dx[idim] += 1.;
            }
            dr[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
            dr[1] = dx[1]*one->ly + dx[2]*one->yz;
            dr[2] = dx[2]*one->lz;
            double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
            int ibin = (r-rmin)*rdr;
            if (ibin >= 0 && ibin<nbin) gr[ibin] += dg;
          }
        }
      }
  
      memory->destroy(insrc);
    }
  
    // normalize the g(r)
    double r = rmin - 0.5*delr;
    for (int i=0; i<nbin; i++){
      r += delr;
      gr[i] /= r*r*nused;
    }
  
    // output the result
    printf("\nPlease input the file to output g(r) [gr.dat]: ");
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr == NULL) strcpy(str, "gr.dat");
    ptr = strtok(str, " \n\t\r\f");
    FILE *fp = fopen(ptr,"w");
    fprintf(fp,"%s", header);
    fprintf(fp,"# r  g(r)\n");
    r = rmin - 0.5*delr;
    for (int i=0; i<nbin; i++){
      r += delr;
      fprintf(fp,"%lg %lg\n", r, gr[i]);
    }
    fclose(fp);
  
    memory->destroy(gr);
    printf("\n%d images were used in the evaluation of g(r), which is written to %s\n", nused, ptr);
    job = 0;
  }

return;
}
