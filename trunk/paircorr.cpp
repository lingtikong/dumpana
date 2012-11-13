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
    printf("  6. g(r) of selected atoms with certain voronoi index;\n");
    printf("  0. Return;\nYour choice [%d]: ", job);
    fgets(str,MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) job = atoi(ptr);
    printf("Your selection : %d\n", job);
    if (job < 1 || job > 6){
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
    } else if (job == 6){ // selected atoms of certain voronoi type

      one = all[0];
  
      // local common variables
      char selcmd[MAXLINE];
  
      // selection commands for atoms
      while (1){
        printf("\nPlease input the selection command for atoms, `h` for help [all]: ");
        if (count_words(fgets(str,MAXLINE,stdin)) > 0){
          strcpy(selcmd, str);
          char *ptr = strtok(str," \n\t\r\f");
          if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
        } else strcpy(selcmd,"all\n");
  
        // check the selection command on the first frame
        one->selection(selcmd); one->SelInfo();
        if (one->nsel < 1){
          printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
          if (count_words(fgets(str,MAXLINE,stdin)) > 0){
            char *ptr = strtok(str," \n\t\r\f");
            if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
          }
        }
        break;
      }

      // voro refinement info
      double voro_mins[3];
      voro_mins[0] = voro_mins[1] = 1.e-4;
      voro_mins[2] = 0.;
      
      printf("\nRefined Voronoi tesselation will be needed to procceed.\n");
      printf("Now please input your criterion for tiny surfaces [%g]: ", voro_mins[0]);
      fgets(str,MAXLINE, stdin); ptr = strtok(str, " \n\t\r\f");
      if (ptr) voro_mins[0] = atof(ptr);
      printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", voro_mins[0]);
    
      printf("Sometimes it might be desirable to keep a minimum # of neighbors when refining\n");
      printf("the Voronoi index, for example, keep at least 14 for a bcc lattice, 12 for hcp\n");
      printf("or fcc. If you prefer to do so, input a positive number now [%d]: ", int(voro_mins[2]));
      if (count_words(fgets(str,MAXLINE, stdin)) > 0){
        double dum = atof(strtok(str, " \n\t\r\f"));
        if (dum > 0.) voro_mins[2] = dum;
        printf("\nA minimum number of %d neighobrs will be kept no matter how tiny the surface is.\n", int(voro_mins[2]));
      }
    
      printf("\nPlease input your criterion for ultra short edge [%g]: ", voro_mins[1]);
      fgets(str,MAXLINE, stdin);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) voro_mins[1] = atof(ptr);
      printf("Edges whose length takes less ratio than %lg will be skipped!\n\n", voro_mins[1]);

      // set clusters to analyse
      set<std::string> voroset; std::string vindex;
      printf("\nPlease input the Voronoi index of the desired clusters, e.g., 0,6,0,8.\n");
      printf("If multiple indices are wanted, separate them by space: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        printf("\n  Clusters centered on Atoms with Voronoi indices: %s  will be analysed.\n", str);
     
        char *ptr = strtok(str," \n\t\r\f");
        while (ptr){
          vindex.assign(ptr);
          voroset.insert(vindex);
     
          ptr = strtok(NULL," \n\t\r\f");
        }
      }
    
      // header of g(r)
      sprintf(header,"# g(r) for atoms selected by: %s# with voronoi indices: ", selcmd);
      for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); it++) sprintf(header,"%s %s", header, it->c_str());
      sprintf(header,"%s\n", header);

      // work space for Voronoi
      int nmax = 24, **neilist, *cenlist;
      int natom = one->natom;
      cenlist = memory->create(cenlist, natom+1, "cenlist");
      neilist = memory->create(neilist, nmax+1, natom+1, "neilist");

      for (int img = istr; img <= iend; img += inc){ // loop over frames
        one = all[img];
        one->selection(selcmd);

        // work space for Voronoi
        if (one->natom > natom){
          cenlist = memory->grow(cenlist, natom+1, "cenlist");

          memory->destroy(neilist);
          neilist = memory->create(neilist, nmax+1, natom+1, "neilist");
        }
        natom = one->natom;
        for (int ii=0; ii<=natom; ii++) cenlist[ii] = 0;
    
        map<int,std::string> voroindex; voroindex.clear();
        // compute the neighbor list and voro info
        FEFF_voro(voroset, voro_mins, nmax, neilist, cenlist, voroindex);
        voroindex.clear(); // they are not needed here
    
        // check the # of selected clusters
        int nc = 0;
        for (int ii=1; ii<=natom; ii++) nc += cenlist[ii];
        if (nc > 1){
          nused++;
          const double dg = one->vol/(2.*tpi*delr*nc*nc);
  
          // need fractional coordinates
          one->car2dir();
  
          // set local variables
          for (int i=1; i<= one->natom; i++){
            if (cenlist[i] == 0) continue;
  
            for (int j=1; j<= one->natom; j++){
              if (cenlist[j] == 0 || i==j) continue;
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
      } // end loop over frames for job == 6
      memory->destroy(cenlist);
      memory->destroy(neilist);
    }
  
    // normalize the g(r)
    double r = rmin - 0.5*delr;
    for (int i=0; i<nbin; i++){
      r += delr;
      gr[i] /= r*r*MAX(1,nused);
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
