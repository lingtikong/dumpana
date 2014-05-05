#include "driver.h"
#include "voro++.hh"
#include <vector>
#include <list>
#include <set>

/*------------------------------------------------------------------------------
 * Method to output selected frames in BGF format with defined properties
 *----------------------------------------------------------------------------*/
void Driver::writebgf()
{
  int job = 0;
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 4; ++i) printf("====");
  printf("   Output in BGF format with defined property   ");
  for (int i = 0; i < 4; ++i) printf("====");
  printf("\nPlease select your way to define property:\n");
  printf("  1. Voronoi volume;\n");
  printf("  2. fraction of N-edged surface;\n");
  printf("  3. Chemical concentration;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 3){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  printf("\n");
  
  // thresholds for surface and edges
  set_cutoffs(1);
  one = all[istr];

  // selection of atoms for each frame
  char selcmd[MAXLINE];
  while (1){
    printf("\nPlease input your selection command of central atoms, `h` for help [all]: ");
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

  // property related
  int nedge = 0;
  char title[MAXLINE];
  double rcut = -1., rcut2;
  double *f5, *chem;
  std::set<int> solute;
  f5 = NULL; chem = NULL; solute.clear();

  if (job == 1){          // Voronoi volume as property
    sprintf(title,"Voronoi volume as atomic property.");

  } else if (job == 2){   // N-edged faces as property
    printf("\nPlease input the # of edge for the desired surfaces: ");
    while (count_words(fgets(str,MAXLINE,stdin)) < 1);
    nedge = atoi(strtok(str," \n\t\r\f"));
    if (nedge < 3 || nedge > 6){
      printf("The # of edges must be within [3, 6].\n");
      for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
      return;
    }
    sprintf(title,"Fraction of %d-edged faces as atomic property.", nedge);

  } else if (job == 3){
    int ijob = 1;
    printf("\nPlease select the way to evaluate local chemical concentration:\n");
    printf(" 1. Based on Voronoi cluster;\n");
    printf(" 2. Within defined radius;\n");
    printf("Your choice [%d]: ", ijob);
    fgets(str,MAXLINE, stdin);
    ijob = atoi(strtok(str, " \n\t\r\f"));
    if (ijob == 2){
      printf("Please input the radius to count neighbors, please note that\n");
      printf("the radius should exceed the 3rd Voronoi shell: ");
      while (count_words(fgets(str,MAXLINE, stdin)) < 1);

      rcut = fabs(atof(strtok(str, " \n\t\r\f")));
      rcut2 = rcut*rcut;
    }

    printf("Please input the atomic type ID(s) of the solutes [1]: ");
    if (count_words(fgets(str,MAXLINE, stdin)) > 1){
      ptr = strtok(str, " \n\t\r\f");
      while (ptr){
        solute.insert(atoi(ptr));
        ptr = strtok(NULL, " \n\t\r\f");
      }

    } else {
      solute.insert(1);
    }

    sprintf(title, "Local concentration of");
    for (std::set<int>::iterator it = solute.begin(); it != solute.end(); ++it)
    sprintf(title, "%s type-%d,", title, *it);
    sprintf(title, "%s as atomic property.", title);
  }

  // output file name
  printf("\nPlease input the output file name [atomcfg.bgf]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"atomcfg.bgf"); ptr = strtok(str, " \n\t\r\f");}
  char *fname = new char[strlen(ptr)+1];
  strcpy(fname, ptr);
  FILE *fp = fopen(fname, "w");
  
  // map atomic type to elements
  if (type2atnum == NULL){
    printf("Mapping of atomic types to actual elements are needed to writen in BGF format.\n");
    MapType2Elem(1, one->ntype); printf("\n");
  }

  // now to do the real job
  for (int img = istr; img <= iend; img += inc){
    one = all[img];
   
    one->ComputeVoro(mins);
    one->selection(selcmd);
 
    if (one->nsel < 1){if (min_mem) one->FreeVoro(); continue;}

    // Vorono volume as properties
    if (job == 1) one->prop = one->volume;

    // Fraction of N-edge surface as properties
    if (job == 2){
      memory->destroy(f5);
      memory->create(f5, one->natom+1, "f5");
      one->prop = f5;

      for (int id = 1; id <= one->natom; ++id){
        if (one->atsel[id] == 0) continue;
        int nv[7]; for (int ii = 0; ii < 7; ++ii) nv[ii] = 0;
        strcpy(str, one->voro[id].c_str());
        ptr = strtok(str, ",");
        for (int ii = 3; ii < 7; ++ii){
          nv[ii] = atoi(ptr);
          ptr = strtok(NULL, ",");
        }
        f5[id] = double(nv[nedge])/double(one->neilist[0][id]);
      }
    }

    // Local chemcial concentration as properties
    if (job == 3){
      memory->destroy(chem);
      memory->create(chem, one->natom+1, "chem");
      one->prop = chem;

      // based on Voronoi cluster
      if (rcut < 0.){
        for (int id = 1; id <= one->natom; ++id){
          if (one->atsel[id] == 0) continue;
          int nhit = 0;
          int ip = one->attyp[id];
          if (solute.find(ip) != solute.end()) ++nhit;
        
          for (int jj = 1; jj <= one->neilist[0][id]; ++jj){
            int jd = one->neilist[jj][id];
            int jp = one->attyp[jd];
            if (solute.find(jd) != solute.end()) ++nhit;
          }
          chem[id] = double(nhit)/double(one->neilist[0][id] + 1);
        }

      // Based on radius, should not exceed 3rd Voronoi shell
      } else {
        // need fractional coordinate
        one->car2dir();
        for (int id = 1; id <= one->natom; ++id){
          if (one->atsel[id] == 0) continue;

          std::list<int> nlist;
          nlist.clear();
          for (int jj = 1; jj <= one->neilist[0][id]; ++jj){
            int jd = one->neilist[jj][id];
            nlist.push_back(jd);
            for (int kk = 1; kk <= one->neilist[0][jd]; ++kk){
              int kd = one->neilist[kk][jd];
              nlist.push_back(kd);
              for (int ll = 1; ll <= one->neilist[0][kd]; ++ll){
                int ld = one->neilist[ll][kd];
                nlist.push_back(ld);
              }
            }
          }
          nlist.sort(); nlist.unique();

          int nnei = 0, nhit = 0;
          while (! nlist.empty() ){
            int jd = nlist.front(); nlist.pop_front();
            if (jd == id) continue;

            double xij[3];
            for (int idim = 0; idim < 3; ++idim){
              xij[idim] = one->atpos[jd][idim] - one->atpos[id][idim];
              while (xij[idim] > 0.5) xij[idim] -= 1.;
              while (xij[idim] <-0.5) xij[idim] += 1.;
            }
            xij[0] = xij[0] * one->lx + xij[1] * one->xy + xij[2] * one->xz;
            xij[1] = xij[1] * one->ly + xij[2] * one->yz;
            xij[2] = xij[2] * one->lz;
            double rij2 = xij[0]*xij[0] + xij[1]*xij[1] + xij[2]*xij[2];
            if (rij2 <= rcut2){
              ++nnei;
              if (solute.find(one->attyp[jd]) != solute.end()) ++nhit;
            }
          }
          if (solute.find(one->attyp[id]) != solute.end()) ++nhit;
          chem[id] = double(nhit)/double(nnei+1);
        }
      }
    }

    // output the results in BGF format
    // need cartesian coordinate
    one->dir2car();
    const double rad2ang = 45./atan(1.);
    double lx = one->lx, ly = sqrt(one->ly*one->ly + one->xy*one->xy);
    double lz = sqrt(one->xz*one->xz + one->yz*one->yz + one->lz*one->lz);
    double alpha = acos((one->xy*one->xz + one->ly*one->yz)/(ly*lz))*rad2ang;
    double beta  = acos((one->lx*one->xz)/(lx*lz))*rad2ang;
    double gamma = acos((one->lx*one->xy)/(lx*ly))*rad2ang;
    fprintf(fp,"BIOGRF 200\nDESCRP %s\n", title);
    fprintf(fp,"REMARK \nFFIELD EAM\n");
    fprintf(fp,"CRYSTX %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", lx, ly, lz, alpha, beta, gamma);
    fprintf(fp,"FORMAT ATOM (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)\n");
    int idum = 0;
    for (int id = 1; id <= one->natom; ++id){
      if (one->atsel[id] == 0) continue;

      char ename[3];
      element->Num2Name(type2atnum[one->attyp[id]], ename);

      ++idum;
      fprintf(fp, "HETATM %5d %5s            %10.5f%10.5f%10.5f %5s%3d%2d %8.5f 0   0  0.\n", idum, ename,
      one->atpos[id][0], one->atpos[id][1], one->atpos[id][2], ename, 0, 0, one->prop[id]);
    }
    fprintf(fp,"END\n");

  } // end of loop over frames
    
  memory->destroy(f5);
  memory->destroy(chem);
  fclose(fp);
  delete []fname;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------*/
