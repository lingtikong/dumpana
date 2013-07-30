#include "driver.h"
#include "voro++.hh"
#include <vector>

/*------------------------------------------------------------------------------
 * Method to output selected Voronoi clusters as an xyz file
 * Atoms forming the Voronoi cluster are written together with the Voronoi
 * vertices, with an atomic type of X.
 *----------------------------------------------------------------------------*/
void Driver::OutputVoroCells( )
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<6; ++i) printf("====");
  printf("   Output  Voronoi  Clusters    ");
  for (int i=0; i<6; ++i) printf("====");

  // thresholds for surface and edges
  set_cutoffs(1);

  // show relevant info if weighted Voronoi is used
  one = all[istr]; ShowRadius4Voro();

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

  printf("\nPlease input the output file name [vorocell.xyz]: ");
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"vorocell.xyz"); ptr = strtok(str, " \n\t\r\f");}
  char *fname = new char[strlen(ptr)+1];
  strcpy(fname, ptr);
  FILE *fp = fopen(fname, "w");
  
  voro::voronoicell_neighbor cell;
  double lmax = MAX(one->lx, MAX(one->ly, one->lz));

  int nclus = 0;
  std::vector<double> vertex, dx, dy, dz;

  // now to do the real job
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    one->ComputeVoro(mins, weighted);
    one->selection(selcmd);

    if (one->nsel < 1) continue;
    one->dir2car();

    for (int id = 1; id <= one->natom; ++id){
      if (one->atsel[id] == 0) continue;
      double xpos = one->atpos[id][0];
      double ypos = one->atpos[id][1];
      double zpos = one->atpos[id][2];
      int nei = one->neilist[0][id];
      
      cell.init(-lmax, lmax, -lmax, lmax, -lmax, lmax);
      dx.clear(); dy.clear(); dz.clear(); vertex.clear();
      dx.resize(nei+1); dy.resize(nei+1); dz.resize(nei+1);

      for (int jj=1; jj<= nei; ++jj){
        int jd = one->neilist[jj][id];
  
        // apply pbc
        double xij = one->atpos[jd][0] - xpos;
        double yij = one->atpos[jd][1] - ypos;
        double zij = one->atpos[jd][2] - zpos;
        
        one->ApplyPBC(xij, yij, zij);
        dx[jj] = xij; dy[jj] = yij; dz[jj] = zij;

        if (weighted){
          double scale = weighted[one->attyp[id]]/(weighted[one->attyp[id]] + weighted[one->attyp[jd]]);
          scale += scale;
          xij *= scale;
          yij *= scale;
          zij *= scale;
        }
  
        cell.nplane(xij,yij,zij,jd);
      }
      cell.vertices(0.,0.,0.,vertex);
      int nv = vertex.size()/3;

      fprintf(fp,"%d\nAtom %d of frame %d: <%s>, [%lg %lg %lg]\n", nv+nei+1, id, one->iframe, one->voro[id].c_str(), xpos, ypos, zpos);
      
      if (type2atnum == NULL){ // no elements assigned, print atomic type num as element
        fprintf(fp,"%d 0. 0.  0.\n", one->attyp[id]);
        for (int jj = 1; jj <= nei; ++jj){
          int jd = one->neilist[jj][id];
          fprintf(fp,"%d %lg %lg %lg\n", one->attyp[jd], dx[jj], dy[jj], dz[jj]);
        }
        for (int kk = 0; kk < 3*nv; kk += 3) fprintf(fp,"%d %lg %lg %lg\n", one->ntype+1, vertex[kk], vertex[kk+1], vertex[kk+2]);

      } else { // in case elements are assigned, print true element names
        char ename[3];
        element->Num2Name(type2atnum[one->attyp[id]], ename);
        fprintf(fp,"%2s 0. 0. 0.\n", ename);
        for (int jj = 1; jj <= nei; ++jj){
          int jd = one->neilist[jj][id];
          element->Num2Name(type2atnum[one->attyp[jd]], ename);
          fprintf(fp,"%2s %lg %lg %lg\n", ename, dx[jj], dy[jj], dz[jj]);
        }
        for (int kk = 0; kk < 3*nv; kk += 3) fprintf(fp,"X %lg %lg %lg\n", vertex[kk], vertex[kk+1], vertex[kk+2]);
      }
      ++nclus;
    }
  }

  fclose(fp);
  vertex.clear(); dx.clear(); dy.clear(); dz.clear();

  printf("\n%d clusters were selected and written to: %s\n", nclus, fname);
  for (int i=0; i<20; ++i) printf("===="); printf("\n");
  delete []fname;
return;
}

/*------------------------------------------------------------------------------*/
