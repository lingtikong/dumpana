#include "driver.h"
#include "voro++.hh"
#include <vector>
#include <list>

/*------------------------------------------------------------------------------
 * Method to output selected Voronoi clusters as an xyz file
 * Atoms forming the Voronoi cluster are written together with the Voronoi
 * vertices, with an atomic type of X.
 *----------------------------------------------------------------------------*/
void Driver::OutputVoroCells( )
{
  int job = 0;
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("   Output  Voronoi  Clusters    ");
  for (int i = 0; i < 6; ++i) printf("====");
  printf("\nPlease select your desired job:\n");
  printf("  1. Output selected Voronoi cluster and vertices;\n");
  printf("  2. Output selected atoms in selected frames;\n");
  printf("  3. Output selected clusters in selected frames;\n");
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
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"vorocell.xyz"); ptr = strtok(str, " \n\t\r\f");}
  char *fname = new char[strlen(ptr)+1];
  strcpy(fname, ptr);
  FILE *fp = fopen(fname, "w");
  
  if (job == 1){ // output selected cluster and the Voronoi vertices
    voro::voronoicell_neighbor cell;
    double lmax = MAX(one->lx, MAX(one->ly, one->lz));
   
    int nclus = 0;
    std::vector<double> vertex, dx, dy, dz;
   
    // now to do the real job
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
   
      one->ComputeVoro(mins, weighted);
      one->selection(selcmd);
   
      if (one->nsel < 1){if (min_mem) one->FreeVoro(); continue;}
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
   
        for (int jj = 1; jj <= nei; ++jj){
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
      if (min_mem) one->FreeVoro();
    }
   
    vertex.clear(); dx.clear(); dy.clear(); dz.clear();
    printf("\n%d clusters were selected and written to: %s\n", nclus, fname);

  } else { // selected atoms or selected clusters
    std::list<int> outlist;
    // now to do the real job
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      one->ComputeVoro(mins, weighted);
      one->selection(selcmd);

      // get the selected list
      outlist.clear();
      for (int id = 1; id <= one->natom; ++id){
        if (one->atsel[id]){
          outlist.push_back(id);
          if (job == 3){
            for (int jj = 1; jj <= one->neilist[0][id]; ++jj) outlist.push_back(one->neilist[jj][id]);
          }
        }
      }
      if (min_mem) one->FreeVoro();

      outlist.sort(); outlist.unique();
      int nlist = outlist.size();
      if (nlist < 1) continue;

      one->dir2car();
      // output as xyz format
      fprintf(fp,"%d\n# Atoms/clusters in frame %d selected by: %s", nlist, one->iframe, selcmd);
      int idx = 0;
      while (! outlist.empty() ){
        int id = outlist.front(); outlist.pop_front();

        if (type2atnum == NULL){ // no elements assigned, print atomic type num as element
          fprintf(fp,"%d %lg %lg %lg", one->attyp[id], one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);
          if (idx < 3) fprintf(fp," crystal_vector %d %lg %lg %lg\n", idx+1, one->axis[idx][0], one->axis[idx][1], one->axis[idx][2]);
          else fprintf(fp,"\n");

        } else { // in case elements are assigned, print true element names
          char ename[3];
          element->Num2Name(type2atnum[one->attyp[id]], ename);
          if (idx < 3){
            fprintf(fp,"%2s %lg %lg %lg crystal_vector %d %lg %lg %lg\n", ename,
            one->atpos[id][0], one->atpos[id][1], one->atpos[id][2], idx+1,
            one->axis[idx][0], one->axis[idx][1], one->axis[idx][2]);
          } else {
            fprintf(fp,"%2s %lg %lg %lg\n", ename, one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);
          }
        }
        ++idx;
      }
    }

    outlist.clear();
  }

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

  fclose(fp);
  delete []fname;
return;
}

/*------------------------------------------------------------------------------*/
