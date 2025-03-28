#include "driver.h"
#include "voro++.hh"
#include <vector>
#include <list>
#include <set>
#include "timer.h"

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
  printf("  4. Output selected atoms and their neighboring atoms/molecules;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  input->read_stdin(str);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 4){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  printf("\n");
  
  // thresholds for surface and edges
  if (job < 4) set_cutoffs(1);
  one = all[istr];

  // selection of atoms for each frame
  char selcmd[MAXLINE];
  while (1){
    printf("\nPlease input your selection command of central atoms, `h` for help [all]: ");
    input->read_stdin(str);
    if (count_words(str) > 0){
      strcpy(selcmd, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(selcmd,"all\n");

    // check the selection command on the first frame
    one->selection(selcmd); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      input->read_stdin(str);
      if (count_words(str) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y") != 0) continue;
      }
    }

    break;
  }

  int flabel = 0, pid = 0;
  // ask whether to label the central atom or not
  if (job == 3){
    printf("\nWould you like to label the central atoms by adding an extra column\n");
    printf("to the output file? (y/n)[n]: ");
    input->read_stdin(str);
    if (count_words(str) > 0){
      ptr = strtok(str, " \n\t\r\f");
      if (strcmp(ptr, "y")==0 || strcmp(ptr, "Y")==0 ) flabel = 1;
    }

  } else if (job == 4 && one->prop_label.size() > 0) { // property id of neighbors

    printf("\nThe available properties are: ");
    for (int i = 0; i < one->prop_label.size(); ++i) printf("%d) %s; ", i+1, one->prop_label[i].c_str());
    printf("\nPlease input the id of the property to be shared by neighbors [%d]: ", pid+1);
    input->read_stdin(str);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) pid = atoi(ptr) - 1;
    pid = MIN(MAX(pid, 0), one->prop_label.size());
    printf("The property chosen is: %s\n", one->prop_label[pid].c_str());
  }

  // ask for output filename
  printf("\nPlease input the output file name [vorocell.xyz]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"vorocell.xyz"); ptr = strtok(str, " \n\t\r\f");}
  char *fname = new char[strlen(ptr)+1];
  strcpy(fname, ptr);
  ConfirmOverwrite(fname);
  FILE *fp = fopen(fname, "w");
  
  if (job == 1){ // output selected cluster and the Voronoi vertices
    voro::voronoicell_neighbor cell;
    double lmax = MAX(one->lx, MAX(one->ly, one->lz));
   
    int nclus = 0;
    std::vector<double> vertex, dx, dy, dz;
   
    // now to do the real job
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
   
      one->ComputeVoro(mins);
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
            double scale = type2radius[one->attyp[id]]/(type2radius[one->attyp[id]] + type2radius[one->attyp[jd]]);
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

  } else if (job == 4) { // selected atoms and their neighboring atoms, and neighbors of same property (e.g., molecule id)
    std::list<int> outList, center;
    std::list<int>::iterator it;
    center.clear(); outList.clear();

    // method to determin the neighbor list
    choose_neighbor_method(1);

    int updateMethod = 1;
    printf("\nHow would you like to update the selection and neighbors:\n");
    printf("  1. Use selection and neighbors identified from Frame-%d only;\n", istr+1);
    printf("  2. Use selection identified from Frame-%d but update neighbors for every frame;\n", istr+1);
    printf("  3. Update selection and neighbors for every frame.\n");
    printf("Your choice [%d]: ", updateMethod);
    input->read_stdin(str);
    if (count_words(str) > 0){
      char *ptr = strtok(str," \n\t\r\f");
      if (ptr) updateMethod = atoi(ptr);
      updateMethod = MAX(0, MIN(3, updateMethod));
    }
    printf("Your selection : %d\n", updateMethod);

    Timer * timer = new Timer();
    int nused = 0;
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      // update selection
      if (updateMethod > 2 || (updateMethod <= 2 && img == istr)){
        one->selection(selcmd);

        center.clear();
        for (int id = 1; id <= one->natom; ++id){
          if (one->atsel[id]) center.push_back(id);
        }
      }
      int nCenter = center.size();

      // update neighbors
      if (nCenter > 0 && (updateMethod > 1 || (updateMethod == 1 && img == istr))){
        std::set<float> propList;
        outList.clear(); propList.clear();
 
        // get neighbor list
        if (neighbor_method == 1) one->ComputeVoro(mins);
        else one->ComputeNeiList(r2cuts);

        // neighbors of center atoms
        for (it = center.begin(); it != center.end(); ++it){
          int id = *it;
          outList.push_back(id);
          for (int jj = 1; jj <= one->neilist[0][id]; ++jj){
            int jd = one->neilist[jj][id];
            outList.push_back(jd);
            if (one->prop_label.size() > 0) propList.insert(one->atprop[jd][pid]);
          }
        }
        // other atoms that share same property (e.g., mol id) as the neighbors
        if (one->prop_label.size() > 0){
          for (int id = 1; id <= one->natom; ++id){
            std::set<float>::iterator jt;
            for (jt = propList.begin(); jt != propList.end(); ++jt){
              if (fabs(one->atprop[id][pid] - *jt) <= ZERO) outList.push_back(id);
            }
          }
          propList.clear();
        }

        // output 
        outList.sort(); outList.unique();
      }
      
      int nlist = outList.size();
      if (nlist < 1) continue;
     
      one->dir2car();
      // output as xyz format
      fprintf(fp,"%d\n# Frame %d, neighbors and atoms sharing same %s for atoms selected by %s", nlist, one->iframe, one->prop_label[pid].c_str(), selcmd);
      int idx = 0;
      for (it = outList.begin(); it != outList.end(); ++it){
        int id = *it;
        if (type2atnum == NULL){ // no elements assigned, print atomic type num as element
          fprintf(fp,"%d %lg %lg %lg %f", one->attyp[id], one->atpos[id][0], one->atpos[id][1], one->atpos[id][2], one->atprop[id][pid]);
          if (flabel) fprintf(fp," %d", id);
          if (idx < 3) fprintf(fp," crystal_vector %d %lg %lg %lg\n", idx+1, one->axis[idx][0], one->axis[idx][1], one->axis[idx][2]);
          else fprintf(fp,"\n");

        } else { // in case elements are assigned, print true element names
          char ename[3];
          element->Num2Name(type2atnum[one->attyp[id]], ename);
          fprintf(fp,"%2s %lg %lg %lg %f", ename, one->atpos[id][0], one->atpos[id][1], one->atpos[id][2], one->atprop[id][pid]);
          if (flabel) fprintf(fp," %d", id);
          if (idx < 3){
            fprintf(fp," crystal_vector %d %lg %lg %lg\n", idx+1, one->axis[idx][0], one->axis[idx][1], one->axis[idx][2]);
          } else {
            fprintf(fp,"\n");
          }
        }
        ++idx;
      }
      ++nused;
    }
    outList.clear(); center.clear();

    timer->stop();
    printf("\n%d frames processed, with total CPU time used: %g seconds.\n", nused, timer->cpu_time());
    delete timer;

  } else { // selected atoms or selected clusters
    std::list<int> outlist;
    std::set<int> center;

    // now to do the real job
    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      one->ComputeVoro(mins);
      one->selection(selcmd);

      // get the selected list
      outlist.clear(); center.clear();
      for (int id = 1; id <= one->natom; ++id){
        if (one->atsel[id]){
          outlist.push_back(id); center.insert(id);
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
        int ic = 0;
        if (center.count(id) > 0) ic = 1;

        if (type2atnum == NULL){ // no elements assigned, print atomic type num as element
          fprintf(fp,"%d %lg %lg %lg", one->attyp[id], one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);
          if (flabel) fprintf(fp," %d-%d", ic, id);
          if (idx < 3) fprintf(fp," crystal_vector %d %lg %lg %lg\n", idx+1, one->axis[idx][0], one->axis[idx][1], one->axis[idx][2]);
          else fprintf(fp,"\n");

        } else { // in case elements are assigned, print true element names
          char ename[3];
          element->Num2Name(type2atnum[one->attyp[id]], ename);
          fprintf(fp,"%2s %lg %lg %lg", ename, one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);
          if (flabel) fprintf(fp," %d-%d", ic, id);
          if (idx < 3){
            fprintf(fp," crystal_vector %d %lg %lg %lg\n", idx+1, one->axis[idx][0], one->axis[idx][1], one->axis[idx][2]);
          } else {
            fprintf(fp,"\n");
          }
        }
        ++idx;
      }
    }

    outlist.clear(); center.clear();
  }

  printf("\nJob done, the results were written to file: %s\n", fname);
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

  fclose(fp);
  delete []fname;
return;
}

/*------------------------------------------------------------------------------*/
