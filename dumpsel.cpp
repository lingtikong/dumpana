#include "driver.h"
#include <list>
#include <set>
#include <map>
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to dump selected frames as a new dump file, with the atom selection
 * information as an extra per atom property.
 *----------------------------------------------------------------------------*/
void Driver::DumpSelection()
{
  char str[MAXLINE];
  int job = 0;
  printf("\n"); for (int i = 0; i < 6; ++i) printf("====");
  printf("   Dump selection as property   ");
  for (int i = 0; i < 6; ++i) printf("====");
  printf("\nPlease select your desired job:\n");
  printf("  1. Select atoms;\n");
  printf("  2. Select atoms and their neighboring atoms/molecules;\n");
  printf("  3. Select atoms and their neighbors within another selection;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  input->read_stdin(str);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 3){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  printf("\n");
  
  // thresholds for surface and edges
  one = all[istr];

  // selection commands
  char srcSel[MAXLINE], desSel[MAXLINE];
  int *insrc; insrc = NULL;

  while (1){
    printf("\nPlease input your selection command, `h` for help [all]: ");
    input->read_stdin(str);
    if (count_words(str) > 0){
      strcpy(srcSel, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(srcSel,"all\n");

    // check the selection command on the first frame
    one->selection(srcSel); one->SelInfo();
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

  if (job == 3){
    while (1){
      printf("\nPlease input your selection command for neighbors, `h` for help [all]: ");
      input->read_stdin(str);
      if (count_words(str) > 0){
        strcpy(desSel, str);
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
      } else strcpy(desSel,"all\n");

      // check the selection command on the first frame
      one->selection(desSel); one->SelInfo();
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
  }

  int pid = 0;
  // whether to output neighbors and other atoms sharing same property
  if (job > 1){
    // method to determin the neighbor list
    choose_neighbor_method(1);

    if (one->prop_label.size() > 0) { // property id of neighbors
      printf("\nThe available properties are: ");
      for (int i = 0; i < one->prop_label.size(); ++i) printf("%d) %s; ", i+1, one->prop_label[i].c_str());
      printf("\nPlease input the id of the property to be shared by neighbors [%d]: ", pid+1);
      input->read_stdin(str);
      ptr = strtok(str, " \n\t\r\f");
      if (ptr) pid = atoi(ptr) - 1;
      if (pid < 0 || pid >= one->prop_label.size()){
        pid = -1;
        printf("The property chosen is: None\n");
      } else printf("The property chosen is: %s\n", one->prop_label[pid].c_str());
    }
  }

  // options to update the selection
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
    updateMethod = MAX(1, MIN(3, updateMethod));
  }
  printf("Your selection : %d\n", updateMethod);

  // options to dump the results
  int dumpStyle = 1;
  printf("\nHow would you like to dump the trajectory information:\n");
  printf("  1. Dump all atoms and append the selection as a property;\n");
  printf("  2. Dump the selected atoms only;\n");
  printf("Your choice [%d]: ", dumpStyle);
  input->read_stdin(str);
  if (count_words(str) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    if (ptr) dumpStyle = atoi(ptr);
    dumpStyle = MAX(1, MIN(2, dumpStyle));
  }
  printf("Your selection : %d\n", dumpStyle);

  // ask for output filename
  printf("\nPlease input the output file name [dump2.lammpstrj]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) {strcpy(str,"dump2.lammpstrj"); ptr = strtok(str, " \n\t\r\f");}
  char *fname = new char[strlen(ptr)+1];
  strcpy(fname, ptr);

  ConfirmOverwrite(fname);
  FILE *fp = fopen(fname, "w");
  
  Timer * timer = new Timer();
  std::set<int> source, outList;
  std::set<int>::iterator it;

  // Loop over frames
  int nused = 0;
  for (int img = istr; img <= iend; img += inc){
    one = all[img];
    // update selection
    if (updateMethod > 2 || (updateMethod <= 2 && img == istr)){
      printf("\nUpdate source, img = %d\n", img);
      one->selection(srcSel);

      source.clear();
      for (int id = 1; id <= one->natom; ++id){
        if (one->atsel[id]) source.insert(id);
      }
    }
    int nSrc = source.size();

    // update neighbors
    if (nSrc > 0 && job > 1 && (updateMethod > 1 || (updateMethod == 1 && img == istr))){
      printf("\nUpdate neighbors, img = %d\n", img);
      // get neighbor list
      if (neighbor_method == 1) one->ComputeVoro(mins);
      else one->ComputeNeiList(r2cuts);

      outList.clear();
      std::set<float> propList;
      propList.clear();
 
      if (job == 3) one->selection(desSel);
      printf("nSrc = %d, job = %d, nsel = %d\n", nSrc, job, one->nsel);

      // neighbors of center atoms
      for (it = source.begin(); it != source.end(); ++it){
        int id = *it;
        outList.insert(id);
        for (int jj = 1; jj <= one->neilist[0][id]; ++jj){
          int jd = one->neilist[jj][id];
          if (job == 2 || (job == 3 && one->atsel[jd])){
            outList.insert(jd);
            if (one->prop_label.size() > 0 && pid >= 0) propList.insert(one->atprop[jd][pid]);
          }
        }
      }
      // other atoms that share same property (e.g., mol id) as the neighbors
      if (propList.size() > 0){
        for (int id = 1; id <= one->natom; ++id){
          std::set<float>::iterator jt;
          for (jt = propList.begin(); jt != propList.end(); ++jt){
            if ((job == 3 && one->atsel[id] == 1) || job == 2){
              if (fabs(one->atprop[id][pid] - *jt) <= ZERO) outList.insert(id);
            }
          }
        }
        propList.clear();
      }
    }
    if (job == 1) outList = source;

    // output 
    int nlist = outList.size();
    if (nlist < 1) continue;
     
    one->car2dir();

    fprintf(fp,"ITEM: TIMESTEP\n%d\n", one->tstep);
    if (dumpStyle == 1) fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n", one->natom);
    else fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n", nlist);

    if (one->triclinic){
      fprintf(fp,"ITEM: BOX BOUNDS pp pp pp xy xz yz\n");
      double xl = one->xlo + MIN(MIN(0., one->xy), MIN(one->xz, one->xy+one->xz));
      double xh = one->xlo + MAX(MAX(0., one->xy), MAX(one->xz, one->xy+one->xz));
      double yl = one->ylo + MIN(0., one->yz);
      double yh = one->ylo + MAX(0., one->yz);
      fprintf(fp,"%lg %lg %lg\n", xl, xh, one->xy);
      fprintf(fp,"%lg %lg %lg\n", yl, yh, one->xz);
      fprintf(fp,"%lg %lg %lg\n", one->zlo, one->zhi, one->yz);
    } else {
      fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
      fprintf(fp,"%lg %lg\n", one->xlo, one->xhi);
      fprintf(fp,"%lg %lg\n", one->ylo, one->yhi);
      fprintf(fp,"%lg %lg\n", one->zlo, one->zhi);
    }
    fprintf(fp,"ITEM: ATOMS id type xs ys zs");
    if (one->image){
      fprintf(fp, " ix iy iz");
    }
    for (int i = 0; i < one->prop_label.size(); ++i) fprintf(fp, " %s", one->prop_label[i].c_str());
    fprintf(fp, " sel\n");
  
    if (dumpStyle == 1){
      for (int id = 1; id <= one->natom; ++id){
        int flag = outList.count(id);
        fprintf(fp,"%d %d %lg %lg %lg", id, one->attyp[id], one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);
        if (one->image){
          fprintf(fp, " %d %d %d", one->image[id][0], one->image[id][1], one->image[id][2]);
        }
        for (int i = 0; i < one->prop_label.size(); ++i) fprintf(fp, " %g", one->atprop[id][i]);
        fprintf(fp, " %d\n", flag);
      }

    } else {
      for (int id = 1; id <= one->natom; ++id){
        int flag = outList.count(id);
        if (flag == 0) continue;
        fprintf(fp,"%d %d %lg %lg %lg", id, one->attyp[id], one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);
        if (one->image){
          fprintf(fp, " %d %d %d", one->image[id][0], one->image[id][1], one->image[id][2]);
        }
        for (int i = 0; i < one->prop_label.size(); ++i) fprintf(fp, " %g", one->atprop[id][i]);
        fprintf(fp, " %d\n", flag);
      }
    }
     ++nused;
  }
  outList.clear(); source.clear(); 
  fclose(fp);

  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  printf("\n%d images were processed and the results are written to %s\n", nused, fname);
  delete []fname;

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------ */
