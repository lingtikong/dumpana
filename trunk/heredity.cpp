#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to compute the heredity of polyhedron
 *----------------------------------------------------------------------------*/
void Driver::heredity()
{
  char str[MAXLINE];

  // menu
  int job = 1;
  printf("\n"); for (int i = 0; i < 7; ++i) printf("====");
  printf("  Heredity of clusters  ");
  for (int i = 0; i < 7; ++i) printf("====");
  printf("\nPlease select your desired job:\n");
  for (int i = 0; i < 20; ++i) printf("----"); printf("\n");
  printf("  1. Heredity of selected cluster between consecutive frames;\n");
  //printf("  2. Sizes of persisting clusters between consecutive frames;\n");
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

  // refinement info
  set_cutoffs(1);
  one = all[istr];
  
  // selection commands
  char selcmd[MAXLINE];
  one = all[istr];
  // selection commands for atoms
  while (1){
    printf("\nPlease input the selection command for central atoms, `h` for help [all]: ");
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

  // output file name
  printf("\nPlease input the file name to output the results [heredity.dat]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str,"heredity.dat");
  ptr = strtok(str," \n\t\r\f");
  FILE *fp = fopen(ptr, "w");

  int nused = 0;
  DumpAtom *now;  // pointer to current frame

  // timer
  Timer * timer = new Timer();

  // prepare for the real calculations
  one = all[istr];
  // compute the neighbor list and voro info
  one->ComputeVoro(mins);

  // now to compute the heredity info
  if (job == 1){  // Heredity of selected Voronoi clusters
    // This is defined as the number of selected atoms in the current frame
    // 1) with the same Voronoi index and same neighbors as in previous frame; perfect heredity
    // 2) with the same Voronoi index but different neighbors as in previous frame; partial heredity
    // 3) different Voronoi index wrt. the previous frame; new-born

    // set clusters to analyse
    set<std::string> voroset; std::string vindex;
    while ( 1 ) {
      printf("\nPlease input the Voronoi index of the desired clusters, e.g., 0,6,0,8.\n");
      printf("If multiple indices are wanted, separate them by space: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        printf("\nClusters centered on Atoms with Voronoi indices: %s will be analysed.\n", str);
   
        char *ptr = strtok(str," \n\t\r\f");
        while (ptr){
          vindex.assign(ptr);
          voroset.insert(vindex);
   
          ptr = strtok(NULL," \n\t\r\f");
        }
        break;

      } else {
        printf("Atoms with any Voronoi index will be examined, procceed? (y/n)[y]: ");
        if (count_words(fgets(str,MAXLINE,stdin)) > 0){
          char *ptr = strtok(str," \n\t\r\f");
          if (strcmp(ptr, "y") == 0 || strcmp(ptr, "Y") == 0) break;
        } else break;
      }
    }

    // write header of output file
    fprintf(fp, "################################################################################\n");
    fprintf(fp, "# Definitions\n#  1) perfect heredity: same Voronoi index and same neighbors as in previous frame;\n");
    fprintf(fp, "#  2) Partial heredity: same Voronoi index but different neighbors as in previous frame;\n");
    fprintf(fp, "#  3) New born: different Voronoi index to that in previous frame.\n");
    fprintf(fp, "################################################################################\n");
    fprintf(fp, "# timestep  n-total  n-perfect  n-partial n-newborn\n");

    for (int img = istr+inc; img <= iend; img += inc){ // loop over frames
      now = all[img];
      // compute the neighbor list and voro info
      now->ComputeVoro(mins);

      if (now->natom != one->natom){
        printf("\nERROR: number of atoms from current (%d) frame differ from previous one!\n", img+1);

        one = now;
        continue;
      }

      now->selection(selcmd);
      if (now->nsel < 1){
        one = now;
        continue;
      }

      int nperf, npart, nnewb;
      nperf = npart = nnewb = 0;
      // loop over all atoms
#pragma omp parallel for default(shared)
      for (int id = 1; id <= now->natom; ++id){
        if (now->atsel[id] == 0) continue;
        if (voroset.size() > 0 && voroset.count(now->voro[id]) == 0) continue;

        if (now->voro[id].compare(one->voro[id]) == 0) {
          int ni = now->neilist[0][id];
          int hit = 0;
          for (int jj = 1; jj <= ni; ++jj){
            int jd = now->neilist[jj][id];
            for (int kk = 1; kk <= ni; ++kk) if (jd == one->neilist[kk][id]) {++hit; break;}
          }
          if (hit == ni){
#pragma omp atomic
            ++nperf;
          } else {
#pragma omp atomic
            ++npart;
          }

        } else {
#pragma omp atomic
          ++nnewb;
        }
         
      }

      if (min_mem) one->FreeVoro();
      one = now;
      ++nused;

      fprintf(fp, "%d %d %d %d %d\n", one->tstep, nperf + npart + nnewb, nperf, npart, nnewb);
    } // end loop over frames

  } else {  // Sizes of persisting clusters between consecutive frames

    int *cid; cid = NULL;
    std::vector<std::list<int> *> ids_all;
    std::list<int> *ids_one;

    for (int img = istr+inc; img <= iend; img += inc){
      now = all[img];
      // compute the neighbor list and voro info
      now->ComputeVoro(mins);

      if (now->natom != one->natom){
        printf("\nERROR: number of atoms from current (%d) frame differ from previous one!\n", img+1);

        one = now;
        continue;
      }

      // select atoms as source
      now->selection(selcmd);
      if (now->nsel < 1){
        one = now;
        continue;
      }

      memory->grow(cid, now->natom+1, "cid");
      for (int id = 1; id <= now->natom; ++id) cid[id] = 0;

      while (! ids_all.empty() ){
        ids_one = ids_all.back();
        ids_one->clear();
        ids_all.pop_back();
      }
      ids_all.clear();

      // loop over all atoms
      for (int id = 1; id <= now->natom; ++id){
        if (now->atsel[id] == 0) continue;

        if (cid[id] == 0){
          ids_one = new std::list<int>;
          ids_one->clear();
          ids_all.push_back(ids_one);

          cid[id] = ids_all.size();
          ids_one->push_back(id);

          std::set<int> nextshell;
          nextshell.clear();
          nextshell.insert(id);
          check_neigh(nextshell, one, now, ids_one, cid);
        }
      }
      fprintf(fp, "# current timestep: %d, previous timestep: %d, number of inheretied clusters: %d\n", now->tstep, one->tstep, ids_all.size());
      if (flag_out&WHereId) fprintf(fp, "# cluster-id nat-in-cluster  atom-IDs\n");
      else fprintf(fp, "# cluster-id nat-in-cluster\n");

      for (int ii = 0; ii < ids_all.size(); ++ii){
        ids_one = ids_all[ii];
        int nat = ids_one->size();
        fprintf(fp, "%d %d", ii+1, nat);
        if (flag_out&WHereId){
          while (! ids_one->empty() ){
            fprintf(fp, " %d", ids_one->front());
            ids_one->pop_front();
          }
        } else ids_one->clear();
        fprintf(fp, "\n");
      }
      ids_all.clear();

      one = now;
      ++nused;
    }
    memory->destroy(cid);
  }

  // print time info
  fclose(fp);
  timer->stop();
  printf("\nTotal CPU time used: %g seconds.\n", timer->cpu_time());
  delete timer;

  printf("\n%d images were used in the evaluation of heredity.\n",  nused);
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------
 * Private recursive method to check if neighbors are inheretied from previous
 * frame.
 *------------------------------------------------------------------------------ */
void Driver::check_neigh(std::set<int> next, DumpAtom *one, DumpAtom *now, std::list<int> * ids, int *cid)
{
  std::set<int> nshell;
  nshell.clear();

  for (std::set<int>::iterator it = next.begin(); it != next.end(); ++it){
    int id = *it;
    int ni = now->neilist[0][id];
    int nj = one->neilist[0][id];
    for (int jj = 1; jj <= ni; ++jj){
      int jd = now->neilist[jj][id];
      if (cid[jd] == cid[id] || cid[jd] == -1) continue;

      int hit = 0;
      for (int kk = 1; kk <= nj; ++kk) if (one->neilist[kk][id] == jd) {++hit; break;}
      if (hit){
        cid[jd] = cid[id];
        ids->push_back(jd);
        nshell.insert(jd);
      } else cid[jd] = -1;
    }
  }
  next.clear();

  if (! nshell.empty()) check_neigh(nshell, one, now, ids, cid);
return;
}
/*------------------------------------------------------------------------------*/
