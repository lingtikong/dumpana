#include "driver.h"

/*------------------------------------------------------------------------------
 * Method to check the connectivity of certain clusters
 *----------------------------------------------------------------------------*/
void Driver::ClusterConnectivity()
{
  int job = 1;
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 8; ++i) printf("===="); printf("  Connectivity  ");
  for (int i = 0; i < 8; ++i) printf("===="); printf("\n");
  printf("  1. Voronoi indices for neighbors of certain clusters;\n");
  printf("  2. Statistics on connectivities of certain clusters;\n");
  printf("  0. return.\nYour choice [%d]: ", job);
  fgets(str,MAXLINE,stdin);

  char *ptr = strtok(str," \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 2){
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  } else for (int i = 0; i < 20; ++i) printf("----"); printf("\n");

  // voro refinement info
  set_cutoffs(1);
  one = all[istr];

  // selection of atoms for each frame
  char selcmd[MAXLINE];
  printf("\nYou can constrain the central atom of the desired clusters by a selection command.\n");
  while (1){
    printf("Please input your selection command, `h` for help [all]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(selcmd, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(selcmd,"all\n");

    // check the selection command on the first frame
    one = all[istr]; one->ComputeVoro(mins);
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

  // set clusters to analyse
  set<std::string> voroset; std::string vindex;
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
  } else {
    printf("\nAre you sure you want to analyse clusters centered on each atom? (y/n)[y]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"n")==0 || strcmp(ptr,"N")==0){
        for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
        return;
      }
    }
  }

  // define output filename
  FILE *fp, *fpx; char *fname;
  fp = fpx = NULL; fname = NULL;

  const int MaxShell = 9;
  int nshell = 1;

  if (job == 1){
    // cluster size related info
    printf("\nWhich shell of the certral atom`s neighbors would you like to analyse? (1-%d)[1]: ",MaxShell);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr = strtok(str," \n\t\r\f");
      nshell = atoi(ptr);
    }
    nshell = MIN(MaxShell,MAX(1,nshell));
    char Nos[][3] = {"","st","nd","rd","th","th","th","th","th","th","th","th","th","th","th"};
    printf("Neighbor atoms of the %d%s shell will be analysed.\n", nshell, Nos[nshell]);

    // output filename
    printf("\nPlease input the output file name [NeighVoro-%d%s.dat]: ", nshell, Nos[nshell]);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ptr = strtok(str," \n\t\r\f");
      fname = new char [strlen(ptr)+1];
      strcpy(fname, ptr);
    } else {
      fname = new char [20];
      sprintf(str,"NeighVoro-%d%s.dat", nshell, Nos[nshell]);
      strcpy(fname, str);
    }
    // open and write header
    fp = fopen(fname,"w");
    fprintf(fp,"#Voronoi refinement info: surf_min = %g, edge_min = %g, nei_min = %d\n", mins[0], mins[2], int(mins[1]));
    fprintf(fp,"#Voronoi indices of the %d%s shell neighbors of atoms selected by: %s", nshell, Nos[nshell], selcmd);
    fprintf(fp,"# centered on clusters: ");
    for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); ++it) fprintf(fp," %s", it->c_str());
    fprintf(fp,"\n# frame atom type voronoi-index\n");

  } else if (job == 2){
    // output filename
    printf("\nPlease input the output file name [Connectivity.dat]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ptr = strtok(str," \n\t\r\f");
      fname = new char [strlen(ptr)+1];
      strcpy(fname, ptr);
    } else {
      fname = new char [18];
      sprintf(str,"Connectivity.dat");
      strcpy(fname, str);
    }
    fp = fopen(fname,"w");

    // detailed info is optional
    printf("If you want to output the cluster info, input the file name now: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ptr = strtok(str," \n\t\r\f");

      // open and write header
      fpx = fopen(ptr,"w");
      fprintf(fp,"#Voronoi refinement info: surf_min = %g, edge_min = %g, nei_min = %d\n", mins[0], mins[2], int(mins[1]));
      fprintf(fpx,"# Cluster connectivity info for atoms selected by: %s", selcmd);
      fprintf(fpx,"# centered on clsuters:" );
      for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); ++it) fprintf(fpx," %s", it->c_str());
      fprintf(fpx,"\n# frame atom SC-id voronoi-index #connected connected-clusters\n");
    }
  }

  // variables for connectivity info
  int nclus, nsuper, niso, counts[5];
  nclus = nsuper = niso = counts[0] = counts[1] = counts[2] = counts[3] = counts[4] = 0;

  // work space for Voronoi
  int *cenlist;
  int natom = one->natom;
  memory->create(cenlist, natom+1, "cenlist");

  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];

    // compute the neighbor list and voro info
    one->ComputeVoro(mins);

    one->selection(selcmd);

    // work space for Voronoi
    if (natom > one->natom) memory->grow(cenlist, natom+1, "cenlist");

    natom = one->natom;
    for (int ii = 0; ii <= natom; ++ii) cenlist[ii] = 0;

    // look for the list of desired centeral atoms
    for (int ii = 1; ii <= natom; ++ii){
      if (one->atsel[ii] == 0) continue;
      if (voroset.size() == 0 || voroset.count(one->voro[ii]) > 0) cenlist[ii] = 1;
    }

    // check the # of selected clusters
    int nc = 0;
    for (int ii = 1; ii <= natom; ++ii) nc += cenlist[ii];
    if (nc < 1) continue;

    // analyse the result
    if (job == 1){ // Voronoi index for neighbors of selected clusters

      for (int id = 1; id <= natom; ++id){
        if (cenlist[id] == 0) continue;

        // find atoms in the the desired shell
        list<int> cluster;
        map<int,int> shell;
        cluster.clear(); shell.clear();
        cluster.push_back(id); shell[id] = 0;
        one->voro_cluster(0, nshell, id, cluster, shell);
          
        cluster.sort(); cluster.unique();

        for (list<int>::iterator it = cluster.begin(); it != cluster.end(); ++it){
          int jd = *it;
          if (shell[jd] == nshell) fprintf(fp,"%d %d %d %s\n", img+1, jd, one->attyp[jd], one->voro[jd].c_str());
        }
      }

    } else if (job == 2){ // Connectivity info

      map<int,int> nconn, SCid;
      map<bigint,int> conns, connt;
      nconn.clear(); SCid.clear(); conns.clear(); connt.clear();

      for (int id = 1; id <= natom; ++id){
        if (cenlist[id] == 0) continue;

        ++nclus;
        if (nconn.count(id) < 1){
          nconn[id] = 0;
        }
          
        for (int jd = id+1; jd <= natom; ++jd){
          if (cenlist[jd] == 0) continue;

          // To determine the num of common neighbors
          int ncomm = 0, paired = 0;
          int ni = one->neilist[0][id];
          int nj = one->neilist[0][jd];
          for (int ii = 1; ii <= ni; ++ii){
            int iid = one->neilist[ii][id];
            if (iid == jd) paired = 1;
            for (int jj = 1; jj <= nj; ++jj){
              int jjd = one->neilist[jj][jd];
              if (iid == jjd) ++ncomm;
            }
          }
          // To define the connection type:
          // 0, not connected; 1, sharing point; 2, sharing edge; 3, sharing face; 4, bonded
          ncomm = MIN(3,ncomm);
          if (paired) ncomm = 4;
          if (ncomm){
            if (nconn.count(jd) < 1) nconn[jd] = 0;
            ++nconn[id]; ++nconn[jd];

            conns[(id-1)*natom+nconn[id]] = jd;
            connt[(id-1)*natom+jd] = ncomm;

            conns[(jd-1)*natom+nconn[jd]] = id;
            connt[(jd-1)*natom+id] = ncomm;

            ++counts[0]; ++counts[ncomm];
          }
        }
        if (nconn[id] == 0) ++niso;
      }

      nsuper = 0;
      for (int id = 1; id <= natom; ++id) SCid[id] = 0;

      for (int id = 1; id <= natom; ++id){
        if (SCid[id] == 0){
           SCid[id] = ++nsuper;
           set<int> checked; checked.clear();
           IterateOverConn(id, nsuper, natom, SCid, nconn, conns, checked);
        }
      }

      // output connectivity info if asked
      if (fpx){
        for (int id = 1; id <= natom; ++id){
          if (cenlist[id] == 0) continue;

          fprintf(fpx,"%d %d %d %s %d", img+1, id, SCid[id], one->voro[id].c_str(), nconn[id]);
          for (int ii = 1; ii <= nconn[id]; ++ii){
            int jd = conns[(id-1)*natom+ii];
            int ct = connt[(id-1)*natom+jd];
            fprintf(fpx," %d-(%s):%d", jd, one->voro[jd].c_str(), ct);
          }
          fprintf(fpx,"\n");
        }
      }
      nconn.clear(); SCid.clear(); conns.clear(); connt.clear();
    }

    if (min_mem) one->FreeVoro();
  }
  memory->destroy(cenlist);

  if (job == 2){
    printf("\n"); for (int i = 0; i < 20; ++i) printf("----");
    printf("\nConnectivity info for clusters:");
    for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); ++it) printf(" %s", it->c_str());
    printf("\n"); for (int i = 0; i < 20; ++i) printf("----");
    printf("\nNClus NSuper Niso NLinked Vertex%% Edge%% Face%% Penetrate%%\n");
    printf("%d %d  %d %d", nclus, nsuper, niso, nsuper-niso); counts[0] = MAX(1,counts[0]);
    for (int ii = 1; ii < 5; ++ii) printf(" %g", double(counts[ii])/double(counts[0])*100.);
    printf("\n"); for (int i = 0; i < 20; ++i) printf("----"); printf("\n");
    
    fprintf(fp,"# Connectivity info for clusters:");
    for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); ++it) fprintf(fp," %s", it->c_str());
    fprintf(fp,"\n#"); for (int i = 0; i < 20; ++i) fprintf(fp,"----");
    fprintf(fp,"\n# Definitions\n#  NClus      : total number of clusters;\n");
    fprintf(fp,"#  NSuper     : total number of super clusters;\n");
    fprintf(fp,"#  Niso       : total number of isolated clusters, ie, not linked;\n");
    fprintf(fp,"#  NLinked    : total number of linked clusters, Niso + NLinked = NSuper;\n");
    fprintf(fp,"#  Vertex%%    : for linked clusters, the percentage of shared via vertex;\n");
    fprintf(fp,"#  Edge%%      : for linked clusters, the percentage of shared via edge;\n");
    fprintf(fp,"#  Face%%      : for linked clusters, the percentage of shared via face;\n");
    fprintf(fp,"#  Penetrate%% : for linked clusters, the percentage of shared via volume,\n");
    fprintf(fp,"#                i.e., nearby cluster centers are bonded to each other;\n#");
    for (int i = 0; i < 20; ++i) fprintf(fp,"----");
    fprintf(fp,"\n#NClus NSuper Niso NLinked Vertex%% Edge%%   Face%%   Penetrate%%\n");
    fprintf(fp,"%d %d  %d %d", nclus, nsuper, niso, nsuper-niso); counts[0] = MAX(1,counts[0]);
    for (int ii = 1; ii < 5; ++ii) fprintf(fp," %g", double(counts[ii])/double(counts[0])*100.); fprintf(fp,"\n");
  }

  fclose(fp);
  if (fpx) fclose(fpx);
  printf("\nJob done, the results are written to file: %s.\n", fname); delete []fname;
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------ */
void Driver::IterateOverConn(int id, int scid, int natom, map<int,int> Sid, map<int,int> nc, map<bigint,int> cns, set<int> cked)
{
  cked.insert(id);
  for (int ii = 1; ii <= nc[id]; ++ii){
    int jd = cns[(id-1)*natom + ii];
    Sid[jd] = scid;

    if (cked.count(jd) < 1) IterateOverConn(jd, scid, natom, Sid, nc, cns, cked);
  }

return;
}
/*------------------------------------------------------------------------------ */
