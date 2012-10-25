#include "driver.h"
#include "voro++.hh"
#include "math.h"
#include <list>
#include "random.h"
#include "time.h"

#define MAXLINE 1024
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to check the connectivity of certain clusters
 *----------------------------------------------------------------------------*/
void Driver::ClusterConnectivity()
{
  int job = 1;
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<8; i++) printf("===="); printf("  Connectivity  ");
  for (int i=0; i<8; i++) printf("===="); printf("\n");
  printf("  1. Voronoi indices for neighbors of certain clusters;\n");
  printf("  2. Statistics on connectivities of certain clusters;\n");
  printf("  0. return.\nYour choice [%d]: ", job);
  fgets(str,MAXLINE,stdin);

  char *ptr = strtok(str," \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 2){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  } else for (int i=0; i<20; i++) printf("----"); printf("\n");

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
    one = all[0];
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
        for (int i=0; i<20; i++) printf("===="); printf("\n");
        return;
      }
    }
  }

  // define output filename
  FILE *fp; char *fname;
  fp = NULL; fname = NULL;

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
    fprintf(fp,"#Voronoi indices of the %d%s shell neighbors of atoms selected by: %s", nshell, Nos[nshell], selcmd);
    fprintf(fp,"# centered on clusters: ");
    for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); it++) fprintf(fp," %s", it->c_str());
    fprintf(fp,"\n# frame atom type voronoi-index\n");

  } else if (job == 2){
    printf("\nIf you want to output the cluster info, input the file name now: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ptr = strtok(str," \n\t\r\f");
      fname = new char [strlen(ptr)+1];
      strcpy(fname, ptr);

      // open and write header
      fp = fopen(fname,"w");
      fprintf(fp,"# Cluster connectivity info for atoms selected by: %s", selcmd);
      fprintf(fp,"# centered on clsuters:" );
      for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); it++) fprintf(fp," %s", it->c_str());
      fprintf(fp,"\n# frame atom SC-id voronoi-index #connected connected-clusters\n");
    }
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
  printf("Edges whose length takes less ratio than %lg will be skiped!\n\n", voro_mins[1]);

  // variables for connectivity info
  int nclus, nsuper, niso, counts[5];
  nclus = nsuper = niso = counts[0] = counts[1] = counts[2] = counts[3] = counts[4] = 0;

  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];
    one->selection(selcmd);

    // work space for Voronoi
    int nmax = 24, **neilist, *cenlist;
    int natom = one->natom;
    cenlist = memory->create(cenlist, natom+1, "cenlist");
    neilist = memory->create(neilist, nmax+1, natom+1, "neilist");
    for (int ii=0; ii<=natom; ii++) cenlist[ii] = 0;

    map<int,std::string> voroindex; voroindex.clear();

    // compute the neighbor list and voro info
    FEFF_voro(voroset, voro_mins, nmax, neilist, cenlist, voroindex);

    // check the # of selected clusters
    int nc = 0;
    for (int ii=1; ii<=natom; ii++) nc += cenlist[ii];
    if (nc < 1){
      voroindex.clear();
      memory->destroy(neilist);
      memory->destroy(cenlist);

      continue;
    }

    // analyse the result
    if (job == 1){ // Voronoi index for neighbors of selected clusters

      for (int id=1; id <= natom; id++){
        if (cenlist[id] == 0) continue;

        // find atoms in the the desired shell
        list<int> cluster;
        map<int,int> shell;
        cluster.clear(); shell.clear();
        cluster.push_back(id); shell[id] = 0;
        FEFF_cluster(0, nshell, neilist, id, cluster, shell);
          
        cluster.sort(); cluster.unique();

        for (list<int>::iterator it = cluster.begin(); it != cluster.end(); it++){
          int jd = *it;
          if (shell[jd] == nshell) fprintf(fp,"%d %d %d %s\n", img+1, jd, one->attyp[jd], voroindex[jd].c_str());
        }
      }

    } else if (job == 2){ // Connectivity info

      map<int,int> nconn, SCid;
      map<bigint,int> conns, connt;
      nconn.clear(); SCid.clear(); conns.clear(); connt.clear();

      for (int id=1; id <= natom; id++){
        if (cenlist[id] == 0) continue;

        nclus++;
        if (nconn.count(id) < 1){
          nconn[id] = 0;
          SCid[id] = ++nsuper;
        }
          
        for (int jd=id+1; jd <= natom; jd++){
          if (cenlist[jd] == 0) continue;

          int ncomm = 0, paired = 0;
          int ni = neilist[0][id];
          int nj = neilist[0][jd];
          for (int ii=1; ii<= ni; ii++){
            int iid = neilist[ii][id];
            if (iid == jd) paired = 1;
            for (int jj=1; jj<= nj; jj++){
              int jjd = neilist[jj][jd];
              if (iid == jjd) ncomm++;
            }
          }
          ncomm = MIN(3,ncomm);
          if (paired) ncomm = 4;
          if (ncomm){
            SCid[jd] = SCid[id];
            if (nconn[jd] < 1) nconn[jd] = 0;
            nconn[id]++; nconn[jd]++;

            conns[(id-1)*natom+nconn[id]] = jd;
            connt[(id-1)*natom+jd] = ncomm;

            conns[(jd-1)*natom+nconn[jd]] = id;
            connt[(jd-1)*natom+id] = ncomm;

            counts[0]++; counts[ncomm]++;
          }
        }
        if (nconn[id] == 0) niso++;
      }

      // output connectivity info if asked
      if (fp){
        for (int id=1; id <= natom; id++){
          if (cenlist[id] == 0) continue;

          fprintf(fp,"%d %d %d %s %d", img+1, id, SCid[id], voroindex[id].c_str(), nconn[id]);
          for (int ii=1; ii<= nconn[id]; ii++){
            int jd = conns[(id-1)*natom+ii];
            int ct = connt[(id-1)*natom+jd];
            fprintf(fp," %d-(%s):%d", jd, voroindex[jd].c_str(), ct);
          }
          fprintf(fp,"\n");
        }
      }
      nconn.clear(); SCid.clear(); conns.clear(); connt.clear();
    }

    voroindex.clear();
    memory->destroy(neilist);
    memory->destroy(cenlist);
  }

  if (job == 2){
    printf("\n"); for (int i=0; i<20; i++) printf("----");
    printf("\nConnectivity info for clusters:");
    for (set<std::string>::iterator it = voroset.begin(); it != voroset.end(); it++) printf(" %s", it->c_str());
    printf("\n"); for (int i=0; i<20; i++) printf("----");
    printf("\nNClus NSuper Niso NLinked Vertex%% Edge%% Face%% Penetrate%%\n");
    printf("%d %d  %d %d", nclus, nsuper, niso, nsuper-niso); counts[0] = MAX(1,counts[0]);
    for (int ii=1; ii<5; ii++) printf(" %g", double(counts[ii])/double(counts[0])*100.); printf("\n");
  }

  if (fp) fclose(fp);
  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}