#include "driver.h"

/*------------------------------------------------------------------------------
 * Method to prepare for the FEFF9 input for XANES or EXAFS calculations.
 *----------------------------------------------------------------------------*/
void Driver::FEFF_main()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<9; i++) printf("===="); printf("  FEFF  ");
  for (int i=0; i<9; i++) printf("===="); printf("\n");
  one = all[istr];

  // map atomic type to elements
  if (type2atnum == NULL){
    printf("Mapping of atomic types to actual elements are needed to generate FEFF input.\n");
    MapType2Elem(1, one->ntype); printf("\n");
  }

  // voro refinement info
  double mins[3];
  mins[0] = 1.e-2; mins[1] = 2.e-3; mins[2] = 0.;

  printf("\nRefined Voronoi tesselation will be needed to procceed.\n");
  printf("Now please input your criterion for tiny surfaces, 0 to keep all [%g]: ", mins[0]);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[0] = atof(ptr);
  printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", mins[0]);

  printf("Sometimes it might be desirable to keep a minimum # of neighbors when refining\n");
  printf("the Voronoi index, for example, keep at least 14 for a bcc lattice, 12 for hcp\n");
  printf("or fcc. If you prefer to do so, input a positive number now [%d]: ", int(mins[2]));
  if (count_words(fgets(str,MAXLINE, stdin)) > 0){
    double dum = atof(strtok(str, " \n\t\r\f"));
    if (dum > 0.) mins[2] = dum;
    printf("\nA minimum number of %d neighobrs will be kept no matter how tiny the surface is.\n", int(mins[2]));
  }

  printf("\nPlease input your criterion for ultra short   [%g]: ", mins[1]);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[1] = atof(ptr);
  printf("Edges whose lengths take less ratio than %lg will be skipped!\n", mins[1]);
  // Show relevant info if Weighted Voronoi is used
  ShowRadius4Voro();
  one->ComputeVoro(mins,weighted);

  int AbsorberType = 0;
  // selection of atoms for each frame
  char selcmd[MAXLINE], workdir[MAXLINE];
  printf("\nPlease define the absorbing atoms, which are usually of a specific type.\n");
  printf("However you can restrict to a fraction of a type by a selection command.\n");
  printf("NOTE: for selection option `voro', if negative MINs are provided, their\n");
  printf("respective default/previous values will be used.\n");
  while (1){
    printf("Please input your selection command, `h` for help [type = 1]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(selcmd, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(selcmd,"type = 1\n");

    // check if `type' is used in the selection command; it is a must unless ntype = 1
    strcpy(str, selcmd);
    char *ptr = strtok(str, " \n\t\r\f");
    while (ptr){
      if (strcmp(ptr,"type") == 0){
        ptr = strtok(NULL," \n\t\r\f");
        ptr = strtok(NULL," \n\t\r\f");
        AbsorberType = atoi(ptr); break;
      }
      ptr = strtok(NULL," \n\t\r\f");
    }
    if (AbsorberType <=0 || AbsorberType > one->ntype){
      if (one->ntype > 1){
        printf("\nNo or wrong atomic type defined in the selection command!\n");
        AbsorberType = 0; continue;
      } else AbsorberType = 1;
    }

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

  // XANES or EXAFS
  int job = 2;
  printf("\nWould you like to do a (1) XANES or (2) EXAFS calculation? (1/2)[%d]: ", job);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    job = 2-atoi(ptr)%2;
  }
  printf("Your selection : %d\n", job);

  // working directory: directory to write all output files
  printf("\nPlease define the working directory, i.e., directory to write all\n");
  printf("output files/directories [.]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    strcpy(workdir, ptr);
  } else strcpy(workdir, ".");
  printf("The working directory will be: %s\n", workdir);

  // some common variables
  char ename[3];
  char dirname[MAXLINE], mkdir[MAXLINE], fname[MAXLINE];
  FILE *fp;

  // file DirList records all directories created
  sprintf(fname,"%s/DirList", workdir);
  FILE *fpx = fopen(fname, "w");
  int ndir = 0;

  // cluster size related info
  const int MaxShell = 12;
  int nshell = 2;
  printf("\nHow many shells would you like to include into the cluster? (1-%d)[%d]: ",MaxShell, nshell);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    nshell = atoi(ptr);
  }
  nshell = MIN(MaxShell,MAX(1,nshell));
  char Nos[][3] = {"","st","nd","rd","th","th","th","th","th","th","th","th","th","th","th"};
  printf("Atoms upto the %d%s shell will be included into the clsuter.\n", nshell, Nos[nshell]);

  // now to do the real job: analyse the selected clusters frame by frame
  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];

    // compute the Voronoi info, so as to get all related info
    one->ComputeVoro(mins, weighted);

    // make selection
    one->selection(selcmd);
    if (one->nsel < 1) continue;

    // loop over all selected absorbing atoms
    for (int id=1; id <= one->natom; id++){
      if (one->atsel[id] == 0) continue;

      // find atoms in the cluster
      list<int> cluster;
      map<int,int> shell;
      cluster.clear(); shell.clear();
      cluster.push_back(id); shell[id] = 0;
      one->voro_cluster(0, nshell, id, cluster, shell);
      
      cluster.sort(); cluster.unique();
      int nclus = cluster.size();

      // check atomic types in the cluster
      map<int,int> attyp2pot;
      attyp2pot.clear(); attyp2pot[0] = 0;
      int npottype = 0;
      for (list<int>::iterator it = cluster.begin(); it != cluster.end(); it++){
        int jd = *it;
        int jp = one->attyp[jd];
        if (jd == id) jp = 0;
        if (attyp2pot.count(jp) == 0) attyp2pot[jp] = ++npottype;
      }

      // open the feff.inp file for current frame
      sprintf(dirname, "%s/F%dA%d", workdir, img+1, id);
      fprintf(fpx, "%s", dirname); ndir++;
      if (flag_out & OutFeff){
        strcpy(mkdir,"mkdir -p "); strcat(mkdir,dirname);
        system(mkdir);
        strcpy(fname, dirname); strcat(fname,"/feff.inp");
        fp = fopen(fname, "w");

        // write the potential part of feff.inp
        element->Num2Name(type2atnum[one->attyp[id]], ename);
        fprintf(fp,"TITLE Cluster centered on atom %d (%s, <%s>) of frame %d (MD steps: %d)\n",
          id, ename, one->voro[id].c_str(), img+1, one->tstep);
        fprintf(fp,"TITLE Total number of atoms in cluster: %d\n\n", nclus);

        fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
        for (map<int,int>::iterator it = attyp2pot.begin(); it != attyp2pot.end(); it++){
          int ip = it->first; if (ip == 0) ip = one->attyp[id];
          int IP = it->second;
          element->Num2Name(type2atnum[ip], ename);
          fprintf(fp,"%4d   %3d  %3s  -1   3\n", IP, type2atnum[ip], ename);
        }

        // write other common info
        FEFF_input(job, fp);
      }

      // storage for coordination number and dist info
      int CN[one->ntype+1], CNtot = 0;
      double nndist[one->ntype+1], nndist2[one->ntype+1];
      for (int ip = 1; ip <= one->ntype; ip++){
        CN[ip] = 0;
        nndist[ip] = nndist2[ip] = 0.;
      }
 
      // atomic positions
      one->car2dir();
 
      if (flag_out & OutFeff) fprintf(fp,"\n* Atomci positions in Angstrom\nATOMS\n* x y z ipot tag ishell dist id\n");
      for (list<int>::iterator it = cluster.begin(); it != cluster.end(); it++){
        int jd = *it;
        int jp = one->attyp[jd];
        if (jd == id) jp = 0;
        double dx[3];
        for (int idim=0; idim<3; idim++){
          dx[idim] = one->atpos[jd][idim] - one->atpos[id][idim];
          while (dx[idim] > 0.5) dx[idim] -= 1.;
          while (dx[idim] <-0.5) dx[idim] += 1.;
        }
        dx[0] = dx[0]*one->lx + dx[1]*one->xy + dx[2]*one->xz;
        dx[1] = dx[1]*one->ly + dx[2]*one->yz;
        dx[2] = dx[2]*one->lz;
 
        double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        double rij = sqrt(r2);
        element->Num2Name(type2atnum[one->attyp[jd]], ename);
        if (flag_out & OutFeff) fprintf(fp,"%15.8f %15.8f %15.8f %d %s %d %g %d\n", dx[0], dx[1], dx[2], attyp2pot[jp], ename, shell[jd], rij, jd);
 
        if (shell[jd] == 1){
          CN[jp]++; CNtot++;
          nndist[jp] += rij;
          nndist2[jp] += r2;
        }
      }
      attyp2pot.clear();
 
      // write coordination number info
      if (flag_out & OutFeff){
        fprintf(fp,"\n* Coordination number and nearest neighbor distance info for atom %d,\n", id);
        fprintf(fp,"* only atoms of the Voronoi neighbors are seen as nearest neighbors.\n* Total: %d\n", CNtot);
      }
      fprintf(fpx," %d :", CNtot);
 
      for (int ip = 1; ip <= one->ntype; ip++){
        double stdv = 0.;
        if (CN[ip] > 0){
          nndist[ip] /= double(CN[ip]);
          stdv = sqrt(nndist2[ip]/double(CN[ip])-nndist[ip]*nndist[ip]);
        }
        element->Num2Name(type2atnum[ip], ename);
        if (flag_out & OutFeff) fprintf(fp,"* %2s  %d  %lg +/- %lg\n", ename, CN[ip], nndist[ip], stdv);
        fprintf(fpx," %2s %d %lg +/- %lg; ", ename, CN[ip], nndist[ip], stdv);
      }
      fprintf(fpx," <%s>\n", one->voro[id].c_str());

      shell.clear(); cluster.clear();
      // clsoe the file
      if (flag_out & OutFeff) fclose(fp);
    } // end of loop over frames
  }
  fclose(fpx);

  // script to run feff
  if (flag_out & OutFeff){
    sprintf(fname,"%s/run_feff", workdir);
    FILE *fp = fopen(fname, "w");
    fprintf(fp,"#!/bin/bash\nfor dir in `ls|grep F`\n");
    fprintf(fp,"do\n  if [ -d %c$dir%c ]; then\n", char(34),char(34));
    fprintf(fp,"    cd $dir; feff\n      echo $dir >> ../jobDone.list\n");
    fprintf(fp,"    cd ..\n  fi\ndone\n\nexit 0\n\n");
    fclose(fp);
    sprintf(fname,"chmod +x %s/run_feff", workdir);
    system(fname);
  }
  printf("\nJob done, %d directorys are created and listed in `%s/DirList`.\n", ndir, workdir);
  for (int i=0; i<20; i++) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------
 * Method to generate the common options for XANES or EXAFS calculations.
 *------------------------------------------------------------------------------
 * job    (in) : 1 for XANES; 2 for EXAFS
 * fp     (in) : FILE to write
 *----------------------------------------------------------------------------*/
void Driver::FEFF_input(int job, FILE *fp)
{
  fprintf(fp, "\n* EDGE label s02\n* EDGE K 1.0\n");
  double rfms1 = pow(90.*one->vol/(12.*one->natom), 1./3.);
  double rfms2 = pow(600.*one->vol/(12.*one->natom), 1./3.);
  double rfms3 = 2.2*pow(50.*one->vol/(12.*one->natom), 1./3.);
  if (job == 1){
    fprintf(fp, "* XANES [xkmax xkstep vixan]\nXANES 5.0\n");
    fprintf(fp, "* SCF rfms1 * rfms1 must be converged for real space calculations\n");
    fprintf(fp, "* SCF %g\n* FMS %g\n\n* RPATH %g\n", rfms1, rfms2, rfms3);
  } else {
    fprintf(fp,"* EXAFS [xkmax]\nEXAFS 20.0\n* RPATH %g\n* NLEG 8\n\n", rfms3);
    fprintf(fp, "* SCF rfms1 * rfms1 must be converged for real space calculations\n");
    fprintf(fp, "* SCF %g\n* FMS %g\n", rfms1, rfms2);
  }
  fprintf(fp, "* COREHOLE type * type = none, RPA, or FSR (default)\n");

  fprintf(fp, "* REAL\n* RMULTIPLIER 1.00\n");
  fprintf(fp, "\n*    emin emax eimag\n* LDOS -20 20 0.2\n");
  fprintf(fp, "*       ipot ixsph ifms ipaths igenfmt iff2x\n");
  fprintf(fp, "CONTROL  1    1     1    1      1       1\n");
  fprintf(fp, "*         ppot pxsph pfms ppaths pgenfmt pff2x\n");
  fprintf(fp, "* PRINT    0    0     0    0      0       0\n");

  fprintf(fp, "\n\n* ***** Other useful or advanced cards *****\n");
  fprintf(fp, "\n* *** pot related ***\n");
  fprintf(fp, "* SO2 s02 * expr. value: 0.8-1.0\n");
  fprintf(fp, "* EXCHANGE ixc vr0 vi0 [ixc0]\n");
  fprintf(fp, "* AFOLP folpx  * folpx = 1.0 - 1.3; default 1.15\n");
  fprintf(fp, "* FOLP ipot folpx\n* RGRID delta * default: 0.05\n");
  fprintf(fp, "* SPIN ispin [x y z] * required for XMCD, SPXAS\n");
  fprintf(fp, "* INTERSTITIAL inters vtot * change the way to get interstitial potential\n");
  fprintf(fp, "* ION ipot ionization * ionizes atoms with unique potential index ipot\n");
  fprintf(fp, "* UNFREEZEF * unfreeze f electrons\n");
  fprintf(fp, "* CONFIG input [configurations] * modify the electron configuration\n");
  fprintf(fp, "* CHSHIFT * calculate initial core state energy level using the final SCF potential\n");
  fprintf(fp, "* CHBROADENING igammach\n* CHWIDTH gamma\n");
  fprintf(fp, "* EGAP * works with MPSE only\n* EPS0 eps0 * works with MPSE only\n");
  fprintf(fp, "* RESTART * restart SCF calculation from existing pot.bin\n");
  fprintf(fp, "* SCREEN parameter value * pass custom parameters to screen\n");
  fprintf(fp, "\n* *** xsph related ***\n");
  fprintf(fp, "* ELLIPTICITY elpty x y z * used with POLARIZATION\n");
  fprintf(fp, "* POLARIZATION  x y z\n* MULTIPOLE le2 [l2lp] * specify which multipole transition to include in the claculations\n");
  fprintf(fp, "* LDEC ld * only for NRIXS\n* LJMAX lj * only for NRIXS\n");
  fprintf(fp, "* MPSE ipl * runs FEFF wiht a many-pole model for the self energy\n");
  fprintf(fp, "* RPHASES  * use real phase shifts rather than complex ones\n");
  fprintf(fp, "* RSIGMA   * neglect the imaginary part of the self-energy\n");
  fprintf(fp, "* TDLDA ixfc * use TDLDA to account for screening of the x-ray field and of the PE-core-hole interaction\n");
  fprintf(fp, "\n* *** fms related ***\n");
  fprintf(fp, "* DEBYE temperature Debye-temperature\n* Used to calculate the Debye-Waller factor.\n");
  fprintf(fp, "* STRFAC eta gmax rmax\n");
  fprintf(fp, "\n* *** path related ***\n");
  fprintf(fp, "* NLEG nleg * limits the # of legs for each scattering path, default 8\n");
  fprintf(fp, "* PCRITERIA pcritk pcrith * limit the # of paths\n");
  fprintf(fp, "* SS index ipot deg rss   * works only with OVERLAP card\n");
  fprintf(fp, "* SYMMETRY isym           * tells FEFF to reduce the list of paths using symmetry\n");
  fprintf(fp, "\n* *** genfmt related ***\n");
  fprintf(fp, "* CRITERIA critcw critpw  * limit the # of paths\n");
  fprintf(fp, "* IORDER iord             * order of approximation used in module genfmt\n");
  fprintf(fp, "* NSTAR                   * to write nstar.dat\n");
  fprintf(fp, "\n* *** ff2x related ***\n");
  fprintf(fp, "* CORRECTIONS vrcorr vicorr\n* ABSOLUTE\n* SIG2 sig2\n");
  fprintf(fp, "\n* *** eels related ***\n");
  fprintf(fp, "* MAGIC emagic\n");
  
  fprintf(fp, "\n* *** other ***\n");
  fprintf(fp, "* DIMS nmax lmax\n\n* EGRID\n* grid_type grid_min grid_max grid_step\n");
  fprintf(fp, "* KMESH nkp usesym symfile\n");

return;
}

/*------------------------------------------------------------------------------*/
