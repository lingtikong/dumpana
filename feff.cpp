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
 * Method to prepare for the FEFF9 input for XANES or EXAFS calculations.
 *----------------------------------------------------------------------------*/
void Driver::FEFF_main()
{
  printf("\n"); for (int i=0; i<9; i++) printf("===="); printf("  FEFF  ");
  for (int i=0; i<9; i++) printf("===="); printf("\n");
  one = all[0];

  // map atomic type to elements
  if (type2atnum == NULL){
    printf("Mapping of atomic types to actual elements are needed to generate FEFF input.\n");
    MapType2Elem(1, one->ntype); printf("\n");
  }

  int AbsorberType = 0;
  // selection of atoms for each frame
  char str[MAXLINE], selcmd[MAXLINE], workdir[MAXLINE];
  printf("Please define the absorbing atoms, which are usually of a specific type.\n");
  printf("However you can restrict to a fraction of a type by a selection command.\n");
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

  // further refining the absorbing atoms by their Voronoi index
  int flag_voro = 0, nmax_voro = 0, seed = 0;
  std::set<std::string> voroset; std::string vindex;
  printf("\nIf you want to futher refine the absorbing atoms by their Voronoi index,\n");
  printf("please input the desired Voronoi index now, e.g., 0,6,0,8. if multiple indices\n");
  printf("are wanted, separate them by space: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    flag_voro = 1;
    printf("\nAtoms selected by:%s with Voronoi indices: %s", selcmd, str);
    printf("will be chosen as absorbing atoms.\n");

    char *ptr = strtok(str," \n\t\r\f");
    while (ptr){
      vindex.assign(ptr);
      voroset.insert(vindex);

      ptr = strtok(NULL," \n\t\r\f");
    }

    if (voroset.size() > 0){
      printf("\nIf you want to limit the # of clusters per frame, input the max # now [0]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0) nmax_voro = atoi(strtok(str," \n\t\r\f"));
      if (nmax_voro > 0){
        printf("Please input the seed of the random generator for refining [%d]: ", seed);
        if (count_words(fgets(str,MAXLINE,stdin)) > 0) seed = atoi(strtok(str," \n\t\r\f"));
        if (seed < 1) seed = time(NULL)%86400+1;
      }
    }
  }

  // job type, CFAVERAGE or Clusters
  int flag_pbc = 0;
  // disable this, since it does not work with FEFF
/*
  printf("\nWould you like to run (1), an auto `average` over all absorbing atoms for\n");
  printf("each frame, or (2), individual calculations for each cluster? (1/2)[2]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    flag_pbc = atoi(ptr)%2;
  }
  printf("Your selection : %d\n", 2-flag_pbc);
*/

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

  // real job
  if (flag_pbc==1 && flag_voro == 0){ // no Voronoi needed

    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      one->selection(selcmd);
  
      // open the feff.inp file for current frame
      sprintf(dirname, "%s/Frame%d", workdir, img+1);
      strcpy(mkdir,"mkdir -p "); strcat(mkdir,dirname);
      system(mkdir); ndir++; fprintf(fpx, "%s\n", dirname); // create direcotry for current frame
      strcpy(fname, dirname); strcat(fname,"/feff.inp");
      fp = fopen(fname, "w");

      // write the potential part of feff.inp
      fprintf(fp,"TITLE Frame %d (MD steps: %d), %d types, %d atoms.\n",
        img+1, one->tstep, one->ntype, one->natom);
      fprintf(fp,"TITLE Averaged over %d type-%d atoms.\n\n", one->nsel, AbsorberType);
      fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
      for (int ip=1; ip <= one->ntype; ip++){
        element->Num2Name(type2atnum[ip], ename);
        fprintf(fp,"%4d   %3d  %3s  -1   3\n", ip, type2atnum[ip], ename);
      }
      element->Num2Name(type2atnum[AbsorberType], ename);
      fprintf(fp,"%4d   %3d  %3s  -1   3\n", one->ntype+1, type2atnum[AbsorberType], ename);

      // reciprocal space calculation
      fprintf(fp,"\n* k-space calculation of crystals\nRECIPROCAL\n");
      fprintf(fp,"\n* Cartesian coordinates, Angstrom units\nCOORDINATES 1\n");
      fprintf(fp,"\n* Core-hole by RPA\nCOREHOLE RPA\n");
      fprintf(fp,"\n* KMesh needed\nKMESH 100 0\n");
      
      // lattice and atomic positions
      fprintf(fp,"\n* the lattice info\nLATTICE P 1.0\n");
      for (int idim=0; idim<3; idim++) fprintf(fp,"%15.8f %15.8f %15.8f\n",
        one->axis[idim][0], one->axis[idim][1], one->axis[idim][2]);

      one->dir2car();
      fprintf(fp,"\n* Atomci positions in Angstrom\nATOMS\n* x y z ipot\n");
      for (int i=1; i <= one->natom; i++){
        int ip = one->attyp[i];
        if (one->atsel[i] == 1 && ip == AbsorberType) ip = one->ntype+1;
        fprintf(fp,"%15.8f %15.8f %15.8f %d\n", one->atpos[i][0], one->atpos[i][1], one->atpos[i][2], ip);
      }
      fprintf(fp,"\n* XAFS is computed by averaging over the following type\n");
      // suitable r for CFAVERAGE
      double rmax = pow(3000.*one->vol/(16.*one->natom), 1./3.);
      fprintf(fp,"CFAVERAGE %d %d %g\n", one->ntype+1, one->nsel, rmax);

      // write other common info
      FEFF_input(job, fp);

      // clsoe the file
      fclose(fp);
    }

  } else { // Voronoi needed

    // voro refinement info
    double voro_mins[3];
    voro_mins[0] = voro_mins[1] = 1.e-4;
    voro_mins[2] = 0.;
  
    printf("\nRefined Voronoi tesselation will be needed to procceed.\n");
    printf("Now please input your criterion for tiny surfaces [%g]: ", voro_mins[0]);
    fgets(str,MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
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
  
    printf("\nPlease input your criterion for ultra short   [%g]: ", voro_mins[1]);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) voro_mins[1] = atof(ptr);
    printf("Edges whose length takes less ratio than %lg will be skiped!\n\n", voro_mins[1]);

    // cluster size related info
    const int MaxShell = 12;
    int nshell = 3;
    if (flag_pbc == 0){
      printf("\nHow many shells would you like to include into the cluster? (1-%d)[3]: ",MaxShell);
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        nshell = atoi(ptr);
      }
      nshell = MIN(MaxShell,MAX(1,nshell));
      char Nos[][3] = {"","st","nd","rd","th","th","th","th","th","th","th","th","th","th","th"};
      printf("Atoms upto the %d%s shell will be included into the clsuter.\n", nshell, Nos[nshell]);
    }

    // random generator for confining the number of clusters when certain voronoi indices are asked
    RanPark * random;
    if (nmax_voro > 0) random = new RanPark(seed);

    for (int img = istr; img <= iend; img += inc){ // loop over frames
      one = all[img];
      // not possible to evaluate voro info for triclinic box
      if (one->triclinic) continue;
      one->selection(selcmd);
      if (one->nsel < 1) continue;

      // work space for Voronoi
      int nmax = 24, **neilist, *cenlist;
      int natom = one->natom;
      cenlist = memory->create(cenlist, natom+1, "cenlist");
      neilist = memory->create(neilist, nmax+1, natom+1, "neilist");
      for (int ii=0; ii<=natom; ii++) cenlist[ii] = 0;

      // compute the neighbor list and voro info
      FEFF_voro(voroset, nmax, neilist, cenlist, voro_mins);

      // analyse the result
      int nc = 0;
      for (int ii=1; ii<=natom; ii++) nc += cenlist[ii];
      if (nc < 1){
        memory->destroy(neilist);
        memory->destroy(cenlist);
        continue;

      } else if (nmax_voro > 0){
        // apply limitation on the total number of clusters when certain voronoi indices are expected
        int ndel = nc - nmax_voro;
        while (ndel > 0){
          int id = MIN(random->uniform()*(natom+1), natom);
          ndel -= cenlist[id];
          cenlist[id] = 0;
        }
      }

      if (flag_pbc){
        // open the feff.inp file for current frame
        sprintf(dirname, "%s/Frame%d", workdir, img+1);
        fprintf(fpx, "%s\n", dirname); ndir++;

        strcpy(mkdir,"mkdir -p "); strcat(mkdir,dirname);
        system(mkdir);
        strcpy(fname, dirname); strcat(fname,"/feff.inp");
        fp = fopen(fname, "w");

        // write the potential part of feff.inp
        fprintf(fp,"TITLE Frame %d (MD steps: %d), %d types, %d atoms.\n",
          img+1, one->tstep, one->ntype, one->natom);
        fprintf(fp,"TITLE Averaged over %d type %d atoms.\n\n", nc, AbsorberType);
        fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
        for (int ip=1; ip <= one->ntype; ip++){
          element->Num2Name(type2atnum[ip], ename);
          fprintf(fp,"%4d   %3d  %3s  -1   3\n", ip, type2atnum[ip], ename);
        }
        element->Num2Name(type2atnum[AbsorberType], ename);
        fprintf(fp,"%4d   %3d  %3s  -1   3\n", one->ntype+1, type2atnum[AbsorberType], ename);
  
        // write some common info
        FEFF_input(job, fp);

        // reciprocal space calculation
        fprintf(fp,"\n* k-space calculation of crystals\nRECIPROCAL\n");
        fprintf(fp,"\n* Cartesian coordinates, Angstrom units\nCOORDINATES 1\n");
      
        // lattice and atomic positions
        fprintf(fp,"\n* the lattice info\nLATTICE P 1.0\n");
        for (int idim=0; idim<3; idim++) fprintf(fp,"%15.8f %15.8f %15.8f\n",
          one->axis[idim][0], one->axis[idim][1], one->axis[idim][2]);
  
        one->dir2car();
        fprintf(fp,"\n* Atomci positions in Angstrom\nATOMS\n* x y z ipot\n");
        for (int i=1; i <= one->natom; i++){
          int ip = one->attyp[i];
          if (cenlist[i] == 1 && ip == AbsorberType) ip = one->ntype+1;
          fprintf(fp,"%15.8f %15.8f %15.8f %d\n", one->atpos[i][0], one->atpos[i][1], one->atpos[i][2], ip);
        }
        fprintf(fp,"\n* XAFS is computed by averaging over the following type\n");
        // suitable r for CFAVERAGE
        double rmax = pow(3000.*one->vol/(16.*one->natom), 1./3.);
        fprintf(fp,"CFAVERAGE %d %d %g\n", one->ntype+1, nc, rmax);
  
        // clsoe the file
        fclose(fp);
  
      } else {

        for (int id=1; id <= one->natom; id++){
          if (cenlist[id] == 0) continue;

          // find atoms in the cluster
          std::list<int> cluster;
          std::map<int,int> shell;
          cluster.clear(); shell.clear(); shell[id] = 0;
          FEFF_cluster(0, nshell, neilist, id, cluster, shell);
          
          cluster.sort(); cluster.unique();
          int nclus = cluster.size();

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
            fprintf(fp,"TITLE Cluster centered on atom %d (%s) of frame %d (MD steps: %d)\n",
              id, ename, img+1, one->tstep);
            fprintf(fp,"TITLE Total number of atoms in cluster: %d\n\n", nclus);
  
            fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
            fprintf(fp,"%4d   %3d  %3s  -1   3\n", 0, type2atnum[one->attyp[id]], ename);
            for (int ip=1; ip <= one->ntype; ip++){
              element->Num2Name(type2atnum[ip], ename);
              fprintf(fp,"%4d   %3d  %3s  -1   3\n", ip, type2atnum[ip], ename);
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
          one->dir2car();
          if (flag_out & OutFeff) fprintf(fp,"\n* Atomci positions in Angstrom\nATOMS\n* x y z ipot tag ishell dist id\n");
          for (std::list<int>::iterator it = cluster.begin(); it != cluster.end(); it++){
            int jd = *it;
            int jp = one->attyp[jd];
            if (jd == id) jp = 0;
            double dx[3], r2; r2 = 0.;
            for (int idim=0; idim<3; idim++){
              dx[idim] = one->atpos[jd][idim] - one->atpos[id][idim];
              while (dx[idim] > one->hbox[idim]) dx[idim] -= one->box[idim];
              while (dx[idim] <-one->hbox[idim]) dx[idim] += one->box[idim];
              r2 += dx[idim]*dx[idim];
            }
            double rij = sqrt(r2);
            element->Num2Name(type2atnum[one->attyp[jd]], ename);
            if (flag_out & OutFeff) fprintf(fp,"%15.8f %15.8f %15.8f %d %s %d %g %d\n", dx[0], dx[1], dx[2], jp, ename, shell[jd], rij, jd);

            if (shell[jd] == 1){
              CN[jp]++; CNtot++;
              nndist[jp] += rij;
              nndist2[jp] += r2;
            }
          }

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
              double stdv = sqrt(nndist2[ip]/double(CN[ip])-nndist[ip]*nndist[ip]);
            }
            element->Num2Name(type2atnum[ip], ename);
            if (flag_out & OutFeff) fprintf(fp,"* %2s  %d  %lg +/- %lg\n", ename, CN[ip], nndist[ip], stdv);
            fprintf(fpx," %2s %d %lg +/- %lg; ", ename, CN[ip], nndist[ip], stdv);
          }
          fprintf(fpx,"\n");

          shell.clear(); cluster.clear();
          // clsoe the file
          if (flag_out & OutFeff) fclose(fp);
        }
      }

      memory->destroy(neilist);
      memory->destroy(cenlist);
    } // end of loop over frames
  }

  fclose(fpx);
  printf("\nJob done, %d directorys are created and listed in `%s/DirList`.\n", ndir, workdir);
  for (int i=0; i<20; i++) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------
 * Method to generate the common options for XANES or EXAFS calculations.
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

/*------------------------------------------------------------------------------
 * Method to find the neighbor list based on refined Voronoi calculation
 *----------------------------------------------------------------------------*/
void Driver::FEFF_voro(std::set<std::string> vlist, int &nmax, int **nlist, int *clist, double *mins)
{
  double surf_min = mins[0];
  double edge_min = mins[1];
  int nminnei = int(mins[2]);
  // set local variables
  double xlo = one->xlo, xhi = one->xhi, lx = one->lx;
  double ylo = one->ylo, yhi = one->yhi, ly = one->ly;
  double zlo = one->zlo, zhi = one->zhi, lz = one->lz;
  double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;
  int natom = one->natom;

  int *attyp = one->attyp;
  double **atpos = one->atpos;

  // need cartesian coordinates
  one->dir2car();

  // compute optimal size for container, and then contrust it
  double l = pow(double(natom)/(5.6*lx*ly*lz), 1./3.);
  int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);
  voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

  // put atoms into the container
  for (int i=1; i<= natom; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);

  voro::voronoicell_neighbor c1, c2, *cell;
  voro::c_loop_all cl(con);
  if (cl.start()) do if (con.compute_cell(c1,cl)){
    int id;
    double x, y, z, vol;
    std::vector<int> ff;       // face_freq
    std::vector<int> neigh;    // neigh list
    std::vector<double> fs;    // face areas
    int index[7];
    for (int i=0; i<7; i++) index[i] = 0;
       
    cl.pos(x,y,z);
    id = cl.pid();
    c1.neighbors(neigh);
    c1.face_areas(fs);
    cell = &c1;

    // refine the voronoi cell by removing tiny surfaces
    c2.init(-lx,lx,-ly,ly,-lz,lz);
  
    int nf = fs.size();
    // sort neighbors by area if asked to keep a minimum # of neighbors
    if (nminnei){
      for (int i=0; i<nf; i++)
      for (int j=i+1; j<nf; j++){
        if (fs[j] > fs[i]){
          double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
          int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
        }
      }
    }

    // add condition on surface
    double fcut = surf_min * cell->surface_area();
    for (int i=0; i<nf; i++){
      if (i < nminnei || fs[i] > fcut){
        int j = neigh[i];
  
        // apply pbc
        double xij = one->atpos[j][0]-x;
        while (xij > hx) xij -= lx;
        while (xij <-hx) xij += lx;

        double yij = one->atpos[j][1]-y;
        while (yij > hy) yij -= ly;
        while (yij <-hy) yij += ly;

        double zij = one->atpos[j][2]-z;
        while (zij > hz) zij -= lz;
        while (zij <-hz) zij += lz;
  
        c2.nplane(xij,yij,zij,j);
      }
    }
    c2.face_areas(fs);
    c2.neighbors(neigh);
    cell = &c2;

    vol = cell->volume();
    cell->face_freq_table(ff);
    int nn = ff.size()-1;
    for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
      
    // refine the voronoi cell if asked by skipping ultra short edges
    std::vector<double> vpos;
    std::vector<int>    vlst;
    double lcut2 = cell->total_edge_distance()*edge_min;
    lcut2 = lcut2*lcut2;
  
    cell->vertices(vpos);
    cell->face_vertices(vlst);
  
    nf = fs.size();
    int ford[nf];
    int k = 0, iface = 0;
    while (k < vlst.size()){
      int ned = vlst[k++];
      int nuc = 0;
      for (int ii=0; ii<ned; ii++){
        int jj = (ii+1)%ned;
        int v1 = vlst[k+ii], v2 = vlst[k+jj];
        double dx = vpos[v1*3]   - vpos[v2*3];
        double dy = vpos[v1*3+1] - vpos[v2*3+1];
        double dz = vpos[v1*3+2] - vpos[v2*3+2];
        double r2 = dx*dx+dy*dy+dz*dz;
        if (r2 <= lcut2) nuc++;
      }
      ford[iface++] = ned - nuc;
      k += ned;
    }
  
    for (int i=3; i<7; i++) index[i] = 0;
    for (int i=0; i<nf; i++){
      if (ford[i] < 7) index[ford[i]] += 1;
    }
  
    if (one->atsel[id] == 1){
      if (vlist.size() > 0){
        char str[MAXLINE];
        sprintf(str,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
        std::string vindex; vindex.assign(str);
        if (vlist.count(vindex) > 0) clist[id] = 1;
      } else clist[id] = 1;
    }
  
    nf = fs.size();
    if (nf > nmax){
      while (nf >= nmax) nmax += 12;
      nlist = memory->grow(nlist, nmax+1, natom+1, "nlist");
    }
    nlist[0][id] = nf;
    for (int i=0; i<nf; i++) nlist[i+1][id] = neigh[i];
  
  } while (cl.inc());

return;
}

/*------------------------------------------------------------------------------
 * Recursive method to find neighbors of id upto max shells
 *----------------------------------------------------------------------------*/
void Driver::FEFF_cluster(int il, const int max, int **nlist, int id,
  std::list<int> &clist, std::map<int,int> &myshell)
{
  if (++il > max) return;

  int nn = nlist[0][id];
  for (int ii=1; ii <= nn; ii++){
    int jd = nlist[ii][id];
    clist.push_back(nlist[ii][id]);
    if (myshell.count(jd)) myshell[jd] = MIN(myshell[jd], il);
    else myshell[jd] = il;

    FEFF_cluster(il, max, nlist, jd, clist, myshell);
  }

return;
}
/*------------------------------------------------------------------------------ */
