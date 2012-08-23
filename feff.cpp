#include "driver.h"
#include "voro++.hh"
#include "math.h"
#include <list>

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
    } else strcpy(selcmd,"type = 1");

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
  int flag_voro = 0, vindex[4]; bigint vprod;
  printf("\nIf you want to futher refine the absorbing atoms by their Voronoi index,\n");
  printf("please input the four index numbers now: ");
  if (count_words(fgets(str,MAXLINE,stdin)) >= 4){
    flag_voro = 1;
    char *ptr = strtok(str," \n\t\r\f");
    for (int i=0; i<3; i++){vindex[i] = atoi(ptr); ptr = strtok(NULL," \n\t\r\f");}
    vindex[3] = atoi(ptr);
    vprod = vindex[0]*1000000 + vindex[1]*10000 + vindex[2]*100 + vindex[3];

    printf("\nAtoms selected by `%s` with index <%d,%d,%d,%d> will be the absorbing atoms.\n",
      selcmd, vindex[0], vindex[1], vindex[2], vindex[3]);
  }

  // job type, CFAVERAGE or Clusters
  int flag_pbc = 0;
  printf("\nWould you like to run (1), an auto `average` over all absorbing atoms for\n");
  printf("each frame, or (2), individual calculations for each cluster? (1/2)[2]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    flag_pbc = atoi(ptr)%2;
  }
  printf("Your selection : %d\n", 2-flag_pbc);

  // XANES or EXAFS
  int job = 2;
  printf("\nWould you like to do a (1) XANES or (2) EXAFS calculation? (1/2)[%d]: ", job);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    job = 2-atoi(ptr)%2;
  }
  printf("Your selection : %d\n", job);

  // working directory: directory to write all output files
  printf("\nPlease define the working directory, i.e., directory to output all\n");
  printf("output files [.]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    strcpy(workdir, ptr);
  } else strcpy(workdir, ".");
  printf("The working directory will be: %s\n", workdir);

  // some common variables
  char dirname[MAXLINE], mkdir[MAXLINE], fname[MAXLINE];
  FILE *fp;

  // file DirList records all directories created
  sprintf(fname,"%s/DirList", workdir);
  FILE *dlist = fopen(fname, "w");
  int ndir = 0;

  // real job
  if (flag_pbc==1 && flag_voro == 0){ // no Voronoi needed

    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      one->selection(selcmd);
  
      // open the feff.inp file for current frame
      sprintf(dirname, "%s/Frame%d", workdir, img+1);
      strcpy(mkdir,"mkdir -p "); strcat(mkdir,dirname);
      system(mkdir); ndir++; fprintf(dlist, "%s\n", dirname); // create direcotry for current frame
      strcpy(fname, dirname); strcat(fname,"/feff.inp");
      fp = fopen(fname, "w");

      // write the potential part of feff.inp
      fprintf(fp,"* Frame %d (MD steps: %d), %d types, %d atoms.\n",
        img+1, one->tstep, one->ntype, one->natom);
      fprintf(fp,"* Averaged over %d type-%d atoms.\n\n", one->nsel, AbsorberType);
      fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
      char ename[3];
      for (int ip=1; ip <= one->ntype; ip++){
        element->Num2Name(type2atnum[ip], ename);
        fprintf(fp,"%4d   %3d  %3s  -1   3\n", ip, type2atnum[ip], ename);
      }
      element->Num2Name(type2atnum[AbsorberType], ename);
      fprintf(fp,"%4d   %3d  %3s  -1   3\n", one->ntype+1, type2atnum[AbsorberType], ename);

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
      FEFF_voro(flag_voro, vprod, nmax, neilist, cenlist);

      // analyse the result
      int nc = 0;
      for (int ii=1; ii<=natom; ii++) nc += cenlist[ii];
      if (nc < 1){
        memory->destroy(neilist);
        memory->destroy(cenlist);
        continue;
      }

      if (flag_pbc){
        // open the feff.inp file for current frame
        sprintf(dirname, "%s/Frame%d", workdir, img+1);
        strcpy(mkdir,"mkdir -p "); strcat(mkdir,dirname);
        system(mkdir); ndir++; fprintf(dlist, "%s\n", dirname); // create direcotry for current frame
        strcpy(fname, dirname); strcat(fname,"/feff.inp");
        fp = fopen(fname, "w");
  
        // write the potential part of feff.inp
        fprintf(fp,"* Frame %d (MD steps: %d), %d types, %d atoms.\n",
          img+1, one->tstep, one->ntype, one->natom);
        fprintf(fp,"* Averaged over %d type %d atoms.\n\n", nc, AbsorberType);
        fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
        char ename[3];
        for (int ip=1; ip <= one->ntype; ip++){
          element->Num2Name(type2atnum[ip], ename);
          fprintf(fp,"%4d   %3d  %3s  -1   3\n", ip, type2atnum[ip], ename);
        }
        element->Num2Name(type2atnum[AbsorberType], ename);
        fprintf(fp,"%4d   %3d  %3s  -1   3\n", one->ntype+1, type2atnum[AbsorberType], ename);
  
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
  
        // write other common info
        FEFF_input(job, fp);
  
        // clsoe the file
        fclose(fp);
  
      } else {

        int const MaxShell = 12;
        int level = 3;
        std::list<int> cluster;

        printf("\nHow many shells would you like to include into the cluster? (1-%d)[3]: ",MaxShell);
        if (count_words(fgets(str,MAXLINE,stdin)) > 0){
          char *ptr = strtok(str," \n\t\r\f");
          level = atoi(ptr);
        }
        level = MIN(MaxShell,MAX(1,level));
        char Nos[][3] = {"","st","nd","rd","th","th","th","th","th","th","th","th","th","th","th"};
        printf("Atoms upto the %d%s shell will be included into the clsuter.\n", level, Nos[level]);
        
        for (int id=1; id <= one->natom; id++){
          if (cenlist[id] == 0) continue;

          // find atoms in the cluster
          int ilevel = 0;
          cluster.clear();
          FEFF_cluster(ilevel, level, neilist, id, cluster);
          
          cluster.sort(); cluster.unique();
          int nclus = cluster.size();

          // open the feff.inp file for current frame
          sprintf(dirname, "%s/F%dA%d", workdir, img+1, id);
          strcpy(mkdir,"mkdir -p "); strcat(mkdir,dirname);
          system(mkdir); ndir++; fprintf(dlist, "%s\n", dirname); // create direcotry for current frame
          strcpy(fname, dirname); strcat(fname,"/feff.inp");
          fp = fopen(fname, "w");
  
          // write the potential part of feff.inp
          char ename[3];
          element->Num2Name(type2atnum[one->attyp[id]], ename);
          fprintf(fp,"* Cluster centered on atom %d (%s) of frame %d (MD steps: %d)\n",
            id, ename, img+1, one->tstep);
          fprintf(fp,"* Total number of atoms in cluster: %d\n\n", nclus);

          fprintf(fp,"POTENTIALS\n*  ipot Z   tag lmax1 lmax2\n");
          fprintf(fp,"%4d   %3d  %3s  -1   3\n", 0, type2atnum[one->attyp[id]], ename);
          for (int ip=1; ip <= one->ntype; ip++){
            element->Num2Name(type2atnum[ip], ename);
            fprintf(fp,"%4d   %3d  %3s  -1   3\n", ip, type2atnum[ip], ename);
          }
  
          // atomic positions
          one->dir2car();
          fprintf(fp,"\n* Atomci positions in Angstrom\nATOMS\n* x y z ipot\n");
          for (std::list<int>::iterator it = cluster.begin(); it != cluster.end(); it++){
            int jd = *it;
            int jp = one->attyp[jd];
            if (jd == id) jp = 0;
            double dx[3];
            for (int idim=0; idim<3; idim++){
              dx[idim] = one->atpos[jd][idim] - one->atpos[id][idim];
              while (dx[idim] > one->hbox[idim]) dx[idim] -= one->box[idim];
              while (dx[idim] <-one->hbox[idim]) dx[idim] += one->box[idim];
            }
            fprintf(fp,"%15.8f %15.8f %15.8f %d\n", dx[0], dx[1], dx[2], jp);
          }

          // write other common info
          FEFF_input(job, fp);
  
          // clsoe the file
          fclose(fp);
        }
      }

      memory->destroy(neilist);
      memory->destroy(cenlist);
    } // end of loop over frames
  }

  printf("\nJob done, %d directorys are created and listed in `%s/DirList`.\n", ndir, workdir);
  for (int i=0; i<20; i++) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------
 * Method to generate the common options for XANES or EXAFS calculations.
 *----------------------------------------------------------------------------*/
void Driver::FEFF_input(int job, FILE *fp)
{
  fprintf(fp, "\n* Other structural information cards\n");
  fprintf(fp, "* CIF cif_file * not needed here\n");
  fprintf(fp, "\n* Tell the code that a real-space calculation is asked\nREAL\n");
  fprintf(fp, "\n* If necessary, provided some comments/references as titles\n* TITLE any_thing\n");
  fprintf(fp, "\n* scale factor of the distances\nRMULTIPLIER 1.00\n");
  fprintf(fp, "\n\n* Spectrum Information cards\n");
  if (job == 1){
    fprintf(fp, "* XANES [xkmax xkstep vixan]\nXANES\n");
    fprintf(fp, "* SCF rfms1\n* FMS rfms\n* RPATH rfms\n");
  } else {
    fprintf(fp,"* EXAFS [xkmax]\nEXAFS\n");
    fprintf(fp,"* RPATH rfms\n* SCF rfms1\n");
  }
  fprintf(fp, "\n* If lDOS needed, uncomment and modify the following line\n");
  fprintf(fp, "*    emin emax eimag\n* LDOS -20 20 0.2");
  fprintf(fp, "\n\n* Control cards\n");
  fprintf(fp, "*         ipot ixsph ifms ipaths igenfmt iff2x\n");
  fprintf(fp, "  CONTROL  1    1     1    1      1       1\n");
  fprintf(fp, "* PRINT \n");
  fprintf(fp, "\n* DIMS nmax lmax\n* nmax: maximum # of atoms in the cluster for FMS matrix inversion\n");
  fprintf(fp, "* lmax: maximal l-value for the potentials and Green`s function\n");
  fprintf(fp, "\n* EGRID\n* This card can be used to customize the energy grid, which is followed by lines\n");
  fprintf(fp, "* specifying the type, interval, incremental of the grid:\n");
  fprintf(fp, "*    grid_type grid_min grid_max grid_step\n");
  fprintf(fp, "* grid_type can be e_grid, k_grid, or exp_grid. In the latter case, Ei = Emin + exp(Estep * i) - 1\n");
  fprintf(fp, "* It can be also user_grid, followed by energy points for the grid, one point in each line.\n");
  fprintf(fp, "\n\n* Some other standard cards\n");
  fprintf(fp, "* AFOLP folpx\n* This automatically overlaps all muffin-tins to a specified maximum value, default 1.15\n");
  fprintf(fp, "* COREHOLE type\n* Determines how the core state is treated: none, equivalent to the NOHOLE card, meaning\n");
  fprintf(fp, "*  there is no core hole; RPA, meaning the screen module calculates an RPA-screened core hole; FSR, a simple\n");
  fprintf(fp, "*  Final state rule core hole (default). It is recommended to use RPA for k-space calculations.\n");
  fprintf(fp, "* EDGE label s02\n* Sets the edge: label is `K`, `L1`...\nEDGE K 1.0\n");
  fprintf(fp, "* SCF rfms1 [lfms1 nscmt ca nmix]\n* Controls FEFF`s automated SCF calculations.\n");
  fprintf(fp, "* rfms1 specifies the radius of the cluster for full multiple scattering during SCF loop.\n");
  fprintf(fp, "* S02 s02\n* Specifies the amplitude reduction factor S02.\n");
  fprintf(fp, "\n* FMS rfms [lfms2 minv toler1 toler2 rdirec]\n");
  fprintf(fp, "* Compute full multiple scattering within a sphere of radius rfms centered on the absorbing atom\n");
  fprintf(fp, "* RPATH rpath\n* Determines the maximum effective (half-path) distance of a given path.\n");
  fprintf(fp, "\n* DEBYE temperature Debye-temperature\n* Used to calculate the Debye-Waller factor.\n");
  fprintf(fp, "\n* END\n");

return;
}

/*------------------------------------------------------------------------------
 * Method to find the neighbor list based on refined Voronoi calculation
 *----------------------------------------------------------------------------*/
void Driver::FEFF_voro(int flag, int vp, int &nmax, int **nlist, int *clist)
{
  char str[MAXLINE];
  double surf_min = 1.e-4, edge_min = 1.e-4;
  int flag_min = 0, nminnei = 14;

  printf("\nRefined Voronoi tesselation will be needed to procceed.\n");
  printf("Now please input your criterion for tiny surfaces [%g]: ", surf_min);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) surf_min = atof(ptr);
  printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", surf_min);

  printf("Would you like to keep a minimum # of neighbors? (y/n)[n]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr != NULL && strcmp(ptr,"y") == 0){
    flag_min = 1;
    printf("Please input your minimum # of neighbors, I would recommend 14 for bcc\n");
    printf("lattice, 12 for hcp and fcc. please input your number [%d]: ", nminnei);
    fgets(str,MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) nminnei = atoi(ptr);
    printf("\nA minimum number of %d neighobrs will be kept no matter how tiny the surface is.\n", nminnei);
  } else nminnei = 0;
  printf("\n");

  printf("Please input your criterion for ultra short   [%g]: ", edge_min);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) edge_min = atof(ptr);
  printf("Edges whose length takes less ratio than %lg will be skiped!\n\n", edge_min);

  // set local variables
  double xlo = one->xlo, xhi = one->xhi, lx = one->lx;
  double ylo = one->ylo, yhi = one->yhi, ly = one->ly;
  double zlo = one->zlo, zhi = one->zhi, lz = one->lz;
  double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;
  int n = one->natom;

  int *attyp = one->attyp;
  double **atpos = one->atpos;

  // need cartesian coordinates
  one->dir2car();

  // compute optimal size for container, and then contrust it
  double l = pow(double(n)/(5.6*lx*ly*lz), 1./3.);
  int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);
  voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

  // put atoms into the container
  for (int i=1; i<= n; i++) con.put(i, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);

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
    if (flag_min){
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
      if (flag){
        int voro = index[3]*1000000 + index[4]*10000 + index[5]*100 + index[6];
        if (voro == vp) clist[id] = 1;
      } else clist[id] = 1;
    }
  
    nf = fs.size();
    if (nf > nmax){
      while (nf >= nmax) nmax += 12;
      nlist = memory->grow(nlist, nmax+1, one->natom+1, "nlist");
    }
    nlist[0][id] = nf;
    for (int i=0; i<nf; i++) nlist[i+1][id] = neigh[i];
  
  } while (cl.inc());

return;
}

/*------------------------------------------------------------------------------
 * Recursive method to find neighbors of id upto max shells
 *----------------------------------------------------------------------------*/
void Driver::FEFF_cluster(int il, int max, int **nlist, int id, std::list<int> &clist)
{
  if (++il > max) return;

  int nn = nlist[0][id];
  for (int ii=1; ii <= nn; ii++){
    int jd = nlist[ii][id];
    clist.push_back(nlist[ii][id]);

    FEFF_cluster(il, max, nlist, jd, clist);
  }

return;
}
/*------------------------------------------------------------------------------ */
