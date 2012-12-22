#include "driver.h"

/*------------------------------------------------------------------------------
 * Constructor of driver, main menu
 *------------------------------------------------------------------------------ */
Driver::Driver(int narg, char** arg)
{
  one = NULL; all.clear();
  dump = NULL;
  nframe = 0;
  type2atnum = NULL; type2radius = weighted = NULL;
  element = NULL;

  memory = new Memory();

  int loop = 1;

  flag_out = 0;
  flag_out |= OutFeff; // by default, feff.inp is written
  flag_out |= WtdVoro; // by default, weighted Voronoi will be performed if element mapping is done

  // analyse command line options
  int iarg = 1;
  while (narg > iarg){
    if (strcmp(arg[iarg],"-h") == 0){
      help();

    } else if (strcmp(arg[iarg], "-1") == 0){ // just do one analysis
      loop = 0;

    } else if (strcmp(arg[iarg], "-os") == 0){ // flag for Voronoi surface ratio outputs
      flag_out |= OutSurf;

    } else if (strcmp(arg[iarg], "-oe") == 0){ // flag for Voronoi edge ratio outputs
      flag_out |= OutEdge;

    } else if (strcmp(arg[iarg], "-ose") == 0){ // flat for Voronoi surface/edge ratio outputs
      flag_out |= OutSurf;
      flag_out |= OutEdge;

    } else if (strcmp(arg[iarg], "-s") == 0){ // to skip writting feff.inp data, while just output the CN info
      flag_out &= ~OutFeff;

    } else if (strcmp(arg[iarg], "-w") == 0){ // to do weighted Voronoi tessellation if possible, default
      flag_out |= WtdVoro;

    } else if (strcmp(arg[iarg], "-x") == 0){ // no weighted Voronoi tessellation even when possible
      flag_out &= ~WtdVoro;

    } else {
      break;
    }

    iarg++;
  }

  // get dump file name if supplied, othewise assume "dump.lammpstrj"
  if (narg > iarg){
    int n = strlen(arg[iarg])+1;
    dump = new char[n];
    strcpy(dump, arg[iarg]);
  } else {
    dump = new char[15];
    strcpy(dump,"dump.lammpstrj");
  }

  // read dump file
  readdump();
  if (nframe < 1) return;

  // main menu
  char str[MAXLINE];
  int job = 1;
  do {
    printf("\n"); for (int i=0; i<20; i++) printf("====");
    printf("\nPlease select your desired task to perform:\n");
    MainMenu();
    printf("  0. Exit.\nYour choice [%d]: ", job);
    fgets(str,MAXLINE,stdin);

    char *ptr = strtok(str," \n\t\r\f");
    if (ptr) job = atoi(ptr);
    printf("Your selection : %d\n", job);
    for (int i=0; i<20; i++) printf("===="); printf("\n");

    // main driver
    switch (job){
    case 1:
      setrange();
      if (nsel > 0) voro();
      break;

    case 2:
      setrange();
      if (nsel > 0) csro();
      break;

    case 3:
      setrange();
      if (nsel > 0) honeycutt_andersen();
      break;

    case 4:
      setrange();
      if (nsel > 0) Compute_CNACNP();
      break;

    case 5:
      setrange();
      if (nsel > 0) FEFF_main();
      break;

    case 6:
      setrange();
      if (nsel > 0) ClusterConnectivity();
      break;

    case 11:
      setrange();
      if (nsel > 0) writesel();
      break;

    case 12:
      setrange();
      if (nsel > 0) avedump();
      break;

    case 13:
      setrange();
      if (nsel > 0) paircorr();
      break;

    case 14:
      setrange();
      if (nsel > 0) strfac();
      break;

    case 15:
      setrange();
      if (nsel > 0) bonds();
      break;

    default:
      loop = 0;
      break;
    }
    job = 0;

  } while (loop);

return;
}

/*------------------------------------------------------------------------------
 * Method to display the main menu of the code
 *------------------------------------------------------------------------------ */
void Driver::MainMenu()
{
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. Voronoi diagram analysis;         |  11. Output selected frames;\n");
  printf("  2. Chemical Short Range Order;       |  12. Average over frames;   \n");
  printf("  3. Honeycutt-Andersen bond index;    |  13. Pair correlation function;\n");
  printf("  4. Common neighbor analysis;         |  14. Static structure factor;\n");
  printf("  5. Prepare for FEFF9;                |  15. Bond length/angles;\n");
  printf("  6. Voronoi cluster connectivity;     | \n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *------------------------------------------------------------------------------ */
Driver::~Driver()
{
  weighted = NULL;
  if (dump) delete []dump;
  if (element) delete element;
  if (type2atnum)  memory->destroy(type2atnum);
  if (type2radius) memory->destroy(type2radius);
  
  for (int img=0; img<nframe; img++){
    one = all[img];
    delete one;
  }
  one = NULL;
  all.clear();

  delete memory;
return;
}

/*------------------------------------------------------------------------------
 * Method to read the dump file from lammps
 *------------------------------------------------------------------------------ */
void Driver::readdump()
{
  // open file
  FILE *fp = fopen(dump,"r");
  if (fp == NULL){
    printf("\nError: file %s not found!\n", dump);
    return;
  }
  
  // set flags
  char flag[4];
  flag[0] = '-'; flag[1] = '\\'; flag[2] = '|'; flag[3] = '/';
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nNow to read atomic configurationss from file: %s...  ", dump);
  int i=0;

  // read file
  while (!feof(fp)){
    one = new DumpAtom(fp);
    if (one->initialized){
      all.push_back(one);
      one->iframe = all.size();
      one = NULL;
    } else {
      delete one;
      break;
    }
    printf("\b%c", flag[++i%4]); fflush(stdout);
  }
  fclose(fp);
  nframe = all.size();
  printf("\bDone!\n");

  // display read dump info
  for (int i=0; i<20; i++) printf("----");
  printf("\n  Total number of  frames  read: %d\n", nframe);
  if (nframe > 0){
    one = all[nframe-1];
    printf("  Number of atoms in last frame: %d\n", one->natom);
    printf("  Number of types in last frame: %d\n", one->ntype);
    printf("  Number of atoms for each type: ");
    for (int i=1; i<=one->ntype; i++) printf("%d, %d; ", i, one->numtype[i]);
    printf("\n");
  }
  MapType2Elem(0, one->ntype);
  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Method to set the frame range for further analysis
 *------------------------------------------------------------------------------ */
void Driver::setrange()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nTotal number of frames read: %d\n", nframe);
  printf("Please input your desired frame, or frame range to analyse [1]: ");
  fgets(str,MAXLINE,stdin);
  int nw = count_words(str);
  if (nw < 1){
    istr = 0; iend = 0; inc = 1;

  } else if (nw == 1){
    istr = iend = atoi(strtok(str, " \n\t\r\f"))-1;
    inc = 1;

  } else if (nw == 2){
    istr = atoi(strtok(str, " \n\t\r\f"))-1;
    iend = atoi(strtok(NULL," \n\t\r\f"))-1;
    inc = 1;

  } else {
    istr = atoi(strtok(str, " \n\t\r\f"))-1;
    iend = atoi(strtok(NULL," \n\t\r\f"))-1;
    inc  = atoi(strtok(NULL," \n\t\r\f"));
  }

  istr = MAX(0,istr);
  iend = MIN(iend,nframe-1);
  inc = MAX(1,inc);

  printf("Frames from No. %d to No. %d with increment of %d will be analysed.\n", istr+1, iend+1, inc);
  for (int i=0; i<20; i++) printf("===="); printf("\n");

  nsel = (iend-istr+1)/inc;

return;
}

/*------------------------------------------------------------------------------
 * Method to output the selected frames
 *------------------------------------------------------------------------------ */
void Driver::writesel()
{
  char str[MAXLINE]; int job = 1;
  printf("\n"); for (int i=0; i<7; i++) printf("====");
  printf(" Output Selected Frames ");
  for (int i=0; i<7; i++) printf("====");
  printf("\nPlease select your desired job:\n");
  printf("  1. Convert into xyz files;\n");
  printf("  2. Write as dump atom;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr) job = atoi(ptr);
  printf("Your selection : %d\n", job);
  if (job < 1 || job > 2){
    for (int i=0; i<20; i++) printf("===="); printf("\n");
    return;
  }
  printf("\n");

  if (job == 1){
    printf("\nPlease input the output file name [%s.xyz]: ", dump);
    if (count_words(fgets(str,MAXLINE,stdin)) < 1) sprintf(str,"%s.xyz", dump);
  } else {
    printf("\nPlease input the output file name [%s]: ", dump);
    if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str, dump);
  }
  ptr = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(ptr, "w");

  int nused = 0;

  if (job == 1){
    for (int img = istr; img<= iend; img += inc){
      one = all[img];
      fprintf(fp,"%d\n", one->natom);
      fprintf(fp,"Frame %d of %s, istep= %d\n", img+1, dump, one->tstep);
  
      one->dir2car();
      if (type2atnum == NULL){ // no elements assigned, print atomic type num as element
        for (int i=1; i<=MIN(3,one->natom); i++){
          fprintf(fp,"%d %lg %lg %lg crystal_vector %d %lg %lg %lg\n",
          one->attyp[i], one->atpos[i][0], one->atpos[i][1], one->atpos[i][2], i,
          one->axis[i-1][0], one->axis[i-1][1], one->axis[i-1][2]);
        }
        for (int i=MIN(3,one->natom)+1; i<= one->natom; i++){
          fprintf(fp,"%d %lg %lg %lg\n", one->attyp[i], one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
        }
      } else { // in case elements are assigned, print true element names
        char ename[3];
        for (int i=1; i<=MIN(3,one->natom); i++){
          element->Num2Name(type2atnum[one->attyp[i]], ename);
          fprintf(fp,"%2s %lg %lg %lg crystal_vector %d %lg %lg %lg\n",
          ename, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2], i,
          one->axis[i-1][0], one->axis[i-1][1], one->axis[i-1][2]);
        }
        for (int i=MIN(3,one->natom)+1; i<= one->natom; i++){
          element->Num2Name(type2atnum[one->attyp[i]], ename);
          fprintf(fp,"%2s %lg %lg %lg\n", ename, one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
        }
      }
      nused++;
    }

  } else if (job == 2){

    for (int img = istr; img <= iend; img += inc){
      one = all[img];
      one->car2dir();

      fprintf(fp,"ITEM: TIMESTEP\n%d\n", one->tstep);
      fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n", one->natom);
      if (one->triclinic){
        fprintf(fp,"ITEM: BOX BOUNDS pp pp pp xy xz yz\n");
        fprintf(fp,"%lg %lg %lg\n", one->xlo, one->xhi, one->xy);
        fprintf(fp,"%lg %lg %lg\n", one->ylo, one->yhi, one->xz);
        fprintf(fp,"%lg %lg %lg\n", one->zlo, one->zhi, one->yz);
      } else {
        fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(fp,"%lg %lg\n", one->xlo, one->xhi);
        fprintf(fp,"%lg %lg\n", one->ylo, one->yhi);
        fprintf(fp,"%lg %lg\n", one->zlo, one->zhi);
      }
      fprintf(fp,"ITEM: ATOMS\n");

      for (int id=1; id<= one->natom; id++) fprintf(fp,"%d %d %lg %lg %lg\n", id, one->attyp[id],
      one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);

      nused++;
    }
  }
  fclose(fp);
  printf("Mission completed, %d frames written to file: %s\n", nused, ptr);
  for (int i=0; i<20; i++) printf("===="); printf("\n"); 
return;
}

/*------------------------------------------------------------------------------
 * Method to write one frame of the dump file to a new file
 *------------------------------------------------------------------------------ */
void Driver::avedump()
{
  DumpAtom *first = all[istr];
  first->car2dir();
  double lx = first->lx, ly = first->ly, lz = first->lz;
  double xy = first->xy, xz = first->xz, yz = first->yz;
  int nfirst = first->natom;
  double **atpos = memory->create(atpos,nfirst+1,3,"avedump:atpos");
  for (int ii=1; ii<=nfirst; ii++)
  for (int idim=0; idim<3; idim++) atpos[ii][idim] = 0.;

  int ncount = 1;
  for (int img = istr+inc; img<= iend; img += inc){
    one = all[img];
    if (one->natom != nfirst) continue;

    one->car2dir();

    lx += one->lx; ly += one->ly; lz += one->lz;
    xy += one->xy; xz += one->xz; yz += one->yz;

    for (int ii=1; ii<= nfirst; ii++){
      for (int idim=0; idim<3; idim++){
        double dx = one->atpos[ii][idim] - first->atpos[ii][idim];
        while (dx >= 0.5) dx -= 1.;
        while (dx < -0.5) dx += 1.;

        atpos[ii][idim] += dx;
      }
    }
    ncount++;
  }
  if (ncount < 1) return;

  for (int ii=1; ii<= nfirst; ii++){
    atpos[ii][0] = atpos[ii][0]/double(ncount) + first->atpos[ii][0];
    atpos[ii][1] = atpos[ii][1]/double(ncount) + first->atpos[ii][1];
    atpos[ii][2] = atpos[ii][2]/double(ncount) + first->atpos[ii][2];
  }

  lx /= double(ncount); ly /= double(ncount); lz /= double(ncount);
  xy /= double(ncount); xz /= double(ncount); yz /= double(ncount);

  // convert fractional into cartesian
  for (int ii=1; ii<= nfirst; ii++){
    atpos[ii][0] = atpos[ii][0]*lx + atpos[ii][1]*xy + atpos[ii][2]*xz;
    atpos[ii][1] = atpos[ii][1]*ly + atpos[ii][2]*yz;
    atpos[ii][2] = atpos[ii][2]*lz;
  }
  
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<20; i++) printf("====");
  printf("\nPlease input the output xyz file name [dumpave.xyz]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str,"dumpave.xyz");

  char *fout = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(fout, "w");

  fprintf(fp,"%d\n", nfirst);
  fprintf(fp,"Averaged over frames from %d to %d with incremental of %d: %lg %lg %lg %lg %lg %lg\n",
  istr+1, iend+1, inc, lx, ly, lz, xy, xz, yz);
  int ii = 1;
  if (nfirst >= 3){
    fprintf(fp,"%d %lg %lg %lg crystal_vector 1 %lg 0.   0.\n", first->attyp[ii], atpos[ii][0], atpos[ii][1], atpos[ii][2], lx); ii++;
    fprintf(fp,"%d %lg %lg %lg crystal_vector 2 %lg %lg  0.\n", first->attyp[ii], atpos[ii][0], atpos[ii][1], atpos[ii][2], xy, ly); ii++;
    fprintf(fp,"%d %lg %lg %lg crystal_vector 3 %lg %lg %lg\n", first->attyp[ii], atpos[ii][0], atpos[ii][1], atpos[ii][2], xz, yz, lz); ii++;
  }
  for (int i=ii; i<= nfirst; i++) fprintf(fp,"%d %lg %lg %lg\n", first->attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
  fclose(fp);

  memory->destroy(atpos);
  for (int i=0; i<20; i++) printf("===="); printf("\n"); 

return;
}

/*------------------------------------------------------------------------------
 * To display help info
 *------------------------------------------------------------------------------ */
void Driver::help()
{
  printf("\n     ######                             #                 \n");
  printf("     #     #  #    #  #    #  #####    # #    #    #    ##  \n");
  printf("     #     #  #    #  ##  ##  #    #  #   #   ##   #   #  # \n");
  printf("     #     #  #    #  # ## #  #    # #     #  # #  #  #    #\n");
  printf("     #     #  #    #  #    #  #####  #######  #  # #  ######\n");
  printf("     #     #  #    #  #    #  #      #     #  #   ##  #    #\n");
  printf("     ######    ####   #    #  #      #     #  #    #  #    #\n");
  printf("\nCode to analyse the atom style dump files of lammps. Functions available:\n");
  MainMenu();
  printf("\nUsage:\n    dumpana [options] [file]\n\nAvailable options:\n");
  printf("    -h       To display this help info;\n");
  printf("    -1       To tell the code to exit once an analysis is done;\n");
  printf("    -os      To output the surface area ratios when analyze Voronoi diagram;\n");
  printf("    -oe      To output the edge length ratios when analyze Voronoi diagram;\n");
  printf("    -ose     To output both the surface area and the edge length ratios\n");
  printf("             when analyze Voronoi diagram;\n");
  printf("    -s       To skip writing feff.inp files when preparing FEFF for desired voronoi\n");
  printf("             clusters; instead, output the CN info only.\n");
  printf("    -w/-x    To or not to perform weighted Voronoi tessellation, if possible;\n");
  printf("             by default, weigthed will be done if element mapping has been done;\n");
  printf("    file     Must be lammps atom style dump file, by default: dump.lammpstrj.\n");
  printf("\n\n");
  exit(0);
return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *------------------------------------------------------------------------------ */
int Driver::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}

/*------------------------------------------------------------------------------
 * Method to get the mapping between atomic type and element names
 *------------------------------------------------------------------------------
 * flag  (in)  : 0, called when knows the # of atom types; 1, otherwise
 * ntype (in)  : total number of atom types
 *------------------------------------------------------------------------------ */
void Driver::MapType2Elem(const int flag, const int ntype)
{
  if (element) delete element; element = NULL;
  if (type2atnum) memory->destroy(type2atnum); type2atnum = NULL;
  if (type2radius)memory->destroy(type2radius);type2radius= NULL;

  char str[MAXLINE];
  if (flag == 0){
    printf("\nIf you want to map the atomic types to elements, input the element\n");
    printf("symbols in sequence now: ");
  } else {
    printf("Total number of atomic types in the system are %d.\n", ntype); 
    printf("Please input the element symbol for each atomic type in sequence: ");
  }

  if (count_words(fgets(str,MAXLINE,stdin)) >= ntype){
    if (element) delete element;
    if (type2atnum)  memory->destroy(type2atnum);
    if (type2radius) memory->destroy(type2radius);
    element = new ChemElements();
    type2atnum = memory->create(type2atnum, ntype+1, "type2atnum");
    type2radius= memory->create(type2radius,ntype+1, "type2radius");

    char *ptr = strtok(str," \n\t\r\f");
    for (int ip=1; ip<= ntype; ip++){
      type2atnum[ip] =  element->Name2Num(ptr);
      type2radius[ip]=  element->Name2Radius(ptr);
      ptr = strtok(NULL, " \n\t\r\f");
    }
    printf("\nThe atomic types are assigned as:");
    for (int ip=1; ip<=ntype; ip++){
      char ename[3];
      int num = type2atnum[ip];
      element->Num2Name(num, ename);
      printf(" %d -> %s(%d);", ip, ename, num);
    } printf("\n");
  }

return;
}

/*------------------------------------------------------------------------------
 * Method to show the atomic radius  for each atomic type
 *------------------------------------------------------------------------------ */
void Driver::ShowRadius4Voro()
{
  weighted = NULL;
  if ((flag_out&WtdVoro)==0 || one==NULL || type2radius==NULL) return;

  printf("\nWeighted Voronoi tesselation will be performed, with atomic radii:\n");
  for (int ip = 1; ip <= one->ntype; ip++){
    char ename[3];
    int num = type2atnum[ip];
    element->Num2Name(num, ename);
    printf("  R(%s) = %g A;", ename, type2radius[ip]);
  } printf("\n\n");

  weighted = type2radius;

return;
}
/*----------------------------------------------------------------------------*/
