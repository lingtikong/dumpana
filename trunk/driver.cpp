#include "driver.h"

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
/*------------------------------------------------------------------------------
 * Constructor of driver, main menu
 *----------------------------------------------------------------------------*/
Driver::Driver(int narg, char** arg)
{
  one = NULL; all.clear();
  dump = NULL;
  nframe = 0;
  memory = new Memory();
  flag_out = 0;
  int loop = 1;

  // analyse command line options
  int iarg = 1;
  while (narg > iarg){
    if (strcmp(arg[iarg],"-h") == 0){
      help();

    } else if (strcmp(arg[iarg], "-1") == 0){ // just do one analysis
      loop = 0;

    } else if (strcmp(arg[iarg], "-o") == 0){ // just do one analysis
      flag_out = atoi(arg[++iarg]);

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
    printf("\nPlease select the job to perform:\n");
    for (int i=0; i<20; i++) printf("----"); printf("\n");
    printf("  1. Voronoi diagram analysis;\n");
    printf("  2. convert to xyz format;\n");
    printf("  3. Average over frames;\n");
    printf("  4. CSRO based on voronoi neighbors;\n");
    printf("  5. Honeycutt-Andersen bond index analysis;\n");
    printf("  6. Common neighbor analysis/parameter;\n");
    printf("  0. Exit;\nYour choice [%d]: ", job);
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
      if (nsel > 0) writexyz();
      break;

    case 3:
      setrange();
      if (nsel > 0) avedump();
      break;

    case 4:
      setrange();
      if (nsel > 0) csro();
      break;

    case 5:
      setrange();
      if (nsel > 0) honeycutt_andersen();
      break;

    case 6:
      setrange();
      if (nsel > 0) Compute_CNACNP();
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
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
Driver::~Driver()
{
  if (dump) delete []dump;
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
 *----------------------------------------------------------------------------*/
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
    for (int i=1; i<=one->ntype; i++) printf("%d: %d  ", i, one->numtype[i]);
    printf("\n");
  }
  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Method to set the frame range for further analysis
 *----------------------------------------------------------------------------*/
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
 * Method to convert the dump file into xyz format
 *----------------------------------------------------------------------------*/
void Driver::writexyz()
{
  char str[MAXLINE];
  char *fout;
  printf("\nPlease input the output xyz file name [%s.xyz]: ", dump);
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL){
    fout = new char[strlen(dump)+5];
    sprintf(fout,"%s.xyz",dump);
  } else {
    fout = new char [strlen(ptr)+1];
    strcpy(fout, ptr);
  }
  FILE *fp = fopen(fout, "w");
  for (int img = istr; img<= iend; img += inc){
    one = all[img];
    fprintf(fp,"%d\n", one->natom);
    fprintf(fp,"Frame %d of %s, istep= %d\n", img+1, dump, one->tstep);

    one->dir2car();
    for (int i=1; i<=MIN(3,one->natom); i++){
      fprintf(fp,"%d %lg %lg %lg crystal_vector %d %lg %lg %lg\n",
      one->attyp[i], one->atpos[i][0], one->atpos[i][1], one->atpos[i][2], i,
      one->axis[i-1][0], one->axis[i-1][1], one->axis[i-1][2]);
    }
    for (int i=MIN(3,one->natom)+1; i<= one->natom; i++){
      fprintf(fp,"%d %lg %lg %lg\n", one->attyp[i], one->atpos[i][0], one->atpos[i][1], one->atpos[i][2]);
    }
  }
  fclose(fp);
  printf("Mission completed, %d frames written to file: %s\n", (iend-istr+1)/inc, fout);
  delete []fout;
return;
}

/*------------------------------------------------------------------------------
 * Method to write one frame of the dump file to a new file
 *----------------------------------------------------------------------------*/
void Driver::avedump()
{
  DumpAtom *first = all[istr];
  first->dir2car();
  double lx = first->lx;
  double ly = first->ly;
  double lz = first->lz;
  int nfirst = first->natom;
  double **atpos = memory->create(atpos,nfirst+1,3,"avedump:atpos");
  for (int ii=1; ii<=nfirst; ii++)
  for (int idim=0; idim<3; idim++) atpos[ii][idim] = 0.;

  int ncount = 1;
  for (int img = istr+inc; img<= iend; img += inc){
    one = all[img];
    if (one->natom != nfirst) continue;

    one->dir2car();

    lx += one->lx;
    ly += one->ly;
    lz += one->lz;

    for (int ii=1; ii<= nfirst; ii++){
      for (int idim=0; idim<3; idim++){
        double dx = one->atpos[ii][idim] - first->atpos[ii][idim];
        while (dx >= one->hbox[idim]) dx -= one->box[idim];
        while (dx < -one->hbox[idim]) dx += one->box[idim];

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

  lx /= double(ncount);
  ly /= double(ncount);
  lz /= double(ncount);
  
  char str[MAXLINE];
  char *fout;
  printf("\nPlease input the output xyz file name [dumpave.xyz]: ");
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL){
    strcpy(str,"dumpave.xyz");
    ptr = strtok(str, " \n\t\r\f");
  }
  fout = new char [strlen(ptr)+1];
  strcpy(fout, ptr);
  FILE *fp = fopen(fout, "w");
  fprintf(fp,"%d\n", nfirst);
  fprintf(fp,"Averaged over frames from %d to %d with incremental of %d: %lg %lg %lg\n",
  istr+1, iend+1, inc, lx, ly, lz);
  int ii = 1;
  if (nfirst >= 3){
    fprintf(fp,"%d %lg %lg %lg crystal_vector 1 %lg 0. 0.\n", first->attyp[ii], atpos[ii][0], atpos[ii][1], atpos[ii][2], lx); ii++;
    fprintf(fp,"%d %lg %lg %lg crystal_vector 2 0. %lg 0.\n", first->attyp[ii], atpos[ii][0], atpos[ii][1], atpos[ii][2], ly); ii++;
    fprintf(fp,"%d %lg %lg %lg crystal_vector 3 0. 0. %lg\n", first->attyp[ii], atpos[ii][0], atpos[ii][1], atpos[ii][2], lz); ii++;
  }
  for (int i=ii; i<= nfirst; i++) fprintf(fp,"%d %lg %lg %lg\n", first->attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
  fclose(fp);

  memory->destroy(atpos);
  delete []fout;

return;
}
/*------------------------------------------------------------------------------
 * To display help info
 *----------------------------------------------------------------------------*/
void Driver::help()
{
  printf("\n  dumpana\nCode to analyse the atom style dump file of lammps.\n");
  printf("\nUsage:\n    dumpana [options] [file]\n");
  printf("\nAvailable options:\n");
  printf("    -h       To display this help info;\n");
  printf("    -1       To tell the code to exit once an analysis is done;\n");
  printf("    -o n     To indicate wether to output the surface area ratio and/or edge length\n");
  printf("             ratio info or not: n = 0, neither; 1, surface; 2, edge; 3, both. default: 0.\n");
  printf("    file     Must be lammps atom style dump file, by default: dump.lammpstrj;\n");
  printf("\n\n");
  exit(0);
return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------------*/