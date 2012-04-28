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

  // analyse command line options
  int iarg = 1;
  while (narg > iarg){
    if (strcmp(arg[iarg],"-h") == 0){
      help();

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
  while (1){
    int job = 0;
    printf("\n"); for (int i=0; i<20; i++) printf("====");
    printf("\nPlease select the job to perform:\n");
    printf("  1. convert to xyz format;\n");
    printf("  2. compute refined voro index info;\n");
    printf("  0. Exit;\nYour choice[0]: ");
    fgets(str,MAXLINE,stdin);

    char *ptr = strtok(str," \n\t\r\f");
    if (ptr) job = atoi(ptr);
    printf("You selected: %d\n", job);
    for (int i=0; i<20; i++) printf("===="); printf("\n");

    if      (job == 1) writexyz();
    else if (job == 2) voro();
    else break;
  }

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
    printf("\b%c", flag[++i%4]);
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
  printf("Please input your desired frame, or frame range to analyse: ");
  fgets(str,MAXLINE,stdin);
  int nw = count_words(str);
  if (nw < 1){
    istr = 0; iend = -1; inc = 1;

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
  inc = MAX(1,inc);

  printf("Frames from No. %d to No. %d with increment of %d will be analysed.\n", istr+1, iend+1, inc);
  for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Method to convert the dump file into xyz format
 *----------------------------------------------------------------------------*/
void Driver::writexyz()
{
  char str[MAXLINE];

  setrange();
  char *fout;
  printf("\nPlease input the output xyz file name [%s.xyz]: ", dump);
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL){
    fout = new char[strlen(dump)+5];
    sprintf(fout,"%s.xyz",dump);
  } else {
    fout = new char [strlen(ptr)]+1;
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
  printf("%d frames written to file %s\n", (iend-istr+1)/inc, fout);
  delete []fout;
return;
}

/*------------------------------------------------------------------------------
 * To display help info
 *----------------------------------------------------------------------------*/
void Driver::help()
{

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
