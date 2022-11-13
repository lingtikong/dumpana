#include "driver.h"
#include "math.h"
#include <vector>

/*------------------------------------------------------------------------------
 * Method to compute the Honeycutt-Andersen index based on the neighbor list
 * for selected frames. Only bonded pairs are analysed.
 *
 * The first index 1 means the selected pair are bonded, the second gives the
 * number of common neighbors between these two, and the third gives the number
 * of bonds formed between the common neighbors; the fourth are generally 1, while
 * for 144?, if all common neighbors are linked together, it is set to 2.
 *----------------------------------------------------------------------------*/
void Driver::honeycutt_andersen()
{
  char str[MAXLINE], *ptr;
  printf("\n"); for (int i = 0; i < 5; ++i) printf("====");
  printf("   Honeycutt-Andersen  Bond  Analysis   ");
  for (int i = 0; i < 5; ++i) printf("===="); printf("\n");

  // choose method to calculate neighbor list
  one = all[istr];
  choose_neighbor_method(0);

  int unbond = 0;
  printf("\nWould you like to analyse un-bonded pairs? (y/n)[n]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr && (strcmp(ptr,"y") == 0 || strcmp(ptr,"Y")==0) )  unbond = 1;

  int outflag = 0;
  printf("\nWould you like to output the common neighbors? (y/n)[n]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr && (strcmp(ptr,"y") == 0 || strcmp(ptr,"Y")==0) )  outflag |= 1;

  printf("\nWould you like to output the atomic coordinates? (y/n)[n]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr && (strcmp(ptr,"y") == 0 || strcmp(ptr,"Y")==0) )  outflag |= 2;

  printf("Please input the file name to output the HA bond index info [ha.dat]: ");
  input->read_stdin(str);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL){
     strcpy(str,"ha.dat\n");
     ptr = strtok(str, " \n\t\r\f");
  }
  char *fname = new char [strlen(ptr)+1];
  strcpy(fname, ptr);
  ConfirmOverwrite(fname);

  FILE *fp = fopen(fname, "w");
  if (outflag & 2) fprintf(fp, "# id  jd index i-x i-y i-z j-x j-y j-z\n");
  else fprintf(fp, "# id  jd index\n");
  fflush(fp);

  // now to loop over all asked images
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // get neighbor list
    if (neighbor_method == 1) one->ComputeVoro(mins);
    else one->ComputeNeiList(r2cuts);

    fprintf(fp,"# frame number: %d\n", img);
    // now to analyse the Honeycutt-Andersen bond type info
    for (int id = 1; id <= one->natom; ++id){
      if (unbond){
        for (int jd = id+1; jd <= one->natom; ++jd){
          count_HA(id, jd, fp, outflag);
        }

      } else {

        int nni = one->neilist[0][id];
        for (int kk = 1; kk <= nni; ++kk){
          int jd  = one->neilist[kk][id];
          if (id > jd) continue;

          count_HA(id, jd, fp, outflag);
        }
      }
    }
    printf("Frame %d done, HA info written to: %s\n", img+1, fname);
    if (min_mem) one->FreeVoro();
  }
  printf("\n"); for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
  fclose(fp);
  delete []fname;

return;
}

/*------------------------------------------------------------------------------
 * Private method to compute the Honeycutt-Andersen index for a pair of atoms
 *------------------------------------------------------------------------------
 * id, jd  (in) : ID of the atom pair
 * fp      (in) : FILE to write the result
 * flag    (in) : 1, write the common neighors; 0, not write
 *----------------------------------------------------------------------------*/
void Driver::count_HA(int id, int jd, FILE * fp, const int flag)
{
  int ibond = 2 - one->bonded(id,jd);
  int **list = one->neilist;
  int nni = list[0][id];
  int nnj = list[0][jd];

  std::vector<int> comms;
  comms.clear();
  for (int ii = 1; ii <= nni; ++ii)
  for (int jj = 1; jj <= nnj; ++jj) if (list[ii][id] == list[jj][jd]) comms.push_back(list[ii][id]);
  int ncomm = comms.size();

  if (ibond == 2 && ncomm < 3) return;

  int nbond = 0;
  for (int mm = 0; mm < ncomm; ++mm)
  for (int nn = mm+1; nn < ncomm; ++nn) nbond += one->bonded(comms[mm],comms[nn]);

  int nconf = 1;
  // needs to distinct same ncomm-nbond for 144, 142
  // See Annals of Physics 324(2):332-342, 2009.
  if (ncomm == 4 && (nbond == 4 || nbond == 2) ){
    int ned[4]; ned[0] = ned[1] = ned[2] = ned[3] = 0;
    for (int mm = 0; mm < ncomm; ++mm)
    for (int nn = mm+1; nn < ncomm; ++nn){
      int md = comms[mm], nd = comms[nn];
      ned[mm] += one->bonded(md, nd);
      ned[nn] += one->bonded(md, nd);
    }
    int nmin = nbond/2;
    for (int mm = 0; mm < ncomm; ++mm) if (ned[mm] < nmin) nconf = 2;
  }

  fprintf(fp,"%d %d %d%d%d%d", id, jd, ibond, ncomm, nbond, nconf);
  if (flag & 2) fprintf(fp," %lg %lg %lg %lg %lg %lg", one->atpos[id][0], 
    one->atpos[id][1], one->atpos[id][2], one->atpos[jd][0], one->atpos[jd][1], one->atpos[jd][2]);
  if (flag & 1) for (int ii = 0; ii < ncomm; ++ii) fprintf(fp," %d", comms[ii]);
  fprintf(fp,"\n");

return;
}

/*------------------------------------------------------------------------------*/
