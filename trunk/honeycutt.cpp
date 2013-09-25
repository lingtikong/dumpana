#include "driver.h"
#include "math.h"
#include <vector>

/*------------------------------------------------------------------------------
 * Method to compute the Honeycutt-Andersen index based on the  voronoi neighbors
 * for selected frames. Only bonded pairs are analysed.
 *
 * The first index 1 means the selected pair are bonded, the second gives the
 * number of common neighbors between these two, and the third gives the number
 * of bonds formed between the common neighbors; the fourth are generally 1, while
 * for 144?, if all common neighbors are linked together, it is set to 2.
 *----------------------------------------------------------------------------*/
void Driver::honeycutt_andersen()
{
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 5; ++i) printf("====");
  printf("   Honeycutt-Andersen  Bond  Analysis   ");
  for (int i = 0; i < 5; ++i) printf("===="); printf("\n");

  // voronoi refinement
  set_cutoffs(0);
  one = all[istr];

  int unbond = 0;
  printf("\nWould you like to analyse un-bonded pairs? (y/n)[n]: ");
  fgets(str, MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr && (strcmp(ptr,"y") == 0 || strcmp(ptr,"Y")==0) )  unbond = 1;

  int outcomm = 0;
  printf("\nWould you like to output the common neighbors? (y/n)[n]: ");
  fgets(str, MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr && (strcmp(ptr,"y") == 0 || strcmp(ptr,"Y")==0) )  outcomm = 1;

  printf("Please input the file name to output the HA bond index info [ha.dat]: ");
  fgets(str, MAXLINE, stdin);
  char *fname = strtok(str, " \n\t\r\f");
  if (fname == NULL){
    strcpy(str,"ha.dat\n");
    fname = strtok(str, " \n\t\r\f");
  }

  FILE *fp = fopen(fname, "w");
  fprintf(fp, "# id  jd index\n"); fflush(fp);

  // now to loop over all asked images
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // Compute Voronoi neighbor info
    one->ComputeVoro(mins);

    fprintf(fp,"# frame number: %d\n", img);
    // now to analyse the Honeycutt-Andersen bond type info
    for (int id = 1; id <= one->natom; ++id){
      if (unbond){
        for (int jd = id+1; jd <= one->natom; ++jd){
          count_HA(id, jd, fp, outcomm);
        }

      } else {

        int nni = one->neilist[0][id];
        for (int kk = 1; kk <= nni; ++kk){
          int jd  = one->neilist[kk][id];
          if (id > jd) continue;

          count_HA(id, jd, fp, outcomm);
        }
      }
    }
    printf("Frame %d done, HA info written to: %s\n", img+1, fname);
    if (min_mem) one->FreeVoro();
  }
  printf("\n"); for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
  fclose(fp);

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
  if (flag) for (int ii = 0; ii < ncomm; ++ii) fprintf(fp," %d", comms[ii]);
  fprintf(fp,"\n");

return;
}

/*------------------------------------------------------------------------------*/
