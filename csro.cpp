#include "driver.h"
#include "math.h"

/*------------------------------------------------------------------------------
 * Method to evaluate the chemical short range order based on the
 * voronoi neighbors for selected frames.
 * The chemical short range order is defined as (Warren-Cowley):
 *  eta = 1. - n_j2i/(N_i*c_j)
 * Where n_j2i is the number of type j neighbors of type i, N_i is the total #
 * of neighbors for type i, c_j is the concentration of type j
 *----------------------------------------------------------------------------*/
void Driver::csro()
{
  char str[MAXLINE];
  int job = 1;
  printf("\n"); for (int i=0; i<6; i++) printf("====");
  printf("   Chemical Short Range Order   ");
  for (int i=0; i<6; i++) printf("====");

  printf("\nPlease select your desired job:\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n");
  printf("  1. Overall CSRO based on Voronoi info;\n");
  printf("  2. Peratom CSRO based on Voronoi info;\n");
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

  // refinement info
  double mins[3]; mins[0] = 1.e-2; mins[1] =-1.; mins[2] =-1.;
  printf("Please input your criterion for tiny surfaces [%g]: ", mins[0]);
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr) mins[0] = atof(ptr);
  printf("Surfaces whose areas take less ratio than %lg will be removed!\n\n", mins[0]);

  // output file name for per atom CSRO
  FILE *fp; fp = NULL;
  printf("Please input the output file name [csro.dat]: ");
  if (count_words(fgets(str,MAXLINE, stdin)) < 1) strcpy(str,"csro.dat");
  ptr = strtok(str, " \n\t\r\f");
  fp = fopen(ptr, "w");

  int ntype = all[istr]->ntype;
  // work spaces
  bigint *NumType, **NumNei, *nnei;
  nnei   = memory->create(nnei, ntype+1, "csro:nnei"); // per atom info
  NumType= memory->create(NumType,ntype+1,"csro:NumType");
  NumNei = memory->create(NumNei,ntype+1,ntype+1,"csro:NumNei");
  for (int i=0; i<= ntype; i++){
    NumType[i] = 0;
    for (int j=0; j<= ntype; j++) NumNei[i][j] = 0;
  }
  // chemical composition for each frame
  double *cc;
  cc = new double [ntype+1];

  // now to do the real job
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // ntype of different frames must be the same
    if (one->ntype != ntype) continue;

    if (job == 2){
      fprintf(fp,"#Voronoi refinement info: surf_min%% = %g, edge_min%% = %g, nei_min = %d\n", mins[0], mins[1], int(mins[2]));
      fprintf(fp,"# Per atom CSRO info for frame %d; istep = %d\n# id type", img+1, one->tstep);
      for (int jp=1; jp<=ntype; jp++) fprintf(fp," csro-%d", jp);
      fprintf(fp,"\n");
    }

    // Compute Vorornoi info, so as to get the neighbor list
    one->ComputeVoro(mins,type2radius);

    // set local variables
    int *attyp = one->attyp;
    double **atpos = one->atpos;
    for (int i=1; i<=ntype; i++) NumType[i] += one->numtype[i];
    for (int i=1; i<=ntype; i++) cc[i] = double(one->numtype[i])/double(one->natom);

    for (int id = 1; id <= one->natom; id++){
      int ip = one->attyp[id];

      // zero per atom info
      for (int jp=1; jp<= ntype; jp++) nnei[jp] = 0;

      for (int ii=1; ii<=one->neilist[0][id]; ii++){
        int jd = one->neilist[ii][id];
        int jp = one->attyp[jd];
        nnei[jp]++;
      }
  
      // accumulate overall info
      for (int jp=1; jp<= ntype; jp++) NumNei[ip][jp] += nnei[jp];

      // output per atom info if needed
      if (job == 2){
        fprintf(fp,"%d %d", id, ip);
        int nn = 0;
        for (int jp=1; jp<= ntype; jp++) nn += nnei[jp];
        for (int jp=1; jp<= ntype; jp++) fprintf(fp," %g", 1.-double(nnei[jp])/(double(nn)*cc[jp]));
        fprintf(fp,"\n");
      }
    }
  }
  
  bigint ntotal = 0;
  double concentration[ntype+1];
  for (int i=1; i<= ntype; i++) ntotal += NumType[i];
  for (int i=1; i<= ntype; i++) concentration[i] = double(NumType[i])/double(ntotal);

  printf("\n"); for (int i=0; i<20; i++) printf("____"); printf("\nType         :");
  fprintf(fp,"# "); for (int i=0; i<20; i++) fprintf(fp,"____"); fprintf(fp,"\n# Type         :");
  if (type2atnum == NULL){
    for (int ip=1; ip<= ntype; ip++) printf("    %2d    ",ip);
    printf("\n"); for (int i=0; i<20; i++) printf("----");
    printf("\nConcentration:"); for (int ip=1; ip<= ntype; ip++) printf("%10.6f", concentration[ip]);
    printf("\n"); for (int i=0; i<20; i++) printf("____");
    printf("\nCSRO         :"); for (int jp=1; jp<= ntype; jp++) printf("    %2d    ",jp);

    for (int ip=1; ip<= ntype; ip++) fprintf(fp,"    %2d    ",ip);
    fprintf(fp,"\n# "); for (int i=0; i<20; i++) fprintf(fp,"----");
    fprintf(fp,"\n# Concentration:"); for (int ip=1; ip<= ntype; ip++) fprintf(fp,"%10.6f", concentration[ip]);
    fprintf(fp,"\n# "); for (int i=0; i<20; i++) fprintf(fp,"____");
    fprintf(fp,"\n# CSRO         :"); for (int jp=1; jp<= ntype; jp++) fprintf(fp,"    %2d    ",jp);

  } else {

    char ename[3];
    for (int ip=1; ip<= ntype; ip++){
      element->Num2Name(type2atnum[ip],ename);
      printf("    %2s    ", ename);
      fprintf(fp,"    %2s    ", ename);
    }
    printf("\n"); for (int i=0; i<20; i++) printf("----");
    printf("\nConcentration:"); for (int ip=1; ip<= ntype; ip++) printf("%10.6f", concentration[ip]);
    printf("\n"); for (int i=0; i<20; i++) printf("____"); printf("\nCSRO         :");
    for (int jp=1; jp<= ntype; jp++){
      element->Num2Name(type2atnum[jp],ename);
      printf("    %2s    ", ename);
    }

    fprintf(fp,"\n# "); for (int i=0; i<20; i++) printf("----");
    fprintf(fp,"\n# Concentration:"); for (int ip=1; ip<= ntype; ip++) fprintf(fp,"%10.6f", concentration[ip]);
    fprintf(fp,"\n"); for (int i=0; i<20; i++) fprintf(fp,"____"); fprintf(fp,"\n# CSRO         :");
    for (int jp=1; jp<= ntype; jp++){
      element->Num2Name(type2atnum[jp],ename);
      fprintf(fp,"    %2s    ", ename);
    }
  }

  printf("\n"); for (int i=0; i<20; i++) printf("----");
  fprintf(fp,"\n# "); for (int i=0; i<20; i++) fprintf(fp,"----");
  for (int ip=1; ip<= ntype; ip++){
    if (type2atnum == NULL){
      printf("\n      %2d     :", ip);
      fprintf(fp,"\n#       %2d     :", ip);
    } else {
      char ename[3]; element->Num2Name(type2atnum[ip], ename);
      printf("\n      %2s     :", ename);
      fprintf(fp,"\n#      %2s     :", ename);
    }
    ntotal = 0;
    for (int jp=1; jp<= ntype; jp++) ntotal += NumNei[ip][jp];
    for (int jp=1; jp<= ntype; jp++){
      printf("%10.6f", 1.-double(NumNei[ip][jp])/double(ntotal)/concentration[jp]);
      fprintf(fp, "%10.6f", 1.-double(NumNei[ip][jp])/double(ntotal)/concentration[jp]);
    }
  }
  fprintf(fp,"\n");
  fclose(fp);

  memory->destroy(nnei);
  memory->destroy(NumType);
  memory->destroy(NumNei);

  printf("\n"); for (int i=0; i<20; i++) printf("===="); printf("\n");

return;
}
