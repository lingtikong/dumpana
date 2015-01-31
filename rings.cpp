#include "driver.h"

/*------------------------------------------------------------------------------
 * Method to prepare for ring statistics by the RINGS code.
 *----------------------------------------------------------------------------*/
void Driver::rings()
{
  char str[MAXLINE];
  printf("\n"); for (int i = 0; i < 7; ++i) printf("===="); printf(" Preparation for  RINGS ");
  for (int i = 0; i < 7; ++i) printf("===="); printf("\n");
  printf("Files needed by RINGS will be generated; please note that this works well only if\n");
  printf("the trajectory is from a constant volume simulation. Besides, the values provided\n");
  printf("for the first valley of g(r) is just a guess, one should measure and replace it\n");
  printf("for serious calculations.\n");

  one = all[istr];
  // map atomic type to elements
  if (type2atnum == NULL){
    printf("Mapping of atomic types to actual elements are required for preparing RINGS.\n");
    MapType2Elem(1, one->ntype); printf("\n");
  }

  // voro refinement info
  set_cutoffs(0);

  // selection of atoms for each frame
  char workdir[MAXLINE], finp[MAXLINE], fcfg[MAXLINE];
  printf("\nPlease input the directory name to output the files [.]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    strcpy(workdir, ptr);
  } else strcpy(workdir, ".");
  printf("The working directory will be: %s\n", workdir);

  printf("\nPlease input the name of the input file for rings [input]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    strcpy(finp, ptr);
  } else strcpy(finp, "input");
  printf("The input file for RINGS will be: %s\n", finp);

  printf("\nPlease input the name of the configuration file [atomcfg.xyz]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");
    strcpy(fcfg, ptr);
  } else strcpy(fcfg, "atomcfg.xyz");
  printf("The configuration file for RINGS will be: %s\n\n", fcfg);

  // open the necessary files
  strcpy(str,"mkdir -p "); strcat(str, workdir); strcat(str, "/data"); system(str);

  char fname[MAXLINE];
  sprintf(fname,"%s/%s", workdir, finp);
  FILE *fp = fopen(fname, "w");

  sprintf(fname,"%s/data/%s", workdir, fcfg);
  FILE *fc = fopen(fname, "w");

  // variables on nearest neighbor distances
  double ***bl;
  int ***nb;
  memory->create(bl, 2, one->ntype+1, one->ntype+1, "bl");
  memory->create(nb, 2, one->ntype+1, one->ntype+1, "nb");
  for (int i = 0; i <= one->ntype; ++i)
  for (int j = 0; j <= one->ntype; ++j){
    bl[0][i][j] = bl[1][i][j] = 0.;
    nb[0][i][j] = nb[1][i][j] = 0;
  }

  int nused = 0;
  char ename[3];
  // Loop over all frames
  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];
    printf("  Now to process frame %4d of %s ...", one->iframe, one->fname);

    // compute the Voronoi info, so as to get all related info
    one->ComputeVoro(mins);

    // Cartesian coordinate needed
    one->dir2car();

    // write header of the configuration file
    fprintf(fc, "%d\nFrame %d from %s.\n", one->natom, one->iframe, one->fname);

    // loop over all atoms
    for (int id = 1; id <= one->natom; ++id){
      int ip = one->attyp[id];
      element->Num2Name(type2atnum[ip], ename);

      fprintf(fc,"%3s %18.12f %18.12f %18.12f\n", ename, one->atpos[id][0], one->atpos[id][1], one->atpos[id][2]);

      list<int> cluster;
      map<int,int> shell;
      cluster.clear(); shell.clear();
      one->voro_cluster(0, 2, id, cluster, shell);

      for (list<int>::iterator it = cluster.begin(); it != cluster.end(); ++it){
        int jd = *it;
        if (id > jd) continue;

        int jp = one->attyp[jd];
        int is = shell[jd]-1;

        double dx = one->atpos[jd][0] - one->atpos[id][0];
        double dy = one->atpos[jd][1] - one->atpos[id][1];
        double dz = one->atpos[jd][2] - one->atpos[id][2];
        one->ApplyPBC(dx, dy, dz);
        double r2 = dx*dx + dy*dy + dz*dz;
        double r  = sqrt(r2);

        bl[is][ip][jp] += r; ++nb[is][ip][jp];
        bl[is][jp][ip] += r; ++nb[is][jp][ip];
        bl[is][0][0]   += r; ++nb[is][0][0];
      }
      cluster.clear(); shell.clear();
    }

    if (min_mem) one->FreeVoro();
    printf("  Done!\n");
    ++nused;
  } // end of loop over frames

  fclose(fc);

  one = all[istr];
  // to write the input file
  fprintf(fp, "#######################################\n");
  fprintf(fp, "#       R.I.N.G.S. input file         #\n");
  fprintf(fp, "#######################################\n");
  fprintf(fp, "Configurations-%d\n", nused);
  fprintf(fp, "%d\n", one->natom);
  fprintf(fp, "%d\n", one->ntype);
  for (int ip = 1; ip <= one->ntype; ++ip){
    element->Num2Name(type2atnum[ip], ename);
    fprintf(fp, "%s ", ename);
  }
  fprintf(fp, "\n%d\n1\n", nused);
  for (int idim = 0; idim < 3; ++idim){
    fprintf(fp, "%18.12f %18.12f %18.12f\n", one->axis[idim][0], one->axis[idim][1], one->axis[idim][2]);
  }
  fprintf(fp, "2.0    # MD timestep in unit of fs\n");
  fprintf(fp, "ANI\n%s\n", fcfg);
  fprintf(fp, "201    # n.p.t for g(r)\n");
  fprintf(fp, "500    # n.p.t for S(q)\n");
  fprintf(fp, "20     # q_max for S(q)\n");
  fprintf(fp, "0.125  # Smoothing factor for S(q)\n");
  fprintf(fp, "180    # Angular discretization\n");
  fprintf(fp, "10     # Real space discretization for voids and ring stat.\n");
  fprintf(fp, "15     # Search depeth for ring stat.\n");
  fprintf(fp, "#######################################\n");
  char ea[3], eb[3];
  for (int ip = 1;  ip <= one->ntype; ++ip)
  for (int jp = ip; jp <= one->ntype; ++jp){
    if (nb[0][ip][jp] > 0) bl[0][ip][jp] /= double(nb[0][ip][jp]);
    if (nb[1][ip][jp] > 0) bl[1][ip][jp] /= double(nb[1][ip][jp]);
    element->Num2Name(type2atnum[ip], ea);
    element->Num2Name(type2atnum[jp], eb);

    fprintf(fp, "%s %s %lg\n", ea, eb, (bl[0][ip][jp] + 4.*bl[1][ip][jp])*0.2);
  }
  if (nb[0][0][0] > 0) bl[0][0][0] /= double(nb[0][0][0]);
  if (nb[1][0][0] > 0) bl[1][0][0] /= double(nb[1][0][0]);
  fprintf(fp, "Grtot   %lg\n", (bl[0][0][0] + 4.*bl[1][0][0])*0.2);
  fprintf(fp, "#######################################\n");

  fclose(fp);
  memory->destroy(bl);
  memory->destroy(nb);

  // write the options file
  sprintf(fname, "%s/options", workdir);
  fp = fopen(fname, "w");
  fprintf(fp, "#######################################\n"); 
  fprintf(fp, "        R.I.N.G.S. options file       #\n");
  fprintf(fp, "#######################################\n");
  fprintf(fp, " PBC             .true.               #\n");
  fprintf(fp, " Frac            .false.              #\n");
  fprintf(fp, " g(r)            .false.              #\n");
  fprintf(fp, " S(q)            .false.              #\n");
  fprintf(fp, " S(k)            .false.              #\n");
  fprintf(fp, " gfft(r)         .false.              #\n");
  fprintf(fp, " MSD             .false.              #\n");
  fprintf(fp, " atMSD           .false.              #\n");
  fprintf(fp, " Bonds           .false.              #\n");
  fprintf(fp, " Angles          .false.              #\n");
  fprintf(fp, " Rings           .true.               #\n");
  fprintf(fp, "---- ! Ring statistics options ! ------\n");
  fprintf(fp, " Species           0                  #\n");
  fprintf(fp, " ABAB            .false.              #\n");
  fprintf(fp, " Rings0          .false.              #\n");
  fprintf(fp, " Rings1          .false.              #\n");
  fprintf(fp, " Rings2          .false.              #\n");
  fprintf(fp, " Rings3          .false.              #\n");
  fprintf(fp, " Rings4          .true.               #\n");
  fprintf(fp, " Prim_Rings      .true.               #\n");
  fprintf(fp, " Str_Rings       .false.              #\n");
  fprintf(fp, " BarycRings      .false.              #\n");
  fprintf(fp, " Prop-1          .false.              #\n");
  fprintf(fp, " Prop-2          .false.              #\n");
  fprintf(fp, " Prop-3          .false.              #\n");
  fprintf(fp, " Prop-4          .false.              #\n");
  fprintf(fp, " Prop-5          .false.              #\n");
  fprintf(fp, "---------------------------------------\n");
  fprintf(fp, " Vacuum          .false.              #\n");
  fprintf(fp, "#######################################\n");
  fprintf(fp, "         Outputting options           #\n");
  fprintf(fp, "#######################################\n");
  fprintf(fp, " Evol            .false.              #\n");
  fprintf(fp, " Dxout           .false.              #\n");
  fprintf(fp, "! OpenDX visualization options !  -----\n");
  fprintf(fp, " RadOut          .false.              #\n");
  fprintf(fp, " RingsOut        .true.               #\n");
  fprintf(fp, " DRngOut         .false.              #\n");
  fprintf(fp, " VoidsOut        .false.              #\n");
  fprintf(fp, " TetraOut        .false.              #\n");
  fprintf(fp, " TrajOut         .false.              #\n");
  fprintf(fp, "---------------------------------------\n");
  fprintf(fp, " Output        my-output.out          #\n");
  fprintf(fp, "#######################################\n");
  fclose(fp);
  printf("\nJob done, the results are witten to directory: %s\n", workdir);
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
return;
}
