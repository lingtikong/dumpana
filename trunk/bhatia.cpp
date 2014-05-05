#include "driver.h"
#include "timer.h"
#include <set>
#include <cmath>

/*------------------------------------------------------------------------------
 * Method to compute the Bhatia-Thornton structure factor for a binary system.
 * Ref:
 *     AB Bhatia and DE Thornton, Phys Rev B 2(8):3004, 1970.
 *----------------------------------------------------------------------------*/
void Driver::bhatia_thornton()
{
  printf("\n"); for (int i = 0; i < 3; ++i) printf("=======");
  printf("  Bhatia-Thornton  Structure  Factor  ");
  for (int i = 0; i < 3; ++i) printf("======="); printf("\n");

  char str[MAXLINE], header[MAXLINE];

  // Check number of atomic types
  one = all[istr];
  if (one->ntype < 2){
    printf("The Bhatia-Thornton structure factor is defined for binary\n");
    printf("systems only, while there is not enough atomic type in your system.\n");
    
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }

  // selection commands
  char srcsel[MAXLINE];
  while (1){
    printf("\nPlease input the selection command for atoms, `h` for help [all]: ");
  
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(srcsel, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(srcsel,"all\n");
    
    // check the selection command on the first frame
    one->selection(srcsel); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }
  sprintf(header,"# Bhatia-Thornton structure factors for frames %d-%d in group: %s", istr+1, iend+1, srcsel);

  std::set<int> A, B, Both;
  A.clear(); B.clear(); Both.clear();
  //  Mapping of atom types
  printf("\nThe number of atom types are %d for the %d frame, while the\n", one->ntype, istr+1);
  printf("Bhatia-Thornton stucture factor works only for binary systems. Now,\n");
  // solvent atoms
  while (B.empty()){
    printf("please input your type list to be seen as solvent [2]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr = strtok(str," \n\t\r\f");
      while (ptr){
        int it = atoi(ptr);
        if (it > 0 && it <= one->ntype){
          B.insert(it); Both.insert(it);
        }

        ptr = strtok(NULL, " \n\t\r\f");
      }
    } else {
      B.insert(2); Both.insert(2);
    }
  }

  // solute atoms
  while (A.empty()){
    printf("Please input your type list to be seen as solute  [1]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr = strtok(str," \n\t\r\f");
      while (ptr){
        int it = atoi(ptr);
        if (it > 0 && it <= one->ntype){
          A.insert(it); Both.insert(it);
        }

        ptr = strtok(NULL, " \n\t\r\f");
      }
    } else {
      A.insert(1); Both.insert(1);
    }
  }
  // make sure no overlap
  if (Both.size() != A.size() + B.size()){
    printf("\nNo element should be seen as both solvent and solute!\n");
    
    for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
    return;
  }
  Both.clear();

  // bound for q
  double qmax = 6.;
  printf("\nNow please input your desired maximum q [%lg]: ", qmax);
  if (count_words(fgets(str,MAXLINE,stdin)) >= 1){
    char *ptr = strtok(str," \n\t\r\f");
    qmax = atof(ptr);

    if (qmax <= 0.){
      printf("\nERROR: q-max must be greater than 0!\n");

      for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
      return;
    }
  }

  int nq[3];
  // ask for the q-mesh size
  nq[0] = nq[1] = nq[2] = 41;

  printf("Now please input the # of q-points along each direction [%d %d %d]: ", nq[0], nq[1], nq[2]);
  if (count_words(fgets(str,MAXLINE,stdin)) >= 3){
    char *ptr = strtok(str," \n\t\r\f");
    for (int i = 0; i < 3; ++i){
      nq[i] = MAX(0, atoi(ptr)-1);

      ptr = strtok(NULL, " \n\t\r\f");
    }
  }
  double dq[3], q[3];
  for (int i = 0; i < 3; ++i) dq[i] = qmax/double(MAX(1, nq[i]));

  // timer
  Timer * timer = new Timer();
  printf("\nComputing, takes time, please be patient\n"); fflush(stdout);

  // Now to initialize the necessary data
  double ***SNN, ***SNC, ***SCC;
  memory->create(SNN, 2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "SNN");
  memory->create(SNC, 2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "SNC");
  memory->create(SCC, 2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "SCC");

  for (int ix = 0; ix < 2*nq[0]+1; ++ix)
  for (int iy = 0; iy < 2*nq[1]+1; ++iy)
  for (int iz = 0; iz <   nq[2]+1; ++iz) SNN[ix][iy][iz] = SNC[ix][iy][iz] = SCC[ix][iy][iz] = 0.;

  int nused = 0;
  const std::complex<double> I0 = std::complex<double>(0.,1.);
  
  std::complex<double> ***N1, ***N2, ***C, ***N;
  memory->create(N1,  2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "N1");
  memory->create(N2,  2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "N2");
  memory->create(C,   2*nq[0]+1, 2*nq[1]+1, nq[2]+1, "C");
  N = N1;

  // now to compute the Fourie transformation
  for (int img = istr; img <= iend; img += inc){ // loop over frames
    one = all[img];
    one->selection(srcsel);
    if (min_mem) one->FreeVoro();
    if (one->nsel < 1) continue;
  
    printf("  Now to process frame %d... ", img+1); fflush(stdout);

    for (int ix = 0; ix < 2*nq[0]+1; ++ix)
    for (int iy = 0; iy < 2*nq[1]+1; ++iy)
    for (int iz = 0; iz <   nq[2]+1; ++iz) N1[ix][iy][iz] = N2[ix][iy][iz] = std::complex<double>(0.,0.);

    // need cartesian coordinates
    one->dir2car();
    
    int na = 0, nb = 0;

      // loop over atoms
    for (int id = 1; id <= one->natom; ++id){
      if (one->atsel[id] == 0) continue;
      int it = one->attyp[id];

      if (A.find(it) != A.end()){
        ++na;

        // loop over q-mesh
        for (int qx = -nq[0]; qx <= nq[0]; ++qx)
        for (int qy = -nq[1]; qy <= nq[1]; ++qy)
        for (int qz =      0; qz <= nq[2]; ++qz){
          int ix = qx + nq[0];
          int iy = qy + nq[1];
          int iz = qz;
          q[0] = qx*dq[0];
          q[1] = qy*dq[1];
          q[2] = qz*dq[2];

          double qr = q[0]*one->atpos[id][0] + q[1]*one->atpos[id][1] + q[2]*one->atpos[id][2];

          N1[ix][iy][iz] += exp(I0*qr);
        }
      }

      if (B.find(it) != B.end()){
        ++nb;

        // loop over q-mesh
        for (int qx = -nq[0]; qx <= nq[0]; ++qx)
        for (int qy = -nq[1]; qy <= nq[1]; ++qy)
        for (int qz =      0; qz <= nq[2]; ++qz){
          int ix = qx + nq[0];
          int iy = qy + nq[1];
          int iz = qz;
          q[0] = qx*dq[0];
          q[1] = qy*dq[1];
          q[2] = qz*dq[2];

          double qr = q[0]*one->atpos[id][0] + q[1]*one->atpos[id][1] + q[2]*one->atpos[id][2];

          N2[ix][iy][iz] += exp(I0*qr);
        }
      }
    }   // end of loop over q-mesh
    ++nused;

    double nt = double(na+nb);
    double inv_n = 1. / nt;
    double c = double(na) * inv_n;

    N1[nq[0]][nq[1]][0] -= double(na);
    N2[nq[0]][nq[1]][0] -= double(nb);


    for (int ix = 0; ix < 2*nq[0]+1; ++ix)
    for (int iy = 0; iy < 2*nq[1]+1; ++iy)
    for (int iz = 0; iz <   nq[2]+1; ++iz) C[ix][iy][iz] = ( (1.-c)*N1[ix][iy][iz] - c*N2[ix][iy][iz] )* inv_n;

    for (int ix = 0; ix < 2*nq[0]+1; ++ix)
    for (int iy = 0; iy < 2*nq[1]+1; ++iy)
    for (int iz = 0; iz <   nq[2]+1; ++iz) N[ix][iy][iz] += N2[ix][iy][iz];

    for (int ix = 0; ix < 2*nq[0]+1; ++ix)
    for (int iy = 0; iy < 2*nq[1]+1; ++iy)
    for (int iz = 0; iz <   nq[2]+1; ++iz) SNN[ix][iy][iz] += real(conj(N[ix][iy][iz])*N[ix][iy][iz]) * inv_n;

    for (int ix = 0; ix < 2*nq[0]+1; ++ix)
    for (int iy = 0; iy < 2*nq[1]+1; ++iy)
    for (int iz = 0; iz <   nq[2]+1; ++iz) SCC[ix][iy][iz] += real(conj(C[ix][iy][iz])*C[ix][iy][iz]) * nt;

    for (int ix = 0; ix < 2*nq[0]+1; ++ix)
    for (int iy = 0; iy < 2*nq[1]+1; ++iy)
    for (int iz = 0; iz <   nq[2]+1; ++iz) SNC[ix][iy][iz] += real(conj(N[ix][iy][iz])*C[ix][iy][iz]);

    printf(" Done! Time used: %g seconds.\n", timer->sincelast());
  } // end loop over frames
  N = NULL;
  memory->destroy(N1);
  memory->destroy(N2);
  memory->destroy(C);

  A.clear(); B.clear();

  // get the average Bhatia-Thornton structure factor info
  double inv_n = 1./double(nused);
  for (int ix = 0; ix < 2*nq[0]+1; ++ix)
  for (int iy = 0; iy < 2*nq[1]+1; ++iy)
  for (int iz = 0; iz <   nq[2]+1; ++iz){
    SNN[ix][iy][iz] *= inv_n;
    SCC[ix][iy][iz] *= inv_n;
    SNC[ix][iy][iz] *= inv_n;
  }
  
  // time info
  timer->stop();
  printf("%d frames were used. Total CPU time cost: %g seconds.\n", nused, timer->cpu_time());
  delete timer;
  
  // output the Bhatia-Thornton structure factor for each q-point
  printf("\nIf you want to output the Bhatia-Thornton structure factor info for\n");
  printf("the whole q-mesh, input the prefix for the files now, enter to skip: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *ptr = strtok(str," \n\t\r\f");

    // SNN
    char fname[3][MAXLINE];
    FILE *fp[3];

    sprintf(fname[0], "%s_SNN.dat", ptr);
    sprintf(fname[1], "%s_SNC.dat", ptr);
    sprintf(fname[2], "%s_SCC.dat", ptr);
    for (int i = 0; i < 3; ++i){
      fp[i] = fopen(fname[i], "w");
      fprintf(fp[i], "%s", header);
      fprintf(fp[i], "# qx qy qz S(q)\n");
    }

    for (int qx = -nq[0]; qx <= nq[0]; ++qx)
    for (int qy = -nq[1]; qy <= nq[1]; ++qy){
      // negative qz
      for (int qz = -nq[2]; qz < 0; ++qz){
        int ix = -qx + nq[0];
        int iy = -qy + nq[1];
        int iz = -qz;
        q[0] = qx*dq[0];
        q[1] = qy*dq[1];
        q[2] = qz*dq[2];

        fprintf(fp[0], "%lg %lg %lg %lg\n", q[0], q[1], q[2], SNN[ix][iy][iz]);
        fprintf(fp[1], "%lg %lg %lg %lg\n", q[0], q[1], q[2], SNC[ix][iy][iz]);
        fprintf(fp[2], "%lg %lg %lg %lg\n", q[0], q[1], q[2], SCC[ix][iy][iz]);
      }
      // positive qz
      for (int qz = 0; qz <= nq[2]; ++qz){
        int ix = qx + nq[0];
        int iy = qy + nq[1];
        int iz = qz;
        q[0] = qx*dq[0];
        q[1] = qy*dq[1];
        q[2] = qz*dq[2];

        fprintf(fp[0], "%lg %lg %lg %lg\n", q[0], q[1], q[2], SNN[ix][iy][iz]);
        fprintf(fp[1], "%lg %lg %lg %lg\n", q[0], q[1], q[2], SNC[ix][iy][iz]);
        fprintf(fp[2], "%lg %lg %lg %lg\n", q[0], q[1], q[2], SCC[ix][iy][iz]);
      }
    }

    for (int i = 0; i < 3; ++i) fclose(fp[i]);
    printf("S info over the whole q-mesh are written to %s_SNN.dat, %s_SNC.dat, %s_SCC.dat, respectively.\n\n", ptr, ptr, ptr);
  }

  // output the radial BT result
  int nbin = 201;
  printf("Please input the # of bins to output S-q [%d]: ", nbin);
  if (count_words(fgets(str,MAXLINE, stdin)) > 0){
    char *ptr = strtok(str, " \n\t\r\f");
    nbin = MAX(1, atoi(ptr));
  }
  double inv_dq = double(MAX(1, nbin-1))/qmax;

  double **Sq;
  memory->create(Sq, 4, nbin, "Sq");
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < nbin; ++j) Sq[i][j] = 0.;
  
  for (int qx = -nq[0]; qx <= nq[0]; ++qx)
  for (int qy = -nq[1]; qy <= nq[1]; ++qy){
    // negative qz
    for (int qz = -nq[2]; qz < 0; ++qz){
      int ix = -qx + nq[0];
      int iy = -qy + nq[1];
      int iz = -qz;
      q[0] = qx*dq[0];
      q[1] = qy*dq[1];
      q[2] = qz*dq[2];

      double qs = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
      int ibin = int(qs * inv_dq + 0.5);
      if (ibin < nbin){
        Sq[0][ibin] += 1.;
        Sq[1][ibin] += SNN[ix][iy][iz];
        Sq[2][ibin] += SNC[ix][iy][iz];
        Sq[3][ibin] += SCC[ix][iy][iz];
      }
    }
    // positive qz
    for (int qz = 0; qz <= nq[2]; ++qz){
      int ix = qx + nq[0];
      int iy = qy + nq[1];
      int iz = qz;
      q[0] = qx*dq[0];
      q[1] = qy*dq[1];
      q[2] = qz*dq[2];

      double qs = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
      int ibin = int(qs * inv_dq + 0.5);
      if (ibin < nbin){
        Sq[0][ibin] += 1.;
        Sq[1][ibin] += SNN[ix][iy][iz];
        Sq[2][ibin] += SNC[ix][iy][iz];
        Sq[3][ibin] += SCC[ix][iy][iz];
      }
    }
  }
  for (int j = 0; j < nbin; ++j){
    if (Sq[0][j] > 0.){
      for (int i = 1; i < 4; ++i) Sq[i][j] /= Sq[0][j];
    }
  }

  printf("Please input the file to output S(q) [sq.dat]: ");
  fgets(str,MAXLINE, stdin);
  char *ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL) strcpy(str, "sq.dat");
  ptr = strtok(str, " \n\t\r\f");
  FILE *fp = fopen(ptr,"w");
  fprintf(fp,"%s", header);
  fprintf(fp,"# q SNN SNC SCC counts\n");
  dq[0] = 1./inv_dq;
  q[0]  = dq[0];
  for (int i = 1; i < nbin; ++i){
    fprintf(fp,"%lg %lg %g %lg %d\n", q[0], Sq[1][i], Sq[2][i], Sq[3][i], int(Sq[0][i]));

    q[0] += dq[0];
  }
  fclose(fp);

  memory->destroy(Sq);
  memory->destroy(SNN);
  memory->destroy(SNC);
  memory->destroy(SCC);
  printf("\n%d images were used in the evaluation of S(q), which is written to %s\n", nused, ptr);

  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");
return;
}

/*------------------------------------------------------------------------------*/
