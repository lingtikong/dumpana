#include "atom.h"
#include "time.h"
#include "random.h"
#include "global.h"

/*------------------------------------------------------------------------------
 * Constructor, to read one image from the atom style dump file of lammps
 *------------------------------------------------------------------------------
 * Parameters:
 *  fp       : (in) File pointer, of dump file
 *  dumpfile : (in) dump file
 *  flag     : (in) 1st bit, least memory or not
 *----------------------------------------------------------------------------*/
DumpAtom::DumpAtom(FILE *fp, const char *dumpfile, const int flag)
{
  iframe = natom = ntype = tstep = 0;
  initialized = triclinic = 0;
  xy = xz = yz = 0.;
  least_memory = cartesian = 0;
  type2radius = NULL;
  flag_smix = 0; smix = 0.;

  int flag_wrap = 0;
  if (flag & 1) least_memory = 1;
  if (flag & 2) flag_wrap = 1;

  realcmd = NULL;
  attyp = atsel = numtype = env = NULL;
  atpos = x = s = NULL;

  atprop = NULL;
  prop_label.clear();

  wted = 0;
  voro.clear();
  neilist = image = NULL;
  prop = volume = NULL;
  MaxNei = 16;
  vmins[0] = vmins[1] = vmins[2] = 0.;

  fname = new char [strlen(dumpfile)+1];
  strcpy(fname, dumpfile);

  // time step info
  char str[MAXLINE];
  if (fgets(str,MAXLINE, fp) == NULL) return;
  fgets(str,MAXLINE, fp);
  tstep = atoi(strtok(str, " \n\t\r\f"));

  // # of atoms info
  fgets(str,MAXLINE, fp);
  fgets(str,MAXLINE, fp);
  natom = atoi(strtok(str, " \n\t\r\f"));
  if (natom < 1) return;

  // box bounds info
  fgets(str,MAXLINE, fp);
  fgets(str,MAXLINE, fp);
  int n = count_words(str);
  xlo = atof(strtok(str, " \n\t\r\f"));
  xhi = atof(strtok(NULL," \n\t\r\f"));
  if (n == 3) xy = atof(strtok(NULL," \n\t\r\f"));

  fgets(str,MAXLINE, fp);
  n = count_words(str);
  ylo = atof(strtok(str, " \n\t\r\f"));
  yhi = atof(strtok(NULL," \n\t\r\f"));
  if (n == 3) xz = atof(strtok(NULL," \n\t\r\f"));

  fgets(str,MAXLINE, fp);
  n = count_words(str);
  zlo = atof(strtok(str, " \n\t\r\f"));
  zhi = atof(strtok(NULL," \n\t\r\f"));
  if (n == 3) yz = atof(strtok(NULL," \n\t\r\f"));
  
  if (xy*xy+xz*xz+yz*yz > ZERO){
    triclinic = 1;

    xlo -= MIN(MIN(0., xy), MIN(xz, xy+xz));
    xhi -= MAX(MAX(0., xy), MAX(xz, xy+xz));
    ylo -= MIN(0., yz);
    yhi -= MAX(0., yz);
  }

  // fields info
  int dcols[9], fcord = 7;
  for (int i = 0; i <= 5; ++i) dcols[i] = i;
  dcols[6] = dcols[7] = dcols[8] = 0;

  map<int,int> col2pid; col2pid.clear();

  fgets(str,MAXLINE, fp);
  char *ptr = strtok(str, " \n\t\r\f");
  for (int i = 0; i < 2; ++i) ptr = strtok(NULL," \n\t\r\f");
  int ic = 1;
  while (ptr){
    int taken = 0;
    if (strcmp(ptr, "id") == 0)  {dcols[1] = ic; taken = 1;}
    if (strcmp(ptr, "type") == 0){dcols[2] = ic; taken = 1;}

    if (strcmp(ptr, "xs") == 0 || strcmp(ptr, "xsu") == 0){dcols[3] = ic; fcord |= 1; taken = 1;}
    if (strcmp(ptr, "ys") == 0 || strcmp(ptr, "ysu") == 0){dcols[4] = ic; fcord |= 2; taken = 1;}
    if (strcmp(ptr, "zs") == 0 || strcmp(ptr, "zsu") == 0){dcols[5] = ic; fcord |= 4; taken = 1;}

    if (strcmp(ptr, "x") == 0 || strcmp(ptr, "xu") == 0){   dcols[3] = ic; fcord &= 6; taken = 1;}
    if (strcmp(ptr, "y") == 0 || strcmp(ptr, "yu") == 0){   dcols[4] = ic; fcord &= 5; taken = 1;}
    if (strcmp(ptr, "z") == 0 || strcmp(ptr, "zu") == 0){   dcols[5] = ic; fcord &= 3; taken = 1;}

    if (strcmp(ptr, "ix") == 0){dcols[6] = ic; taken = 1;}
    if (strcmp(ptr, "iy") == 0){dcols[7] = ic; taken = 1;}
    if (strcmp(ptr, "iz") == 0){dcols[8] = ic; taken = 1;}

    if (taken == 0){
      prop_label.push_back(ptr);
      col2pid[ic] = prop_label.size() - 1;
    }

    ++ic;
    ptr = strtok(NULL," \n\t\r\f");
  }
  fcord &= 7;
  if (fcord != 7 && fcord != 0) return;
  int flag_img = 0;
  if (dcols[6] > 0 && dcols[7] > 0 && dcols[8] > 0){
    memory->create(image, natom+1, 3, "image");
    flag_img = 1;
  }

  memory = new Memory();
  memory->create(attyp, natom+1, "attyp");
  memory->create(atsel, natom+1, "atsel");

  for (int i = 1; i <= natom; ++i) atsel[i] = 1; atsel[0] = 0;
  realcmd = new char [MAXLINE]; strcpy(realcmd, "all\n");
  nsel = natom;

  // always assumes fractional coordinate
  memory->create(s, natom+1, 3, "s");
  atpos = s;

  // property related
  int nprop = prop_label.size();
  if (nprop >= 1) memory->create(atprop, natom+1, nprop, "atprop");
  
  // read coordinate
  int id, ip, ix, iy, iz;
  double xp, yp, zp;
  for (int i = 0; i < natom; ++i){
    fgets(str,MAXLINE, fp);

    int ic = 1, frd = 0, fimg = 0;
    ptr = strtok(str, " \n\t\r\f");
    vector<string> prop; prop.clear();
    while (ptr){
      if (ic == dcols[1]){ id = atoi(ptr); frd |=  1; }
      if (ic == dcols[2]){ ip = atoi(ptr); frd |=  2; }
      if (ic == dcols[3]){ xp = atof(ptr); frd |=  4; }
      if (ic == dcols[4]){ yp = atof(ptr); frd |=  8; }
      if (ic == dcols[5]){ zp = atof(ptr); frd |= 16; }
      if (flag_img){
        if (ic == dcols[6]){ ix = atoi(ptr); fimg |= 1; }
        if (ic == dcols[7]){ iy = atoi(ptr); fimg |= 2; }
        if (ic == dcols[8]){ iz = atoi(ptr); fimg |= 4; }
      }
      if (col2pid.count(ic) > 0) prop.push_back(ptr);

      ptr = strtok(NULL," \n\t\r\f"); ++ic;
    }
    if (frd == 31){
      attyp[id] = ip; ntype = MAX(ip, ntype);
      atpos[id][0] = xp;
      atpos[id][1] = yp;
      atpos[id][2] = zp;
      if (fimg == 7){
        image[id][0] = ix;
        image[id][1] = iy;
        image[id][2] = iz;
      }
      for (int ii = 0; ii < nprop; ++ii)
        atprop[id][ii] = atof(prop[ii].c_str());

    } else { return; } // insufficient info, return
  }

  lx = box[0] = axis[0][0] = xhi - xlo;
  ly = box[1] = axis[1][1] = yhi - ylo;
  lz = box[2] = axis[2][2] = zhi - zlo;
  hx = 0.5*lx;  hy = 0.5*ly;  hz = 0.5*lz;

  vol = lx*ly*lz;
  box[3] = axis[1][0] = xy;
  box[4] = axis[2][0] = xz;
  box[5] = axis[2][1] = yz;
  axis[0][1] = axis[0][2] = axis[1][2] = 0.;
  for (int idim = 0; idim < 3; ++idim) hbox[idim] = 0.5*box[idim];

  h_inv[0] = 1./box[0];
  h_inv[1] = 1./box[1];
  h_inv[2] = 1./box[2];
  h_inv[3] = -box[3] / (box[1]*box[2]);
  h_inv[4] = (box[3]*box[5] - box[1]*box[4]) / (box[0]*box[1]*box[2]);
  h_inv[5] = -box[5] / (box[0]*box[1]);

  if ((fcord&7) == 0){ // in case cartesian coordinate read, convert to fractional
    if (least_memory){
      double x0[3];
      for (int id = 1; id <= natom; ++id){
        x0[0] = s[id][0] - xlo;
        x0[1] = s[id][1] - ylo;
        x0[2] = s[id][2] - zlo;
        s[id][0] = h_inv[0]*x0[0] + h_inv[5]*x0[1] + h_inv[4]*x0[2];
        s[id][1] = h_inv[1]*x0[1] + h_inv[3]*x0[2];
        s[id][2] = h_inv[2]*x0[2];
      }

    } else {
      memory->create(x, natom+1, 3, "x");
      for (int id = 1; id <= natom; ++id)
      for (int idim = 0; idim < 3; ++idim) x[id][idim] = s[id][idim];
   
      double x0[3];
      for (int id = 1; id <= natom; ++id){
        x0[0] = x[id][0] - xlo;
        x0[1] = x[id][1] - ylo;
        x0[2] = x[id][2] - zlo;
        s[id][0] = h_inv[0]*x0[0] + h_inv[5]*x0[1] + h_inv[4]*x0[2];
        s[id][1] = h_inv[1]*x0[1] + h_inv[3]*x0[2];
        s[id][2] = h_inv[2]*x0[2];
      }
    }
  }

  // wrap all atoms into the central box
  if (flag_wrap){
    if (image){
      for (int id = 1; id <= natom; ++id){
        while (s[id][0] >= 1.){ s[id][0] -= 1.; image[id][0]++;}
        while (s[id][0] < 0. ){ s[id][0] += 1.; image[id][0]--;}
        while (s[id][1] >= 1.){ s[id][1] -= 1.; image[id][1]++;}
        while (s[id][1] < 0. ){ s[id][1] += 1.; image[id][1]--;}
        while (s[id][2] >= 1.){ s[id][2] -= 1.; image[id][2]++;}
        while (s[id][2] < 0. ){ s[id][2] += 1.; image[id][2]--;}
      }
    } else {
      for (int id = 1; id <= natom; ++id){
        while (s[id][0] >= 1.) s[id][0] -= 1.;
        while (s[id][0] <  0.) s[id][0] += 1.;
        while (s[id][1] >= 1.) s[id][1] -= 1.;
        while (s[id][1] <  0.) s[id][1] += 1.;
        while (s[id][2] >= 1.) s[id][2] -= 1.;
        while (s[id][2] <  0.) s[id][2] += 1.;
      }
    }

    if (least_memory == 0) dir2car();
  }

  memory->create(numtype,ntype+1,"numtype");
  for (int i = 0; i <=ntype; ++i) numtype[i] = 0;
  for (int i = 1; i <=natom; ++i) ++numtype[attyp[i]];

  col2pid.clear();
  initialized = 1;

return;
}

/*------------------------------------------------------------------------------
 * Deconstructor, to free allocated memory
 *----------------------------------------------------------------------------*/
DumpAtom::~DumpAtom()
{
  if (fname)   delete []fname;
  if (realcmd) delete []realcmd;
  memory->destroy(attyp);
  memory->destroy(atsel);
  memory->destroy(atprop);
  memory->destroy(image);
  memory->destroy(numtype);

  atpos = NULL;
  memory->destroy(x);
  memory->destroy(s);

  memory->destroy(neilist);
  memory->destroy(volume);

  voro.clear();
  prop_label.clear();
  type2radius = NULL;

  memory->destroy(env);
  memory->destroy(prop);
  delete memory;
}


/*------------------------------------------------------------------------------
 * Method to convert fractional coordinates into cartesian
 *----------------------------------------------------------------------------*/
void DumpAtom::dir2car()
{
  if (cartesian) return;

  if (least_memory == 0 && x){
    atpos = x; cartesian = 1;
    return;
  }

  if (least_memory){
    double s0[3];
    if (triclinic){
      for (int id = 1; id <= natom; ++id){
        for (int idim = 0; idim < 3; ++idim) s0[idim] = s[id][idim];
        s[id][0] = s0[0]*lx + s0[1]*xy + s0[2]*xz + xlo;
        s[id][1] = s0[1]*ly + s0[2]*yz + ylo;
        s[id][2] = s0[2]*lz + zlo;
      }
    } else {
      for (int id = 1; id <= natom; ++id){
        for (int idim = 0; idim < 3; ++idim) s0[idim] = s[id][idim];
        s[id][0] = s0[0]*lx + xlo;
        s[id][1] = s0[1]*ly + ylo;
        s[id][2] = s0[2]*lz + zlo;
      }
    }
    atpos = s;

  } else {

    if (x == NULL){
      memory->create(x, natom+1, 3,"x");
      if (triclinic){
        for (int id = 1; id <= natom; ++id){
          x[id][0] = s[id][0]*lx + s[id][1]*xy + s[id][2]*xz + xlo;
          x[id][1] = s[id][1]*ly + s[id][2]*yz + ylo;
          x[id][2] = s[id][2]*lz + zlo;
        }
      } else {
        for (int id = 1; id <= natom; ++id){
          x[id][0] = s[id][0]*lx + xlo;
          x[id][1] = s[id][1]*ly + ylo;
          x[id][2] = s[id][2]*lz + zlo;
        }
      }
    }
    atpos = x;
  }
  cartesian = 1;
return;
}

/*------------------------------------------------------------------------------
 * Method to get the pbc distance for a coordinator vector
 *----------------------------------------------------------------------------*/
double DumpAtom::get_dist2(double xij, double yij, double zij)
{
  if (cartesian) return xij*xij + yij*yij + zij*zij;
  double rij[3];
  if (triclinic){
     rij[0] = xij*lx + yij*xy + zij*xz;
     rij[1] = yij*ly + zij*yz;
     rij[2] = zij*lz;

  } else {
     rij[0] = xij*lx;
     rij[1] = yij*ly;
     rij[2] = zij*lz;
  }
  return rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinates into fractional
 *----------------------------------------------------------------------------*/
void DumpAtom::car2dir()
{
  if (cartesian == 0) return;

  double x0[3];
  if (least_memory){
    if (triclinic){
      for (int id = 1; id <= natom; ++id){
        x0[0] = s[id][0] - xlo;
        x0[1] = s[id][1] - ylo;
        x0[2] = s[id][2] - zlo;
        s[id][0] = h_inv[0]*x0[0] + h_inv[5]*x0[1] + h_inv[4]*x0[2];
        s[id][1] = h_inv[1]*x0[1] + h_inv[3]*x0[2];
        s[id][2] = h_inv[2]*x0[2];
      }

    } else {

      for (int id = 1; id <= natom; ++id){
        s[id][0] = h_inv[0]*(s[id][0] - xlo);
        s[id][1] = h_inv[1]*(s[id][1] - ylo);
        s[id][2] = h_inv[2]*(s[id][2] - zlo);
      }
    }
  }

  atpos = s;
  cartesian = 0;
return;
}

/*------------------------------------------------------------------------------
 * Method to select some atoms, not implemented yet
 *----------------------------------------------------------------------------*/
void DumpAtom::selection(const char *line)
{
  int n = strlen(line) + 1;
  char *selcmd;
  memory->create(selcmd, n, "selcmd");
  strcpy(selcmd,line);

  char *key, *oper, *ptr;
  int ilow, ihigh;
  double rlow, rhigh;
  int logand = 1;

  char onecmd[MAXLINE];
  strcpy(realcmd," ");

  // by default, all are selected
  for (int id = 1; id <= natom; ++id) atsel[id] = 1;
  key = strtok(selcmd, " \n\t\r\f");

  while (key != NULL){
    strcpy(onecmd,key);

    if (strcmp(key, "&") == 0){
      logand = 1; strcat(onecmd," ");

    } else if (strcmp(key, "|") == 0){
      logand = 0; strcat(onecmd," ");

    } else if (strcmp(key, "type") == 0){ // selection by type
      oper = strtok(NULL, " \n\t\r\f");
      if (oper == NULL) break;
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      ilow = atoi(ptr);

      strcat(onecmd," "); strcat(onecmd,oper);
      strcat(onecmd," "); strcat(onecmd,ptr);

      if (strcmp(oper, "=") == 0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id] != ilow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id] == ilow) atsel[id] = 1;
        }

      } else if (strcmp(oper,">") == 0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id] <= ilow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id] > ilow) atsel[id] = 1;
        }

      } else if (strcmp(oper,">=") == 0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id] < ilow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id] >= ilow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<") == 0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id] >= ilow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id] < ilow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<=") == 0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id] > ilow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id] <= ilow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id]<ilow || attyp[id]>ihigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id]>=ilow && attyp[id]<=ihigh) atsel[id] = 1;
        }

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id]>ilow && attyp[id]<ihigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id]<=ilow || attyp[id]>=ihigh) atsel[id] = 1;
        }

      } else if (strcmp(oper,"%") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);
        if (logand){ for (int id = 1; id <= natom; ++id) if (attyp[id]%ilow != ihigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (attyp[id]%ilow == ihigh) atsel[id] = 1;
        }

      } else break;

    } else if (strcmp(key,"x")==0 || strcmp(key,"y")==0 || strcmp(key,"z")==0 ){ // selection by fractional position
      car2dir();

      int dir = key[0]-'x';
      oper = strtok(NULL, " \n\t\r\f");
      if (oper == NULL) break;
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      rlow = atof(ptr);

      strcat(onecmd," "); strcat(onecmd,oper);
      strcat(onecmd," "); strcat(onecmd,ptr);

      if (strcmp(oper,">")==0 || strcmp(oper,">=")==0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir] < rlow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir] >= rlow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<")==0 || strcmp(oper,"<=")==0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir] > rlow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir] <= rlow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir]<rlow || atpos[id][dir]>rhigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir]>=rlow && atpos[id][dir]<=rhigh) atsel[id] = 1;
        }

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir]>rlow && atpos[id][dir]<rhigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir]<=rlow || atpos[id][dir]>=rhigh) atsel[id] = 1;
        }

      } else break;

    } else if (strcmp(key,"X")==0 || strcmp(key,"Y")==0 || strcmp(key,"Z")==0 ){ // selection by cartesian position
      dir2car();

      int dir = key[0]-'X';
      oper = strtok(NULL, " \n\t\r\f");
      if (oper == NULL) break;
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      rlow = atof(ptr);

      strcat(onecmd," "); strcat(onecmd,oper);
      strcat(onecmd," "); strcat(onecmd,ptr);

      if (strcmp(oper,">")==0 || strcmp(oper,">=")==0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir] < rlow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir] >= rlow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<")==0 || strcmp(oper,"<=")==0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir] > rlow) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir] <= rlow) atsel[id] = 1;
        }

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir]<rlow || atpos[id][dir]>rhigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir]>=rlow && atpos[id][dir]<=rhigh) atsel[id] = 1;
        }

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (logand){ for (int id = 1; id <= natom; ++id) if (atpos[id][dir]>rlow && atpos[id][dir]<rhigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (atpos[id][dir]<=rlow || atpos[id][dir]>=rhigh) atsel[id] = 1;
        }

      } else break;

    } else if (strcmp(key,"id") == 0) { // selection by atomic id
      oper = strtok(NULL, " \n\t\r\f");
      if (oper == NULL) break;
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      ilow = atoi(ptr);

      strcat(onecmd," "); strcat(onecmd,oper);
      strcat(onecmd," "); strcat(onecmd,ptr);

      if (strcmp(oper,"=") == 0){
        if (logand){ for (int id = 1; id <= natom; ++id) if (id != ilow) atsel[id] = 0;
        } else atsel[ilow] = 1;

      } else if (strcmp(oper,">") == 0){
        if (logand) for (int id = 1; id <= MIN(ilow,natom); ++id) atsel[id] = 0;
        else for (int id = ilow+1; id <= natom; ++id) atsel[id] = 1;

      } else if (strcmp(oper,">=") == 0){
        if (logand) for (int id = 1; id < MIN(ilow,natom); ++id) atsel[id] = 0;
        else for (int id = ilow; id <= natom; ++id) atsel[id] = 1;

      } else if (strcmp(oper,"<") == 0){
        if (logand) for (int id = ilow; id <= natom; ++id) atsel[id] = 0;
        else for (int id = 1; id < MIN(ilow,natom); ++id) atsel[id] = 1;

      } else if (strcmp(oper,"<=") == 0){
        if (logand) for (int id = ilow+1; id <= natom; ++id) atsel[id] = 0;
        else for (int id = 1; id <= MIN(ilow,natom); ++id) atsel[id] = 1;

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);

        if (logand){
          for (int id = 1; id < MIN(ilow,natom); ++id) atsel[id] = 0;
          for (int id = ihigh+1; id <= natom; ++id) atsel[id] = 0;
        } else for (int id = ilow; id <= MIN(ihigh,natom); ++id) atsel[id] = 1;

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);

        if (logand) for (int id = ilow; id <= MIN(ihigh,natom); ++id) atsel[id] = 0;
        else {
          for (int id = 1; id < MIN(ilow,natom); ++id) atsel[id] = 1;
          for (int id = ihigh+1; id <= natom; ++id) atsel[id] = 1;
        }

      } else if (strcmp(oper,"%") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);

        if (logand){ for (int id = 1; id <= natom; ++id) if (id%ilow != ihigh) atsel[id] = 0;
        } else { for (int id = 1; id <= natom; ++id) if (id%ilow == ihigh) atsel[id] = 1;
        }

      } else break;

    } else if (strcmp(key,"ran") == 0){ // random selection from current selection
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      ilow = atoi(ptr);

      strcat(onecmd," "); strcat(onecmd,ptr);

      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      int seed = atoi(ptr);
      if (seed < 1) seed = time(NULL)%86400+1;

      sprintf(onecmd,"%s %d", onecmd, seed);

      seed += iframe;
      RanPark * random = new RanPark(seed);

      nsel = 0;
      for (int id = 1; id <= natom; ++id) nsel += atsel[id];
      int ndel = nsel - ilow;
      while (ndel > 0){
        int id = MIN(random->uniform()*(natom+1), natom);
        if (atsel[id] == 1){
          atsel[id] = 0;
          ndel--;
        }
      }
      delete random;

    } else if (strcmp(key,"voro") == 0){ // selected by Voronoi indices
      set<string> voroset; string vindex;
      double mins[3];
      for (int i = 0; i < 3; ++i){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        mins[i] = atof(ptr);
      }
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      strcat(onecmd," "); strcat(onecmd,ptr);
      int flag_ex = 0;
      int nv = atoi(ptr);
      if (nv < 0){
        flag_ex = 1;
        nv = -nv;
      }

      voroset.clear();
      for (int i = 0; i < nv; ++i){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        vindex.assign(ptr); voroset.insert(vindex);
      }
      ComputeVoro(mins);

      // make the selection
      std::set<int> vatset; vatset.clear();
      if (voroset.size() > 0){
        for (int id = 1; id <= natom; ++id){
          vindex = voro[id];
          if (voroset.count(vindex) > 0 && flag_ex == 0) vatset.insert(id);
          else if (voroset.count(vindex) < 1 && flag_ex == 1) vatset.insert(id);
        }
      }

      if (logand){
        for (int id = 1; id <= natom; ++id){
          if (vatset.count(id) < 1) atsel[id] = 0;
        }

      } else {
        for (int id = 1; id <= natom; ++id){
          if (vatset.count(id) > 0) atsel[id] = 1;
        }
      }
      vatset.clear();
      voroset.clear();

    } else if (strcmp(key,"VORO") == 0){ // selected by Voronoi indices and their neighbors
      set<string> voroset; string vindex;
      double mins[3];
      for (int i = 0; i < 3; ++i){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        mins[i] = atof(ptr);
      }
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      strcat(onecmd," "); strcat(onecmd,ptr);
      int flag_ex = 0;
      int nv = atoi(ptr);
      if (nv < 0){
        flag_ex = 1;
        nv = -nv;
      }

      voroset.clear();
      for (int i = 0; i < nv; ++i){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        vindex.assign(ptr); voroset.insert(vindex);
      }
      ComputeVoro(mins);

      // make the selection
      std::set<int> vatset; vatset.clear();
      if (voroset.size() > 0){
        for (int id = 1; id <= natom; ++id){
          vindex = voro[id];
          if (voroset.count(vindex) > 0 && flag_ex == 0){
            vatset.insert(id);
            for (int jj = 1; jj <= neilist[0][id]; ++jj){
              int jd = neilist[jj][id];
              vatset.insert(jd);
            }

          } else if (voroset.count(vindex) < 1 && flag_ex == 1){
            vatset.insert(id);
            for (int jj = 1; jj <= neilist[0][id]; ++jj){
              int jd = neilist[jj][id];
              vatset.insert(jd);
            }
          }
        }
      }

      if (logand){
        for (int id = 1; id <= natom; ++id){
          if (vatset.count(id) < 1) atsel[id] = 0;
        }

      } else {
        for (int id = 1; id <= natom; ++id){
          if (vatset.count(id) > 0) atsel[id] = 1;
        }
      }
      vatset.clear();
      voroset.clear();

    } else if (strcmp(key,"vol")==0){ // selection by Voronoi volume

      oper = strtok(NULL, " \n\t\r\f");
      if (oper == NULL) break;
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      rlow = atof(ptr);

      strcat(onecmd," "); strcat(onecmd,oper);
      strcat(onecmd," "); strcat(onecmd,ptr);

      if (strcmp(oper,">")==0 || strcmp(oper,">=")==0){
        if (volume){
          if (logand){ for (int id = 1; id <= natom; ++id) if (volume[id] < rlow) atsel[id] = 0;
          } else { for (int id = 1; id <= natom; ++id) if (volume[id] >= rlow) atsel[id] = 1;
          }
        }

      } else if (strcmp(oper,"<")==0 || strcmp(oper,"<=")==0){
        if (volume){
          if (logand){ for (int id = 1; id <= natom; ++id) if (volume[id] > rlow) atsel[id] = 0;
          } else { for (int id = 1; id <= natom; ++id) if (volume[id] <= rlow) atsel[id] = 1;
          }
        }

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (volume){
          if (logand){ for (int id = 1; id <= natom; ++id) if (volume[id]<rlow || volume[id]>rhigh) atsel[id] = 0;
          } else { for (int id = 1; id <= natom; ++id) if (volume[id]>=rlow && volume[id]<=rhigh) atsel[id] = 1;
          }
        }

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (volume){
          if (logand){ for (int id = 1; id <= natom; ++id) if (volume[id]>rlow && volume[id]<rhigh) atsel[id] = 0;
          } else { for (int id = 1; id <= natom; ++id) if (volume[id]<=rlow || volume[id]>=rhigh) atsel[id] = 1;
          }
        }

      } else break;

    } else if (strcmp(key,"R")==0 || strcmp(key, "r")==0){ // selection by region
      double inside = 1.;                     // R, inside region
      if (strcmp(key, "r")==0) inside = -1.;  // r, outside region

      int dir = 3;
      double c0[4], radius;
      oper = strtok(NULL, " \n\t\r\f");
      if (strcmp(oper, "C") == 0 || strcmp(oper, "c") == 0){ // Cylindrical: dir c1 c2 Radius
        strcat(onecmd," "); strcat(onecmd,oper);

        ptr = strtok(NULL, " \n\t\r\f");
        if (strcmp(ptr, "x") == 0 || strcmp(ptr, "X") == 0){
          dir = 0; c0[0] = 0.;

        } else if (strcmp(ptr, "y") == 0 || strcmp(ptr, "Y") == 0){
          dir = 1;
        } else {
          dir = 2;
        }

      } else {                                               // Spherical: cx cy cz Radius
        strcat(onecmd," "); strcat(onecmd,"S");
        ptr = strtok(NULL, " \n\t\r\f");
        c0[0] = atof(ptr);
      }
      strcat(onecmd," "); strcat(onecmd,ptr);

      for (int ii = 1; ii < 4; ++ii){
        ptr = strtok(NULL, " \n\t\r\f");
        c0[ii]  = atof(ptr);
        strcat(onecmd," "); strcat(onecmd,ptr);
      }
      radius = c0[3];
      if (dir == 1) c0[0] = c0[1];
      if (dir == 2){ c0[0] = c0[1]; c0[1] = c0[2];}

      // need cartesian coordinate system
      dir2car();
      double dr[4], r2 = radius * radius;
      for (int id = 1; id <= natom; ++id){
        for (int idim = 0; idim < 3; ++idim) dr[idim] = atpos[id][idim] - c0[idim];
        ApplyPBC(dr[0], dr[1], dr[2]);
        dr[dir] = 0.;
        double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        double side = (r2 - dr2)*inside;
        if (logand  && side < 0.) atsel[id] = 0;
        else if (!logand && side > 0.) atsel[id] = 1;
      }

    } else if (strcmp(key,"all") == 0){ // select all; it will discard all previous selections
      strcpy(realcmd,""); strcpy(onecmd,key);
      for (int id = 1; id <= natom; ++id) atsel[id] = 1;

    } else if (strcmp(key,"!") == 0 || strcmp(key,"~") == 0){  // reverse the previous selection
      strcpy(onecmd,key);
      for (int id = 1; id <= natom; ++id) atsel[id] = 1 - atsel[id];

    } else break;

    strcat(realcmd," "); strcat(realcmd,onecmd);
    key = strtok(NULL, " \n\t\r\f");
  }

  nsel = 0;
  for (int id = 1; id <= natom; ++id) nsel += atsel[id];
  memory->destroy(selcmd);

return;
}

/*------------------------------------------------------------------------------
 * Method to print out selection info
 *----------------------------------------------------------------------------*/
void DumpAtom::SelInfo()
{
  printf("\nThe realized selection command is: %s,\n", realcmd);
  printf("and %d of %d atoms are selected for frame %d.\n", nsel, natom, iframe);

return;
}

/*------------------------------------------------------------------------------
 * Help info related to the selection function
 *----------------------------------------------------------------------------*/
void DumpAtom::SelHelp()
{
  printf("\n"); for (int i = 0; i < 20; ++i) printf("----");
  printf("\nThe grammar for the selection command is:\n\n");
  printf("  key op values [& key2 op2 values2 [| key3 op3 values3]]\n");
  printf("\nwhere `key` is either `type`, `x`,`X', `y`, 'Y', `z`, 'Z', 'vol',\n");
  printf("or `id`. Lower case indicates fractional while upper Cartesian.\n");
  printf("It can also be `all`, which takes no argument and selects all atoms;\n");
  printf("or `!`, which takes no argument and reverses the previous selection;\n");
  printf("or `ran num seed`, which takes no other argument and selects `num` atoms\n");
  printf("from the current selection randomly; `seed` is the seed for the uniform\n");
  printf("random number generator, if a non-positive number is provided, it will be\n");
  printf("set automatically based on current time. The logical operation before `ran`\n");
  printf("is always assumed to be AND, even when `|' is specified;\n");
  printf("or `voro MinSurf NMinNei MinEdge NVoroIndex VoroIndices', which will\n");
  printf("select atoms with certain Voronoi indices. `MinSurf' and `MinEdge' defines\n");
  printf("the thresholds for surface and edge, respectively; they can be zero; `NMinNei'\n");
  printf("defines the Minimum # of neighbors wanted, zero means no limitations; `NVoroIndex'\n");
  printf("defines the # of Voronoi indices that will follow, it can be zero; in case\n");
  printf("`NVoroIndex` is negative, atoms whose Voronoi index are not in the following\n");
  printf("will be selected instead. The following `NVoroIndex' arguments should be\n");
  printf("the Voronoi indices, for example: 0,6,0,8.\n");
  printf("If the key is `VORO` instead, the central atom and its neighbors will also\n");
  printf("be included or excluded, depending on if `NVoroIndex' is positive or negative.\n");
  printf("\n`op` is either `=`, `>`, `>=`, `<`, `<=`, `<>`, `><`, or `%%`,\n");
  printf("while neithor `=` nor `%%` is available for key=`x-z, X-Z`.\n");
  printf("\n`values` can be one or two numbers, depend on the type of `op`:\n");
  printf("for `<>`, `><`, and `%%` two numbers are need, otherwise one.\n");
  printf("`key <> num1 num2` selects `num1 <= kye <= num2`, `key >< num1 num2`\n");
  printf("selects `key <= num1 or key >= num2`, while `key %% num1 num2` selects\n");
  printf("atoms satisify `key%%num1 == num2`.\n");
  printf("\n`key` could also be `R` or `r`, corresponding to space within (R) or beyond (r)\n");
  printf("a cylindrical (when `op` is `C`) or a spherical (when `op` is `S`) region.\n");
  printf("In this case, the grammar is:\n");
  printf("  `R/r C dir c1 c1 radius`  or `R/r S cx cy cz radius`\nwhere radius is always in length unit.\n");
  printf("\nMultiple `key op values` could be combined together, by either `&`\n");
  printf("or `|`, which means logical `and` or `or`, respectively. In this case\n");
  printf("the selections take effects sequentially.\n");
  printf("Selection according to `vol` is only enabled if Voronoi has been done before hand.\n");
  printf("\nFor example, `type = 1 & x <> 0.1 0.5 & id >< 200 800` will select atoms\n");
  printf("of type 1 within 0.1 < x < 0.5 (fractional) and has id <= 200 or id >= 800;\n");
  printf("`type = 1 | type = 2` will select atoms of type 1 or 2; while `type = 1 & type = 2`\n");
  printf("will select nothing. `type = 1 & ran 100 0` will randomly select 100 atoms\n");
  printf("from all of type 1. `X <> 0 10 & voro 0 0 0 1 0,0,12,0' will select all atoms\n");
  printf("that have a Voronoi index of 0,0,12,0 within [0,10] along the x direction.\n");
  for (int i = 0; i < 20; ++i) printf("----"); printf("\n\n");
    
return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int DumpAtom::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy, n, "copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) ++n;

  memory->sfree(copy);
  return n;
}

/*------------------------------------------------------------------------------
 * Private method to compute the direct Voronoi index, Voronoi neighbor list info
 * mins     (in)  : surf_min  Min#Nei edge_min; For either surf_min or edge_min,
 *                : if a positive number is provided, it is the absolute surface
 *                : area or edge length threshold; otherwise, it is the minimum
 *                : ratio over total surface area or circumference.
 * fp       (in)  : file pointer to write full voro info; if NULL, write nothing
 * fpsurf   (in)  : file pointer to write surf ratio info; if NULL, write nothing
 * fpedge   (in)  : file pointer to write edge ratio info; if NULL, write nothing
 *----------------------------------------------------------------------------*/
void DumpAtom::Direct_Voro(double *mins, FILE *fp, FILE *fpsurf, FILE *fpedge)
{
  double diff = 0.;
  for (int i = 0; i < 3; ++i) diff += (mins[i] - vmins[i])*(mins[i] - vmins[i]);
  if (diff <= ZERO && voro.size() == natom && wted == 0) return;

  wted = 0;
  int flag_surf_ratio = 0;
  if (mins[0] < 0.) flag_surf_ratio = 1;

  for (int i = 0; i < 3; ++i) vmins[i] = mins[i];

  double surf_min = fabs(vmins[0]);
  int nminnei = abs(int(vmins[1]));
  double edge_min = vmins[2];

  voro.clear();
  if (neilist) memory->destroy(neilist);
  if (volume)  memory->destroy(volume);
  memory->create(neilist, MaxNei+1, natom+1, "neilist");
  memory->create(volume,  natom+1, "volume");
  
  // need cartesian coordinates
  dir2car();

  // compute optimal size for container, and then contrust it
  double l = pow(double(natom)/(5.6*lx*ly*lz), 1./3.);
  int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);

  // common variables
  int id, index[7];
  double xpos, ypos, zpos, vol;
  std::vector<int> ff, neigh; // face_freq, neigh list
  std::vector<double> fs;     // face areas

  // real job
  if (triclinic){

    voro::container_periodic con(lx,xy,ly,xz,yz,lz,nx,ny,nz,8);
    // put atoms into the container
    for (int i = 1; i <= natom; ++i) con.put(i, atpos[i][0], atpos[i][1], atpos[i][2]);

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor c1, c2, *cell;
    voro::c_loop_all_periodic cl(con);
    if (cl.start()) do if (con.compute_cell(c1,cl)){
      for (int i = 0; i < 7; ++i) index[i] = 0;
       
      cl.pos(xpos,ypos,zpos);
      id = cl.pid();
      c1.neighbors(neigh);
      c1.face_areas(fs);
      cell = &c1;

      // print surface ratios if required
      if (fpsurf){
        int nf = fs.size();
        double wt = 1./cell->surface_area();
        for (int i = 0; i < nf; ++i) fprintf(fpsurf, "%lg %lg\n", fs[i]*wt, fs[i]);
      }

      // refine the voronoi cell if asked by removing tiny surfaces
      if (surf_min > ZERO){
        c2.init(-lx,lx,-ly,ly,-lz,lz);
        double fs_scale = 1.;
        if (flag_surf_ratio) fs_scale = 1./cell->surface_area();
  
        int nf = fs.size();
        // sort neighbors by area if asked to keep a minimum # of neighbors
        if (nminnei > 0){
          for (int i = 0; i < nf; ++i)
          for (int j = i+1; j < nf; ++j){
            if (fs[j] > fs[i]){
              double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
              int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
            }
          }
        }

        // add condition on surface
        for (int i = 0; i < nf; ++i){
          if (i < nminnei || (fs[i]*fs_scale) > surf_min){
            int jd = neigh[i];
            if (jd <= 0) continue;
  
            // apply pbc
            double xij = atpos[jd][0]-xpos;
            double yij = atpos[jd][1]-ypos;
            double zij = atpos[jd][2]-zpos;

            ApplyPBC(xij, yij, zij);

            c2.nplane(xij,yij,zij,jd);
          }
        }
        c2.face_areas(fs);
        c2.neighbors(neigh);
        cell = &c2;
      }

      vol = cell->volume();
      cell->face_freq_table(ff);
      int nn = ff.size()-1;
      for (int i = 3; i <= MIN(6,nn); ++i) index[i] = ff[i];
    
      // refine the voronoi cell if asked by skipping ultra short edges
      if (fabs(edge_min) > ZERO || fpedge) RefineEdge(fs.size(), cell, index, edge_min, fpedge);

      // assign voronoi index and neighbor list info
      int nf = fs.size();
      char vstr[MAXLINE];
      sprintf(vstr,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
      voro[id].assign(vstr);
      if (nf > MaxNei){
        MaxNei = nf + 2;
        memory->grow(neilist, MaxNei+1, natom+1, "neilist");
      }
      int fid = 0;
      for (int i = 0; i < nf; ++i){
        if (neigh[i] > 0) neilist[++fid][id] = neigh[i];;
      }
      neilist[0][id] = fid;
      volume[id] = vol;

      // output voro index info
      if (fp){
        fprintf(fp,"%d %d %lg %lg %lg %lg %s %d %d", id, attyp[id], xpos, ypos, zpos, vol, vstr, index[5], nf);
        for (int i = 0; i < nf; ++i) fprintf(fp," %d", neigh[i]);
        for (int i = 0; i < nf; ++i) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
      }

    } while (cl.inc());
    cell = NULL;

  } else {  // orthogonal box

    voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

    // put atoms into the container
    for (int i = 1; i <= natom; ++i) con.put(i, atpos[i][0], atpos[i][1], atpos[i][2]);

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor c1, c2, *cell;
    voro::c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(c1,cl)){
      for (int i = 0; i < 7; ++i) index[i] = 0;
       
      cl.pos(xpos,ypos,zpos);
      id = cl.pid();
      c1.neighbors(neigh);
      c1.face_areas(fs);
      cell = &c1;

      // print surface ratios if required
      if (fpsurf){
        int nf = fs.size();
        double wt = 1./cell->surface_area();
        for (int i = 0; i < nf; ++i) fprintf(fpsurf, "%lg %lg\n", fs[i]*wt, fs[i]);
      }

      // refine the voronoi cell if asked by removing tiny surfaces
      if (surf_min > ZERO){
        c2.init(-lx,lx,-ly,ly,-lz,lz);
        double fs_scale = 1.;
        if (flag_surf_ratio) fs_scale = 1./cell->surface_area();
  
        int nf = fs.size();
        // sort neighbors by area if asked to keep a minimum # of neighbors
        if (nminnei > 0){
          for (int i = 0; i < nf; ++i)
          for (int j = i+1; j < nf; ++j){
            if (fs[j] > fs[i]){
              double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
              int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
            }
          }
        }

        // add condition on surface
        for (int i = 0; i < nf; ++i){
          if (i < nminnei || (fs[i]*fs_scale) > surf_min){
            int jd = neigh[i];
            if (jd <= 0) continue;
  
            // apply pbc
            double xij = atpos[jd][0]-xpos;
            double yij = atpos[jd][1]-ypos;
            double zij = atpos[jd][2]-zpos;

            ApplyPBC(xij, yij, zij);

            c2.nplane(xij,yij,zij,jd);
          }
        }
        c2.face_areas(fs);
        c2.neighbors(neigh);
        cell = &c2;
      }

      vol = cell->volume();
      cell->face_freq_table(ff);
      int nn = ff.size()-1;
      for (int i = 3; i <= MIN(6,nn); ++i) index[i] = ff[i];
    
      // refine the voronoi cell if asked by skipping ultra short edges
      if (fabs(edge_min) > ZERO || fpedge) RefineEdge(fs.size(), cell, index, edge_min, fpedge);

      // assign voronoi index and neighbor list info
      int nf = fs.size();
      char vstr[MAXLINE];
      sprintf(vstr,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
      voro[id].assign(vstr);
      if (nf > MaxNei){
        MaxNei = nf + 2;
        memory->grow(neilist, MaxNei+1, natom+1, "neilist");
      }
      int fid = 0;
      for (int i = 0; i < nf; ++i){
        if (neigh[i] > 0) neilist[++fid][id] = neigh[i];;
      }
      neilist[0][id] = fid;
      volume[id] = vol;

      // output voro index info
      if (fp){
        fprintf(fp,"%d %d %lg %lg %lg %lg %s %d %d", id, attyp[id], xpos, ypos, zpos, vol, vstr, index[5], nf);
        for (int i = 0; i < nf; ++i) fprintf(fp," %d", neigh[i]);
        for (int i = 0; i < nf; ++i) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
      }

    } while (cl.inc());
    cell = NULL;
  }
return;
}

/*------------------------------------------------------------------------------
 * Recursive method to find neighbors of id upto max shells
 *------------------------------------------------------------------------------
 * il    (in)  : current shell level
 * max   (in)  : max shell #
 * clist (out) : list stores the IDs of atoms in the cluster
 * myshell (out) : shell # of each atom in the cluster
 *------------------------------------------------------------------------------
 * This subroutine cannot be called before the Voronoi info is constructed
 *----------------------------------------------------------------------------*/
void DumpAtom::voro_cluster(int il, const int max, int id, list<int> &clist, map<int,int> &myshell)
{
  if (++il > max) return;

  int nn = neilist[0][id];
  for (int ii = 1; ii <= nn; ++ii){
    int jd = neilist[ii][id];

    clist.push_back(neilist[ii][id]);
    if (myshell.count(jd)) myshell[jd] = MIN(myshell[jd], il);
    else myshell[jd] = il;

    voro_cluster(il, max, jd, clist, myshell);
  }

return;
}

/*------------------------------------------------------------------------------
 * Private method to check if a pair of atoms are bonded or not
 * Voronoi neighbors are seen as bonded
 * id  (in)  : atom ID of i
 * jd  (in)  : atom ID of j
 *----------------------------------------------------------------------------*/
int DumpAtom::bonded(int id, int jd)
{
  int ni = neilist[0][id];
  for (int jj = 1; jj <= ni; ++jj){
    if (neilist[jj][id] == jd) return 1;
  }

return 0;
}

/*------------------------------------------------------------------------------
 * Public method to compute the Voronoi index, Voronoi neighbor list info
 * mins     (in)  : surf_min Min#Nei edge_min ; if negative, default/previous
 *                : values will be taken.
 * fp       (in)  : file pointer to write full voro info; if NULL, write nothing
 * fpsurf   (in)  : file pointer to write surf ratio info; if NULL, write nothing
 * fpedge   (in)  : file pointer to write edge ratio info; if NULL, write nothing
 *----------------------------------------------------------------------------*/
void DumpAtom::ComputeVoro(double *mins){ ComputeVoro(mins, NULL, NULL, NULL);}
void DumpAtom::ComputeVoro(double *mins, FILE *fp, FILE *fpsurf, FILE *fpedge)
{
  if (type2radius) Radica_Voro(mins, fp, fpsurf, fpedge);
  else Direct_Voro(mins, fp, fpsurf, fpedge);

return;
}

/*------------------------------------------------------------------------------
 * Private method to compute the Radical Voronoi index, Voronoi neighbor list info
 * mins     (in)  : surf_min  Min#Nei edge_min; For either surf_min or edge_min,
 *                : if a positive number is provided, it is the absolute surface
 *                : area or edge length threshold; otherwise, it is the minimum
 *                : ratio over total surface area or circumference.
 * fp       (in)  : file pointer to write full voro info; if NULL, write nothing
 * fpsurf   (in)  : file pointer to write surf ratio info; if NULL, write nothing
 * fpedge   (in)  : file pointer to write edge ratio info; if NULL, write nothing
 *----------------------------------------------------------------------------*/
void DumpAtom::Radica_Voro(double *mins, FILE *fp, FILE *fpsurf, FILE *fpedge)
{
  // now to compute weighted voronoi info
  double diff = 0.;
  for (int i = 0; i < 3; ++i) diff += (mins[i] - vmins[i])*(mins[i] - vmins[i]);
  if (diff <= ZERO && voro.size() == natom && wted) return;

  wted = 1;
  int flag_surf_ratio = 0;
  if (mins[0] < 0.) flag_surf_ratio = 1;

  for (int i = 0; i < 3; ++i) vmins[i] = mins[i];

  double surf_min = fabs(vmins[0]);
  int nminnei = abs(int(vmins[1]));
  double edge_min = vmins[2];

  voro.clear();
  if (neilist) memory->destroy(neilist);
  if (volume)  memory->destroy(volume);
  memory->create(neilist, MaxNei+1, natom+1, "neilist");
  memory->create(volume,  natom+1, "volume");
  
  // need cartesian coordinates
  dir2car();

  // compute optimal size for container, and then contrust it
  double l = pow(double(natom)/(5.6*lx*ly*lz), 1./3.);
  int nx = int(lx*l+1), ny = int(ly*l+1), nz = int(lz*l+1);

  // common variables
  int id, index[7];
  double xpos, ypos, zpos, vol;
  std::vector<int> ff, neigh; // face_freq, neigh list
  std::vector<double> fs;     // face areas

  // real job
  if (triclinic){

    voro::container_periodic_poly con(lx,xy,ly,xz,yz,lz,nx,ny,nz,8);
    // put atoms into the container
    for (int i = 1; i <= natom; ++i){
      double ri = type2radius[attyp[i]];
      con.put(i, atpos[i][0], atpos[i][1], atpos[i][2], ri);
    }

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor c1, c2, *cell;
    voro::c_loop_all_periodic cl(con);
    if (cl.start()) do if (con.compute_cell(c1,cl)){
      for (int i = 0; i < 7; ++i) index[i] = 0;
       
      cl.pos(xpos,ypos,zpos);
      id = cl.pid();
      c1.neighbors(neigh);
      c1.face_areas(fs);
      cell = &c1;

      // print surface ratios if required
      if (fpsurf){
        int nf = fs.size();
        double wt = 1./cell->surface_area();
        for (int i = 0; i < nf; ++i) fprintf(fpsurf, "%lg %lg\n", fs[i]*wt, fs[i]);
      }

      // refine the voronoi cell if asked by removing tiny surfaces
      if (surf_min > ZERO){
        c2.init(-lx,lx,-ly,ly,-lz,lz);
        double fs_scale = 1.;
        if (flag_surf_ratio) fs_scale = 1./cell->surface_area();
  
        int nf = fs.size();
        // sort neighbors by area if asked to keep a minimum # of neighbors
        if (nminnei > 0){
          for (int i = 0; i < nf; ++i)
          for (int j = i+1; j < nf; ++j){
            if (fs[j] > fs[i]){
              double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
              int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
            }
          }
        }

        // add condition on surface
        for (int i = 0; i < nf; ++i){
          if (i < nminnei || (fs[i]*fs_scale) > surf_min){
            int jd = neigh[i];
            if (jd <= 0) continue;
  
            // apply pbc
            double xij = atpos[jd][0]-xpos;
            double yij = atpos[jd][1]-ypos;
            double zij = atpos[jd][2]-zpos;

            ApplyPBC(xij, yij, zij);

            double scale = type2radius[attyp[id]]/(type2radius[attyp[id]] + type2radius[attyp[jd]]);
            scale += scale;
            xij *= scale;
            yij *= scale;
            zij *= scale;
  
            c2.nplane(xij,yij,zij,jd);
          }
        }
        c2.face_areas(fs);
        c2.neighbors(neigh);
        cell = &c2;
      }

      vol = cell->volume();
      cell->face_freq_table(ff);
      int nn = ff.size()-1;
      for (int i = 3; i <= MIN(6,nn); ++i) index[i] = ff[i];
    
      // refine the voronoi cell if asked by skipping ultra short edges
      if (fabs(edge_min) > ZERO || fpedge) RefineEdge(fs.size(), cell, index, edge_min, fpedge);

      // assign voronoi index and neighbor list info
      int nf = fs.size();
      char vstr[MAXLINE];
      sprintf(vstr,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
      voro[id].assign(vstr);
      if (nf > MaxNei){
        MaxNei = nf + 2;
        memory->grow(neilist, MaxNei+1, natom+1, "neilist");
      }
      int fid = 0;
      for (int i = 0; i < nf; ++i){
        if (neigh[i] > 0) neilist[++fid][id] = neigh[i];;
      }
      neilist[0][id] = fid;
      volume[id] = vol;

      // output voro index info
      if (fp){
        fprintf(fp,"%d %d %lg %lg %lg %lg %s %d %d", id, attyp[id], xpos, ypos, zpos, vol, vstr, index[5], nf);
        for (int i = 0; i < nf; ++i) fprintf(fp," %d", neigh[i]);
        for (int i = 0; i < nf; ++i) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
      }

    } while (cl.inc());
    cell = NULL;

  } else {  // orthogonal box

    voro::container_poly con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

    // put atoms into the container
    for (int i = 1; i <= natom; ++i){
      double ri = type2radius[attyp[i]];
      con.put(i, atpos[i][0], atpos[i][1], atpos[i][2], ri);
    }

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor c1, c2, *cell;
    voro::c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(c1,cl)){
      for (int i = 0; i < 7; ++i) index[i] = 0;
       
      cl.pos(xpos,ypos,zpos);
      id = cl.pid();
      c1.neighbors(neigh);
      c1.face_areas(fs);
      cell = &c1;

      // print surface ratios if required
      if (fpsurf){
        int nf = fs.size();
        double wt = 1./cell->surface_area();
        for (int i = 0; i < nf; ++i) fprintf(fpsurf, "%lg %lg\n", fs[i]*wt, fs[i]);
      }

      // refine the voronoi cell if asked by removing tiny surfaces
      if (surf_min > ZERO){
  
        int nf = fs.size();
        double fs_scale = 1.;
        if (flag_surf_ratio) 1./cell->surface_area();

        // sort neighbors by area if asked to keep a minimum # of neighbors
        if (nminnei > 0){
          for (int i = 0; i < nf; ++i)
          for (int j = i+1; j < nf; ++j){
            if (fs[j] > fs[i]){
              double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
              int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
            }
          }
        }

        c2.init(-lx,lx,-ly,ly,-lz,lz);
        // add condition on surface
        for (int i = 0; i < nf; ++i){
          if (i < nminnei || (fs[i]*fs_scale) > surf_min){
            int jd = neigh[i];
            if (jd <= 0) continue;
  
  
            // apply pbc
            double xij = atpos[jd][0]-xpos;
            double yij = atpos[jd][1]-ypos;
            double zij = atpos[jd][2]-zpos;

            ApplyPBC(xij, yij, zij);

            double scale = type2radius[attyp[id]]/(type2radius[attyp[id]] + type2radius[attyp[jd]]);
            scale += scale;
            xij *= scale;
            yij *= scale;
            zij *= scale;
  
            c2.nplane(xij,yij,zij,jd);
          }
        }
        c2.face_areas(fs);
        c2.neighbors(neigh);
        cell = &c2;
      }

      vol = cell->volume();
      cell->face_freq_table(ff);
      int nn = ff.size()-1;
      for (int i = 3; i <= MIN(6,nn); ++i) index[i] = ff[i];
    
      // refine the voronoi cell if asked by skipping ultra short edges
      if (fabs(edge_min) > ZERO || fpedge) RefineEdge(fs.size(), cell, index, edge_min, fpedge);

      // assign voronoi index and neighbor list info
      int nf = fs.size();
      char vstr[MAXLINE];
      sprintf(vstr,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
      voro[id].assign(vstr);
      if (nf > MaxNei){
        MaxNei = nf + 2;
        memory->grow(neilist, MaxNei+1, natom+1, "neilist");
      }
      int fid = 0;
      for (int i = 0; i < nf; ++i){
        if (neigh[i] > 0) neilist[++fid][id] = neigh[i];;
      }
      neilist[0][id] = fid;
      volume[id] = vol;

      // output voro index info
      if (fp){
        fprintf(fp,"%d %d %lg %lg %lg %lg %s %d %d", id, attyp[id], xpos, ypos, zpos, vol, vstr, index[5], nf);
        for (int i = 0; i < nf; ++i) fprintf(fp," %d", neigh[i]);
        for (int i = 0; i < nf; ++i) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
      }

    } while (cl.inc());
    cell = NULL;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Private method to refine the Voronoi edge info
 * -----------------------------------------------------------------------------
 * cell (in)  : voro::voronoicell_neighbor
 * idx  (out) : carries the refined edge counts
 * ---------------------------------------------------------------------------*/
void DumpAtom::RefineEdge(int nf, voro::voronoicell_neighbor *cell, int *idx, double e_min_in, FILE *fpedge)
{
  std::vector<double> vpos;
  std::vector<int>    vlst;

  cell->vertices(vpos);
  cell->face_vertices(vlst);

  double edge_min = fabs(e_min_in);
  int flag_edge_ratio = 0;
  if (e_min_in < 0.) flag_edge_ratio = 1;

  int nv = vlst.size();
  double perim = 0., edges[nv];
  int k = 0, is = 0, ie = 0;
  while (k < nv){
    int ned = vlst[k++];
    for (int ii = 0; ii < ned; ++ii){
      int jj = (ii+1)%ned;
      int v1 = vlst[k+ii], v2 = vlst[k+jj];
      double dx = vpos[v1*3]   - vpos[v2*3];
      double dy = vpos[v1*3+1] - vpos[v2*3+1];
      double dz = vpos[v1*3+2] - vpos[v2*3+2];
      double r = sqrt(dx*dx+dy*dy+dz*dz);
      edges[ie++] = r;
      perim += r;
    }
    ++is; k += ned;
  }
  
  int ford[nf];
  k = is = ie = 0;
  double inv_perim = 2./perim;
  double ed_scale = 1.;
  if (flag_edge_ratio) ed_scale = inv_perim;

  while (k < nv){
    int ned = vlst[k++];

    int nuc = 0;
    for (int ii = 0; ii < ned; ++ii){
      double rwt = edges[ie] * inv_perim;
      if (fpedge) fprintf(fpedge, "%lg %g\n", rwt, edges[ie]);
      if ((edges[ie]*ed_scale) <= edge_min) ++nuc;
      ++ie;
    }
    ford[is++] = ned - nuc;
    k += ned;
  }

  if (edge_min > ZERO){
    for (int i = 3; i < 7; ++i) idx[i] = 0;
    for (int i = 0; i < nf; ++i){
      if (ford[i] < 7) idx[ford[i]] += 1;
    }
  }

return;
}

/* -----------------------------------------------------------------------------
 * Private method to apply PBC
 * -----------------------------------------------------------------------------
 * xij (inout) : x
 * yij (inout) : y
 * zij (inout) : z
 * ---------------------------------------------------------------------------*/
void DumpAtom::ApplyPBC(double &xij, double &yij, double &zij)
{
  if (cartesian) {
     if (triclinic) {  // non-orthogonal box
       while (zij > hz){
         xij -= xz;
         yij -= yz;
         zij -= lz;
       }
       while (zij <-hz){
         xij += xz;
         yij += yz;
         zij += lz;
       }
       
       while (yij > hy){
         xij -= xy;
         yij -= ly;
       }
       while (yij <-hy){
         xij += xy;
         yij += ly;
       }
       
       while (xij > hx) xij -= lx;
       while (xij <-hx) xij += lx;
   
     } else { // orthogonal box
   
       while (xij > hx) xij -= lx;
       while (xij <-hx) xij += lx;
     
       while (yij > hy) yij -= ly;
       while (yij <-hy) yij += ly;
     
       while (zij > hz) zij -= lz;
       while (zij <-hz) zij += lz;
     }

  } else {
    while (xij >= 0.5) xij -= 1.;
    while (xij < -0.5) xij += 1.;
    while (yij >= 0.5) yij -= 1.;
    while (yij < -0.5) yij += 1.;
    while (zij >= 0.5) zij -= 1.;
    while (zij < -0.5) zij += 1.;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Private method to free memory of Voro info
 * ----------------------------------------------------------------------------- */
void DumpAtom::FreeVoro()
{
  voro.clear();
  if (neilist) memory->destroy(neilist); neilist = NULL;
  if (volume)  memory->destroy(volume);  volume  = NULL;
  if (prop) memory->destroy(prop); prop = NULL;
  if (env)  memory->destroy(env);  env = NULL;

return;
}

/*------------------------------------------------------------------------------
 * Private method to identify the local env based on the CSP
 * parameter computed and the threshold given.
 *------------------------------------------------------------------------------ */
void DumpAtom::identify_env(const double thr)
{
  const int NMaxIt = 2000;
  int *swap, *otyp;
  memory->create(otyp, natom+1, "identify_env:otyp");

  if (env) memory->destroy(env);
  memory->create(env,  natom+1, "identify_env:env");

  for (int i = 1; i <= natom; ++i){
    if (prop[i] <= thr) env[i] = 0;   // crystal
    else env[i] = 2;                  // liquid or amorphous
  }

  // if nearly all neighbors of a crystal are liquid/amorphous, set it to be liquid/amorphous; and vice versa
  int nit = 0;
  while (nit < NMaxIt){
    nit++;
    swap = otyp; otyp = env; env = swap; swap = NULL;

    int nreset = 0;
    for (int id = 1; id <= natom; ++id){
      int neisum = 0;
      int ni = neilist[0][id];
      for (int in = 1; in <= ni; ++in){
        int jd = neilist[in][id];
        neisum += otyp[jd];
      }

      env[id] = otyp[id];
      if ( (otyp[id] == 2) && (neisum <= 2      ) ){ env[id] = 0; nreset++; }
      if ( (otyp[id] == 0) && (neisum >= ni+ni-2) ){ env[id] = 2; nreset++; }
    }
    if (nreset == 0) break;
  }
  if (nit >= NMaxIt) printf("\nWarning from identify_env: your cutoff might be improper!\n\n");

  swap = otyp; otyp = env; env = swap; swap = NULL;
  // assign type: 0, crystal; 1, crystal at interface; 2, liquid/amorphous; 3, liquid/amorhpous at interface
  for (int id = 1; id <= natom; ++id){
    int it = otyp[id];
    int ni = neilist[0][id];
    int neisum = 0;
    for (int in = 1; in <= ni; in++){
      int jd = neilist[in][id];
      neisum += otyp[jd];
    }
    if (neisum > 4 && neisum < ni+ni-4) ++it;
    else it = (neisum+4)/ni;

    env[id] = it;
  }

  nit = 0;
  while (nit < NMaxIt){
    ++nit;
    swap = otyp; otyp = env; env = swap; swap = NULL;

    int nreset = 0;
    for (int id = 1; id <= natom; ++id){
      int has[4];
      has[0] = has[1] = has[2] = has[3] = 0;
      int ni = neilist[0][id];
      for (int in = 1; in <= ni; ++in){
        int jd = neilist[in][id];
        int jt = otyp[jd];
        has[jt] = 1;
      }

      int it = otyp[id];
      env[id] = it;

      if (it == 0 && has[2] == 1) { env[id] = 1; nreset++; }
      if (it == 2 && has[0] == 1) { env[id] = 3; nreset++; }
      if (it == 1){
        if(has[0] == 0){
          nreset++;
          if (has[3]) env[id] = 3;
          else env[id] = 2;
        }
        if (has[2] == 0 && has[3] == 0){
          nreset++;
          env[id] = 0;
        }
      }
      if (it == 3){
        if (has[2] == 0){
          nreset++;
          if (has[1]) env[id] = 1;
          else env[id] = 0;
        }
        if (has[0] == 0 && has[1] == 0){
          nreset++;
          env[id] = 2;
        }
      }

    }

    if (nreset == 0) break;
  }
  if (nit >= NMaxIt) printf("\nWarning from identify_env: your cutoff might be inappropriate!\n\n");

  memory->destroy(otyp);

return;
}
/*------------------------------------------------------------------------------ */
