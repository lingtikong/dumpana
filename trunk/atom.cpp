#include "atom.h"
#include "time.h"
#include "random.h"

#define MAXLINE 512
#define ZERO 1.e-8
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Constructor, to read one image from the atom style dump file of lammps
 *----------------------------------------------------------------------------*/
DumpAtom::DumpAtom(FILE *fp)
{
  natom = ntype = tstep = 0;
  initialized = cartesian = triclinic = 0;
  xy = xz = yz = 0.;

  realcmd = NULL;
  attyp = atsel = numtype = NULL; atpos = NULL;

  char str[MAXLINE];
  if (fgets(str,MAXLINE, fp) == NULL) return;
  fgets(str,MAXLINE, fp);
  tstep = atoi(strtok(str, " \n\t\r\f"));

  fgets(str,MAXLINE, fp);
  fgets(str,MAXLINE, fp);
  natom = atoi(strtok(str, " \n\t\r\f"));
  if (natom < 1) return;

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
  
  if (xy*xy+xz*xz+yz*yz > ZERO) triclinic = 1;

  fgets(str,MAXLINE, fp);

  memory = new Memory();
  attyp = memory->create(attyp, natom+1, "attyp");
  atsel = memory->create(atsel, natom+1, "atsel");
  atpos = memory->create(atpos, natom+1, 3, "atpos");
  for (int i=1; i<=natom; i++) atsel[i] = 1;
  realcmd = new char [MAXLINE]; strcpy(realcmd, "all");
  nsel = natom;

  for (int i=0; i<natom; i++){
    fgets(str,MAXLINE, fp);
    int id = atoi(strtok(str, " \n\t\r\f"));
    int ip = atoi(strtok(NULL," \n\t\r\f"));
    attyp[id] = ip; ntype = MAX(ip,ntype);
    atpos[id][0] = atof(strtok(NULL," \n\t\r\f"));
    atpos[id][1] = atof(strtok(NULL," \n\t\r\f"));
    atpos[id][2] = atof(strtok(NULL," \n\t\r\f"));
  }

  lx = box[0] = axis[0][0] = xhi - xlo;
  ly = box[1] = axis[1][1] = yhi - ylo;
  lz = box[2] = axis[2][2] = zhi - zlo;
  vol = lx*ly*lz;
  axis[1][0] = xy;
  axis[2][0] = xz;
  axis[2][1] = yz;
  axis[0][1] = axis[0][2] = axis[1][2] = 0.;
  for (int idim=0; idim<3; idim++) hbox[idim] = 0.5*box[idim];

  numtype = memory->create(numtype,ntype+1,"numtype");
  for (int i=0; i<=ntype; i++) numtype[i] = 0;
  for (int i=1; i<=natom; i++) numtype[attyp[i]]++;

  initialized = 1;

return;
}

/*------------------------------------------------------------------------------
 * Deconstructor, to free allocated memory
 *----------------------------------------------------------------------------*/
DumpAtom::~DumpAtom()
{
  if (realcmd) delete []realcmd;
  memory->destroy(attyp);
  memory->destroy(atsel);
  memory->destroy(atpos);
  memory->destroy(numtype);

  delete memory;
}


/*------------------------------------------------------------------------------
 * Method to convert fractional coordinates into cartesian
 *----------------------------------------------------------------------------*/
void DumpAtom::dir2car()
{
  if (cartesian == 1) return;
  double s[3];

  if (triclinic){
    for (int i=1; i<=natom; i++){
      for (int idim=0; idim<3; idim++){
        s[idim] = atpos[i][idim];
        while (s[idim] >= 1.) s[idim] -= 1.;
        while (s[idim] <  0.) s[idim] += 1.;
      }
      atpos[i][0] = s[0]*lx + s[1]*xy + s[2]*xz + xlo;
      atpos[i][1] = s[1]*ly + s[2]*yz + ylo;
      atpos[i][2] = s[2]*lz + zlo;
    }
  } else {
    for (int i=1; i<=natom; i++){
      for (int idim=0; idim<3; idim++) s[idim] = atpos[i][idim];
      atpos[i][0] = s[0]*lx + xlo;
      atpos[i][1] = s[1]*ly + ylo;
      atpos[i][2] = s[2]*lz + zlo;
    }
  }
  cartesian = 1;
return;
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinates into fractional
 *----------------------------------------------------------------------------*/
void DumpAtom::car2dir()
{
  if (cartesian == 0) return;
  double h[6], h_inv[6], delta[3];

  h[0] = lx; h[1] = ly; h[2] = lz;
  h[3] = yz; h[4] = xz; h[5] = xy;
  h_inv[0] = 1./h[0];
  h_inv[1] = 1./h[1];
  h_inv[2] = 1./h[2];

  if (triclinic){
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);

    for (int i = 1; i <= natom; i++) {
      delta[0] = atpos[i][0] - xlo;
      delta[1] = atpos[i][1] - ylo;
      delta[2] = atpos[i][2] - zlo;
  
      atpos[i][0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
      atpos[i][1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
      atpos[i][2] = h_inv[2]*delta[2];
    }
  } else {
    for (int i = 1; i <= natom; i++) {
      delta[0] = atpos[i][0] - xlo;
      delta[1] = atpos[i][1] - ylo;
      delta[2] = atpos[i][2] - zlo;
  
      atpos[i][0] = h_inv[0]*delta[0];
      atpos[i][1] = h_inv[1]*delta[1];
      atpos[i][2] = h_inv[2]*delta[2];
    }
  }

  cartesian = 0;
return;
}

/*------------------------------------------------------------------------------
 * Method to select some atoms, not implemented yet
 *----------------------------------------------------------------------------*/
void DumpAtom::selection(const char *line)
{
  int n = strlen(line) + 1;
  char *selcmd = (char *) memory->smalloc(n*sizeof(char),"selcmd");
  strcpy(selcmd,line);

  char *key, *oper, *ptr;
  int ilow, ihigh;
  double rlow, rhigh;
  int logand = 1;

  char onecmd[MAXLINE];
  strcpy(realcmd," ");

  // by default, all are selected
  for (int i = 1; i <= natom; i++) atsel[i] = 1;

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
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i] != ilow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i] == ilow) atsel[i] = 1;
        }

      } else if (strcmp(oper,">") == 0){
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i] <= ilow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i] > ilow) atsel[i] = 1;
        }

      } else if (strcmp(oper,">=") == 0){
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i] < ilow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i] >= ilow) atsel[i] = 1;
        }

      } else if (strcmp(oper,"<") == 0){
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i] >= ilow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i] < ilow) atsel[i] = 1;
        }

      } else if (strcmp(oper,"<=") == 0){
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i] > ilow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i] <= ilow) atsel[i] = 1;
        }

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i]<ilow || attyp[i]>ihigh) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i]>=ilow && attyp[i]<=ihigh) atsel[i] = 1;
        }

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i]>ilow && attyp[i]<ihigh) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i]<=ilow || attyp[i]>=ihigh) atsel[i] = 1;
        }

      } else break;

    } else if (strcmp(key,"x")==0 || strcmp(key,"y")==0 || strcmp(key,"z")==0 ){
    // selection by fractional position
      int dir = key[0]-'x';
      oper = strtok(NULL, " \n\t\r\f");
      if (oper == NULL) break;
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      rlow = atof(ptr);

      strcat(onecmd," "); strcat(onecmd,oper);
      strcat(onecmd," "); strcat(onecmd,ptr);

      if (strcmp(oper,">")==0 || strcmp(oper,">=")==0){
        if (logand){ for (int i = 1; i <= natom; i++) if (atpos[i][dir] < rlow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (atpos[i][dir] >= rlow) atsel[i] = 1;
        }

      } else if (strcmp(oper,"<")==0 || strcmp(oper,"<=")==0){
        if (logand){ for (int i = 1; i <= natom; i++) if (atpos[i][dir] > rlow) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (atpos[i][dir] <= rlow) atsel[i] = 1;
        }

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (logand){ for (int i = 1; i <= natom; i++) if (atpos[i][dir]<rlow || atpos[i][dir]>rhigh) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (atpos[i][dir]>=rlow && atpos[i][dir]<=rhigh) atsel[i] = 1;
        }

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        rhigh = atof(ptr);

        if (logand){ for (int i = 1; i <= natom; i++) if (atpos[i][dir]>rlow && atpos[i][dir]<rhigh) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (atpos[i][dir]<=rlow || atpos[i][dir]>=rhigh) atsel[i] = 1;
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
        if (logand){ for (int i = 1; i <= natom; i++) if (i != ilow) atsel[i] = 0;
        } else atsel[ilow] = 1;

      } else if (strcmp(oper,">") == 0){
        if (logand) for (int i = 1; i <= MIN(ilow,natom); i++) atsel[i] = 0;
        else for (int i = ilow+1; i <= natom; i++) atsel[i] = 1;

      } else if (strcmp(oper,">=") == 0){
        if (logand) for (int i = 1; i < MIN(ilow,natom); i++) atsel[i] = 0;
        else for (int i = ilow; i <= natom; i++) atsel[i] = 1;

      } else if (strcmp(oper,"<") == 0){
        if (logand) for (int i = ilow; i <= natom; i++) atsel[i] = 0;
        else for (int i = 1; i < MIN(ilow,natom); i++) atsel[i] = 1;

      } else if (strcmp(oper,"<=") == 0){
        if (logand) for (int i = ilow+1; i <= natom; i++) atsel[i] = 0;
        else for (int i = 1; i <= MIN(ilow,natom); i++) atsel[i] = 1;

      } else if (strcmp(oper,"<>") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);

        if (logand){
          for (int i = 1; i < MIN(ilow,natom); i++) atsel[i] = 0;
          for (int i = ihigh+1; i <= natom; i++) atsel[i] = 0;
        } else for (int i = ilow; i <= MIN(ihigh,natom); i++) atsel[i] = 1;

      } else if (strcmp(oper,"><") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);

        if (logand) for (int i = ilow; i <= MIN(ihigh,natom); i++) atsel[i] = 0;
        else {
          for (int i = 1; i < MIN(ilow,natom); i++) atsel[i] = 1;
          for (int i = ihigh+1; i <= natom; i++) atsel[i] = 1;
        }

      } else break;

    } else if (strcmp(key,"ran") == 0){ // random selection from current selection
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      ilow = atoi(ptr);

      strcat(onecmd," "); strcat(onecmd,ptr);

      nsel = 0;
      for (int i = 1; i <= natom; i++) nsel += atsel[i];
      int ndel = nsel - ilow;
      RanPark * random = new RanPark(time(NULL));
      while (ndel > 0){
        int id = random->uniform()*natom;
        if (atsel[id] == 1){
          atsel[id] = 0;
          ndel--;
        }
      }
      delete random;

    } else if (strcmp(key,"all") == 0){ // select all; it will discard all previous selections
      strcpy(realcmd,"");
      strcpy(onecmd,key);
      for (int i = 1; i <= natom; i++) atsel[i] = 1;

    } else break;

    strcat(realcmd," "); strcat(realcmd,onecmd);
    key = strtok(NULL, " \n\t\r\f");
  }

  nsel = 0;
  for (int i = 1; i <= natom; i++) nsel += atsel[i];

return;
}

/*------------------------------------------------------------------------------
 * Method to print out selection info
 *----------------------------------------------------------------------------*/
void DumpAtom::SelInfo()
{
  printf("\nThe realized selection command is: %s,\n", realcmd);
  printf("and %d of %d atoms are in selection.\n", nsel, natom);

return;
}

/*------------------------------------------------------------------------------
 * Help info related to the selection function
 *----------------------------------------------------------------------------*/
void DumpAtom::SelHelp()
{
  printf("\n"); for (int i=0; i<20; i++) printf("----");
  printf("\nThe grammar for the selection command is:\n\n");
  printf("  key op values [& key2 op2 values2 [| key3 op3 values3]]\n");
  printf("\nwhere %ckey%c is either %ctype%c, %cx%c,%cy%c,%cz%c, or %cid%c.\n",
    char(96),char(39),char(96),char(39),char(96),char(39),char(96),char(39),
    char(96),char(39),char(96),char(39));
  printf("It can also be %call%c, which takes no argument and selects all atoms;\n",
    char(96),char(39));
  printf("or %cran num%c, which takes no other argument and selects %cnum%c atoms\n",
    char(96),char(39), char(96),char(39));
  printf("from the current selection randomly. The logical operation before %cran%c\n",
    char(96),char(39));
  printf("is always assumed to be AND, no matter what is defined.\n");
  printf("\n%cop%c is either %c=%c (N.A. for %cx-z%c), %c>%c, %c>=%c, %c<%c, %c<=%c, %c<>%c, or %c><%c.\n",
    char(96),char(39),char(96),char(39),char(96),char(39),char(96),char(39),
    char(96),char(39),char(96),char(39),char(96),char(39),char(96),char(39),
    char(96),char(39));
  printf("\n%cvalues%c can be one or two numbers, depend on the type of %cop%c:\n",
    char(96),char(39),char(96),char(39));
  printf("for %c<>%c and %c><%c two numbers are need, otherwise one.\n",
    char(96),char(39),char(96),char(39));
  printf("\nMultiple %ckey op values%c could be combined together, by either %c&%c\n",
    char(96),char(39),char(96),char(39),char(96),char(39));
  printf("or %c|%c, which means logical %cand%c or %cor%c, respectively. In this case\n",
    char(96),char(39),char(96),char(39),char(96),char(39));
  printf("the selections take effect sequentially.\n");
  printf("\nFor example, %ctype = 1 & x <> 0.1 0.5 & id >< 200 800%c will select atoms\n",
    char(96),char(39));
  printf("of type 1 within 0.1 < x < 0.5 (fractional) and has id <= 200 or id >= 800;\n");
  printf("%ctype = 1 | type = 2%c will select atoms of type 1 or 2; while %ctype = 1 & type = 2%c\n",
    char(96),char(39),char(96),char(39));
  printf("will select nothing. %ctype = 1 & ran 100%c will randomly select 100 atoms\n",
    char(96),char(39));
  printf("from all of type 1.\n");
  for (int i=0; i<20; i++) printf("----"); printf("\n\n");
    
return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int DumpAtom::count_words(const char *line)
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
/*------------------------------------------------------------------------------ */
