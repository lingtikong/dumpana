#include "atom.h"

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
      for (int idim=0; idim<3; idim++){
        s[idim] = atpos[i][idim];
        //while (s[idim] >= 1.) s[idim] -= 1.;
        //while (s[idim] <  0.) s[idim] += 1.;
      }
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
 * Method to convert fractional coordinates into cartesian
 *----------------------------------------------------------------------------*/
void DumpAtom::selection(const char *selcmd)
{

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
