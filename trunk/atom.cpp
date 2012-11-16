#include "atom.h"
#include "time.h"
#include "random.h"
#include "global.h"
#include "voro++.hh"

/*------------------------------------------------------------------------------
 * Constructor, to read one image from the atom style dump file of lammps
 *----------------------------------------------------------------------------*/
DumpAtom::DumpAtom(FILE *fp)
{
  iframe = natom = ntype = tstep = 0;
  initialized = triclinic = 0;
  xy = xz = yz = 0.;

  realcmd = NULL;
  attyp = atsel = numtype = NULL; atpos = x = s = NULL;

  voro.clear();
  neilist = NULL;
  MaxNei = 20;
  vmins[0] = vmins[1] = vmins[2] = 0.;

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

  for (int i=1; i<=natom; i++) atsel[i] = 1; atsel[0] = 0;
  realcmd = new char [MAXLINE]; strcpy(realcmd, "all");
  nsel = natom;

  s = memory->create(s, natom+1, 3, "s");
  atpos = s;

  for (int i=0; i<natom; i++){
    fgets(str,MAXLINE, fp);
    int id = atoi(strtok(str, " \n\t\r\f"));
    int ip = atoi(strtok(NULL," \n\t\r\f"));
    attyp[id] = ip; ntype = MAX(ip,ntype);
    s[id][0] = atof(strtok(NULL," \n\t\r\f"));
    s[id][1] = atof(strtok(NULL," \n\t\r\f"));
    s[id][2] = atof(strtok(NULL," \n\t\r\f"));
  }

  lx = box[0] = axis[0][0] = xhi - xlo;
  ly = box[1] = axis[1][1] = yhi - ylo;
  lz = box[2] = axis[2][2] = zhi - zlo;
  vol = lx*ly*lz;
  box[3] = axis[1][0] = xy;
  box[4] = axis[2][0] = xz;
  box[5] = axis[2][1] = yz;
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
  memory->destroy(numtype);

  atpos = NULL;
  memory->destroy(x);
  memory->destroy(s);

  memory->destroy(neilist);

  voro.clear();

  delete memory;
}


/*------------------------------------------------------------------------------
 * Method to convert fractional coordinates into cartesian
 *----------------------------------------------------------------------------*/
void DumpAtom::dir2car()
{
  if (x){
    atpos = x;
    return;
  }

  x = memory->create(x,natom+1, 3,"x");

  if (triclinic){
    for (int i=1; i<=natom; i++){
      x[i][0] = s[i][0]*lx + s[i][1]*xy + s[i][2]*xz + xlo;
      x[i][1] = s[i][1]*ly + s[i][2]*yz + ylo;
      x[i][2] = s[i][2]*lz + zlo;
    }
  } else {
    for (int i=1; i<=natom; i++){
      x[i][0] = s[i][0]*lx + xlo;
      x[i][1] = s[i][1]*ly + ylo;
      x[i][2] = s[i][2]*lz + zlo;
    }
  }

  atpos = x;
return;
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinates into fractional
 *----------------------------------------------------------------------------*/
void DumpAtom::car2dir()
{
  atpos = s;
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

      } else if (strcmp(oper,"%") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);
        if (logand){ for (int i = 1; i <= natom; i++) if (attyp[i]%ilow != ihigh) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (attyp[i]%ilow == ihigh) atsel[i] = 1;
        }

      } else break;

    } else if (strcmp(key,"x")==0 || strcmp(key,"y")==0 || strcmp(key,"z")==0 ){
    // selection by fractional position
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

      } else if (strcmp(oper,"%") == 0){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        ihigh = atoi(ptr);

        if (logand){ for (int i = 1; i <= natom; i++) if (i%ilow != ihigh) atsel[i] = 0;
        } else { for (int i = 1; i <= natom; i++) if (i%ilow == ihigh) atsel[i] = 1;
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
      for (int i = 1; i <= natom; i++) nsel += atsel[i];
      int ndel = nsel - ilow;
      while (ndel > 0){
        int id = MIN(random->uniform()*(natom+1), natom);
        if (atsel[id] == 1){
          atsel[id] = 0;
          ndel--;
        }
      }
      delete random;

    } else if (strcmp(key,"voro") == 0){ // selected by Voronoi indices; it must be the last selection option
      set<string> voroset; string vindex;
      double mins[3];
      for (int i=0; i<3; i++){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        mins[i] = atof(ptr);
      }
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) break;
      strcat(onecmd," "); strcat(onecmd,ptr);
      int nv = atoi(ptr);

      for (int i=0; i<nv; i++){
        ptr = strtok(NULL, " \n\t\r\f");
        if (ptr == NULL) break;
        strcat(onecmd," "); strcat(onecmd,ptr);
        vindex.assign(ptr); voroset.insert(vindex);
      }
      if (voroset.size() > 0) ComputeVoro(mins);

      // make the selection
      if (logand){
        for (int i = 1; i <= natom; i++){
          vindex = voro[i];
          if (voroset.count(vindex) < 1) atsel[i] = 0;
        }
      } else {
        for (int i = 1; i <= natom; i++){
          vindex = voro[i];
          if (voroset.count(vindex) > 0) atsel[i] = 1;
        }
      }

    } else if (strcmp(key,"all") == 0){ // select all; it will discard all previous selections
      strcpy(realcmd,""); strcpy(onecmd,key);
      for (int i = 1; i <= natom; i++) atsel[i] = 1;

    } else break;

    strcat(realcmd," "); strcat(realcmd,onecmd);
    key = strtok(NULL, " \n\t\r\f");
  }

  nsel = 0;
  for (int i = 1; i <= natom; i++) nsel += atsel[i];
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
  printf("\n"); for (int i=0; i<20; i++) printf("----");
  printf("\nThe grammar for the selection command is:\n\n");
  printf("  key op values [& key2 op2 values2 [| key3 op3 values3]]\n");
  printf("\nwhere `key` is either `type`, `x`,`y`,`z`, or `id`.\n");
  printf("It can also be `all`, which takes no argument and selects all atoms;\n");
  printf("or `ran num seed`, which takes no other argument and selects `num` atoms\n");
  printf("from the current selection randomly; `seed` is the seed for the uniform\n");
  printf("random number generator, if a non-positive number is provided, it will be\n");
  printf("set automatically based on current time. The logical operation before `ran`\n");
  printf("is always assumed to be AND, no matter what is defined;\n");
  printf("or `voro MinSurf%% MinEdge%% NMinNei NVoroIndex VoroIndices', which will\n");
  printf("select atoms with certain Voronoi indices. `MinSurf%%' and `MinEdge%%' defines\n");
  printf("the thresholds for surface and edge, respectively; they can be zero; `NMinNei'\n");
  printf("defines the Minimum # of neighbors wanted, zero means no limitations; `NVoroIndex'\n");
  printf("defines the # of Voronoi indices that will follow; and the following `NVoroIndex'\n");
  printf("arguments should be the Voronoi indices, for example: 0,6,0,8\n");
  printf("\n`op` is either `=`, `>`, `>=`, `<`, `<=`, `<>`, `><`, or `%%`,\n");
  printf("while neithor `=` nor `%%` is available for key=`x`, `y`, or `z`.\n");
  printf("\n`values` can be one or two numbers, depend on the type of `op`:\n");
  printf("for `<>`, `><`, and `%%` two numbers are need, otherwise one.\n");
  printf("`key <> num1 num2` selects `num1 <= kye <= num2`, `key >< num1 num2`\n");
  printf("selects `key <= num1 or key >= num2`, while `key %% num1 num2` selects\n");
  printf("atoms satisify `key%%num1 == num2`.\n");
  printf("\nMultiple `key op values` could be combined together, by either `&`\n");
  printf("or `|`, which means logical `and` or `or`, respectively. In this case\n");
  printf("the selections take effect sequentially.\n");
  printf("\nFor example, `type = 1 & x <> 0.1 0.5 & id >< 200 800` will select atoms\n");
  printf("of type 1 within 0.1 < x < 0.5 (fractional) and has id <= 200 or id >= 800;\n");
  printf("`type = 1 | type = 2` will select atoms of type 1 or 2; while `type = 1 & type = 2`\n");
  printf("will select nothing. `type = 1 & ran 100 0` will randomly select 100 atoms\n");
  printf("from all of type 1. `voro 0 0 0 1 0,0,12,0' will select all atoms that have\n");
  printf("a Voronoi index of 0,0,12,0.\n");
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

/*------------------------------------------------------------------------------
 * Method to compute the Voronoi index, Voronoi neighbor list info
 * mins     (in)  : surf_min% edge_min% Min#Nei; if negative, default/previous
 *                : values will be taken.
 * fp       (in)  : file pointer to write full voro info; if NULL, write nothing
 * fpsurf   (in)  : file pointer to write surf ratio info; if NULL, write nothing
 * fpedge   (in)  : file pointer to write edge ratio info; if NULL, write nothing
 *----------------------------------------------------------------------------*/
void DumpAtom::ComputeVoro(double *mins){ ComputeVoro(mins, NULL, NULL, NULL); }
void DumpAtom::ComputeVoro(double *mins, FILE *fp, FILE *fpsurf, FILE *fpedge)
{
  for (int i=0; i<3; i++) if (mins[i] < 0.) mins[i] = vmins[i];

  double diff = 0.;
  for (int i=0; i<3; i++) diff += (mins[i] - vmins[i])*(mins[i] - vmins[i]);
  if (diff <= ZERO && voro.size() == natom) return;

  for (int i=0; i<3; i++) vmins[i] = fabs(mins[i]);

  double surf_min = fabs(mins[0]);
  double edge_min = fabs(mins[1]);
  int nminnei = int(mins[2]);

  voro.clear();
  if (neilist) memory->destroy(neilist);
  neilist = memory->create(neilist, MaxNei+1, natom+1, "neilist");
  
  // set local variables
  double hx = 0.5*lx, hy = 0.5*ly, hz = 0.5*lz;

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
    for (int i=1; i<= natom; i++) con.put(i, atpos[i][0], atpos[i][1], atpos[i][2]);

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor c1, c2, *cell;
    voro::c_loop_all_periodic cl(con);
    if (cl.start()) do if (con.compute_cell(c1,cl)){
      for (int i=0; i<7; i++) index[i] = 0;
       
      cl.pos(xpos,ypos,zpos);
      id = cl.pid();
      c1.neighbors(neigh);
      c1.face_areas(fs);
      cell = &c1;

      // refine the voronoi cell if asked by removing tiny surfaces
      if (surf_min > ZERO){
        c2.init(-lx,lx,-ly,ly,-lz,lz);
  
        int nf = fs.size();
        // sort neighbors by area if asked to keep a minimum # of neighbors
        if (nminnei > 0){
          for (int i=0; i<nf; i++)
          for (int j=i+1; j<nf; j++){
            if (fs[j] > fs[i]){
              double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
              int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
            }
          }
        }

        // add condition on surface
        double fcut = surf_min * cell->surface_area();
        for (int i=0; i<nf; i++){
          if (fpsurf) fprintf(fpsurf, "%lg\n", fs[i]/cell->surface_area());
          if (i < nminnei || fs[i] > fcut){
            int j = neigh[i];
  
            // apply pbc
            double xij = atpos[j][0]-xpos;
            double yij = atpos[j][1]-ypos;
            double zij = atpos[j][2]-zpos;

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
  
            c2.nplane(xij,yij,zij,j);
          }
        }
        c2.face_areas(fs);
        c2.neighbors(neigh);
        cell = &c2;
      }

      vol = cell->volume();
      cell->face_freq_table(ff);
      int nn = ff.size()-1;
      for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
    
      // refine the voronoi cell if asked by skipping ultra short edges
      if (edge_min > ZERO){
        std::vector<double> vpos;
        std::vector<int>    vlst;
        double lcut2 = cell->total_edge_distance()*edge_min;
        lcut2 = lcut2*lcut2;

        cell->vertices(vpos);
        cell->face_vertices(vlst);

        int nf = fs.size();
        int ford[nf];
        int k = 0, iface = 0;
        while (k < vlst.size()){
          int ned = vlst[k++];
          int nuc = 0;
          for (int ii=0; ii<ned; ii++){
            int jj = (ii+1)%ned;
            int v1 = vlst[k+ii], v2 = vlst[k+jj];
            double dx = vpos[v1*3]   - vpos[v2*3];
            double dy = vpos[v1*3+1] - vpos[v2*3+1];
            double dz = vpos[v1*3+2] - vpos[v2*3+2];
            double r2 = dx*dx+dy*dy+dz*dz;
            if ((fpedge) && (v1 > v2)) fprintf(fpedge, "%lg\n", sqrt(r2)/cell->total_edge_distance());
            if (r2 <= lcut2) nuc++;
          }
          ford[iface++] = ned - nuc;
          k += ned;
        }

        for (int i=3; i<7; i++) index[i] = 0;
        for (int i=0; i<nf; i++){
          if (ford[i] < 7) index[ford[i]] += 1;
        }
      }

      // assign voronoi index and neighbor list info
      int nf = fs.size();
      char vstr[MAXLINE];
      sprintf(vstr,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
      voro[id].assign(vstr);
      if (nf > MaxNei){
        MaxNei = nf + 2;
        neilist = memory->grow(neilist, MaxNei+1, natom+1, "neilist");
      }
      neilist[0][id] = nf;
      for (int i=0; i<nf; i++) neilist[i+1][id] = neigh[i];

      // output voro index info
      if (fp){
        double wf = double(index[5])/double(nf)*100.;
        fprintf(fp,"%d %d %lg %lg %lg %lg %s %g %d", id, attyp[id], xpos, ypos, zpos, vol, vstr, wf, nf);
        for (int i=0; i<nf; i++) fprintf(fp," %d", neigh[i]);
        for (int i=0; i<nf; i++) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
      }

    } while (cl.inc());

  } else {  // orthogonal box

    voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,true,true,true,8);

    // put atoms into the container
    for (int i=1; i<= natom; i++) con.put(i, atpos[i][0], atpos[i][1], atpos[i][2]);

    // loop over all particles and compute their voronoi cell
    voro::voronoicell_neighbor c1, c2, *cell;
    voro::c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(c1,cl)){
      for (int i=0; i<7; i++) index[i] = 0;
       
      cl.pos(xpos,ypos,zpos);
      id = cl.pid();
      c1.neighbors(neigh);
      c1.face_areas(fs);
      cell = &c1;

      // refine the voronoi cell if asked by removing tiny surfaces
      if (surf_min > ZERO){
        c2.init(-lx,lx,-ly,ly,-lz,lz);
  
        int nf = fs.size();
        // sort neighbors by area if asked to keep a minimum # of neighbors
        if (nminnei > 0){
          for (int i=0; i<nf; i++)
          for (int j=i+1; j<nf; j++){
            if (fs[j] > fs[i]){
              double dswap = fs[i]; fs[i] = fs[j]; fs[j] = dswap;
              int ik = neigh[i]; neigh[i] = neigh[j]; neigh[j] = ik;
            }
          }
        }

        // add condition on surface
        double fcut = surf_min * cell->surface_area();
        for (int i=0; i<nf; i++){
          if (fpsurf) fprintf(fpsurf, "%lg\n", fs[i]/cell->surface_area());
          if (i < nminnei || fs[i] > fcut){
            int j = neigh[i];
  
            // apply pbc
            double xij = atpos[j][0]-xpos;
            while (xij > hx) xij -= lx;
            while (xij <-hx) xij += lx;

            double yij = atpos[j][1]-ypos;
            while (yij > hy) yij -= ly;
            while (yij <-hy) yij += ly;

            double zij = atpos[j][2]-zpos;
            while (zij > hz) zij -= lz;
            while (zij <-hz) zij += lz;
  
            c2.nplane(xij,yij,zij,j);
          }
        }
        c2.face_areas(fs);
        c2.neighbors(neigh);
        cell = &c2;
      }

      vol = cell->volume();
      cell->face_freq_table(ff);
      int nn = ff.size()-1;
      for (int i=3; i<= MIN(6,nn); i++) index[i] = ff[i];
    
      // refine the voronoi cell if asked by skipping ultra short edges
      if (edge_min > ZERO){
        std::vector<double> vpos;
        std::vector<int>    vlst;
        double lcut2 = cell->total_edge_distance()*edge_min;
        lcut2 = lcut2*lcut2;

        cell->vertices(vpos);
        cell->face_vertices(vlst);

        int nf = fs.size();
        int ford[nf];
        int k = 0, iface = 0;
        while (k < vlst.size()){
          int ned = vlst[k++];
          int nuc = 0;
          for (int ii=0; ii<ned; ii++){
            int jj = (ii+1)%ned;
            int v1 = vlst[k+ii], v2 = vlst[k+jj];
            double dx = vpos[v1*3]   - vpos[v2*3];
            double dy = vpos[v1*3+1] - vpos[v2*3+1];
            double dz = vpos[v1*3+2] - vpos[v2*3+2];
            double r2 = dx*dx+dy*dy+dz*dz;
            if ((fpedge) && (v1 > v2)) fprintf(fpedge, "%lg\n", sqrt(r2)/cell->total_edge_distance());
            if (r2 <= lcut2) nuc++;
          }
          ford[iface++] = ned - nuc;
          k += ned;
        }

        for (int i=3; i<7; i++) index[i] = 0;
        for (int i=0; i<nf; i++){
          if (ford[i] < 7) index[ford[i]] += 1;
        }
      }

      // assign voronoi index and neighbor list info
      int nf = fs.size();
      char vstr[MAXLINE];
      sprintf(vstr,"%d,%d,%d,%d", index[3], index[4], index[5], index[6]);
      voro[id].assign(vstr);
      if (nf > MaxNei){
        MaxNei = nf + 2;
        neilist = memory->grow(neilist, MaxNei+1, natom+1, "neilist");
      }
      neilist[0][id] = nf;
      for (int i=0; i<nf; i++) neilist[i+1][id] = neigh[i];

      // output voro index info
      if (fp){
        double wf = double(index[5])/double(nf)*100.;
        fprintf(fp,"%d %d %lg %lg %lg %lg %s %g %d", id, attyp[id], xpos, ypos, zpos, vol, vstr, wf, nf);
        for (int i=0; i<nf; i++) fprintf(fp," %d", neigh[i]);
        for (int i=0; i<nf; i++) fprintf(fp," %lg", fs[i]);
        fprintf(fp,"\n");
      }

    } while (cl.inc());
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
  for (int ii=1; ii <= nn; ii++){
    int jd = neilist[ii][id];

    clist.push_back(neilist[ii][id]);
    if (myshell.count(jd)) myshell[jd] = MIN(myshell[jd], il);
    else myshell[jd] = il;

    voro_cluster(il, max, jd, clist, myshell);
  }

return;
}

/*------------------------------------------------------------------------------ */
