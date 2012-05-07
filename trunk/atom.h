#ifndef DUMP_ATOM_H
#define DUMP_ATOM_H

/* -----------------------------------------------------------------------------
 * A class which contains one frame of the atom style dump of lammps
 * -------------------------------------------------------------------------- */

#include "memory.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

using namespace std;

class DumpAtom {
public:
  DumpAtom(FILE *fp);
  ~DumpAtom();

  int natom, ntype, tstep;
  int initialized, cartesian;
  int triclinic;

  Memory *memory;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz;
  double lx, ly, lz, box[3], hbox[3];
  int *attyp, *atsel;  // note: atom IDs go from 1 to natom
  int *numtype;
  double **atpos;
  double axis[3][3];

  void selection(const char *);
  void car2dir();
  void dir2car();

private:
  int count_words(const char *);
};
#endif
