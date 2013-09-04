#ifndef DUMP_ATOM_H
#define DUMP_ATOM_H

/* -----------------------------------------------------------------------------
 * A class which contains one frame of the atom style dump of lammps
 * -------------------------------------------------------------------------- */

#include "memory.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <map>
#include <set>
#include <list>
#include <string>
#include "voro++.hh"

using namespace std;

class DumpAtom {
public:
  DumpAtom(FILE *fp, const char *, const int);
  ~DumpAtom();

  int iframe;
  int natom, ntype, tstep, nsel;
  int initialized;
  int triclinic;

  char *fname;

  Memory *memory;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz, vol;
  double lx, ly, lz, box[6], hbox[3], h_inv[6];
  double hx, hy, hz;
  int *attyp, *atsel;  // note: atom IDs go from 1 to natom; type ID from 1 to ntype
  int *numtype;
  double **atpos;
  double axis[3][3];

  int weighted;
  // Compute the Voronoi info, equal-distance method
  void ComputeVoro(double *);
  void ComputeVoro(double *, FILE *, FILE *, FILE *);

  // Compute the Voronoi info, weighted-distance method
  void ComputeVoro(double *, double *); 
  void ComputeVoro(double *, FILE *, FILE *, FILE *, double *);

  // the following four variables/functions cannot be called before ComputeVoro
  int **neilist;         // Voronoi neighbor list, only available if voro is computed
  map<int,string> voro;  // Voronoi index for each atom
  int bonded(int,int);   // check if two atoms are neighbors to each other
  void voro_cluster(int, const int, int, list<int> &, map<int,int> &); // find neighbors of an atom upto certain Voronoi shells
  double *volume;        // Voronoi volume
  void FreeVoro();       // free memory of neilist, voro, and volume

  void selection(const char *);
  void SelInfo();
  void SelHelp();

  void car2dir();
  void dir2car();

  void ApplyPBC(double &, double &, double &);

private:
  int MaxNei;
  double vmins[3];

  int least_memory, cartesian;
  double **x, **s;

  char *realcmd;
  int count_words(const char *);

  // used by ComputeVoro to refine edges
  void RefineEdge(int,voro::voronoicell_neighbor *,int *,double, FILE *);
};
#endif
