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
#include <vector>
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
  int **image;
  int *numtype;
  double **atpos;
  double axis[3][3];

  // Atomic properties, anything in the dump file besides id, type, x, y, z, ix, iy, iz
  vector<string> prop_label;
  double **atprop;

  // Public interface to compute the Voronoi tessellation
  double *type2radius;
  void ComputeVoro(double *);
  void ComputeVoro(double *, FILE *, FILE *, FILE *);

  // the following four variables/functions cannot be called before ComputeVoro
  int neiMethod;         // methoed to generate the neighbor list
  int **neilist;         // Voronoi neighbor list, only available if voro is computed
  map<int,string> voro;  // Voronoi index for each atom
  int bonded(int,int);   // check if two atoms are neighbors to each other
  void voro_cluster(int, const int, int, list<int> &, map<int,int> &); // find neighbors of an atom upto certain Voronoi shells
  double *volume;        // Voronoi volume
  void FreeVoro();       // free memory of neilist, voro, and volume

  // Neibhbor List calculation based on distance instead of voro
  void ComputeNeiList(double ***);

  void selection(const char *);
  void SelInfo();
  void SelHelp();

  void car2dir();
  void dir2car();

  void ApplyPBC(double &, double &, double &);
  double get_dist2(double, double, double);

  int *env;
  double *prop;          // Property, can by anything
  void identify_env(const double);

  double smix;           // Mixing entropy
  int flag_smix;         // 1, computed; 0, not yet

private:
  // private method to compute the Voronoi tessellation; the real engine
  int wted;
  void Direct_Voro(double *, FILE *, FILE *, FILE *);
  void Radica_Voro(double *, FILE *, FILE *, FILE *);

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
