#ifndef COMPUTE_CNA_ATOM_H
#define COMPUTE_CNA_ATOM_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"

class ComputeCNAAtom {
public:
  ComputeCNAAtom(const int, const int, int **, double **, double *, FILE *);
  ~ComputeCNAAtom();

private:
  void init();
  void compute_cna();
  void compute_cnp();
  void centro_atom(const int);

  void apply_pbc(double &, double &, double &);

  void output(FILE *);

  Memory *memory;

  int natom;
  int **nearest;
  int non_ortho_box;
  double **x, L[3], hL[3];
  double xy, xz, yz;
  double *pattern;

  int bonded(int,int);

  // copied from LAMMPS compute_centro_atom
  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

#endif
