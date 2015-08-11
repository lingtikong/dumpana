#ifndef COMPUTE_CNA_ATOM_H
#define COMPUTE_CNA_ATOM_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "atom.h"

class ComputeCNAAtom {
public:
  ComputeCNAAtom(const int, DumpAtom *, FILE *, double);
  ~ComputeCNAAtom();

private:
  void init();
  void compute_cna();
  void compute_cnp();
  void centro_atom(const int);

  void output(FILE *);

  Memory *memory;

  DumpAtom *one;
  int *attyp, **neilist;
  double *pattern, **x;
  double *lop_sum;

  // determine local environment
  int flag_env;
  double thr_env;
  
  // copied from LAMMPS compute_centro_atom
  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

#endif
