#ifndef COMPUTE_CNA_ATOM_H
#define COMPUTE_CNA_ATOM_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"

class ComputeCNAAtom {
public:
  ComputeCNAAtom(const int, const int, int **, int **, double **, double *, FILE *);
  ~ComputeCNAAtom();

private:
  void init();
  void compute_cna();
  void compute_cnp();

  void output(FILE *);

  Memory *memory;

  int natom;
  int **nearest;
  int **bonded;
  double **x, L[3], hL[3];
  double *pattern;

};

#endif
