#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

/* -----------------------------------------------------------------------------
 * A class that calcualtes the Spherical Harmonics function Y(l, m, theta, phi)
 * Instead of "theta" and "phi", the code is fed with the vector r.
 * -------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"

using namespace std;

class SphericalHarmonics {
public:
  SphericalHarmonics();
  double Y(int, int, double*);
  double w3j(int,int,int,int);

private:
  double K(int, int);
  double P(int, int, double);
  double dfact(int);

  double sign[2];
};
#endif
