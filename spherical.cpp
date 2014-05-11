#include <cmath>
#include "spherical.h"
#include "global.h"
#include "gsl/gsl_sf_coupling.h"

/*------------------------------------------------------------------------------
 * Constructor, initialize some global variables
 *----------------------------------------------------------------------------*/
SphericalHarmonics::SphericalHarmonics()
{
  sign[0] =  1.;
  sign[1] = -1.;

return;
}

/*------------------------------------------------------------------------------
 * Function to get the double Factorial of a positive number; for any non-positive
 * number, 1 is returned
 *
 * (2k-1)!! = prod_i=1^k (2i-1)
 *----------------------------------------------------------------------------*/
double SphericalHarmonics::dfact(int k)
{
  if (k <= 1) return 1.;
  return double(k+k-1)*dfact(k-1);
}
  
/*------------------------------------------------------------------------------
 * K_l^m, l must be positive, and l >= abs(m), this won't be checked
 *----------------------------------------------------------------------------*/
double SphericalHarmonics::K(int l, int m)
{
  const double inv_fpi = 1./(16.*atan(1.));

  double tmp = double(l+l+1)*inv_fpi;
  double fac = 1.;
  for (int i = l-abs(m)+1; i <= l+abs(m); ++i) fac *= double(i);

return sqrt(tmp/fac);
}
  

/*------------------------------------------------------------------------------
 * P_l^m(z), both l and m must be positive, and l >= abs(m), this won't be checked
 *----------------------------------------------------------------------------*/
double SphericalHarmonics::P(int l, int m, double z)
{
  m = abs(m);
  double res;
  if (l == m && l == 0) res = 1.;
  else if (l == m)      res = dfact(m)*pow(1.-z*z, 0.5*double(m));
  else if (l == m+1)    res = z*double(m+m+1)*P(m, m, z);
  else                  res = z*double(l+l-1)/double(l-m)*P(l-1, m, z) - double(l+m-1)/double(l-m)*P(l-2, m, z);

return res;
}

/*------------------------------------------------------------------------------
 * Y_l^m(r), l must be positive, and l >= abs(m), this won't be checked
 * r is the vector
 *----------------------------------------------------------------------------*/
double SphericalHarmonics::Y(int l, int m, double *r)
{
  const double PI = M_PI;
  const double rq2 = sqrt(2.);
  double r0   = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  double cosT = r[2]/r0;
  double phi;

  if (fabs(r[0]) < ZERO) phi = PI*0.5;
  else phi = atan(r[1]/r[0]);

  if (r[0] < 0.) phi += PI;
  else if (r[1] < 0.) phi += PI+PI;

  double res = P(l, m, cosT);
  if (m == 0)     res *= K(l, m);
  else if (m > 0) res *= K(l, m) * rq2 * cos( double(m)*phi);
  else            res *= K(l, m) * rq2 * sin(-double(m)*phi) *sign[m%2];

return res;
}

/*------------------------------------------------------------------------------
 * Wigner 3j
 *----------------------------------------------------------------------------*/
double SphericalHarmonics::w3j(int l, int m1, int m2, int m3)
{
  return gsl_sf_coupling_3j(l+l, l+l, l+l, m1+m1, m2+m2, m3+m3);
}

/*----------------------------------------------------------------------------*/
