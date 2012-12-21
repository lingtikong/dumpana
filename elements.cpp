#include "elements.h"
#include "string.h"

/* -----------------------------------------------------------------------------
 * Function to return the atomic number of an input element symbol; 0 if not found
 * ---------------------------------------------------------------------------*/
int ChemElements::Name2Num(const char *ename)
{
  int AtNum;
  for (AtNum = 1; AtNum <= NumMax; AtNum++) if (strcmp(ename, symbol[AtNum]) == 0) break;
  if (AtNum > NumMax) AtNum = 0;

return AtNum;
}

/* -----------------------------------------------------------------------------
 * Function to return the element symbol for an input atomic number
 * ---------------------------------------------------------------------------*/
void ChemElements::Num2Name(const int AtNum, char * ename)
{
  int num = AtNum;
  if (num < 0||num > NumMax) num = 0;
  strcpy(ename, symbol[num]);

return;
}

/* -----------------------------------------------------------------------------
 * Return the molar mass of an input element symbol
 * ---------------------------------------------------------------------------*/
double ChemElements::Name2Mass(const char * ename)
{
  int AtNum = Name2Num(ename);

return weight[AtNum];
}

/* -----------------------------------------------------------------------------
 * Return the molar mass of an input atomic number
 * ---------------------------------------------------------------------------*/
double ChemElements::Num2Mass(const int AtNum)
{
  int num = AtNum;
  if (num < 0||num > NumMax) num = 0;

return weight[num];
}

/* -----------------------------------------------------------------------------
 * Return the atomic radius of an input element symbol
 * ---------------------------------------------------------------------------*/
double ChemElements::Name2Radius(const char * ename)
{
  int AtNum = Name2Num(ename);

return radius[AtNum];
}

/* -----------------------------------------------------------------------------
 * Return the atomic radius of an input atomic number
 * ---------------------------------------------------------------------------*/
double ChemElements::Num2Radius(const int AtNum)
{
  int num = AtNum;
  if (num < 0||num > NumMax) num = 0;

return radius[num];
}
/* -----------------------------------------------------------------------------
 * Atomic weight data taken from:
 * Pure Appl. Chem., Vol. 83, No. 2, pp. 359–396, 2011.
 * Atomic weights of the elements 2009 (IUPAC Technical Report)
 * ---------------------------------------------------------------------------*/
const int ChemElements::NumMax = 112;
const double ChemElements::weight[] = {0.,                                  // Null
    1.00797,    4.0026,     6.939,      9.012182,    10.811,                // H - B
   12.01115,   14.0067,    15.9994,    18.9984032,   20.17976,   22.989769, // C - Na
   24.30506,   26.9815386, 28.086,     30.973762,    32.064,                // Mg - S
   35.453,     39.948,     39.0983,    40.078,       44.955912,             // Cl - Sc
   47.867,     50.9415,    51.9961,    54.938045,    55.845,                // Ti - Fe
   58.933195,  58.6934,    63.546,     65.38,        69.723,                // Co - Ga
   72.63,      74.9216,    78.96,      79.904,       83.798,                // Ge - Kr
   85.4678,    87.62,      88.90585,   91.224,       92.90638,              // Rb - Nb
   95.96,      98.9062,   101.07,     102.9055,     106.42,                 // Mo - Pd
  107.8682,   112.411,    114.818,    118.71,       121.76,                 // Ag - Sb
  127.6,      126.90447,  131.293,    132.9054519,  137.327,                // Te - Ba
  138.90547,  140.116,    140.90765,  144.242,      147.,                   // La - Pm
  150.36,     151.964,    157.25,     158.92535,    162.5,                  // Sm - Dy
  164.93032,  167.259,    168.93421,  173.054,      174.9668,               // Ho - Lu
  178.49,     180.94788,  183.84,     186.207,      190.23,                 // Hf - Os
  192.217,    195.084,    196.966569, 200.59,       204.383,                // Ir - Tl
  207.2,      208.9804,   209.,       210.,         222.,                   // Pb - Rn
  223.,       226.025,    227.028,    232.03806,    231.03588,              // Fr - Pa
  238.02891,  237.048,    244.,       243.,         247.,                   // U - Cm
  247.,       251.,       252.,       257.,         258.,                   // Bk - Md
  259.,       260.,       261.11,     262.11,       263.12,                 // No - Sg
  262.12,     265.,       266.,       269.,         272.,      285.         // Bh - Cn
};

/* -----------------------------------------------------------------------------
 * Atomic radius data taken from:
 * Wikipage of Atomic radius of elements
 *
 * As atomic radius are not uniquely determined, the priority in selecting
 * the value for each element is as follows:
 *   Metals    : Metallic radius -> empirical -> van der Waals -> Calculated
 *   Non-metal : Covalent -> empirical -> van der Waals -> Calculated
 *   RareGas   : van der Waals -> empirical -> Calculated -> Covalent
 *   NoData    : 1.111
 * ---------------------------------------------------------------------------*/
const double ChemElements::radius[] = {1.111,// Null
  0.380, 0.320, 1.520, 1.120, 0.820, // H  - B   1  - 5
  0.770, 0.750, 0.730, 0.710, 1.540, // C  - Ne  6  - 10
  1.860, 1.600, 1.430, 1.110, 1.060, // Na - P   11 - 15
  1.020, 0.990, 1.880, 2.270, 1.970, // S  - Ca  16 - 20
  1.620, 1.470, 1.340, 1.280, 1.270, // Sc - Mn  21 - 25
  1.260, 1.250, 1.240, 1.280, 1.340, // Fe - Zn  26 - 30
  1.350, 1.220, 1.190, 1.160, 1.140, // Ga - Br  31 - 35
  2.020, 2.480, 2.150, 1.800, 1.600, // Kr - Zr  36 - 40
  1.460, 1.390, 1.360, 1.340, 1.340, // Nb - Rh  41 - 45
  1.370, 1.440, 1.510, 1.670, 1.450, // Pd - Sn  46 - 50
  1.450, 1.400, 1.330, 2.160, 2.650, // Sb - Cs  51 - 55
  2.220, 1.870, 1.818, 1.824, 1.814, // Ba - Nd  56 - 60
  1.834, 1.804, 1.804, 1.804, 1.773, // Pm - Tb  61 - 65
  1.781, 1.762, 1.761, 1.759, 1.760, // Dy - Yb  66 - 70
  1.738, 1.590, 1.460, 1.390, 1.370, // Lu - Re  71 - 75
  1.350, 1.355, 1.385, 1.440, 1.510, // Os - Hg  76 - 80
  1.700, 1.800, 0.160, 1.900, 2.020, // Tl - At  81 - 85
  2.200, 3.480, 2.150, 1.950, 1.790, // Rn - Th  86 - 90
  1.630, 1.560, 1.550, 1.590, 1.730, // Pa - Am  91 - 95
  1.740, 1.700, 1.860, 1.860, 1.111, // Cm - Fm  95 - 100
  1.111, 1.111, 1.111, 1.310, 1.260, // Md - Db  101- 105
  1.210, 1.190, 1.180, 1.130, 1.120, // Sg - Ds  106- 110
  1.180, 1.300                       // Rg - Cn  111- 112
};

/* -----------------------------------------------------------------------------
 * Element symbols
 * ---------------------------------------------------------------------------*/
const char ChemElements::symbol[][3] = { "X",//0
    "H",  "He", "Li", "Be", "B",           // 1  - 5
    "C",  "N",  "O",  "F",  "Ne", "Na",    // 6  - 11
    "Mg", "Al", "Si", "P",  "S",           // 12 - 16
    "Cl", "Ar", "K",  "Ca", "Sc",          // 17 - 21
    "Ti", "V",  "Cr", "Mn", "Fe",          // 22 - 26
    "Co", "Ni", "Cu", "Zn", "Ga",          // 27 - 31
    "Ge", "As", "Se", "Br", "Kr",          // 32 - 36
    "Rb", "Sr", "Y",  "Zr", "Nb",          // 37 - 41
    "Mo", "Tc", "Ru", "Rh", "Pd",          // 42 - 46
    "Ag", "Cd", "In", "Sn", "Sb",          // 47 - 51
    "Te", "I",  "Xe", "Cs", "Ba",          // 52 - 56
    "La", "Ce", "Pr", "Nd", "Pm",          // 57 - 61
    "Sm", "Eu", "Gd", "Tb", "Dy",          // 62 - 66
    "Ho", "Er", "Tm", "Yb", "Lu",          // 67 - 71
    "Hf", "Ta", "W",  "Re", "Os",          // 72 - 76
    "Ir", "Pt", "Au", "Hg", "Tl",          // 77 - 81
    "Pb", "Bi", "Po", "At", "Rn",          // 82 - 86
    "Fr", "Ra", "Ac", "Th", "Pa",          // 87 - 91
    "U",  "Np", "Pu", "Am", "Cm",          // 92 - 96
    "Bk", "Cf", "Es", "Fm", "Md",          // 97 - 101
    "No", "Lr", "Rf", "Db", "Sg",          // 102- 106 
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"     // 107- 112
};
