#ifndef DRIVER_H
#define DRIVER_H

#include "atom.h"
#include "memory.h"
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include "elements.h"

#define OutSurf 1
#define OutEdge 2
#define OutFeff 4

using namespace std;

class Driver {
public:
  Driver(int, char**);
  ~Driver();

private:
  char *dump;                         // input file name
  int nframe;                         // total # of frames
  int istr, iend, inc;                // frame range
  DumpAtom *one;                      // pointer to one frame
  vector<DumpAtom *> all;             // all frames from lammps atom style dump
  void readdump();                    // to read in the dump file

  // map atomic type to elements
  ChemElements *element;
  int *type2atnum;
  void MapType2Elem(const int, const int);

  // to define the range of frames to be analysed
  void setrange();
  int nsel;

  // to output selected frames as xyz file
  void writexyz();

  // to do voronoi diagram analysis
  int flag_out;
  void voro();

  // Chemical short range order based on voronoi neighbors
  void csro();

  // to average up all selected frames
  void avedump();

  // Hondeycutt-Andersen bond type analysis
  void honeycutt_andersen();
  void count_HA(const int, const int, int **, FILE *, const int);
  int bonded(int, int, int **);

  // Common neighber analysis or common neighbor parameter
  void Compute_CNACNP();

  // Pair Correlation and static structure factor
  void paircorr();
  void strfac();

  // Prepare for FEFF9
  void FEFF_main();
  void FEFF_input(int, FILE *);
  void FEFF_voro(set<string>, double *, int &, int **, int *, map<int,string>&); // get the neighbor list, list of atoms with selected voronoi index, and the voronoi index for all atoms
  void FEFF_cluster(int, const int, int **, int, list<int> &, map<int,int>&);    // get the list of atoms upto required shell of an atom, and their shell number

  // Connectivity info for certain clusters; needs FEFF_voro
  void ClusterConnectivity();

  // help info
  void help();

  Memory *memory;
  int count_words(const char *);
};
#endif
