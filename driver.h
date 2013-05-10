#ifndef DRIVER_H
#define DRIVER_H

#include "atom.h"
#include "global.h"
#include "memory.h"
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include "elements.h"


using namespace std;

class Driver {
public:
  Driver(int, char**);
  ~Driver();

private:
  char *dump;                         // input file name
  int spk;                            // SPPARKS trajectory, cartesian coordinate
  int nframe;                         // total # of frames
  int istr, iend, inc;                // frame range
  DumpAtom *one;                      // pointer to one frame
  vector<DumpAtom *> all;             // all frames from lammps atom style dump
  void readdump();                    // to read in the dump file
  void MainMenu();

  // private method to set the cutoff for face or edge, needed by most others
  double mins[3];
  void set_cutoffs(int);

  // map atomic type to elements
  ChemElements *element;
  int *type2atnum;
  double *type2radius, *weighted;
  void MapType2Elem(const int, const int);
  void ShowRadius4Voro();

  // to define the range of frames to be analysed
  void setrange();
  int nsel;

  // to output selected frames
  void writesel();

  // to do voronoi diagram analysis
  int flag_out;
  void voro();

  // to output selected Voronoi cells as xyz file
  void OutputVoroCells();

  // Chemical short range order based on voronoi neighbors
  void csro();

  // to average up all selected frames
  void avedump();

  // Hondeycutt-Andersen bond type analysis
  void honeycutt_andersen();
  void count_HA(const int, const int, FILE *, const int);

  // Common neighber analysis or common neighbor parameter
  void Compute_CNACNP();

  // Pair Correlation and static structure factor
  void paircorr();
  void strfac();

  // Prepare for FEFF9
  void FEFF_main();
  void FEFF_input(int, FILE *);

  // Connectivity info for certain clusters
  void ClusterConnectivity();

  // bond length or bond angle distribution; Voronoi neighbors are seen as bonded
  void bonds();

  // spatial distribution of atoms
  void spatial();

  // help info
  void help();
  void ShowVersion();

  Memory *memory;
  int count_words(const char *);
};
#endif
