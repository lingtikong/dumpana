#ifndef DRIVER_H
#define DRIVER_H

#include "atom.h"
#include "memory.h"
#include <list>
#include <vector>
#include "elements.h"

class Driver {
public:
  Driver(int, char**);
  ~Driver();

private:
  char *dump;                         // input file name
  int nframe;                         // total # of frames
  int istr, iend, inc;                // frame range
  DumpAtom *one;                      // pointer to one frame
  std::vector<DumpAtom *> all;        // all frames from lammps atom style dump
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
  void count_HA(const int, const int, int **, int**, FILE *, const int);

  // Common neighber analysis or common neighbor parameter
  void Compute_CNACNP();

  // Pair Correlation
  void paircorr();

  // Prepare for FEFF9
  void FEFF_main();
  void FEFF_input(int, FILE *);
  void FEFF_voro(int, int &, int **, int *, double *);
  void FEFF_cluster(int, int, int **, int, std::list<int> &);

  // help info
  void help();

  Memory *memory;
  int count_words(const char *);
};
#endif
