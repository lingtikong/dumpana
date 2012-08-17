#ifndef DRIVER_H
#define DRIVER_H

#include "atom.h"
#include "memory.h"
#include <vector>

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

  // help info
  void help();

  Memory *memory;
  int count_words(const char *);
};
#endif
