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
  void voro();

  // to average up all selected frames
  void avedump();

  // help info
  void help();

  Memory *memory;
  int count_words(const char *);
};
#endif
