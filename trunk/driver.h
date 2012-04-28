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
  int once; // indicate that the only one analysis is needed
  char *dump;
  int nframe, img;
  int istr, iend, inc;
  DumpAtom *one;
  std::vector<DumpAtom *> all;
  void readdump();

  void setrange();
  void writexyz();
  void voro(const int);

  void help();

  Memory *memory;
  int count_words(const char *);
};
#endif
