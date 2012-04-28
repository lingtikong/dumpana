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
  char *dump;
  int nframe, img;
  int istr, iend, inc;
  DumpAtom *one;
  std::vector<DumpAtom *> all;
  void readdump();

  void setrange();
  void writexyz();
  void voro();

  void help();

  Memory *memory;
  int count_words(const char *);
};
#endif
