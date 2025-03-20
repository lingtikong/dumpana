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
#include "input.h"


using namespace std;

class Driver {
public:
  Driver(int, char**);
  ~Driver();

private:
  char *dump;                             // dump file name; should carries the last file name
  int flag_dump;                          // flag related to the reading of dump
  int min_mem;                            // Flag, whether to minimize memory usage
  int nframe;                             // total # of frames
  int istr, iend, inc;                    // frame range
  DumpAtom *one;                          // pointer to one frame
  UserInput *input;                       // pointer to user input
  vector<DumpAtom *> all;                 // all frames from lammps atom style dump
  void readdump(const int,int,char **);   // to read in the dump file
  void guess_image();                     // to guess the image info, if not supplied.
  void MainMenu();

  // private method to set the cutoff for face or edge, needed by most others
  double mins[3];
  void set_cutoffs(int);

  // private method to set the cutoff for bond distances for each pair
  double ***r2cuts;
  void set_r2cuts();

  // choose the method to define nearest neighbors
  int neighbor_method;
  void choose_neighbor_method(int);

  // map atomic type to elements
  ChemElements *element;
  int *type2atnum, weighted;
  double *type2radius;
  void MapType2Elem(const int, const int);
  void ShowRadius4Voro();

  // to define the range of frames to be analysed
  void setrange();
  int nsel;

  // to output selected frames
  void writesel();

  // to unwrap so bonded atoms will be in correct positions
  void unwrap();
  void unwrap_neighbors(int, int *);

  // to output selected frames in bgf format with properties
  void writebgf();

  // to do voronoi diagram analysis
  int flag_out;
  void voro();

  // to output selected Voronoi cells as xyz file
  void OutputVoroCells();

  // Chemical short range order based on voronoi neighbors
  void csro();

  // to average up all selected frames
  void avedump();
  // dump with selection info as extra per atom property
  void DumpSelection();

  // Hondeycutt-Andersen bond type analysis
  void honeycutt_andersen();
  void count_HA(const int, const int, FILE *, const int);

  // Common neighber analysis or common neighbor parameter
  void Compute_CNACNP();

  // Bond order parameter based on Spherical Harmonics
  void compute_sh();

  // Pair Correlation and static structure factor
  void paircorr();
  void write_gr(double, double, int, double **, int, char *, char *);
  void strfac();

  // Property Pair Correlation
  void property_pc();
  void write_prop_pc(double, double, int, double **, int, char *, char *);

  // Average distance for pais of same Property
  void ave_dist_same_property();

  // Orientation for pais of same Property
  void orient_same_property();

  // Prepare for FEFF9
  void FEFF_main();
  void FEFF_input(int, FILE *);

  // prepare for RINGS
  void rings();

  // Connectivity info for certain clusters
  void ClusterConnectivity();
  void IterateOverConn(int, const int, int);
  map<int,int> nconn, SCid;
  map<bigint,int> conns;
  set<int> checked;

  // bond length or bond angle distribution; Voronoi neighbors are seen as bonded
  void bonds();

  // spatial distribution of atoms
  void spatial();

  // radial distribution of atoms
  void radial();

  // compare the rmsd between frames
  void compare_rmsd();

  // To compute the Bhatia-Thornton structure factor
  void bhatia_thornton();

  // Compute MSD for atoms in selection; image info from dump file is needed;
  // besides, the selection will be intersected among all frames
  void compute_msd();

  // Compute the configurational mixing entropy
  void compute_smix();

  // Evaluate the heredity of cluster between consecutive frames
  void heredity();
  void check_neigh(std::set<int>, DumpAtom *, DumpAtom *, std::list<int> *, int *);

  // Count atoms within selection wrt tstep
  void count_selected();

  // help info
  void help();
  void ShowVersion();

  Memory *memory;
  int count_words(const char *);
  int flag_overwrite;
  void ConfirmOverwrite(char *);
};
#endif
