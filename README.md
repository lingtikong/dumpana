**DumpAna**

Code to analyse the atom style dump file from LAMMPS. It relies heavily on the
voro++ code to determine the nearest neighbors, but for most functions, neighbors
can also be defined based on cutoff distances. For a complete function of this
code, please check by:

    dumpana -h

The dump file from LAMMPS must contain atomic id, type, and positions (in either
cartesian or fractional). The code is driven by menu, and therefore no manual is
provided.

For compilation, the following libraries are needed:
1) voro++ (https://math.lbl.gov/voro++/)
2) gsl    (https://www.gnu.org/software/gsl/)

One should have these libraries installed before compilation and set the correct
path to them in the Makefile.

Once compiled, one can run this code as:

    dumpana [options] file

where `options` will change some of the performance of the code, and `file` should
be the file name of the LAMMPS dump file. For available options or help, please run:

    dumpana -h

Any bug found, please report to konglt@sjtu.edu.cn, thanks.

(R) LT Kong 2022
konglt@sjtu.edu.cn
