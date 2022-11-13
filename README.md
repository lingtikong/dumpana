*DumpAna*

Code to analyse the atom style dump file from LAMMPS.

The dump file from LAMMPS must contain atomic id, type, and positions (in either
cartesian or fractional). The code is driven by menu, and therefore no manual is
provided.

For compilation, the following libraries are needed:
1) voro++
2) gsl
One should have these libraries installed before compilation and set the correct
path to them in the Makefile.

Once compiled, one can run this code as:
    dumpana [options] file
where `options` will change some of the performance of the code, and `file` should
be the file name of the LAMMPS dump file.

For help, one can invoke:
    dumpana -h

(R) LT Kong 2022
konglt@sjtu.edu.cn
