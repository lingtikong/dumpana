.SUFFIXES : .o .cpp
# compiler and flags
#CC     = /opt/intel/bin/icc
CC     = g++ -Wno-unused-result
LINK   = $(CC) ${PARALIB}
CFLAGS = -O3 $(UFLAG) $(DEBUG)
#
OFLAGS = -O3 $(DEBUG)
INC    = $(FFTINC) $(LPKINC) $(USRINC) $(VoroINC) $(GslINC)
LIB    = $(FFTLIB) $(LPKLIB) $(USRLIB) $(VoroLIB) $(GslLIB)
#
# fftw 3 library; not needed by this code
#FFTINC    = -I/opt/fftw/fftw3/include
#FFTLIB    = -L/opt/fftw/fftw3/lib -lfftw3

# Lapack library; not needed by this code
#LPKINC = -I/opt/clapack/3.2.1/include
#LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c -lm

# Voro++, Needed.
VoroINC = -I/opt/voro/include/voro++
VoroLIB = -L/opt/voro/lib -lvoro++

# GSL, needed
GslINC  = -I/opt/gsl/include
GslLIB  = -L/opt/gsl/lib -lgsl -lgslcblas

# Parallization related, can be switch off
PARAINC = -DOMP -fopenmp
PARALIB = -fopenmp

# User flag
#UFLAG =
# Debug flags
#DEBUG = -g -O1
#====================================================================
# executable name
BASE   = dumpana
EXE    = ${BASE}

#================= Do not modify the following ======================
# source and rules
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
#====================================================================
all: ver ${EXE}

${EXE}:  $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${BASE}.tar; tar -czvf ${BASE}.tar.gz *.cpp  *.h Makefile README

ver:
	@echo "#define VERSION `git log|grep commit|wc -l`" > version.h; cat version.h

.f.o:
	$(FC) $(FFLAGS) $(FREE) $(PARAINC) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(PARAINC) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) $(PARAINC) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(PARAINC) $(INC) -c $<
