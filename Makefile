.SUFFIXES : .o .cpp
# compiler and flags
CC     = g++ -Wno-unused-result
LINK   = $(CC)
CFLAGS = -O3 $(UFLAG) $(DEBUG)
#
OFLAGS = -O3 $(DEBUG)
INC    = $(FFTINC) $(LPKINC) $(USRINC) $(VoroINC)
LIB    = $(FFTLIB) $(LPKLIB) $(USRLIB) $(VoroLIB)
#
# fftw 3 library; not needed by this code
#FFTINC    = -I/opt/fftw/fftw3/include
#FFTLIB    = -L/opt/fftw/fftw3/lib -lfftw3

# Lapack library; not needed by this code
#LPKINC = -I/opt/clapack/3.2.1/include
#LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c -lm

# Voro++
VoroINC = -I/opt/libs/voro_svn/src
VoroLIB = -L/opt/libs/voro_svn/src -lvoro++

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
all:  ${EXE}

${EXE}:  $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${BASE}.tar; tar -czvf ${BASE}.tar.gz *.cpp  *.h Makefile README

.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<
