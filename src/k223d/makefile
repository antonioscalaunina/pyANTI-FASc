# Makefile
# #
#FC = ifort     # Fortran 90 compiler
FC = gfortran     # Fortran 90 compiler


OBJ = forparse.o utils.o lateration.o typedef.o makepdf.o k223d.o

#FCFLAGS= -O3 -r8 -assume byterecl
#FCFLAGS= -O2
#FCFLAGS= -O -fdefault-real-8 -g
FCFLAGS= -O -g -fbounds-check -fcheck=all -fbacktrace  #-ffpe-trap=invalid,zero,overflow -g -Wall -Wextra -Werror -pedantic
#FCFLAGS= -O -fdefault-real-8 -g -Warray-temporaries
#FCFLAGS = -O -check all -g -heap-arrays -traceback -assume byterecl
#FCFLAGS += -I/usr/local/include
#FCFLAGS+=-I/soft/intel/Compiler/2013.1.117/mkl/include/fftw

LIB = #-lfftw3

#############################
.SUFFIXES: .f90

all: k223d.x

############################
# Objects
############################

forparse.o: forparse.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

utils.o: utils.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

typedef.o: typedef.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

lateration.o: lateration.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

makepdf.o: makepdf.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

k223d.o : k223d.f90
	$(FC) $(FCFLAGS) -c $<

k223d.x :$(OBJ)
	$(FC)  $(FFLAGS) -o $@ $(OBJ)


clean:
	rm *.o *.mod
