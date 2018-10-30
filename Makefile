# Makefile for NFLP
# Make sure NLFP_SYSTEM is a defined environment variable

ifndef NLFP_SYSTEM
$(error Set the environment variable NLFP_SYSTEM)
endif
include Makefiles/Makefile.$(NLFP_SYSTEM)

OPTS = -g -cpp -ffree-line-length-none

.DEFAULT_GOAL := nlfp
nlfp: nlfp.o io.o mp.o
	$(FLINKER) -o nlfp nlfp.o mp.o io.o -I${PETSC_DIR}/include $(PETSC_LIB) $(NETCDF_INC) $(NETCDF_LIB)

nlfp.o: nlfp.f90  io.o mp.o
	$(FC) -c nlfp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o nlfp.o

io.o: io.f90  mp.o
	$(FC) -c io.f90 $(OPTS) $(PETSC_FC_INCLUDES) $(NETCDF_INC) $(NETCDF_LIB) -o io.o

mp.o: mp.f90 
	$(FC) -c mp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o mp.o

clean::
	rm -f nlfp *.mod *.o 
