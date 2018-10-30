# Makefile for NFLP
# Make sure NLFP_SYSTEM is a defined environment variable

ifndef NLFP_SYSTEM
$(error Set the environment variable NLFP_SYSTEM)
endif
include Makefiles/Makefile.$(NLFP_SYSTEM)

OPTS = -g -cpp -ffree-line-length-none

.DEFAULT_GOAL := nlfp
nlfp: nlfp.o inputoutput.o mp.o
	$(FLINKER) -o nlfp nlfp.o mp.o inputoutput.o -I${PETSC_DIR}/include $(PETSC_LIB) $(NETCDF_INC) $(NETCDF_LIB)

nlfp.o: nlfp.f90  inputoutput.o mp.o
	$(FC) -c nlfp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o nlfp.o

inputoutput.o: inputoutput.f90  mp.o
	$(FC) -c inputoutput.f90 $(OPTS) $(PETSC_FC_INCLUDES) $(NETCDF_INC) $(NETCDF_LIB) -o inputoutput.o

mp.o: mp.f90 
	$(FC) -c mp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o mp.o

clean::
	rm -f nlfp *.mod *.o 
