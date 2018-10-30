# Makefile for NFLP
# Make sure NLFP_SYSTEM is a defined environment variable

ifndef NLFP_SYSTEM
$(error Set the environment variable NLFP_SYSTEM)
endif
include Makefiles/Makefile.$(NLFP_SYSTEM)

GIT_MOD='"$(shell git diff-index HEAD)"'
ifneq ($(GIT_MOD),'""')
	GIT_MOD=Modified_from_
else
	GIT_MOD=
endif
GIT_HASH='"$(GIT_MOD)$(shell git rev-list HEAD -n 1)"'

OPTS = -g -cpp -ffree-line-length-none -DGIT_HASH=$(GIT_HASH)

.DEFAULT_GOAL := nlfp
nlfp: nlfp.o io.o mp.o constants.o
	$(FLINKER) -o nlfp nlfp.o mp.o io.o constants.o -I${PETSC_DIR}/include $(PETSC_LIB) $(NETCDF_INC) $(NETCDF_LIB)

nlfp.o: nlfp.f90  io.o mp.o
	$(FC) -c nlfp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o nlfp.o

io.o: io.f90 mp.o constants.o
	$(FC) -c io.f90 $(OPTS) $(PETSC_FC_INCLUDES) $(NETCDF_INC) $(NETCDF_LIB) -o io.o

mp.o: mp.f90 
	$(FC) -c mp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o mp.o

constants.o: constants.f90 
	$(FC) -c constants.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o constants.o

diffusion.o: diffusion.f90 
	$(FC) -c diffusion.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o diffusion.o

source.o: source.f90 
	$(FC) -c source.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o source.o

distribution.o: distribution.f90 
	$(FC) -c distribution.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o distribution.o

clean::
	rm -f nlfp *.mod *.o 
