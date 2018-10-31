# Makefile for NFLP
# Make sure NLFP_SYSTEM is a defined environment variable
ifndef NLFP_SYSTEM
$(error Set the environment variable NLFP_SYSTEM)
endif
include Makefiles/Makefile.$(NLFP_SYSTEM)

OPTS = -g -Wall -cpp -fdefault-real-8 -ffree-line-length-none -DGIT_HASH=$(GIT_HASH)

GIT_MOD='"$(shell git diff-index HEAD)"'
ifneq ($(GIT_MOD),'""')
	GIT_MOD=Modified_from_
else
	GIT_MOD=
endif
GIT_HASH='"$(GIT_MOD)$(shell git rev-list HEAD -n 1)"'


.DEFAULT_GOAL := nlfp
nlfp: nlfp.o input.o output.o mp.o constants.o source.o diffusion.o distribution.o grids.o geometry.o matrix.o
	$(FLINKER) -o nlfp nlfp.o mp.o input.o output.o source.o diffusion.o distribution.o grids.o constants.o geometry.o matrix.o -I${PETSC_DIR}/include $(PETSC_LIB) $(NETCDF_INC) $(NETCDF_LIB)

nlfp.o: nlfp.f90  input.o mp.o grids.o output.o matrix.o source.o
	$(FC) -c nlfp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o nlfp.o

input.o: input.f90 mp.o constants.o
	$(FC) -c input.f90 $(OPTS) $(PETSC_FC_INCLUDES) $(NETCDF_INC) $(NETCDF_LIB) -o input.o

output.o: output.f90 mp.o constants.o input.o grids.o
	$(FC) -c output.f90 $(OPTS) $(PETSC_FC_INCLUDES) $(NETCDF_INC) $(NETCDF_LIB) -o output.o

mp.o: mp.f90 
	$(FC) -c mp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o mp.o

matrix.o: matrix.f90 diffusion.o geometry.o source.o mp.o
	$(FC) -c matrix.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o matrix.o

constants.o: constants.f90 
	$(FC) -c constants.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o constants.o

geometry.o: geometry.f90 input.o
	$(FC) -c geometry.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o geometry.o

grids.o: grids.f90 constants.o input.o
	$(FC) -c grids.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o grids.o

diffusion.o: diffusion.f90 
	$(FC) -c diffusion.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o diffusion.o

source.o: source.f90 input.o
	$(FC) -c source.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o source.o

distribution.o: distribution.f90 
	$(FC) -c distribution.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o distribution.o

tests::
	make -C tests/test_index -B
	exec tests/test_index/./test_index	

clean::
	rm -f nlfp *.mod *.o 
