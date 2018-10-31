# Makefile for NFLP
# Make sure NLFP_SYSTEM is a defined environment variable
ifndef NLFP_SYSTEM
$(error Set the environment variable NLFP_SYSTEM)
endif
include ../../Makefiles/Makefile.$(NLFP_SYSTEM)

OPTS = -g -Wall -cpp -ffree-line-length-none -DGIT_HASH=$(GIT_HASH)

GIT_MOD='"$(shell git diff-index HEAD)"'
ifneq ($(GIT_MOD),'""')
	GIT_MOD=Modified_from_
else
	GIT_MOD=
endif
GIT_HASH='"$(GIT_MOD)$(shell git rev-list HEAD -n 1)"'


.DEFAULT_GOAL := test_index
test_index: test_index.o grids.o input.o 
	$(FLINKER) -o test_index test_index.o grids.o input.o mp.o -I${PETSC_DIR}/include $(PETSC_LIB) $(NETCDF_INC) $(NETCDF_LIB)

test_index.o: test_index.f90 input.o grids.o
	$(FC) -c test_index.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o test_index.o

input.o: ../../input.f90 mp.o constants.o
	$(FC) -c ../../input.f90 $(OPTS) $(PETSC_FC_INCLUDES) $(NETCDF_INC) $(NETCDF_LIB) -o input.o

grids.o: ../../grids.f90 constants.o input.o
	$(FC) -c ../../grids.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o grids.o

mp.o: ../../mp.f90 
	$(FC) -c ../../mp.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o mp.o

constants.o: ../../constants.f90 
	$(FC) -c ../../constants.f90 $(OPTS) $(PETSC_FC_INCLUDES) -o constants.o

clean::
	rm -f test_index *.mod *.o 
