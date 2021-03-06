!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NLFP (Nonlinear Fokker-Planck Solver) !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Written by: George Wilkie (2018)
! Solves the 3D Fokker-Planck equation with a time-derivative nonlinearity on an advection term. Utilizes the PETSc library for
! parallelization of solution and matrix population.

!> The main wrapper for a single run of NLFP.
program NLFP
#include <petsc/finclude/petscsnes.h>
use petscsnes
!use mpi
use input, only: init_input
use output, only: init_output, write_initial_data, finish_output
use source, only: init_source
use diffusion, only: init_diffusion
use geometry, only: init_geometry
use matrix, only: init_matrix, finish_matrix, set_matrix_size, matrix_size
use contexts, only: init_precomputes
use mp, only: mp_end, mp_init, iproc, set_ownership
use grids, only: init_grids
implicit none
integer::  nargs, l
character(len=64):: runname

   ! Call basic initializations
   call mp_init()

   ! Get the run name from command line.
   nargs= command_argument_count() 
   if (nargs == 0) then
      print*, "ERROR: Must pass the runname as an argument."
      stop 
   end if
   call get_command_argument(1,runname)
   runname = trim(runname)
   l = len_trim(runname)
   if (l > 3) then
      if (runname(l-2:l) == ".in") then
         runname = runname(1:l-3)
      end if
   end if

   ! Read input file and set up output file
   call init_input(runname)
   if (iproc == 0) then
      call init_output(runname)
   end if

   call init_geometry()

   call init_grids()

   call set_matrix_size()

   call set_ownership(matrix_size)

   call init_source()

   call init_diffusion()

   call init_precomputes()

   call init_matrix()

   if (iproc == 0) then
      call write_initial_data()
   end if

   ! Finish output
   if (iproc == 0) then
      call finish_output()
   end if

   call finish_matrix()

   ! Finish mpi and PETSc
   call mp_end()

end program

