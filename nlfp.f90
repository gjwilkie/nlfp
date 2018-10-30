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
use mpi
use io, only: init_io, finish_io
use mp, only: mp_end, mp_init, iproc, nproc
implicit none
integer:: i, nargs, l
character(len=64):: arg

   ! Call basic initializations
   call mp_init()

   ! Get the run name from command line.
   nargs= command_argument_count() 
   if (nargs == 0) then
      print*, "ERROR: Must pass the runname as an argument."
      stop 
   end if
   call get_command_argument(1,arg)
   arg = trim(arg)
   l = len_trim(arg)
   if (l > 3) then
      if (arg(l-2:l) == ".in") then
         arg = arg(1:l-3)
      end if
   end if

   call init_io(arg)

   if (iproc == 0) then
      call finish_io()
   end if

   call mp_end()

end program

