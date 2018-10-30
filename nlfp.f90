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
!use nlfp_main
use io, only: read_input, close_output, init_output
use mp, only: mp_end, mp_init, iproc, nproc
implicit none
integer:: i, nargs, l
character(len=100)::arg, runname

   call mp_init()

   ! Get the run name from command line.
   nargs= command_argument_count() 
   if (nargs == 0) then
      write(*,*) "ERROR: Must pass the runname as an argument."
      stop
   end if
   call get_command_argument(1,arg)
   runname = trim(arg)
   l = len_trim(runname)
   if (l > 3) then
      if (runname(l-2:l) == ".in") then
         runname = runname(1:l-3)
      end if
   end if

   call read_input(runname)

   if (iproc == 0) then
      call init_output()
   end if

   if (iproc == 0) then
      call close_output()
   end if

   call mp_end()
end program

