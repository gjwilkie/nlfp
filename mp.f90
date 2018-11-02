!> Includes routines and wrappers needed globally to handle PETSc and MPI
module mp
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscsys.h>
use petscsnes
use mpi
   implicit none

   private

   public :: mp_init, mp_end, snes, mpicomm

   !> Integer representing the process index
   public :: iproc

   !> Integer representing the total number of processes
   public :: nproc

   !> The PETSc error code
   public :: ierr

   PetscMPIInt:: iproc,nproc
   PetscErrorCode:: ierr
   SNES:: snes
   MPI_Comm:: mpicomm

   contains

   subroutine mp_init()

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
!      call PetscInitialize("options",ierr)
      if (ierr .ne. 0) then
         print*, "Unable to initialize PETSc"
         stop
      endif

      mpicomm = PETSC_COMM_WORLD

      call SNESCreate(PETSC_COMM_WORLD,snes,ierr)
      if (ierr .ne. 0) then
         print*, "Unable to initialize SNES"
         stop
      endif

      call MPI_Comm_size(PETSC_COMM_WORLD,nproc,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,iproc,ierr)

   end subroutine mp_init

   subroutine mp_end()
      implicit none
      call SNESDestroy(snes,ierr)
      call PetscFinalize(ierr)
   end subroutine mp_end

end module mp
