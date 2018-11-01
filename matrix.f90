!> Where the matrix gets built
module matrix
#include <petsc/finclude/petscsnes.h>
use petscsnes
implicit none

public :: init_matrix, matrix_size, build_matrix

PetscInt :: matrix_size

private

Vec:: stateVector, rhs, residualVector

contains
   subroutine init_matrix()
      use mp
      use input
      use distribution
      implicit none

      ! Determine size of matrix
      matrix_size = Nr*Np*Nx
      if (efield_option == "diffusive") then
         matrix_size = matrix_size + Nr
      end if

      ! Create state vectors
      call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrix_size, stateVector, ierr)
      call VecDuplicate(stateVector, residualVector, ierr)
      call VecDuplicate(stateVector, rhs, ierr)

      call init_distribution(stateVector)

   end subroutine init_matrix

   subroutine build_matrix()
      use input
      implicit none

   end subroutine

end module matrix
