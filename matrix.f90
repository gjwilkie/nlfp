!> Where the matrix gets built
module matrix
#include <petsc/finclude/petscsnes.h>
use petscsnes
implicit none

public :: init_matrix, matrix_size, build_matrix

PetscInt :: matrix_size

private

!Vec:: stateVector, rhs, residualVector
!Mat:: global_matrix

contains
   subroutine init_matrix()
      use mp
      use input
      use distribution
      implicit none
      integer,dimension(:),allocatable:: nNonZeros_offdiag, nNonZeros_ondiag
      integer:: firstLocalRow,lastLocalRow, localNrows, jproc, i, j, cumRows_temp!, cumRows
      PetscInt,dimension(:),allocatable:: nRowsAll
      PetscInt:: cumRows

      ! Determine size of matrix
      matrix_size = Nr*Np*Nx
      if (efield_option == "diffusive") then
         matrix_size = matrix_size + Nr
      end if

      ! Predict the number of nonzeros for preallocation
      allocate(nNonZeros_offdiag(matrix_size))
      allocate(nNonZeros_ondiag(matrix_size))

      nNonZeros_offdiag(:) = 0
      nNonZeros_ondiag(:) = 1
      if (Nr > 1) then
         nNonZeros_ondiag(:) = nNonZeros_ondiag(:) + 2
      end if
      if (Nx > 1) then
         nNonZeros_ondiag(:) = nNonZeros_ondiag(:) + 2
      end if
      if (Np > 1) then
         nNonZeros_ondiag(:) = nNonZeros_ondiag(:) + 2
      end if

      ! Extra rows in electric field calculation are dense
      if (efield_option == "diffusive") then
         nNonZeros_ondiag(:) = nNonZeros_ondiag(:) + 5
         nNonZeros_offdiag(matrix_size-Nr+1:matrix_size) = matrix_size
      end if

      localNrows = PETSC_DECIDE
      call PetscSplitOwnership(PETSC_COMM_WORLD, localNrows, matrix_size,ierr)

      allocate(nRowsAll(nproc))
      nRowsAll(:) = 0.0;
      nRowsAll(iproc+1) = localNrows;
      
      call MPI_Bcast(nRowsAll,nproc,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
      call MPI_Barrier(PETSC_COMM_WORLD,ierr)
      call MPI_Scan(nRowsAll,cumRows,nproc,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
      if (ierr /=0 ) then
         print*,"MPI_Scan returned error code ", ierr
      end if


      ! Create state vectors
!      call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrix_size, stateVector, ierr)
!      call VecDuplicate(stateVector, residualVector, ierr)
!      call VecDuplicate(stateVector, rhs, ierr)

!      call init_distribution(stateVector)

      deallocate(nRowsAll)
      deallocate(nNonZeros_offdiag)
      deallocate(nNonZeros_ondiag)

   end subroutine init_matrix

   subroutine build_matrix()
      use input
      implicit none

   end subroutine

end module matrix
