!> Where the matrix gets built
module matrix
#include <petsc/finclude/petscsnes.h>
use petscsnes
implicit none

public :: init_matrix, matrix_size, build_matrix

PetscInt :: matrix_size

private

Vec:: stateVector, rhs, residualVector
Mat:: global_matrix

contains
   subroutine init_matrix()
      use mp
      use input
      use distribution
      implicit none
      integer,dimension(:),allocatable:: nNonZeros_offdiag, nNonZeros_ondiag
      integer:: firstLocalRow,lastLocalRow, localNrows, globalNrows,i, j, cumRows_temp,localguess!, cumRows
      PetscMPIInt:: jproc
      PetscInt:: cumRows

      ! Determine size of matrix
      matrix_size = Nr*Np*Nx
      if (efield_option == "diffusive") then
         matrix_size = matrix_size + Nr
      end if

      ! Predict the number of nonzeros for preallocation
      allocate(nNonZeros_offdiag(matrix_size))
      allocate(nNonZeros_ondiag(matrix_size))

      ! First guess, to be adjusted below
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

      ! If it's less than the resolution, off-diagonal will pick up more elements still
      ! and parallelization won't be as effective.
      if (nproc < N_ordered(1)) then
         nNonZeros_offdiag(:) = nNonZeros_ondiag(:) 
      end if
      ! If nproc > twice the resoultion of the slowest-changing coordinate,
      ! then the local matrix will have no nonzero elements outsize
      ! the diagonal block. 
      if (nproc < 2*N_ordered(1)) then
         nNonZeros_ondiag(:) = nNonZeros_ondiag(:) - 2
         nNonZeros_offdiag(:) = nNonZeros_offdiag(:) + 2
      end if


      localNrows = PETSC_DECIDE
      call PetscSplitOwnership(PETSC_COMM_WORLD, localNrows, matrix_size,ierr)

      cumRows = 0
      
      call MPI_Scan(localNrows,cumRows,1,MPI_INTEGER,MPI_SUM,mpicomm,ierr)
      if (ierr /=0 ) then
         print*,"MPI_Scan returned error code ", ierr
      end if

      firstLocalRow = cumRows-localNrows+1
      lastLocalRow = cumRows

      call MatCreate(PETSC_COMM_WORLD,global_matrix,ierr)
      call MatSetType(global_matrix, MATAIJ, ierr)
      call PetscSplitOwnership(PETSC_COMM_WORLD, localNrows, matrix_size, ierr)
      call MatSetSizes(global_matrix, localNrows, localNrows, matrix_size, matrix_size,ierr)


      call MatMPIAIJSetPreallocation(global_matrix, 0, nNonZeros_ondiag(firstLocalRow:lastLocalRow), &
            0, nNonZeros_offdiag(firstLocalRow:lastLocalRow),ierr)
         

      ! Create state vectors
!      call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrix_size, stateVector, ierr)
!      call VecDuplicate(stateVector, residualVector, ierr)
!      call VecDuplicate(stateVector, rhs, ierr)

!      call init_distribution(stateVector)


      deallocate(nNonZeros_offdiag)
      deallocate(nNonZeros_ondiag)

   end subroutine init_matrix

   subroutine build_matrix()
      use input
      implicit none

   end subroutine

end module matrix
