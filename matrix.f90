!> Defines some things that don't really belong anywhere else
module matrix_submodule
   public:: resContext
   type resContext
      integer:: opt
   end type resContext
end module matrix_submodule

!> Where the matrix gets built
module matrix
#include <petsc/finclude/petscsnes.h>
use petscsnes
implicit none

public :: init_matrix, matrix_size, finish_matrix, init_precomputes

PetscInt :: matrix_size

private

Vec:: stateVector, rhs, residualVector
Mat:: global_matrix

interface SNESSetApplicationContext
   subroutine SNESSetApplicationContext(snes_in,ctx,localerr)
      use petscsnes
      use matrix_submodule
      SNES:: snes_in
      type(resContext),intent(in):: ctx
      PetscErrorCode:: localerr
   end subroutine SNESSetApplicationContext
end interface

contains
   subroutine init_matrix()
      use mp
      use input
      use distribution
      use matrix_submodule
      implicit none
      integer,dimension(:),allocatable:: nNonZeros_offdiag, nNonZeros_ondiag
      integer:: firstLocalRow,lastLocalRow, i, j
      real:: rtol,atol,stol
      integer:: maxit, maxf
      PetscMPIInt:: jproc
      PetscInt:: cumRows
      type(resContext):: ctx
      KSP:: ksp
      PC:: pc
      PetscViewerAndFormat:: viewer

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

      ! If nproc < the resolution of the slowest-changing coordinate, off-diagonal
      ! will pick up more elements and parallelization won't be as effective.
      if (nproc < N_ordered(1)) then
         nNonZeros_offdiag(:) = nNonZeros_ondiag(:) 
      end if
      ! If nproc > twice the resoultion of the slowest-changing coordinate,
      ! then the local matrix will have no nonzero elements outsize
      ! the diagonal block. Otherwise, there will two additional elements.
      if (nproc < 2*N_ordered(1)) then
         nNonZeros_ondiag(:) = nNonZeros_ondiag(:) - 2
         nNonZeros_offdiag(:) = nNonZeros_offdiag(:) + 2
      end if

      ! Get how many rows this processor is responsible for
      localNrows = PETSC_DECIDE
      call PetscSplitOwnership(PETSC_COMM_WORLD, localNrows, matrix_size,ierr)

      ! Calculate the total number of rows up to and including this processor
      ! (i.e., from proc0 to iproc)
      cumRows = 0
      call MPI_Scan(localNrows,cumRows,1,MPI_INTEGER,MPI_SUM,mpicomm,ierr)
      if (ierr /=0 ) then
         print*,"MPI_Scan returned error code ", ierr
         stop
      end if

      ! Use the prefix sum above to calculate the row indices this processor is responsible for
      firstLocalRow = cumRows-localNrows+1
      lastLocalRow = cumRows

      ! Create the parallel matrix and share the rows
      call MatCreate(PETSC_COMM_WORLD,global_matrix,ierr)
      call MatSetType(global_matrix, MATAIJ, ierr)
      call PetscSplitOwnership(PETSC_COMM_WORLD, localNrows, matrix_size, ierr)
      ! Passing both localNrows and matrix_size is redundant, but doing so checks for inconsistencies
      call MatSetSizes(global_matrix, localNrows, localNrows, matrix_size, matrix_size,ierr)
      if (ierr /=0 ) then
         print*,"MatSetSizes returned error code ", ierr
         stop
      end if

      ! Preallocate the matrix based on the row-wise estimate of the total number of nonzero elements
      call MatMPIAIJSetPreallocation(global_matrix, 0, nNonZeros_ondiag(firstLocalRow:lastLocalRow), &
            0, nNonZeros_offdiag(firstLocalRow:lastLocalRow),ierr)
         
      ! Create state vector, RHS, and residual vector
      call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrix_size, stateVector, ierr)
      call VecDuplicate(stateVector, residualVector, ierr)
      call VecDuplicate(stateVector, rhs, ierr)

      ! Initialize the state vector
      call set_initDistribution(stateVector)

      ! Define context for residual functions
      ctx%opt=0

      call SNESSetFunction(snes, residualVector, residual_func, ctx, ierr)

      call SNESSetJacobian(snes, global_matrix, global_matrix, jacobian_func, ctx, ierr)

      call SNESGetKSP(snes,ksp,ierr)

      call KSPGetPC(ksp,pc,ierr)
      call PCSetType(pc,PCLU,ierr)
      if (solver_tol < 0.0) then

         print*, "ERROR: Direct solver not yet implemented"
         stop

      else

         call KSPSetType(ksp,KSPGMRES,ierr)
         call KSPSetTolerances(ksp,solver_tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
         call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,viewer,ierr)
         call KSPMonitorSet(ksp,KSPMonitorDefault,viewer,PetscViewerAndFormatDestroy,ierr)

      end if
      call KSPSetFromOptions(ksp,ierr)

      if (nonlinear) then
         call SNESSetType(snes, SNESNEWTONLS, ierr)
      else
         call SNESSetType(snes, SNESKSPONLY, ierr)
      end if

      deallocate(nNonZeros_offdiag)
      deallocate(nNonZeros_ondiag)

   end subroutine init_matrix

   subroutine init_precomputes()
      use input
      implicit none

   end subroutine

   subroutine residual_func(snes_in, stateVector, residualVector, context, localerr)
      use matrix_submodule
      implicit none
      SNES:: snes_in
      Vec,intent(in):: stateVector
      Vec,intent(inout):: residualVector
      type(resContext):: context
      PetscErrorCode:: localerr

   end subroutine residual_func

   subroutine jacobian_func(snes_in, stateVector, matrix_local, matrix_pc, context, localerr)
      use matrix_submodule
      implicit none
      SNES:: snes_in
      Vec,intent(in):: stateVector
      Mat:: matrix_local, matrix_pc
      type(resContext):: context
      PetscErrorCode:: localerr

   end subroutine jacobian_func

   subroutine finish_matrix()
      use mp
      implicit none

      call VecDestroy(stateVector,ierr)
      call VecDestroy(rhs,ierr)
      call VecDestroy(residualVector,ierr)
      call MatDestroy(global_matrix,ierr)

   end subroutine finish_matrix

end module matrix
