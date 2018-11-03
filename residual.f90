subroutine residual_func(snes_in, stateVector, residualVector, context, localerr)
#include <petsc/finclude/petscsnes.h>
   use petscsnes
   use contexts
   implicit none
   SNES:: snes_in
   Vec,intent(in):: stateVector
   Vec,intent(inout):: residualVector
   type(resContext):: context
   PetscErrorCode:: localerr

end subroutine residual_func

subroutine jacobian_func(snes_in, stateVector, matrix_local, matrix_pc, context, localerr)
#include <petsc/finclude/petscsnes.h>
   use petscsnes
   use contexts
   implicit none
   SNES:: snes_in
   Vec,intent(in):: stateVector
   Mat:: matrix_local, matrix_pc
   type(resContext):: context
   PetscErrorCode:: localerr

end subroutine jacobian_func

