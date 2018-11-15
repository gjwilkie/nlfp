module contexts
#include <petsc/finclude/petscsnes.h>
use mp

   public:: resContext, init_precomputes, precomp

   private 

   type resContext
      real:: time
      real:: delta_t
   end type resContext
   
   type(resContext):: precomp

   contains

   subroutine init_precomputes()
      use input
      use grids
      implicit none
      integer:: idx
      

      precomp%time = 0.0
      precomp%delta_t = 0.0
      ! Define context for residual functions

      ! For now this is empty. Premature optimization is the root of all evil.
      

   end subroutine

end module contexts


