module distribution
#include <petsc/finclude/petscsnes.h>
use petscsnes
   implicit none

   public:: set_initDistribution

   private

   contains 

   subroutine set_initDistribution(f0)
      use input
      use mp
      implicit none
      integer:: firstLocalRow, lastLocalRow
      Vec,intent(inout)::f0

      ! Set the initial distribution
      select case (initial_condition)
      case ("zero")
         call VecSet(f0,0.0,ierr)
      case ("test1d")
         call VecSet(f0,1.0,ierr)
      case default
         print*,"ERROR: initial_condition=",initial_condition," not recognized."
         stop
      end select

      call VecAssemblyBegin(f0,ierr)
      call VecAssemblyEnd(f0,ierr)

   end subroutine set_initDistribution

end module distribution
