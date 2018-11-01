module distribution
#include <petsc/finclude/petscsnes.h>
use petscsnes
   implicit none

   public:: initialdistribution_test1d, init_distribution

   private

   contains 

   subroutine init_distribution(f0)
      use input
      use mp
      implicit none
      integer:: firstLocalRow, lastLocalRow

      Vec,intent(inout)::f0
      select case (initial_condition)
      case ("test1d")
         call VecSet(f0,1.0,ierr)
      case default
         print*,"ERROR: initial_condition=",initial_condition," not recognized."
         stop
      end select

      call VecAssemblyBegin(f0,ierr)
      call VecAssemblyEnd(f0,ierr)

   end subroutine init_distribution

   real function initialdistribution_test1d(r,v,x)
      implicit none
      real,intent(in)::r,v,x

      initialdistribution_test1d = 1.0

   end function initialdistribution_test1d

end module distribution
