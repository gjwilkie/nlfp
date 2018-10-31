!> The diffusion (and advection) coefficient module for nlfp.
!! Diffusion coefficients are defined at cell centers (half-grid-points, same as f)
!! Advection coefficients are defined at the center of the cell face in the relevant dimension. 
!! For example, advection coefficient into the cell centered at [x_(i+1/2),y_(j+1/2),z_(k+1/2)] 
!!   is defined at [x_i, y_(j+1/2), z_(k+1/2)].
!! This module manages all the coefficients in an abstract sense; includes collision operator
!! and Efield acceleration operators.
module diffusion
   implicit none


   public:: init_diffusion, Drr, Dpp, Dxx, Ar, Ap, Ax

   private

   abstract interface
      function coeff_iface(r,p,x,t)
         real,intent(in):: r,p,x,t
         real:: coeff_iface
      end function coeff_iface
   end interface

   procedure (coeff_iface), pointer:: Drr => null()
   procedure (coeff_iface), pointer:: Dpp => null()
   procedure (coeff_iface), pointer:: Dxx => null()
   procedure (coeff_iface), pointer:: Ar => null()
   procedure (coeff_iface), pointer:: Ap => null()
   procedure (coeff_iface), pointer:: Ax => null()

   contains 

   subroutine init_diffusion()
      use input
      implicit none

      select case (diffusion_type)
      case ("test1d")
         Drr => diffcoeff_test1d
         Ar => advcoeff_test1d
      case default
         print*, "ERROR: source_type=",source_type," not recognized."
         stop
      end select

      ! Fill these in with collision operator and electric field advection
         Ap => coeff_zero
         Ax => coeff_zero
         Dpp => coeff_zero
         Dxx => coeff_zero

   end subroutine init_diffusion

   real function diffcoeff_test1d(r,p,x,t)
      implicit none
      real,intent(in)::r,p,x,t
      diffcoeff_test1d = exp(-r)
   end function diffcoeff_test1d

   real function advcoeff_test1d(r,p,x,t)
      implicit none
      real,intent(in)::r,p,x,t
      advcoeff_test1d = r**2
   end function advcoeff_test1d

   real function coeff_zero(r,p,x,t)
      implicit none
      real,intent(in)::r,p,x,t
      coeff_zero = 0.0
   end function coeff_zero

end module diffusion
