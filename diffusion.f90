! The diffusion (and advection) coefficient module for nlfp.
! Diffusion coefficients are defined at cell centers (half-grid-points, same as f)
! Advection coefficients are defined at the center of the cell face in the relevant dimension. 
! For example, flux_x into the cell centered at [x_(i+1/2),y_(j+1/2),z_(k+1/2)] 
!   is defined at [x_i, y_(j+1/2), z_(k+1/2)].
module diffusion
   implicit none

   private

   public:: diffcoeff_test1d
   public:: advcoeff_test1d

   contains 

   real function diffcoeff_test1d(r,v,x)
      implicit none
      real,intent(in)::r,v,x
      diffcoeff_test1d = exp(-r)
   end function diffcoeff_test1d

   real function advcoeff_test1d(r,v,x)
      implicit none
      real,intent(in)::r,v,x
      advcoeff_test1d = r**2
   end function advcoeff_test1d

end module diffusion
