module geometry
   implicit none

   public::Vprime

   private

   contains

   function Vprime(r)
      use input, only: geometry_type, rmin, rmax, rmaj
      use constants
      implicit none
      real:: r
      real:: Vprime

      ! For cartesian geometry, treat rmaj as a cross-sectional radius
      if (geometry_type == "flat") then
         Vprime = pi*rmaj**2
      else if (geometry_type == "cylindrical") then
         Vprime = (4*pi*r)*(2*pi*rmaj)
      end if

      return
   end function Vprime


end module geometry
