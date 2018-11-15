module geometry
   implicit none

   public::Vprime, init_geometry, jacob

   private

   abstract interface

      function Vprime_iface(r)
         real,intent(in):: r
         real:: Vprime_iface
      end function Vprime_iface

   end interface

   procedure (Vprime_iface), pointer:: Vprime => null()

   contains

   subroutine init_geometry()
      use input
      implicit none

      select case (geometry_type)
      case("flat")
         Vprime => Vprime_flat
      case("cylindrical")
         Vprime => Vprime_cylindrical
      case default
         print*,"ERROR. ", geometry_type, " is an invalid geometry_type."
      end select

   end subroutine

   !> The Jacobian of our coordinates. ***Not to be confused with the Jacobian of the residual for the linear solver.***
   !! This is the phase space volume element for our chosen coordinates. 
   !! Leaving in terms of x ini case that gets changed.
   function jacob(r,p,x)
      use constants
      implicit none
      real, intent(in):: r,p,x
      real:: jacob
      jacob =  2*pi*Vprime(r)*p**2
      return
   end function

   function Vprime_flat(r)
      use input
      use constants
      implicit none
      real,intent(in):: r
      real:: Vprime_flat

      ! For cartesian geometry, treat rmaj as a cross-sectional radius
      Vprime_flat = pi*rmaj**2
      return
   end function Vprime_flat


   function Vprime_cylindrical(r)
      use input
      use constants
      implicit none
      real,intent(in):: r
      real:: Vprime_cylindrical

      Vprime_cylindrical = (4*pi*r)*(2*pi*rmaj)
      return
   end function Vprime_cylindrical

end module geometry
