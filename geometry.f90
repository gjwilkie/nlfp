module geometry
   implicit none

   public::Vprime

   private

   contains

   function Vprime(r)
      use input, only: geometry_type
      implicit none
      real,dimension(:):: r
      real,dimension(size(r)):: Vprime

      if (geometry_type == "flat") then
         Vprime(:) = 1.0   
      else if (geometry_type == "cylindrical") then
         Vprime(:) = r
      end if

      return
   end function Vprime

end module geometry
