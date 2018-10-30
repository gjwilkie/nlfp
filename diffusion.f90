module diffusion
   implicit none

   private

   public:: diffcoeff_test1d
   public:: advcoeff_test1d

   contains 

   real function diffcoeff_test1d(r,v,xi)
      implicit none
      real,intent(in)::r,v,xi
      diffcoeff_test1d = exp(-r)
   end function diffcoeff_test1d

   real function advcoeff_test1d(r,v,xi)
      implicit none
      real,intent(in)::r,v,xi
      advcoeff_test1d = r**2
   end function advcoeff_test1d

end module diffusion
