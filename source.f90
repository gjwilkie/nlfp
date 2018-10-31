module source
   implicit none

   private

   public:: source_test1d

   contains 

   real function source_test1d(r,t)
      implicit none
      real,intent(in)::r,t
      real:: fac
      fac = 0.75*(1.0-r**2) - 1.5*r**3 + exp(-r**2)*(r-1.0)
      source_test1d = 0.25*(1.0-r**2)*cos(t/3.0) + 3.0 + 3.0*sin(t/3.0)*fac
   end function source_test1d

end module source
