module source
   implicit none

   public:: init_source, sourcefunc

   private

   abstract interface
      function sourcefunc_iface(r,p,x,t)
         real,intent(in):: r,p,x,t
         real:: sourcefunc_iface
      end function sourcefunc_iface
   end interface

   procedure (sourcefunc_iface), pointer:: sourcefunc => null()

   contains 

   subroutine init_source()
      use input
      implicit none

      select case (source_type)
      case ("test1d")
         sourcefunc => source_test1d
      case default
         print*, "ERROR: source_type=",source_type," not recognized."
         stop
      end select

   end subroutine init_source

   real function source_test1d(r,p,x,t)
      implicit none
      real,intent(in)::r,p,x,t
      real:: fac
      fac = 0.75*(1.0-r**2) - 1.5*r**3 + exp(-r**2)*(r-1.0)
      source_test1d = 0.25*(1.0-r**2)*cos(t/3.0) + 3.0 + 3.0*sin(t/3.0)*fac
   end function source_test1d

end module source
