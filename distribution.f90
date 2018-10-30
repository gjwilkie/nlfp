module distribution
   implicit none

   private

   public:: initialdistribution_test1d

   contains 

   real function initialdistribution_test1d(r,v,xi)
      implicit none
      real,intent(in)::r,v,xi

      initialdistribution_test1d = 1.0

   end function initialdistribution_test1d

end module distribution
