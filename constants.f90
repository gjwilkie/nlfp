module constants
   implicit none

   public:: me, c, eps0, el

   real,parameter:: c = 2.99792458e8
   real,parameter:: me = 9.10938356e-31
   real,parameter:: eps0 = 8.85418782e-12
   real,parameter:: el = 1.60217662e-19
   real,parameter:: pi = acos(-1.0)
end module constants
