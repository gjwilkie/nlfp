module grids
#include <petsc/finclude/petscsnes.h>
   use petscsnes
   implicit none

   public :: init_grids,rgrid,pgrid,xigrid,rgrid_edge,pgrid_edge,xigrid_edge

   real, dimension(:), allocatable:: rgrid, pgrid, xigrid
   real, dimension(:), allocatable:: rgrid_edge, pgrid_edge, xigrid_edge

   private

   contains

   subroutine init_grids()
      use input
      implicit none

      allocate(rgrid(Nr))
      allocate(pgrid(Np))
      allocate(xigrid(Nxi))
      allocate(rgrid_edge(Nr+1))
      allocate(pgrid_edge(Np+1))
      allocate(xigrid_edge(Nxi+1))

      if (rgrid_opt == "uniform") then
         call get_grid_uniform(Nr,rmin,rmax,rgrid,rgrid_edge)
      else
         print*, "ERROR: rgrid_opt ",rgrid_opt," not valid"
         stop
      end if

      if (pgrid_opt == "uniform") then
         call get_grid_uniform(Np,0.0,pmax,pgrid,pgrid_edge)
      else
         print*, "ERROR: pgrid_opt ",pgrid_opt," not valid"
         stop
      end if

      if (xigrid_opt == "uniform") then
         call get_grid_uniform(Nxi,-1.0,1.0,xigrid,xigrid_edge)
      else
         print*, "ERROR: xigrid_opt ",xigrid_opt," not valid"
         stop
      end if

   end subroutine init_grids

   subroutine get_grid_uniform(Nx,xmin,xmax,xgrid,xgrid_edge)
      implicit none
      integer,intent(in):: Nx
      real,intent(in):: xmin,xmax
      real,dimension(:),intent(inout):: xgrid, xgrid_edge
      real:: dx
      integer:: i

      dx = (xmax-xmin)/Nx

      do i=1,Nx
         xgrid(i) = xmin + dx*(i-0.5)
         xgrid_edge(i) = xmin + dx*(i-1.0)
      end do
      xgrid_edge(Nx+1) = xmax
      
   end subroutine get_grid_uniform

end module grids
   
