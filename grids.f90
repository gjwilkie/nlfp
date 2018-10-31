module grids
#include <petsc/finclude/petscsnes.h>
   use petscsnes
   implicit none

   public :: init_grids,rgrid,pgrid,xgrid,rgrid_edge,pgrid_edge,xgrid_edge
   public :: get_idx, get_idx_r, get_idx_p, get_idx_x

   real, dimension(:), allocatable:: rgrid, pgrid, xgrid
   real, dimension(:), allocatable:: rgrid_edge, pgrid_edge, xgrid_edge

   private
   abstract interface
      function get_global_index(i,j,k)
         integer,intent(in):: i,j,k
         integer:: get_global_index
      end function get_global_index
      function get_coord_index(idx)
         integer,intent(in):: idx
         integer:: get_coord_index
      end function get_coord_index
   end interface

   procedure (get_global_index), pointer:: get_idx => null()
   procedure (get_coord_index), pointer:: get_idx_r => null()
   procedure (get_coord_index), pointer:: get_idx_x => null()
   procedure (get_coord_index), pointer:: get_idx_p => null()
   contains

   subroutine init_grids()
      use input
      implicit none

      allocate(rgrid(Nr))
      allocate(pgrid(Np))
      allocate(xgrid(Nx))
      allocate(rgrid_edge(Nr+1))
      allocate(pgrid_edge(Np+1))
      allocate(xgrid_edge(Nx+1))

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

      if (xgrid_opt == "uniform") then
         call get_grid_uniform(Nx,-1.0,1.0,xgrid,xgrid_edge)
      else
         print*, "ERROR: xgrid_opt ",xgrid_opt," not valid"
         stop
      end if

      select case (layout)
      case ("xrp")
         get_idx => get_idx_xrp
         get_idx_r => idx_r_xrp
         get_idx_p => idx_p_xrp
         get_idx_x => idx_x_xrp
      case ("xpr")
         get_idx => get_idx_xpr
         get_idx_r => idx_r_xpr
         get_idx_p => idx_p_xpr
         get_idx_x => idx_x_xpr
      case ("rxp")
         get_idx => get_idx_rxp
         get_idx_r => idx_r_rxp
         get_idx_p => idx_p_rxp
         get_idx_x => idx_x_rxp
      case ("rpx")
         get_idx => get_idx_rpx
         get_idx_r => idx_r_rpx
         get_idx_p => idx_p_rpx
         get_idx_x => idx_x_rpx
      case ("pxr")
         get_idx => get_idx_pxr
         get_idx_r => idx_r_pxr
         get_idx_p => idx_p_pxr
         get_idx_x => idx_x_pxr
      case ("prx")
         get_idx => get_idx_prx
         get_idx_r => idx_r_prx
         get_idx_p => idx_p_prx
         get_idx_x => idx_x_prx
      case default
         print*, "ERROR: variable order ",layout," not recognized."
         stop
      end select

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

   !> Obtain global vector index from the coordinate indices.
   !! Arguments will always be in the order (r,p,x) regardless of layout.
   !! Only one from  the following 6 will be used depending on the
   !! user's chosen layout.
   integer function get_idx_xrp(ir,ip,ix)
      use input, only: Nr, Nx, Np
      implicit none
      integer, intent(in)::ir,ip,ix
      get_idx_xrp = (ix-1)*(Np*Nr) + (ir-1)*Np + ip
   end function get_idx_xrp
   integer function get_idx_xpr(ir,ip,ix)
      use input, only: Nr, Nx, Np
      implicit none
      integer, intent(in)::ir,ip,ix
      get_idx_xpr = (ix-1)*(Nr*Np) + (ip-1)*Nr + ir
   end function get_idx_xpr
   integer function get_idx_pxr(ir,ip,ix)
      use input, only: Nr, Nx, Np
      implicit none
      integer, intent(in)::ir,ip,ix
      get_idx_pxr = (ip-1)*(Nr*Nx) + (ix-1)*Nr + ir
   end function get_idx_pxr
   integer function get_idx_prx(ir,ip,ix)
      use input, only: Nr, Nx, Np
      implicit none
      integer, intent(in)::ir,ip,ix
      get_idx_prx = (ip-1)*(Nx*Nr) + (ir-1)*Nx + ix
   end function get_idx_prx
   integer function get_idx_rpx(ir,ip,ix)
      use input, only: Nr, Nx, Np
      implicit none
      integer, intent(in)::ir,ip,ix
      get_idx_rpx = (ir-1)*(Nx*Np) + (ip-1)*Nx + ix
   end function get_idx_rpx
   integer function get_idx_rxp(ir,ip,ix)
      use input, only: Nr, Nx, Np
      implicit none
      integer, intent(in)::ir,ip,ix
      get_idx_rxp = (ir-1)*(Np*Nx) + (ix-1)*Np + ip
   end function get_idx_rxp

   !> Obtain coordinate indices from the global vector index.
   !! Long list of unpleasant index-algebra functions. 
   !! only one from each set of 6 will be used depending on the
   !! user's chosen layout
   integer function idx_x_xrp(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_x_xrp= (idx-1)/(Nr*Np) + 1
   end function idx_x_xrp
   integer function idx_x_xpr(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_x_xpr=(idx-1)/(Nr*Np) + 1
   end function idx_x_xpr
   integer function idx_x_pxr(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_x_pxr=(mod(idx-1,Nr*Nx))/Nr + 1
   end function idx_x_pxr
   integer function idx_x_prx(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_x_prx= mod(mod(idx-1,Nx*Nr),Nx) + 1
   end function idx_x_prx
   integer function idx_x_rpx(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_x_rpx= mod(mod(idx-1,Np*Nx),Nx) + 1
   end function idx_x_rpx
   integer function idx_x_rxp(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_x_rxp=(mod(idx-1,Np*Nx))/Np + 1
   end function idx_x_rxp

   integer function idx_p_pxr(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_p_pxr= (idx-1)/(Nx*Nr) + 1
   end function idx_p_pxr
   integer function idx_p_prx(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_p_prx= (idx-1)/(Nx*Nr) + 1
   end function idx_p_prx
   integer function idx_p_rpx(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_p_rpx=(mod(idx-1,Nx*Np))/Nx + 1
   end function idx_p_rpx
   integer function idx_p_rxp(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_p_rxp=mod(mod(idx-1,Nx*Np),Np) + 1
   end function idx_p_rxp
   integer function idx_p_xrp(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_p_xrp=mod(mod(idx-1,Nr*Np),Np) + 1
   end function idx_p_xrp
   integer function idx_p_xpr(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_p_xpr=(mod(idx-1,Nr*Np))/Nr + 1
   end function idx_p_xpr

   integer function idx_r_rxp(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_r_rxp=(idx-1)/(Nx*Np) + 1
   end function idx_r_rxp
   integer function idx_r_rpx(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_r_rpx=(idx-1)/(Nx*Np) + 1
   end function idx_r_rpx
   integer function idx_r_prx(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_r_prx=(mod(idx-1,Nx*Nr))/Nx + 1
   end function idx_r_prx
   integer function idx_r_pxr(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_r_pxr=mod(mod(idx-1,Nx*Nr),Nr) + 1
   end function idx_r_pxr
   integer function idx_r_xpr(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_r_xpr=mod(mod(idx-1,Np*Nr),Nr) + 1
   end function idx_r_xpr
   integer function idx_r_xrp(idx)
      use input, only: Nr, Nx, Np
      implicit none
      integer,intent(in):: idx
      idx_r_xrp=(mod(idx-1,Np*Nr))/Np + 1
   end function idx_r_xrp




end module grids
   
