module input
use mp
use netcdf
use constants
implicit none

public :: init_input, finish_input, check, ncid, runname
public :: p_dim, r_dim, xi_dim
public :: Nr, Nxi, Np, rmin, rmax, rgrid_opt,pgrid_opt,xigrid_opt, pmax
public :: initial_conditinputn,diffusinputn_type,source_type,geometry_type
public :: dt,Nt


!!!!! Input optinputns
integer :: Nr !< Number of spatial grid points
integer :: Np !< Number of momentum magnitude grid points
integer :: Nxi !< Number of pitch angle grid points
real :: rmin, rmax !< Boundaries of spatial grid
real :: pmax !< Maximum momentum to simulate. Input as pmax_mc - in units of m_e*c
character(len=16) :: rgrid_opt, pgrid_opt, xigrid_opt !< Optinputns specifying how spatial grid points are distributed
character(len=16) :: initial_conditinputn !< Optinputn specifying the initial conditinputn
character(len=16) :: diffusinputn_type !< Optinputn specifying the radial diffusinputn operator
character(len=16) :: source_type !< Optinputn specifying the source
character(len=16) :: geometry_type !< Optinputn specifying the geometry
real :: dt !< Time step size in seconds. Input may be otherwise normalized.

integer :: ncid !> NetCDF control id for output file
integer :: p_dim, r_dim, xi_dim !> NetCDF dimensinputn ids
character(len=64) runname !< The user-specified string which refers to the input and output files

private


contains

   !> Reads the namelist input file runname.in and sets the optinputns accordingly.
   !! Also creates the output file and immediately writes all input parameters.
   subroutine init_input(runname)
      implicit none
      character(len=64), intent(in):: runname
      character(len=64) :: infile
      namelist /grid_params/ Nxi,Np,Nr,pmax_mc,rgrid_opt,xigrid_opt,pgrid_opt
      namelist /init_params/ initial_conditinputn
      namelist /diffusinputn_params/ diffusinputn_type
      namelist /source_params/ source_type
      namelist /geometry_params/ geometry_type,rmin,rmax
      namelist /time_params/ dt,Nt
      integer:: inunit = 11, l, nargs, inputerr
      real:: pmax_mc

      ! Open input file
      infile = trim(runname)//".in"
      open(unit=inunit, file=infile, actinputn="read", status="old", inputstat=inputerr)
      if (inputerr /= 0) then
         write(*,*) "ERROR: could not open file ",infile
         stop
      end if

      ! Read namelists. 
      read(inunit,nml=grid_params); rewind inunit
      read(inunit,nml=init_params); rewind inunit
      read(inunit,nml=source_params); rewind inunit
      read(inunit,nml=diffusinputn_params); rewind inunit
      read(inunit,nml=geometry_params); rewind inunit
      read(inunit,nml=time_params); rewind inunit
      close(inunit)

      ! Rescale normalized input
      pmax = pmax_mc * me*c


   end subroutine init_input

   subroutine finish_input()
         implicit none
         call check(nf90_close(ncid))
   end subroutine finish_input

   subroutine check(flag)
      integer, intent(in):: flag
      if (flag /= nf90_noerr) then
              print*, trim(nf90_strerror(flag))
              stop 
      end if
   end subroutine check

end module input
