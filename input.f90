module input
use mp
use netcdf
use constants
implicit none

public :: init_input
public :: Nr, Nxi, Np, rmin, rmax, rgrid_opt,pgrid_opt,xigrid_opt, pmax
public :: initial_condition,diffusion_type,source_type,geometry_type
public :: dt,Nt


!!!!! Input options
integer :: Nr !< Number of spatial grid points
integer :: Np !< Number of momentum magnitude grid points
integer :: Nxi !< Number of pitch angle grid points
real :: rmin, rmax !< Boundaries of spatial grid
real :: pmax !< Maximum momentum to simulate. Input as pmax_mc - in units of m_e*c
character(len=16) :: rgrid_opt, pgrid_opt, xigrid_opt !< Options specifying how spatial grid points are distributed
character(len=16) :: initial_condition !< Option specifying the initial condition
character(len=16) :: diffusion_type !< Option specifying the radial diffusion operator
character(len=16) :: source_type !< Option specifying the source
character(len=16) :: geometry_type !< Option specifying the geometry
real :: dt !< Time step size in seconds. Input may be otherwise normalized.
integer :: Nt !< Total number of computational timesteps.


private

real :: pmax_mc

contains

   !> Reads the namelist input file runname.in and sets the options accordingly.
   !! Also creates the output file and immediately writes all input parameters.
   subroutine init_input(runname)
      implicit none
      character(len=64), intent(in):: runname
      character(len=64) :: infile
      integer:: inunit = 11,  inputerr
      namelist /grid_params/ Nxi,Np,Nr,pmax_mc,rgrid_opt,xigrid_opt,pgrid_opt
      namelist /init_params/ initial_condition
      namelist /diffusion_params/ diffusion_type
      namelist /source_params/ source_type
      namelist /geometry_params/ geometry_type,rmin,rmax
      namelist /time_params/ dt,Nt

      ! Open input file
      infile = trim(runname)//".in"
      open(unit=inunit, file=infile, action="read", status="old", iostat=inputerr)
      if (inputerr /= 0) then
         write(*,*) "ERROR: could not open file ",infile
         stop
      end if

      call set_defaults()

      ! Read namelists. 
      read(inunit,nml=grid_params); rewind inunit
      read(inunit,nml=init_params); rewind inunit
      read(inunit,nml=source_params); rewind inunit
      read(inunit,nml=diffusion_params); rewind inunit
      read(inunit,nml=geometry_params); rewind inunit
      read(inunit,nml=time_params); rewind inunit
      close(inunit)

      call parse_inputs()


   end subroutine init_input

   subroutine set_defaults()
      implicit none

      Nr = 1; Np = 1; Nxi = 1
      pmax = 2.0
      pmax_mc = -1.0
      rgrid_opt = "uniform"; pgrid_opt = "uniform"; xigrid_opt = "uniform"
      initial_condition = "maxwellian"
      diffusion_type = "none"
      source_type = "none"
      geometry_type = "flat"
      rmin = 0.0
      rmax = 1.0
      dt = -1.0
      Nt = 0

   end subroutine set_defaults

   subroutine parse_inputs()
      implicit none

      ! Rescale normalized input
      if (pmax_mc > 0.0) then
         pmax = pmax_mc * me * c
      end if

      ! Sanity checks
      if ( (Nr == 1) .AND. (Np == 1) .AND. (Nxi == 1) ) then
         print*, "ERROR: No resolutions set. Nothing to do."
         stop
      end if
      if ( (Nr == 1) .AND. (diffusion_type == "none") ) then
         print*, "WARNING: Running with spatial dependence, but no diffusion"
      end if
      if ( (dt < 0.0) .AND. (Nt > 0) ) then
         print*, "ERROR: Running a finite number of timesteps, but no step size option has been defined."
         stop
      end if
      if ( (dt > 0.0) .AND. (Nt == 0) ) then
         print*, "WARNING: Time step has been set, but Nt=0, indicating a steady-state simulation. Are you sure?"
         stop
      end if
         
   end subroutine parse_inputs 

end module input
