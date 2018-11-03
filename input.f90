module input
use mp
use netcdf
use constants
implicit none

public :: init_input
public :: Nr, Nx, Np, rmin, rmax, rgrid_opt,pgrid_opt,xgrid_opt, pmax, rmaj
public :: initial_condition,diffusion_type,source_type,geometry_type
public :: dt,Nt, layout, solver_tol
public :: set_resolutions, set_layout, set_gridopts_uniform
public :: efield_option, nonlinear, N_ordered

!!!!! Input options
integer :: Nr !< Number of spatial grid points
integer :: Np !< Number of momentum magnitude grid points
integer :: Nx !< Number of pitch angle grid points
real :: rmin, rmax !< Boundaries of spatial grid
real :: pmax !< Maxmum momentum to simulate. Input as pmax_mc - in units of m_e*c
character(len=16) :: rgrid_opt, pgrid_opt, xgrid_opt !< Options specifying how spatial grid points are distributed
character(len=16) :: initial_condition !< Option specifying the initial condition
character(len=16) :: diffusion_type !< Option specifying the radial diffusion operator
character(len=16) :: source_type !< Option specifying the source
character(len=16) :: geometry_type !< Option specifying the geometry
real :: dt !< Time step size in seconds. Input may be otherwise normalized.
integer :: Nt !< Total number of computational timesteps.
character(len=3) :: layout !< 3-character string specifying the order in which the coordinates are laid out in memory.
character(len=16) :: efield_option !< Option determining the behavior of the electric field: "none", "diffusive", "constant", "time"
integer,dimension(3) :: N_ordered !< The resolutions in each dimension according to the chosen layout
real :: rmaj !< Major radius or other normalizing quantity. Used in calculating volume integral.
real :: solver_tol !< Relative tolerance for the iterative solver, negative for direct solve

! Input-derived paramters
logical:: nonlinear !< Determines weather a nonlinear solve is needed

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
      namelist /grid_params/ Nx,Np,Nr,pmax_mc,rgrid_opt,xgrid_opt,pgrid_opt
      namelist /init_params/ initial_condition
      namelist /diffusion_params/ diffusion_type
      namelist /source_params/ source_type
      namelist /geometry_params/ geometry_type,rmin,rmax, rmaj
      namelist /time_params/ dt,Nt
      namelist /control_params/ layout, solver_tol
      namelist /physics_params/ efield_option

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
      read(inunit,nml=control_params); rewind inunit
      read(inunit,nml=physics_params); rewind inunit
      close(inunit)

      call process_inputs()


   end subroutine init_input

   subroutine set_defaults()
      implicit none

      Nr = 1; Np = 1; Nx = 1
      pmax = 2.0
      pmax_mc = -1.0
      rgrid_opt = "uniform"; pgrid_opt = "uniform"; xgrid_opt = "uniform"
      initial_condition = "maxwellian"
      diffusion_type = "none"
      source_type = "none"
      geometry_type = "flat"
      layout = "rxp"
      solver_tol = 1.0e-4
      rmin = 0.0
      rmax = 1.0
      dt = -1.0
      Nt = 0

   end subroutine set_defaults

   subroutine process_inputs()
      implicit none
      integer:: i

      ! Rescale normalized input
      if (pmax_mc > 0.0) then
         pmax = pmax_mc * me * c
      end if

      ! Set N_ordered according to layout
      do i =1,3
         select case (layout(i:i))
         case ("r")
            N_ordered(i) = Nr
         case ("x")
            N_ordered(i) = Nx
         case ("p")
            N_ordered(i) = Np
         case default
            print*,"ERROR: layout character ",layout(i:i)," not recognized."
            stop
         end select
      end do

      ! Sanity checks
      if ( (Nr == 1) .AND. (Np == 1) .AND. (Nx == 1) ) then
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
         
   end subroutine process_inputs 

   subroutine set_layout(layout_in)
      implicit none
      character(len=3), intent(in):: layout_in
      layout = layout_in
   end subroutine set_layout

   subroutine set_gridopts_uniform()
      implicit none
      rgrid_opt="uniform"; pgrid_opt="uniform"; xgrid_opt="uniform"
   end subroutine set_gridopts_uniform

   subroutine set_resolutions(Nr_in,Np_in,Nx_in)
      implicit none
      integer, intent(in):: Nr_in, Np_in, Nx_in
      Nr=Nr_in; Np=Np_in; Nx=Nx_in
   end subroutine set_resolutions

end module input
