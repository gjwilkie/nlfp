module io
use mp
use netcdf
use constants
implicit none

public :: init_io, finish_io
public :: Nr, Nxi, Np, r_0, r_end, rgrid_opt

private

!!!!! Input options
integer :: Nr !< Number of spatial grid points
integer :: Np !< Number of momentum magnitude grid points
integer :: Nxi !< Number of pitch angle grid points
real :: r_0, r_end !< Boundaries of spatial grid
real :: pmax !< Maximum momentum to simulate. Input as pmax_mc - in units of m_e*c
character(len=16) :: rgrid_opt, pgrid_opt, xigrid_opt !< Options specifying how spatial grid points are distributed
character(len=16) :: initial_condition !< Option specifying the initial condition
character(len=16) :: diffusion_type !< Option specifying the radial diffusion operator
character(len=16) :: source_type !< Option specifying the source


!!!!! Module-internal options
integer :: ncid !> NetCDF control id for output file
character(len=64) runname !< The user-specified string which refers to the input and output files

contains

   !> Sets the runname 
   subroutine set_runname(arg)
      implicit none
      character(len=64), intent(in):: arg

      runname = arg
     
   end subroutine set_runname

   !> Reads the namelist input file runname.in and sets the options accordingly.
   !! Also creates the output file and immediately writes all input parameters.
   subroutine init_io(arg)
      implicit none
      character(len=64), intent(in):: arg
      character(len=64) :: infile
      namelist /grid_params/ Nxi,Np,Nr,pmax_mc,r_0,r_end,rgrid_opt,xigrid_opt,pgrid_opt
      namelist /init_params/ initial_condition
      namelist /diffusion_params/ diffusion_type
      namelist /source_params/ source_type
      integer:: inunit = 11, l, nargs, ioerr
      real:: pmax_mc

      call set_runname(arg)

      ! Open input file
      infile = trim(runname)//".in"
      open(unit=inunit, file=infile, action="read", status="old", iostat=ioerr)
      if (ioerr /= 0) then
         write(*,*) "ERROR: could not open file ",infile
         stop
      end if

      ! Read namelists. 
      read(inunit,nml=grid_params); rewind inunit
      read(inunit,nml=init_params); rewind inunit
      read(inunit,nml=source_params); rewind inunit
      read(inunit,nml=diffusion_params); rewind inunit
      close(inunit)

      ! Rescale normalized input
      pmax = pmax_mc * me*c

      ! Initialize output
      if (iproc == 0) then
         call init_output()
      end if

   end subroutine init_io

   subroutine init_output()
      implicit none
      character(len=64) :: outfile, infile, git_hash
      character(:), allocatable:: inputfile_text
      character(len=8):: date, localtime
      character(len=5):: timezone
      character(len=23):: datetime

      integer:: inunit = 11, filesize, infiletext_id, infiletext_dim, hash_id, hash_dim, ioerr

      ! Define output file name from runname
      outfile = trim(runname)//".nc"

      git_hash = GIT_HASH

      call date_and_time(date,localtime,timezone)
      datetime = date // "_"//localtime(1:4)//"_"//timezone

      ! Create output file
      call check( nf90_create(outfile, NF90_CLOBBER, ncid) )

      ! Read the input file as a string
      infile = trim(runname)//".in"
      open(unit=inunit, file=infile, action="read", status="old", form="unformatted",access="stream",iostat=ioerr)
      inquire(unit=inunit,size=filesize)
      allocate(character(filesize)::inputfile_text)
      read(unit=inunit)  inputfile_text
      close(inunit)

      ! Define dimensions
!      call check( nf90_def_dim(ncid,"filestring_len",len(trim(inputfile_text)),infiletext_dim))
!      call check( nf90_def_dim(ncid,"hash_len",len(git_hash),hash_dim))

      ! Define variables
!      call check( nf90_def_var(ncid,"inputfile_text",NF90_CHAR,infiletext_dim,infiletext_id) )
!      call check( nf90_def_var(ncid,"git_version",NF90_CHAR,hash_dim,hash_id) )

      ! Write some global attributes
      call check( nf90_put_att(ncid,NF90_GLOBAL,"inputfile_text",inputfile_text))
      call check( nf90_put_att(ncid,NF90_GLOBAL,"git_version",git_hash))
      call check( nf90_put_att(ncid,NF90_GLOBAL,"date_created",datetime))

      ! End define mode
      call check( nf90_enddef(ncid) )

!      call check( nf90_put_var(ncid,infiletext_id,inputfile_text))
!      call check( nf90_put_var(ncid,hash_id,git_hash))
   
   end subroutine init_output

   subroutine finish_io()
         implicit none
         call check(nf90_close(ncid))
   end subroutine finish_io

   subroutine check(flag)
      integer, intent(in):: flag
      if (flag /= nf90_noerr) then
              print*, trim(nf90_strerror(flag))
              stop 
      end if
   end subroutine check

end module io
