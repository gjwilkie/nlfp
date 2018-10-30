module inputoutput
use mp
use netcdf
implicit none

public :: Nr, Nxi, Np, ncid, ncerr



!> Number of spatial grid points
integer :: Nr

!> Number of momentum magnitude grid points
integer :: Np

!> Number of pitch angle grid points
integer :: Nxi

!> NetCDF control ids
integer :: ncid, ncerr

!> The user-specified string which refers to the input and output files
character(len=100) runname

contains

   !> Sets the runname 
   subroutine set_runname(arg)
      implicit none
      character(len=100), intent(in):: arg

      runname = arg
     
   end subroutine

   !> Reads the namelist input file runname.in and sets the options accordingly.
   !! Also creates the output file and immediately writes all input parameters.
   subroutine read_input(runname)
      implicit none
      character(len=100), intent(in):: runname
      character(len=100) :: infile
      namelist /grid_knobs/ Nxi, Np, Nr
      integer:: inunit = 11, ioerr

      call set_runname(runname)

      infile = trim(runname)//".in"
      
      open(unit=inunit, file=infile, action="read", status="old", iostat=ioerr)
      if (ioerr /= 0) then
         write(*,*) "ERROR: could not open file ",infile
         stop
      end if

      read(inunit,nml=grid_knobs)

      close(inunit)

   end subroutine

   subroutine init_output()
      implicit none
      character(len=100) :: outfile, infile
      character(:), allocatable:: inputfile_text
      integer:: ncerr, ncid, inunit = 11, filesize, ioerr, infile_id, version_id

      ! Define output file name from runname
      outfile = trim(runname)//".nc"

      ! Read the input file as a string
      infile = trim(runname)//".in"
      open(unit=inunit, file=infile, action="read", status="old", form="unformatted",access="stream",iostat=ioerr)
      inquire(unit=inunit,size=filesize)
      allocate(character(filesize)::inputfile_text)
      read(unit=inunit)  inputfile_text
      close(inunit)


      ! Create output file
      ncerr = nf90_create(outfile, NF90_CLOBBER, ncid)
      if (ncerr /= 0) then
         write(*,*) "ERROR: nf90_create returned error code ",ncerr, " when opening file ",outfile
         stop
      end if

      ! Write the input file used as a global attribute
      ncerr = nf90_put_att(ncid,NF90_GLOBAL,"inputfile_text",inputfile_text)
      if (ncerr /= 0) then
         print*, "ERROR: NetCDF could not write attribute inputfile_txt"
         print*, ncerr
      end if
   
   end subroutine

   subroutine close_output()
         implicit none
         nf90_close(ncid)
   end subroutine

end module inputoutput
