module output
use netcdf
implicit none

public:: init_output, write_initial_data, finish_output

private

integer :: ncid !> NetCDF control id for output file
integer :: r_dim,p_dim,xi_dim
integer :: pgrid_id,xigrid_id,rgrid_id

contains
   subroutine init_output(runname)
      use input
      implicit none
      character(len=64), intent(in):: runname
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
      call check( nf90_def_dim(ncid,"r",Nr,r_dim))
      call check( nf90_def_dim(ncid,"p",Np,p_dim))
      call check( nf90_def_dim(ncid,"xi",Nxi,xi_dim))

      ! Define variables
      call check( nf90_def_var(ncid,"rgrid",NF90_DOUBLE,r_dim,rgrid_id) )
      call check( nf90_def_var(ncid,"pgrid",NF90_DOUBLE,p_dim,pgrid_id) )
      call check( nf90_def_var(ncid,"xigrid",NF90_DOUBLE,xi_dim,xigrid_id) )

      ! Write some global attributes
      call check( nf90_put_att(ncid,NF90_GLOBAL,"inputfile_text",inputfile_text))
      call check( nf90_put_att(ncid,NF90_GLOBAL,"git_version",git_hash))
      call check( nf90_put_att(ncid,NF90_GLOBAL,"date_created",datetime))

      ! End define mode
      call check( nf90_enddef(ncid) )

   end subroutine init_output

   subroutine write_initial_data()
      implicit none
   end subroutine

   subroutine finish_output()
      implicit none
      call check(nf90_close(ncid))
   end subroutine finish_output

   subroutine check(flag)
      integer, intent(in):: flag
      if (flag /= nf90_noerr) then
              print*, trim(nf90_strerror(flag))
              stop 
      end if
   end subroutine check


end module output
