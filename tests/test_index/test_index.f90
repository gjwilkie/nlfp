program test_index
   use grids
   use input
   implicit none
   integer:: ncase = 12, i, ir, ip, ix, idx_return, idx
   integer,dimension(:),allocatable:: idx_case
   logical:: passed

   allocate(idx_case(ncase))

   idx_case = (/5, 7, 1, 420, 60, 61, 41, 42, 43, 138, 270, 388 /)

   call set_resolutions(6,10,7)
   call set_layout("xpr")
   call set_gridopts_uniform()
   call init_grids()

   passed = .true.

   ! test_index: Test return of global index
   idx = get_idx(4,4,6)
   if (idx /= 322) then
      print*, "test_index: Test failed. idx=",idx,", but should =388"
      passed = .false.
   end if

   ! test_index: Test return of coord indices
   ir = get_idx_r(idx_case(12))
   ip = get_idx_p(idx_case(12))
   ix = get_idx_x(idx_case(12))
   if (ix /= 7) then
      print*, "test_index: Test failed. ix=",ix,", but should =7"
      passed = .false.
   end if
   if (ip /= 5) then
      print*, "test_index: Test failed. ip=",ip,", but should =5"
      passed = .false.
   end if
   if (ir /= 4) then
      print*, "test_index: Test failed. ir=",ir,", but should =4"
      passed = .false.
   end if

   do i = 1,ncase
      idx = idx_case(i)
      ir = get_idx_r(idx)
      ip = get_idx_p(idx)
      ix = get_idx_x(idx)
      idx_return = get_idx(ir,ip,ix)
      if (idx /= idx_return) then
         print*, "test_index: Test failed. idx=",idx," idx_return=",idx_return,". ir = ",ir,", ip = ",ip,", ix = ",ix
         passed = .false.
      end if
   end do

   if (passed) then
      print*, "test_index: Tests passed."
   end if

end program test_index
