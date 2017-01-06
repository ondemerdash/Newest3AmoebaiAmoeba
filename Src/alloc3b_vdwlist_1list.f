      subroutine alloc3b_vdwlist_1list
      use sizes
      use atoms
      use deriv3b
      use limits
      use aprx
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep3bmut=0.0d0
      ep3b=0.0d0
      ntriples=0
      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
         if (.not. allocated(dep3bmut)) allocate (dep3bmut(3,n))
      end if
        do i=1,3
          do j=1,3
            virep3b(j,i)=0.0d0
            virep3bmut(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
            dep3bmut(j,i) = 0.0d0
         end do
      end do

      if(use_vlist) call vlist
      if (approxmode.ne.'1BODYMODE') call mollist2bodyOO6_8_par_1list
      return
      end
