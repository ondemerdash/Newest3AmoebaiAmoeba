      subroutine alloc3b_uind
      use sizes
      use atoms
      use deriv3b
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      if (first2) then
         first2 = .false.
         if (.not. allocated(uind3b)) allocate (uind3b(3,n))
      end if
      do i = 1, n
         do j=1,3
            uind3b(j,i) = 0.0d0
         end do
      end do

      return
      end

