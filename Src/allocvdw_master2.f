      subroutine allocvdw_master2 
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use iounit
      use limits
      use potent
      use virial
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time

      save first
      data first  / .true. /


      ev = 0.0d0

      if (first) then
         first = .false.
         if (.not. allocated(dev))  allocate (dev(3,n))
      end if


      do i = 1, n
         do j=1,3
            dev(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            virev(j,i)=0.0d0
         end do
      end do

      return
      end

