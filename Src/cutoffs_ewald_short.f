c
c
      subroutine cutoffs_ewald_short
      use sizes
      use atoms
      use bound
      use hescut
      use keys
      use neigh
      use limits
      use limits2
      use neigh2
      use polpot
      use tarray
      use sizes2
      implicit none
      integer i,next
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set defaults for spherical energy cutoff distances
c
      ewaldcutshort = 3.5d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffers
c
      lbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      use_mlist = .false.
      domlst2 = .true.
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:6) .eq. 'EWALD ') then
            use_ewald = .true.
         else if (keyword(1:18) .eq. 'EWALD-CUTOFFSHORT ') then
            read (string,*,err=10,end=10)  ewaldcutshort
c
c
c
c     get values for the tapering style and neighbor method
c
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.
         end if
   10    continue
      end do
c
c     set buffer region limits for pairwise neighbor lists
c
      lbuf2 = (0.5d0*lbuffer)**2
      mbuf2_short = (ewaldcutshort+lbuffer)**2
      mbufx_short = (ewaldcutshort+2.0d0*lbuffer)**2
c
c     perform dynamic allocation of some global arrays
c
      if (use_clist .or. use_mlist) then
         if (.not.allocated(nelst2))  allocate (nelst2(n))
         if (.not.allocated(elst2))  allocate (elst2(maxelst2,n))
      end if
      if (use_mlist) then
         if (.not.allocated(xmold2))  allocate (xmold2(n))
         if (.not.allocated(ymold2))  allocate (ymold2(n))
         if (.not.allocated(zmold2))  allocate (zmold2(n))
      end if
      return
      end

