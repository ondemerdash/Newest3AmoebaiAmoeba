c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2010 by T. Darden, D. Gohara & Jay W. Ponder  ##
c     ##                      All Rights Reserved                     ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below implement various B-spline and coordinate  ##
c     ##  manipulations for particle mesh Ewald summation; spatial  ##
c     ##  grid assignment by David Gohara; modified from original   ##
c     ##  PME code by Thomas Darden, NIEHS, Research Triangle, NC   ##
c     ##                                                            ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine bspline_fill  --  get PME B-spline coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "bspline_fill" finds B-spline coefficients and derivatives
c     for PME atomic sites along the fractional coordinate axes
c
c
      subroutine bspline_fill_omp1bclust(cnt)
      use sizes
      use atoms
      use boxes
      use boxes1bclust
      use pme
      implicit none
      integer i,ifr
      real*8 xi,yi,zi
      real*8 w,fr,eps
      logical first
      integer cnt
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c

      if (first) then
         first = .false.
         if (.not. allocated(igrid))  allocate (igrid(3,n))
      end if

c
c     offset used to shift sites off exact lattice bounds
c
      eps = 1.0d-8
      !print*,"OMP version of bspline_fill"
c
c     get the B-spline coefficients for each atomic site
c
!$OMP PARALLEL default(shared) private(xi,yi,zi,fr,w,ifr,
!$OMP& i
!$OMP& )
!$OMP DO

      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip1b(1,1,cnt) + yi*recip1b(2,1,cnt) 
     &        + zi*recip1b(3,1,cnt)
         fr = dble(nfft1) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(1,i) = ifr - bsorder
         call bsplgen_omp (w,thetai1(1,1,i))
         w = xi*recip1b(1,2,cnt) + yi*recip1b(2,2,cnt) 
     &         + zi*recip1b(3,2,cnt)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(2,i) = ifr - bsorder
         call bsplgen_omp (w,thetai2(1,1,i))
         w = xi*recip1b(1,3,cnt) + yi*recip1b(2,3,cnt) 
     &        + zi*recip1b(3,3,cnt)
         fr = dble(nfft3) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(3,i) = ifr - bsorder
         call bsplgen_omp (w,thetai3(1,1,i))
      end do
!$OMP END DO
!$OMP END PARALLEL

      
c      print*,"Successful end of bspline_fill"
      return
      end
c
c
