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
      subroutine bspline_fill_omp
      use sizes
      use atoms
      use boxes
      use pme
      implicit none
      integer i,ifr
      real*8 xi,yi,zi
      real*8 w,fr,eps
      logical first
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
         w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr = dble(nfft1) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(1,i) = ifr - bsorder
         call bsplgen_omp (w,thetai1(1,1,i))
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(2,i) = ifr - bsorder
         call bsplgen_omp (w,thetai2(1,1,i))
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bsplgen_omp  --  B-spline coefficients for an atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bsplgen_omp" gets B-spline coefficients and derivatives for
c     a single PME atomic site along a particular direction
c
c
      subroutine bsplgen_omp (w,thetai)
      use sizes
      use pme
      use potent
      implicit none
      integer i,j,k
      integer level
      real*8 w,denom
      real*8 thetai(4,*)
      real*8 temp(maxorder,maxorder)
c
c
c     set B-spline depth for partial charges or multipoles
c
      level = 2
      if (use_mpole .or. use_polar)  level = 4
c
c     initialization to get to 2nd order recursion
c
      temp(2,2) = w
      temp(2,1) = 1.0d0 - w
c
c     perform one pass to get to 3rd order recursion
c
      temp(3,3) = 0.5d0 * w * temp(2,2)
      temp(3,2) = 0.5d0 * ((1.0d0+w)*temp(2,1)+(2.0d0-w)*temp(2,2))
      temp(3,1) = 0.5d0 * (1.0d0-w) * temp(2,1)
c
c     compute standard B-spline recursion to desired order
c
      do i = 4, bsorder
         k = i - 1
         denom = 1.0d0 / dble(k)
         temp(i,i) = denom * w * temp(k,k)
         do j = 1, i-2
            temp(i,i-j) = denom * ((w+dble(j))*temp(k,i-j-1)
     &                             +(dble(i-j)-w)*temp(k,i-j))
         end do
         temp(i,1) = denom * (1.0d0-w) * temp(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      temp(k,bsorder) = temp(k,bsorder-1)
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
c
c     get coefficients for the B-spline second derivative
c
      if (level .eq. 4) then
         k = bsorder - 2
         temp(k,bsorder-1) = temp(k,bsorder-2)
         do i = bsorder-2, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
         temp(k,bsorder) = temp(k,bsorder-1)
         do i = bsorder-1, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
c
c     get coefficients for the B-spline third derivative
c
         k = bsorder - 3
         temp(k,bsorder-2) = temp(k,bsorder-3)
         do i = bsorder-3, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
         temp(k,bsorder-1) = temp(k,bsorder-2)
         do i = bsorder-2, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
         temp(k,bsorder) = temp(k,bsorder-1)
         do i = bsorder-1, 2, -1
            temp(k,i) = temp(k,i-1) - temp(k,i)
         end do
         temp(k,1) = -temp(k,1)
      end if
c
c     copy coefficients from temporary to permanent storage
c
      do i = 1, bsorder
         do j = 1, level
            thetai(j,i) = temp(bsorder-j+1,i)
         end do
      end do
      return
      end
c
c
