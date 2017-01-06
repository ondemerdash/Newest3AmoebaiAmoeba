c
      subroutine emrecip1_3b_Perm2_noqpole
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use virial
      implicit none
      integer i,j,k,ii
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 cphim(4),cphid(4)
      real*8 cphip(4)
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
c
c     derivative indices into the fphi and fphidp arrays
c
      print*,"NO Quadrupole Version"
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,n))
      allocate (trq(3,npole))
c      allocate (cmp(10,npole))
c      allocate (fmp(10,npole))
c      allocate (fphi(20,npole))
c      allocate (cphi(10,npole))

      allocate (cmp(4,npole))
      allocate (fmp(4,npole))
      allocate (fphi(10,npole))
      allocate (cphi(4,npole))

c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
c         cmp(5,i) = rpole(5,i)
c         cmp(6,i) = rpole(9,i)
c         cmp(7,i) = rpole(13,i)
c         cmp(8,i) = 2.0d0 * rpole(6,i)
c         cmp(9,i) = 2.0d0 * rpole(7,i)
c         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     compute the arrays of B-spline coefficients
c
c      if (.not. use_polar) then
         call bspline_fill
         call table_fill
c      end if
c
c     perform dynamic allocation of some local arrays
c
c      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c

c      if (use_polar) then
c         end do
c      else
         call cmp_to_fmp_noqpole (cmp,fmp)
         call grid_mpole_noqpole (fmp)
         call fftfront
      ! qgrip does not appear to be necessary for just multipoles  
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
c                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
c               end do
c            end do
c         end do
c      end if

c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)*qgrid(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrid(2,k1,k2,k3)
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vyx = vyx + h2*h1*vterm
            vzx = vzx + h3*h1*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vzy = vzy + h3*h2*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
       ! Condense step below w/ current loop
       !  qfac(k1,k2,k3) = expterm
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
      end do
c
c     assign just the induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c

c      if (use_polar .and. poltyp.eq.'DIRECT') then
c      end if
c
c     perform deallocation of some local arrays
c
c      deallocate (qgrip)
c
c     transform permanent multipoles without induced dipoles
c

c      if (use_polar) then
c      end if

c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = 0.5d0 * expterm * struc2
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c

c      do k = 1, nfft3
c         do j = 1, nfft2
c            do i = 1, nfft1
c               term = qfac(i,j,k)
c               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
c               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
c            end do
c         end do
c      end do


c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole_noqpole (fphi)
      do i = 1, npole
c         do j = 1, 20
         do j = 1, 10
            fphi(j,i) = electric * fphi(j,i)
         end do
      end do
      call fphi_to_cphi_noqpole (fphi,cphi)
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
c         do k = 1, 10
         do k = 1, 4
            e = e + fmp(k,i)*fphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
        ! Condensed from below
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + frc(1,i)
         dem(2,ii) = dem(2,ii) + frc(2,i)
         dem(3,ii) = dem(3,ii) + frc(3,i)

c         trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
c     &                 + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
c     &                 + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
c     &                 - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
c         trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
c     &                 + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
c     &                 + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
c     &                 - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
c         trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
c     &                 + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
c     &                 + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
c     &                 - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)

         trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
         trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
         trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)

         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0

c         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
c     &             - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
c         vyx = vyx - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
c     &             - (cmp(5,i)+cmp(6,i))*cphi(8,i)
c     &             - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
c     &             - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
c         vzx = vzx - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
c     &             - (cmp(5,i)+cmp(7,i))*cphi(9,i)
c     &             - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
c     &             - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
c         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
c     &             - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
c         vzy = vzy - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
c     &             - (cmp(6,i)+cmp(7,i))*cphi(10,i)
c     &             - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
c     &             - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
c         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
c     &             - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)


         vxx = vxx - cmp(2,i)*cphi(2,i) 
         vyx = vyx - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
         vzx = vzx - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
         vyy = vyy - cmp(3,i)*cphi(3,i) 
         vzy = vzy - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
         vzz = vzz - cmp(4,i)*cphi(4,i) 

      end do
      e = 0.5d0 * e
      em = em + e

      !do i = 1, npole
      !   ii = ipole(i)
      !   dem(1,ii) = dem(1,ii) + frc(1,i)
      !   dem(2,ii) = dem(2,ii) + frc(2,i)
      !   dem(3,ii) = dem(3,ii) + frc(3,i)
      !end do
c
c     distribute torques into the permanent multipole gradient
c
      !do i = 1, npole
      !   trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
c     &                 + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
c     &                 + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
c     &                 - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
c         trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
c     &                 + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
c     &                 + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
c     &                 - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
c         trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
c     &                 + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
c     &                 + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
c     &                 - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
c      end do
      !do i = 1, n
      !   frc(1,i) = 0.0d0
      !   frc(2,i) = 0.0d0
      !   frc(3,i) = 0.0d0
      !end do
      call torque2 (trq,frc)
      do i = 1, n
         dem(1,i) = dem(1,i) + frc(1,i)
         dem(2,i) = dem(2,i) + frc(2,i)
         dem(3,i) = dem(3,i) + frc(3,i)
      end do
c
c     permanent multipole contribution to the internal virial
c
c      do i = 1, npole
c         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
c     &             - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
c         vyx = vyx - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
c     &             - (cmp(5,i)+cmp(6,i))*cphi(8,i)
c     &             - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
c     &             - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
c         vzx = vzx - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
c     &             - (cmp(5,i)+cmp(7,i))*cphi(9,i)
c     &             - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
c     &             - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
c         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
c     &             - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
c         vzy = vzy - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
c     &             - (cmp(6,i)+cmp(7,i))*cphi(10,i)
c     &             - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
c     &             - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
c         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
c     &             - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)
c      end do
c
c     convert Cartesian induced dipoles to fractional coordinates
c

c      if (use_polar) then
c      end if

c
c     increment the internal virial tensor components
c
      viremrecip(1,1) = viremrecip(1,1) + vxx
      viremrecip(2,1) = viremrecip(2,1) + vyx
      viremrecip(3,1) = viremrecip(3,1) + vzx
      viremrecip(1,2) = viremrecip(1,2) + vyx
      viremrecip(2,2) = viremrecip(2,2) + vyy
      viremrecip(3,2) = viremrecip(3,2) + vzy
      viremrecip(1,3) = viremrecip(1,3) + vzx
      viremrecip(2,3) = viremrecip(2,3) + vzy
      viremrecip(3,3) = viremrecip(3,3) + vzz
c
c     perform deallocation of some local arrays
c
c      print*,"End of emrecip1_3b_Perm",em
      deallocate (frc)
      deallocate (trq)
c      deallocate (fuind)
c      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
c      deallocate (fphid)
c      deallocate (fphip)
c      deallocate (fphidp)
      deallocate (cphi)
      return
      end
