c
      subroutine emrecip1_3b_Polartest(npole3b,pnum,uind,uinp,
     & eptemp,deptemp,virtemp,emtemp,demtemp)
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
      use polar, only: polarity, thole, pdamp
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
      integer npole3b,pnum(*),cnt,l1
      real*8 uind(3,*),eptemp,deptemp(3,*)
      real*8 virtemp(3,3),uinp(3,*)
      real*8 emtemp,demtemp(3,*),viremtemp(3,3)
c
c     derivative indices into the fphi and fphidp arrays
c
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
      allocate (trq(3,n))
      allocate (fuind(3,npole3b))
      allocate (fuinp(3,npole3b))
      allocate (cmp(10,npole3b))
      allocate (fmp(10,npole3b))
      allocate (fphi(20,npole3b))
      allocate (fphid(10,npole3b))
      allocate (fphip(10,npole3b))
      allocate (fphidp(20,npole3b))
      allocate (cphi(10,npole3b))
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
      !do i = 1, npole
      do l1 = 1, npole3b
         i=pnum(l1)
         cmp(1,l1) = rpole(1,i)
         cmp(2,l1) = rpole(2,i)
         cmp(3,l1) = rpole(3,i)
         cmp(4,l1) = rpole(4,i)
         cmp(5,l1) = rpole(5,i)
         cmp(6,l1) = rpole(9,i)
         cmp(7,l1) = rpole(13,i)
         cmp(8,l1) = 2.0d0 * rpole(6,i)
         cmp(9,l1) = 2.0d0 * rpole(7,i)
         cmp(10,l1) = 2.0d0 * rpole(10,i)
      end do
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     compute the arrays of B-spline coefficients
c
      if (.not. use_polar) then
         call bspline_fill
         call table_fill
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
      if (use_polar) then
c         do i = 1, npole
         do l1=1,npole3b
            i=pnum(l1)
            do j = 2, 4
               cmp(j,l1) = cmp(j,l1) + uinp(j-1,l1)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
c         do i = 1, npole
         do l1=1,npole3b
            i=pnum(l1)
            do j = 2, 4
               cmp(j,l1) = cmp(j,l1) + uind(j-1,l1) - uinp(j-1,l1)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         !do i = 1, npole
         do l1 = 1, npole3b
            do j = 2, 4
               cmp(j,l1) = cmp(j,l1) - uind(j-1,l1)
            end do
         end do
      else
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
      end if
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
            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vyx = vyx + h2*h1*vterm
            vzx = vzx + h3*h1*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vzy = vzy + h3*h2*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     assign just the induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
      if (use_polar .and. poltyp.eq.'DIRECT') then
         do i = 1, npole3b
            do j = 1, 10
               cmp(j,i) = 0.0d0
            end do
            do j = 2, 4
               cmp(j,i) = uinp(j-1,i)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
         do i = 1, npole3b
            do j = 2, 4
               cmp(j,i) = uind(j-1,i)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do l1 = 1, npole3b
            i=pnum(l1)
            cmp(1,l1) = rpole(1,i)
            cmp(2,l1) = rpole(2,i)
            cmp(3,l1) = rpole(3,i)
            cmp(4,l1) = rpole(4,i)
            cmp(5,l1) = rpole(5,i)
            cmp(6,l1) = rpole(9,i)
            cmp(7,l1) = rpole(13,i)
            cmp(8,l1) = 2.0d0 * rpole(6,i)
            cmp(9,l1) = 2.0d0 * rpole(7,i)
            cmp(10,l1) = 2.0d0 * rpole(10,i)
         end do
c
c     make the scalar summation over reciprocal lattice
c
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
               struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                     + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
               eterm = 0.5d0 * electric * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx - h1*h1*vterm + eterm
               vyx = vyx - h2*h1*vterm
               vzx = vzx - h3*h1*vterm
               vyy = vyy - h2*h2*vterm + eterm
               vzy = vzy - h3*h2*vterm
               vzz = vzz - h3*h3*vterm + eterm
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)
c
c     transform permanent multipoles without induced dipoles
c
c      if (use_polar) then
c         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
c         call grid_mpole3b_new (npole3b,pnum,fmp)
c         call fftfront
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
c      call fftback
c      call fphi_mpole3b_new (npole3b,pnum,fphi)
c      do i = 1, npole3b
c         do j = 1, 20
c            fphi(j,i) = electric * fphi(j,i)
c         end do
c      end do
c      call fphi_to_cphi3b (npole3b,pnum,fphi,cphi)
c
c     increment the permanent multipole energy and gradient
c
c      e = 0.0d0
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         f1 = 0.0d0
c         f2 = 0.0d0
c         f3 = 0.0d0
c         do k = 1, 10
c            e = e + fmp(k,l1)*fphi(k,l1)
c            f1 = f1 + fmp(k,l1)*fphi(deriv1(k),l1)
c            f2 = f2 + fmp(k,l1)*fphi(deriv2(k),l1)
c            f3 = f3 + fmp(k,l1)*fphi(deriv3(k),l1)
c         end do
c         f1 = dble(nfft1) * f1
c         f2 = dble(nfft2) * f2
c         f3 = dble(nfft3) * f3
c         frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
c         frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
c         frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
c      end do
c      e = 0.5d0 * e
c      emtemp = emtemp + e
c      do l1 = 1, npole3b
c         i=pnum(l1)
c      !   ii = ipole(i)
cc         demtemp(1,l1) = demtemp(1,l1) + frc(1,i)
cc         demtemp(2,l1) = demtemp(2,l1) + frc(2,i)
cc         demtemp(3,l1) = demtemp(3,l1) + frc(3,i)
c
c         demtemp(1,i) = demtemp(1,i) + frc(1,i)
c         demtemp(2,i) = demtemp(2,i) + frc(2,i)
c         demtemp(3,i) = demtemp(3,i) + frc(3,i)
c      end do
c
c     distribute torques into the permanent multipole gradient
c
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         trq(1,i) = cmp(4,l1)*cphi(3,l1) - cmp(3,l1)*cphi(4,l1)
c     &                 + 2.0d0*(cmp(7,l1)-cmp(6,l1))*cphi(10,l1)
c     &                 + cmp(9,l1)*cphi(8,l1) + cmp(10,l1)*cphi(6,l1)
c     &                 - cmp(8,l1)*cphi(9,l1) - cmp(10,l1)*cphi(7,l1)
c         trq(2,i) = cmp(2,l1)*cphi(4,l1) - cmp(4,l1)*cphi(2,l1)
c     &                 + 2.0d0*(cmp(5,l1)-cmp(7,l1))*cphi(9,l1)
c     &                 + cmp(8,l1)*cphi(10,l1) + cmp(9,l1)*cphi(7,l1)
c     &                 - cmp(9,l1)*cphi(5,l1) - cmp(10,l1)*cphi(8,l1)
c         trq(3,i) = cmp(3,l1)*cphi(2,l1) - cmp(2,l1)*cphi(3,l1)
c     &                 + 2.0d0*(cmp(6,l1)-cmp(5,l1))*cphi(8,l1)
c     &                 + cmp(8,l1)*cphi(5,l1) + cmp(10,l1)*cphi(9,l1)
c     &                 - cmp(8,l1)*cphi(6,l1) - cmp(9,l1)*cphi(10,l1)
c      end do
c      do i = 1, n
c         frc(1,i) = 0.0d0
c         frc(2,i) = 0.0d0
c         frc(3,i) = 0.0d0
c      end do
c      call torque2_3bnew (npole3b,pnum,trq,frc)
c      do l1 = 1, npole3b
c         i=pnum(l1)
cc         demtemp(1,l1) = demtemp(1,l1) + frc(1,i)
cc         demtemp(2,l1) = demtemp(2,l1) + frc(2,i)
cc         demtemp(3,l1) = demtemp(3,l1) + frc(3,i)
c
c         demtemp(1,i) = demtemp(1,i) + frc(1,i)
c         demtemp(2,i) = demtemp(2,i) + frc(2,i)
c         demtemp(3,i) = demtemp(3,i) + frc(3,i)
c      end do
c
c     permanent multipole contribution to the internal virial
c
c      do i = 1, npole3b
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
      if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         do i = 1, npole3b
            do j = 1, 3
               fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
     &                          + a(j,3)*uind(3,i)
               fuinp(j,i) = a(j,1)*uinp(1,i) + a(j,2)*uinp(2,i)
     &                          + a(j,3)*uinp(3,i)
            end do
         end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
         call grid_uind3b_new (npole3b,pnum,fuind,fuinp)
         call fftfront
c
c     account for the zeroth grid point for a finite system
c
         if (.not. use_bounds) then
            expterm = 0.5d0 * pi / xbox
            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
            e = 0.5d0 * expterm * struc2
         end if
c
c     complete the transformation of the PME grid
c
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  term = qfac(i,j,k)
                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
               end do
            end do
         end do
c
c     perform 3-D FFT backward transform and get potential
c
         call fftback
         call fphi_uind3b_new (npole3b,pnum,fphid,fphip,fphidp)
         do i = 1, npole3b
            do j = 1, 10
               fphid(j,i) = electric * fphid(j,i)
               fphip(j,i) = electric * fphip(j,i)
            end do
            do j = 1, 20
               fphidp(j,i) = electric * fphidp(j,i)
            end do
         end do
c
c     increment the induced dipole energy and gradient
c
         e = 0.0d0
         do l1 = 1, npole3b
            i=pnum(l1)
            f1 = 0.0d0
            f2 = 0.0d0
            f3 = 0.0d0
            do k = 1, 3
               j1 = deriv1(k+1)
               j2 = deriv2(k+1)
               j3 = deriv3(k+1)
               e = e + fuind(k,l1)*fphi(k+1,l1)
               f1 = f1 + (fuind(k,l1)+fuinp(k,l1))*fphi(j1,l1)
     &                 + fuind(k,l1)*fphip(j1,l1)
     &                 + fuinp(k,l1)*fphid(j1,l1)
               f2 = f2 + (fuind(k,l1)+fuinp(k,l1))*fphi(j2,l1)
     &                 + fuind(k,l1)*fphip(j2,l1)
     &                 + fuinp(k,l1)*fphid(j2,l1)
               f3 = f3 + (fuind(k,l1)+fuinp(k,l1))*fphi(j3,l1)
     &                 + fuind(k,l1)*fphip(j3,l1)
     &                 + fuinp(k,l1)*fphid(j3,l1)
               if (poltyp .eq. 'DIRECT') then
                  f1 = f1 - fuind(k,l1)*fphip(j1,l1)
     &                    - fuinp(k,l1)*fphid(j1,l1)
                  f2 = f2 - fuind(k,l1)*fphip(j2,l1)
     &                    - fuinp(k,l1)*fphid(j2,l1)
                  f3 = f3 - fuind(k,l1)*fphip(j3,l1)
     &                    - fuinp(k,l1)*fphid(j3,l1)
               end if
            end do
            do k = 1, 10
               f1 = f1 + fmp(k,l1)*fphidp(deriv1(k),l1)
               f2 = f2 + fmp(k,l1)*fphidp(deriv2(k),l1)
               f3 = f3 + fmp(k,l1)*fphidp(deriv3(k),l1)
            end do
            f1 = 0.5d0 * dble(nfft1) * f1
            f2 = 0.5d0 * dble(nfft2) * f2
            f3 = 0.5d0 * dble(nfft3) * f3
            frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
            frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
            frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         end do
         e = 0.5d0 * e
         eptemp = eptemp + e
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
c            deptemp(1,l1) = deptemp(1,l1) + frc(1,i)
c            deptemp(2,l1) = deptemp(2,l1) + frc(2,i)
c            deptemp(3,l1) = deptemp(3,l1) + frc(3,i)

            deptemp(1,i) = deptemp(1,i) + frc(1,i)
            deptemp(2,i) = deptemp(2,i) + frc(2,i)
            deptemp(3,i) = deptemp(3,i) + frc(3,i)
         end do
c
c     set the potential to be the induced dipole average
c
         do i = 1, npole3b
            do k = 1, 10
               fphidp(k,i) = 0.5d0 * fphidp(k,i)
            end do
         end do
         call fphi_to_cphi3b (npole3b,pnum,fphidp,cphi)
c
c     distribute torques into the induced dipole gradient
c
         do l1 = 1, npole3b
            i=pnum(l1)
            trq(1,i) = cmp(4,l1)*cphi(3,l1) - cmp(3,l1)*cphi(4,l1)
     &                    + 2.0d0*(cmp(7,l1)-cmp(6,l1))*cphi(10,l1)
     &                    + cmp(9,l1)*cphi(8,l1) + cmp(10,l1)*cphi(6,l1)
     &                    - cmp(8,l1)*cphi(9,l1) - cmp(10,l1)*cphi(7,l1)
            trq(2,i) = cmp(2,l1)*cphi(4,l1) - cmp(4,l1)*cphi(2,l1)
     &                    + 2.0d0*(cmp(5,l1)-cmp(7,l1))*cphi(9,l1)
     &                    + cmp(8,l1)*cphi(10,l1) + cmp(9,l1)*cphi(7,l1)
     &                    - cmp(9,l1)*cphi(5,l1) - cmp(10,l1)*cphi(8,l1)
            trq(3,i) = cmp(3,l1)*cphi(2,l1) - cmp(2,l1)*cphi(3,l1)
     &                    + 2.0d0*(cmp(6,l1)-cmp(5,l1))*cphi(8,l1)
     &                    + cmp(8,l1)*cphi(5,l1) + cmp(10,l1)*cphi(9,l1)
     &                    - cmp(8,l1)*cphi(6,l1) - cmp(9,l1)*cphi(10,l1)
         end do
         do i = 1, n
            frc(1,i) = 0.0d0
            frc(2,i) = 0.0d0
            frc(3,i) = 0.0d0
         end do
         call torque2_3bnew (npole3b,pnum,trq,frc)
         do l1 = 1, npole3b
            i=pnum(l1)
c            deptemp(1,l1) = deptemp(1,l1) + frc(1,i)
c            deptemp(2,l1) = deptemp(2,l1) + frc(2,i)
c            deptemp(3,l1) = deptemp(3,l1) + frc(3,i)

            deptemp(1,i) = deptemp(1,i) + frc(1,i)
            deptemp(2,i) = deptemp(2,i) + frc(2,i)
            deptemp(3,i) = deptemp(3,i) + frc(3,i)
         end do
c
c     induced dipole contribution to the internal virial
c
         do i = 1, npole3b
            l1=i
            do j = 2, 4
               cphim(j) = 0.0d0
               cphid(j) = 0.0d0
               cphip(j) = 0.0d0
               do k = 2, 4
                  cphim(j) = cphim(j) + ftc(j,k)*fphi(k,i)
                  cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
                  cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
               end do
            end do
            vxx = vxx - cphi(2,i)*cmp(2,i)
     &                - 0.5d0*(cphim(2)*(uind(1,l1)+uinp(1,l1))
     &                        +cphid(2)*uinp(1,l1)+cphip(2)*uind(1,l1))
            vyx = vyx - 0.5d0*(cphi(2,i)*cmp(3,i)+cphi(3,i)*cmp(2,i))
     &                - 0.25d0*(cphim(2)*(uind(2,l1)+uinp(2,l1))
     &                         +cphim(3)*(uind(1,l1)+uinp(1,l1))
     &                         +cphid(2)*uinp(2,l1)+cphip(2)*uind(2,l1)
     &                         +cphid(3)*uinp(1,l1)+cphip(3)*uind(1,l1))
            vzx = vzx - 0.5d0*(cphi(2,i)*cmp(4,i)+cphi(4,i)*cmp(2,i))
     &                - 0.25d0*(cphim(2)*(uind(3,l1)+uinp(3,l1))
     &                         +cphim(4)*(uind(1,l1)+uinp(1,l1))
     &                         +cphid(2)*uinp(3,l1)+cphip(2)*uind(3,l1)
     &                         +cphid(4)*uinp(1,l1)+cphip(4)*uind(1,l1))
            vyy = vyy - cphi(3,i)*cmp(3,i)
     &                - 0.5d0*(cphim(3)*(uind(2,l1)+uinp(2,l1))
     &                        +cphid(3)*uinp(2,l1)+cphip(3)*uind(2,l1))
            vzy = vzy - 0.5d0*(cphi(3,i)*cmp(4,i)+cphi(4,i)*cmp(3,i))
     &                - 0.25d0*(cphim(3)*(uind(3,l1)+uinp(3,l1))
     &                         +cphim(4)*(uind(2,l1)+uinp(2,l1))
     &                         +cphid(3)*uinp(3,l1)+cphip(3)*uind(3,l1)
     &                         +cphid(4)*uinp(2,l1)+cphip(4)*uind(2,l1))
            vzz = vzz - cphi(4,i)*cmp(4,i)
     &                - 0.5d0*(cphim(4)*(uind(3,l1)+uinp(3,l1))
     &                        +cphid(4)*uinp(3,l1)+cphip(4)*uind(3,l1))
            vxx = vxx - 2.0d0*cmp(5,i)*cphi(5,i) - cmp(8,i)*cphi(8,i)
     &                - cmp(9,i)*cphi(9,i)
            vyx = vyx - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &                - 0.5d0*(cmp(8,i)*(cphi(6,i)+cphi(5,i))
     &                     +cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
            vzx = vzx - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &                - 0.5d0*(cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &                     +cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
            vyy = vyy - 2.0d0*cmp(6,i)*cphi(6,i) - cmp(8,i)*cphi(8,i)
     &                - cmp(10,i)*cphi(10,i)
            vzy = vzy - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &                - 0.5d0*(cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &                     +cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
            vzz = vzz - 2.0d0*cmp(7,i)*cphi(7,i) - cmp(9,i)*cphi(9,i)
     &                - cmp(10,i)*cphi(10,i)
            if (poltyp .eq. 'DIRECT') then
             vxx = vxx + 0.5d0*(cphid(2)*uinp(1,l1)+cphip(2)*uind(1,l1))
             vyx = vyx + 0.25d0*(cphid(2)*uinp(2,l1)+cphip(2)*uind(2,l1)
     &                         +cphid(3)*uinp(1,l1)+cphip(3)*uind(1,l1))
             vzx = vzx + 0.25d0*(cphid(2)*uinp(3,l1)+cphip(2)*uind(3,l1)
     &                         +cphid(4)*uinp(1,l1)+cphip(4)*uind(1,l1))
             vyy = vyy + 0.5d0*(cphid(3)*uinp(2,l1)+cphip(3)*uind(2,l1))
             vzy = vzy + 0.25d0*(cphid(3)*uinp(3,l1)+cphip(3)*uind(3,l1)
     &                         +cphid(4)*uinp(2,l1)+cphip(4)*uind(2,l1))
             vzz = vzz + 0.5d0*(cphid(4)*uinp(3,l1)+cphip(4)*uind(3,l1))
            end if
         end do
      end if
c
c     increment the internal virial tensor components
c
      virtemp(1,1) = virtemp(1,1) + vxx
      virtemp(2,1) = virtemp(2,1) + vyx
      virtemp(3,1) = virtemp(3,1) + vzx
      virtemp(1,2) = virtemp(1,2) + vyx
      virtemp(2,2) = virtemp(2,2) + vyy
      virtemp(3,2) = virtemp(3,2) + vzy
      virtemp(1,3) = virtemp(1,3) + vzx
      virtemp(2,3) = virtemp(2,3) + vzy
      virtemp(3,3) = virtemp(3,3) + vzz
c
c     perform deallocation of some local arrays
c
      deallocate (frc)
      deallocate (trq)
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
      deallocate (cphi)
      return
      end
