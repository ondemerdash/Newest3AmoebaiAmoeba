c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emrecip1  --  mpole Ewald recip energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to atomic multipole
c     interactions and dipole polarizability
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip1_3b_Polar_usemlist_direct 
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
c      use pme, only: nfft1,nfft2,nfft3,maxorder,bsorder,igrid,thetai1,
c     & thetai2,thetai3,bsmod1,bsmod2,bsmod3
      use pme
c      use polar, only: polarity, thole, pdamp
      use polar
      use polpot
      use potent
      use virial
      use chunks
      use fft
      use totfield
      use deriv3b
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
c      integer npole3b,pnum(*),l1,l3
c      real*8 eptemp,dep3b_recip(3,npole),virep3b_recip(3,3)
c      real*8 uind(3,npole3b),uinp(3,npole3b)
c      real*8 qgrid(2,nfft1,nfft2,nfft3)
c      real*8 qfac3b(nfft1,nfft2,nfft3)
c      real*8 thetai1_3b(4,bsorder,n)
c      real*8 thetai2_3b(4,bsorder,n)
c      real*8 thetai3_3b(4,bsorder,n)
c      integer pmetable3b(n,nchunk)
c      integer igrid3b(3,n)
c      integer*8 planf3b
c      integer*8 planb3b
c      integer iprime3b(maxprime,3)
c      real*8 ffttable3b(maxtable,3)
c      real*8, allocatable :: cphi3b(:,:)
c      logical doi
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

      allocate (frc(3,npole))
      allocate (trq(3,npole))
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (fphi(20,npole))
      allocate (fphid(10,npole))
      allocate (fphip(10,npole))
      allocate (fphidp(20,npole))

      allocate (cphi(10,npole))
      if(.not.allocated(dep3b_recip)) allocate (dep3b_recip(3,n))

c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0

       ep3b_recip=0.0d0

       do i = 1, npole
          dep3b_recip(1,i) = 0.0d0 
          dep3b_recip(2,i) = 0.0d0 
          dep3b_recip(3,i) = 0.0d0 
       end do

       do i=1,3
         do j=1,3
          virep3b_recip(j,i) = 0.0d0
         end do
       end do
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     compute the arrays of B-spline coefficients
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
c
c     make the scalar summation over reciprocal lattice
c
c
C    VIRIAL CONTRIBUTION FROM PERMANENT MULTIPOLES SAVED FROM
C    BEFORE.
          call cmp_to_fmp(cmp,fmp)

        vxx = vxx - vxx_perm
        vyx = vyx - vyx_perm
        vzx = vzx - vzx_perm
        vyy = vyy - vyy_perm
        vzy = vzy - vzy_perm
        vzz = vzz - vzz_perm
c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)

      !if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         do i=1,npole
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
         call grid_uind(fuind,fuinp)
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
                 ! term = qfac3b(i,j,k)
                  term = qfac(i,j,k)
                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
               end do
            end do
         end do
c
c     perform 3-D FFT backward transform and get potential
c

c         call fftback3b(planb3b,qgrid,iprime3b,ffttable3b)
c        call fphi_uind3b_totfield_513 (npole3b,pnum,fphid,fphip,fphidp,
c     & qgrid)
        call fftback
        call fphi_uind (fphid,fphip,fphidp)

c         print*,"After call to fphi_uind"
         do i = 1, npole
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
         do i = 1, npole
            f1 = 0.0d0
            f2 = 0.0d0
            f3 = 0.0d0
            do k = 1, 3
               j1 = deriv1(k+1)
               j2 = deriv2(k+1)
               j3 = deriv3(k+1)
c               e = e + fuind(k,i)*fphi(k+1,i)
c               f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi(j1,i)
c     &                 + fuind(k,i)*fphip(j1,i)
c     &                 + fuinp(k,i)*fphid(j1,i)
c               f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi(j2,i)
c     &                 + fuind(k,i)*fphip(j2,i)
c     &                 + fuinp(k,i)*fphid(j2,i)
c               f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi(j3,i)
c     &                 + fuind(k,i)*fphip(j3,i)
c     &                 + fuinp(k,i)*fphid(j3,i)

               !e = e + fuind(k,i)*fphi(k+1,l1)
               !f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi(j1,l1)
     &          !       + fuind(k,i)*fphip(j1,l1)
     &          !       + fuinp(k,i)*fphid(j1,l1)
               !f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi(j2,l1)
     &         !        + fuind(k,i)*fphip(j2,l1)
     &         !        + fuinp(k,i)*fphid(j2,l1)
               !f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi(j3,l1)
     &          !       + fuind(k,i)*fphip(j3,l1)
     &          !       + fuinp(k,i)*fphid(j3,l1)

c               if(doi) then
               e = e + fuind(k,i)*fphi_totfield(k+1,i)
               f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi_totfield(j1,i)
c     &                 + fuind(k,i)*fphip(j1,i)
c     &                 + fuinp(k,i)*fphid(j1,i)
               f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi_totfield(j2,i)
c     &                 + fuind(k,i)*fphip(j2,i)
c     &                 + fuinp(k,i)*fphid(j2,i)
               f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi_totfield(j3,i)
c     &                 + fuind(k,i)*fphip(j3,i)
c     &                 + fuinp(k,i)*fphid(j3,i)
c               end if

c               if ((poltyp .eq. 'DIRECT').and.(doi)) then
c                  f1 = f1 - fuind(k,i)*fphip(j1,i)
c     &                    - fuinp(k,i)*fphid(j1,i)
c                  f2 = f2 - fuind(k,i)*fphip(j2,i)
c     &                    - fuinp(k,i)*fphid(j2,i)
c                  f3 = f3 - fuind(k,i)*fphip(j3,i)
c     &                    - fuinp(k,i)*fphid(j3,i)
c               end if
            end do
            do k = 1, 10
               f1 = f1 + fmp(k,i)*fphidp(deriv1(k),i)
               f2 = f2 + fmp(k,i)*fphidp(deriv2(k),i)
               f3 = f3 + fmp(k,i)*fphidp(deriv3(k),i)
             !  f1 = f1 + fmp(k,l1)*fphidp(deriv1(k),l1)
             !  f2 = f2 + fmp(k,l1)*fphidp(deriv2(k),l1)
             !  f3 = f3 + fmp(k,l1)*fphidp(deriv3(k),l1)

             !  f1 = f1 + fmp(k,i)*fphidp(deriv1(k),l1)
             !  f2 = f2 + fmp(k,i)*fphidp(deriv2(k),l1)
             !  f3 = f3 + fmp(k,i)*fphidp(deriv3(k),l1)
            end do
            f1 = 0.5d0 * dble(nfft1) * f1
            f2 = 0.5d0 * dble(nfft2) * f2
            f3 = 0.5d0 * dble(nfft3) * f3
            frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
            frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
            frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         end do
         e = 0.5d0 * e
         ep3b_recip=ep3b_recip+e

         do i = 1, npole
            ii = ipole(i)
            dep3b_recip(1,ii) = dep3b_recip(1,ii) + frc(1,i)
            dep3b_recip(2,ii) = dep3b_recip(2,ii) + frc(2,i)
            dep3b_recip(3,ii) = dep3b_recip(3,ii) + frc(3,i)
         end do
c
c     set the potential to be the induced dipole average
c
         do i = 1, npole
            do k = 1, 10
               fphidp(k,i) = 0.5d0 * fphidp(k,i)
            end do
         end do
          call fphi_to_cphi (fphidp,cphi)
c
c     distribute torques into the induced dipole gradient
c
         do i = 1, npole
            trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &                    + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &                    + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &                    - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
            trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &                    + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &                    + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &                    - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
            trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &                    + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &                    + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &                    - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
         end do

         do i = 1, n
            frc(1,i) = 0.0d0
            frc(2,i) = 0.0d0
            frc(3,i) = 0.0d0
         end do
         call torque2 (trq,frc)
         do i = 1, n
            dep3b_recip(1,i) = dep3b_recip(1,i) + frc(1,i)
            dep3b_recip(2,i) = dep3b_recip(2,i) + frc(2,i)
            dep3b_recip(3,i) = dep3b_recip(3,i) + frc(3,i)
         end do
c
c     induced dipole contribution to the internal virial
c
         do i = 1, npole
            do j = 2, 4
               cphim(j) = 0.0d0
               cphid(j) = 0.0d0
               cphip(j) = 0.0d0
               do k = 2, 4
                  cphim(j) = cphim(j) + ftc(j,k)*fphi_totfield(k,i)
                  cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
                  cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
               end do
            end do

            vxx = vxx - cphi(2,i)*cmp(2,i)
     &                - 0.5d0*(cphim(2)*(uind(1,i)+uinp(1,i)))
     &                      !  +cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
           vyx = vyx-0.5d0*(cphi(2,i)*cmp(3,i)+cphi(3,i)*cmp(2,i))
     &                - 0.25d0*(cphim(2)*(uind(2,i)+uinp(2,i))
     &                         +cphim(3)*(uind(1,i)+uinp(1,i)))
     &                      !   +cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
     &                      !   +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
          vzx = vzx -0.5d0*(cphi(2,i)*cmp(4,i)+cphi(4,i)*cmp(2,i))
     &                - 0.25d0*(cphim(2)*(uind(3,i)+uinp(3,i))
     &                         +cphim(4)*(uind(1,i)+uinp(1,i)))
     &                      !   +cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
     &                      !   +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
            vyy = vyy - cphi(3,i)*cmp(3,i)
     &                - 0.5d0*(cphim(3)*(uind(2,i)+uinp(2,i)))
     &                     !   +cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
          vzy = vzy -0.5d0*(cphi(3,i)*cmp(4,i)+cphi(4,i)*cmp(3,i))
     &                - 0.25d0*(cphim(3)*(uind(3,i)+uinp(3,i))
     &                         +cphim(4)*(uind(2,i)+uinp(2,i)))
     &                     !    +cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
     &                     !    +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
            vzz = vzz - cphi(4,i)*cmp(4,i)
     &                - 0.5d0*(cphim(4)*(uind(3,i)+uinp(3,i)))
     &                      !  +cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))

          vxx = vxx - 2.0d0*cmp(5,i)*cphi(5,i) -cmp(8,i)*cphi(8,i)
     &                - cmp(9,i)*cphi(9,i)
            vyx = vyx - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &                - 0.5d0*(cmp(8,i)*(cphi(6,i)+cphi(5,i))
     &                   +cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
            vzx = vzx - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &                - 0.5d0*(cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &                   +cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
           vyy = vyy - 2.0d0*cmp(6,i)*cphi(6,i)-cmp(8,i)*cphi(8,i)
     &                - cmp(10,i)*cphi(10,i)
            vzy = vzy - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &                - 0.5d0*(cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &                     +cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
          vzz = vzz - 2.0d0*cmp(7,i)*cphi(7,i) -cmp(9,i)*cphi(9,i)
     &                - cmp(10,i)*cphi(10,i)

c            if ((poltyp .eq. 'DIRECT').and.doi) then
c             vxx = vxx + 0.5d0*(cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
c             vyx = vyx + 0.25d0*(cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
c     &                         +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
c             vzx = vzx + 0.25d0*(cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
c     &                         +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
c             vyy = vyy + 0.5d0*(cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
c             vzy = vzy + 0.25d0*(cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
c     &                         +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
c               vzz =vzz +0.5d0*(cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))
c            end if
         end do
      !end if
c
c     increment the internal virial tensor components
c

      virep3b_recip(1,1) = virep3b_recip(1,1) + vxx
      virep3b_recip(2,1) = virep3b_recip(2,1) + vyx
      virep3b_recip(3,1) = virep3b_recip(3,1) + vzx
      virep3b_recip(1,2) = virep3b_recip(1,2) + vyx
      virep3b_recip(2,2) = virep3b_recip(2,2) + vyy
      virep3b_recip(3,2) = virep3b_recip(3,2) + vzy
      virep3b_recip(1,3) = virep3b_recip(1,3) + vzx
      virep3b_recip(2,3) = virep3b_recip(2,3) + vzy
      virep3b_recip(3,3) = virep3b_recip(3,3) + vzz

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
      !deallocate (cphi3b)
      return
      end
