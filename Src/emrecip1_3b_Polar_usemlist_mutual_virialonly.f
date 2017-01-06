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
      subroutine emrecip1_3b_Polar_usemlist_mutual_virialonly
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

c      allocate (frc(3,npole))
c      allocate (trq(3,npole))
c      allocate (fuind(3,npole))
c      allocate (fuinp(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
c      allocate (fphi(20,npole))
c      allocate (fphid(10,npole))
c      allocate (fphip(10,npole))
c      allocate (fphidp(20,npole))

c      allocate (cphi(10,npole))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0

c       ep3b_recip=0.0d0

c       do i = 1, npole
c          dep3b_recip(1,i) = 0.0d0 
c          dep3b_recip(2,i) = 0.0d0 
c          dep3b_recip(3,i) = 0.0d0 
c       end do

       do i=1,3
         do j=1,3
          virep3b_recip_pt1(j,i) = 0.0d0
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
      call bspline_fill_omp
      print*,"After bsplin emrecip1_3b_Polar_usemlist_mutual_virialonly"

      call table_fill_omp
c      call table_fill
      print*,"After table emrecip1_3b_Polar_usemlist_mutual_virialonly"
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,nfft1,nfft2,nfft3))
      print*,"use_polar=",use_polar
c      print*,"uind in emrecip1_3b_Polar_mutual_virialonly",uind 
      if (use_polar) then
         do i = 1, npole
            do j = 2, 4
               cmp(j,i) = cmp(j,i) + uinp(j-1,i)
            end do
         end do
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
      print*,"After cmp_to_fmp,grid_mpole,fftfront" 
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
         do i = 1, npole
            do j = 2, 4
               cmp(j,i) = cmp(j,i) + uind(j-1,i) - uinp(j-1,i)
            end do
         end do
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
         do i = 1, npole
            do j = 2, 4
               cmp(j,i) = cmp(j,i) - uind(j-1,i)
            end do
         end do
      else
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
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
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
c
c     make the scalar summation over reciprocal lattice
c

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
C    VIRIAL CONTRIBUTION FROM PERMANENT MULTIPOLES SAVED FROM
C    BEFORE.
c          call cmp_to_fmp(cmp,fmp)

c        vxx = vxx - vxx_perm
c        vyx = vyx - vyx_perm
c        vzx = vzx - vzx_perm
c        vyy = vyy - vyy_perm
c        vzy = vzy - vzy_perm
c        vzz = vzz - vzz_perm
      virep3b_recip_pt1(1,1) = virep3b_recip_pt1(1,1) + vxx
      virep3b_recip_pt1(2,1) = virep3b_recip_pt1(2,1) + vyx
      virep3b_recip_pt1(3,1) = virep3b_recip_pt1(3,1) + vzx
      virep3b_recip_pt1(1,2) = virep3b_recip_pt1(1,2) + vyx
      virep3b_recip_pt1(2,2) = virep3b_recip_pt1(2,2) + vyy
      virep3b_recip_pt1(3,2) = virep3b_recip_pt1(3,2) + vzy
      virep3b_recip_pt1(1,3) = virep3b_recip_pt1(1,3) + vzx
      virep3b_recip_pt1(2,3) = virep3b_recip_pt1(2,3) + vzy
      virep3b_recip_pt1(3,3) = virep3b_recip_pt1(3,3) + vzz


c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)


c
c     perform deallocation of some local arrays
c
c      deallocate (frc)
c      deallocate (trq)
c      deallocate (fuind)
c      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
c      deallocate (fphi)
c      deallocate (fphid)
c      deallocate (fphip)
c      deallocate (fphidp)
c      deallocate (cphi)
      !deallocate (cphi3b)
      return
      end
