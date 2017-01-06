      subroutine ufield0c_3b_vacompbox1bclust(npole3b,pnum,field,fieldp,
     &   uind,uinp,cnt)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar, only: polarity, thole, pdamp
      implicit none
      integer i,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 uind(3,*),uinp(3,*)
      integer pnum(*),npole3b,l1,cnt
c
c
c     zero out the electrostatic field at each site
c
      do i = 1, npole3b
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     get the reciprocal space part of the electrostatic field
c
      
      call umutual1_3b_newbox1bclust(npole3b,pnum,field,fieldp,
     &      uind,uinp,cnt)
c
c     get the real space portion of the electrostatic field
c
      call umutual2b_3b_vacomp(npole3b,pnum,field,fieldp,
     &   uind,uinp)

c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      !do i = 1, npole
      do l1 = 1, npole3b
         do j = 1, 3
            field(j,l1) = field(j,l1) + term*uind(j,l1)
            fieldp(j,l1) = fieldp(j,l1) + term*uinp(j,l1)
         end do
      end do

      return
      end

c     #################################################################
c     ##                                                             ##
c
c
c
c
      subroutine umutual1_3b_newbox1bclust(npole3b,pnum,field,fieldp,
     &           uind,uinp,cnt)
      use sizes
      use boxes
      use boxes1bclust
      use ewald
      use math
      use mpole
      use pme
      use polar, only: polarity, thole, pdamp
      implicit none
      integer i,j,k
      real*8 term
      real*8 a(3,3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
      real*8 uind(3,*),uinp(3,*)
      integer pnum(*),npole3b,l1,cnt
c
c
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole3b))
      allocate (fuinp(3,npole3b))
      allocate (fdip_phi1(10,npole3b))
      allocate (fdip_phi2(10,npole3b))
      allocate (fdip_sum_phi(20,npole3b))
      allocate (dipfield1(3,npole3b))
      allocate (dipfield2(3,npole3b))

c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip1b(i,1,cnt)
         a(2,i) = dble(nfft2) * recip1b(i,2,cnt)
         a(3,i) = dble(nfft3) * recip1b(i,3,cnt)
      end do
      do i = 1, npole3b
         do k = 1, 3
            fuind(k,i) = a(k,1)*uind(1,i) + a(k,2)*uind(2,i)
     &                      + a(k,3)*uind(3,i)
            fuinp(k,i) = a(k,1)*uinp(1,i) + a(k,2)*uinp(2,i)
     &                      + a(k,3)*uinp(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind3b_new (npole3b,pnum,fuind,fuinp)
      call fftfront
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
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_uind3b_new (npole3b,pnum,
     &    fdip_phi1,fdip_phi2,fdip_sum_phi)
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip1b(i,1,cnt)
         a(i,2) = dble(nfft2) * recip1b(i,2,cnt)
         a(i,3) = dble(nfft3) * recip1b(i,3,cnt)
      end do

      do i = 1, npole3b
         do k = 1, 3
            dipfield1(k,i) = a(k,1)*fdip_phi1(2,i)
     &                          + a(k,2)*fdip_phi1(3,i)
     &                          + a(k,3)*fdip_phi1(4,i)
            dipfield2(k,i) = a(k,1)*fdip_phi2(2,i)
     &                          + a(k,2)*fdip_phi2(3,i)
     &                          + a(k,3)*fdip_phi2(4,i)
         end do
      end do
c
c     increment the field at each multipole site
c
      do i = 1, npole3b
         do k = 1, 3
            field(k,i) = field(k,i) - dipfield1(k,i)
            fieldp(k,i) = fieldp(k,i) - dipfield2(k,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end

