      subroutine ufield0c_3b_vacomp(npole3b,pnum,field,fieldp,
     &   uind,uinp)
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
      integer pnum(*),npole3b,l1
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
      
      call umutual1_3b_new (npole3b,pnum,field,fieldp,uind,uinp)
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
c     ##  subroutine ufield0a_3b  --  mutual induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ufield0a_3b" computes the mutual electrostatic field due to
c     induced dipole moments via a double loop
c
c
      subroutine umutual2b_3b_vacomp(npole3b,pnum,field,fieldp,
     &   uind,uinp)
      use sizes
      use atoms
      use bound
      use cell
      use group
      use mpole
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use shunt
      use tarray
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8, allocatable :: dscale(:)
      logical proceed
      character*6 mode
      integer npole3b,pnum(*),l1,l3
      real*8 field(3,npole3b)
      real*8 fieldp(3,npole3b)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)

c
c
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole3b))
      allocate (fieldtp(3,npole3b))

c
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,uind,uinp,ntpair,tindex,
!$OMP& tdipdip,field,fieldp,fieldt,fieldtp,npole3b)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole3b
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
      do m = 1, ntpair
         i = tindex(1,m)
         k = tindex(2,m)
         fimd(1) = tdipdip(1,m)*uind(1,k) + tdipdip(2,m)*uind(2,k)
     &                + tdipdip(3,m)*uind(3,k)
         fimd(2) = tdipdip(2,m)*uind(1,k) + tdipdip(4,m)*uind(2,k)
     &                + tdipdip(5,m)*uind(3,k)
         fimd(3) = tdipdip(3,m)*uind(1,k) + tdipdip(5,m)*uind(2,k)
     &                + tdipdip(6,m)*uind(3,k)
         fkmd(1) = tdipdip(1,m)*uind(1,i) + tdipdip(2,m)*uind(2,i)
     &                + tdipdip(3,m)*uind(3,i)
         fkmd(2) = tdipdip(2,m)*uind(1,i) + tdipdip(4,m)*uind(2,i)
     &                + tdipdip(5,m)*uind(3,i)
         fkmd(3) = tdipdip(3,m)*uind(1,i) + tdipdip(5,m)*uind(2,i)
     &                + tdipdip(6,m)*uind(3,i)
         fimp(1) = tdipdip(1,m)*uinp(1,k) + tdipdip(2,m)*uinp(2,k)
     &                + tdipdip(3,m)*uinp(3,k)
         fimp(2) = tdipdip(2,m)*uinp(1,k) + tdipdip(4,m)*uinp(2,k)
     &                + tdipdip(5,m)*uinp(3,k)
         fimp(3) = tdipdip(3,m)*uinp(1,k) + tdipdip(5,m)*uinp(2,k)
     &                + tdipdip(6,m)*uinp(3,k)
         fkmp(1) = tdipdip(1,m)*uinp(1,i) + tdipdip(2,m)*uinp(2,i)
     &                + tdipdip(3,m)*uinp(3,i)
         fkmp(2) = tdipdip(2,m)*uinp(1,i) + tdipdip(4,m)*uinp(2,i)
     &                + tdipdip(5,m)*uinp(3,i)
         fkmp(3) = tdipdip(3,m)*uinp(1,i) + tdipdip(5,m)*uinp(2,i)
     &                + tdipdip(6,m)*uinp(3,i)
c
c     increment the field at each site due to this interaction
c
         do j = 1, 3
            fieldt(j,i) = fieldt(j,i) + fimd(j)
            fieldt(j,k) = fieldt(j,k) + fkmd(j)
            fieldtp(j,i) = fieldtp(j,i) + fimp(j)
            fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
         end do
      end do
!$OMP END DO
c
c     end OpenMP directives for the major loop structure
c
!$OMP DO
      do i = 1, npole3b
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c
      subroutine umutual1_3b_new(npole3b,pnum,field,fieldp,uind,uinp)
      use sizes
      use boxes
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
      integer pnum(*),npole3b,l1
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
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
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
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
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

