c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns particle mesh Ewald parameters and options
c     for a periodic system
c
c
      subroutine kewald2bclust(cnt)
      use sizes
      use atoms
      use bound
      use boxes
      use boxes2bclust
      use chunks
      use ewald
      use fft
      use inform
      use iounit
      use keys
      use limits
      use openmp
      use pme
      implicit none
      integer maxpower
      parameter (maxpower=54)
      integer i,k,next
      integer ifft1,ifft2,ifft3
      integer multi(maxpower)
      real*8 delta,rmax,dens
      character*20 keyword
      character*120 record
      character*120 string
      integer cnt
c
c     PME grid size must be even with factors of only 2, 3 and 5
c
      data multi  /   2,   4,   6,   8,  10,  12,  16,  18,  20,
     &               24,  30,  32,  36,  40,  48,  50,  54,  60,
     &               64,  72,  80,  90,  96, 100, 108, 120, 128,
     &              144, 150, 160, 162, 180, 192, 200, 216, 240,
     &              250, 256, 270, 288, 300, 320, 324, 360, 384,
     &              400, 432, 450, 480, 486, 500, 512, 540, 576 /
c
c
c     default boundary treatment, B-spline order and grid density
c
c      ffttyp = 'FFTPACK'
c      if (nthread .gt. 1)  ffttyp = 'FFTW'
c      ffttyp = 'FFTW'
c      boundary = 'TINFOIL'
c      bsorder = 5
      dens = 1.2d0
c
c     estimate an optimal value for the Ewald coefficient
c
c      call ewaldcof (aewald,ewaldcut)
c
c     set the system extent for nonperiodic Ewald summation
c
c      if (.not. use_bounds) then
c         call extent (rmax)
c         xbox = 2.0d0 * (rmax+ewaldcut)
c         ybox = xbox
c         zbox = xbox
c         alpha = 90.0d0
c         beta = 90.0d0
c         gamma = 90.0d0
c         orthogonal = .true.
c         call lattice
c         boundary = 'NONE'
c         dens = 0.75d0
c      end if
c
c     get default grid size from periodic system dimensions
c
      delta = 1.0d-8
      ifft1 = int(xbox2b(cnt)*dens-delta) + 1
      ifft2 = int(ybox2b(cnt)*dens-delta) + 1
      ifft3 = int(zbox2b(cnt)*dens-delta) + 1
c
c     search keywords for Ewald summation commands
c
c      do i = 1, nkey
c         record = keyline(i)
c         next = 1
c         call upcase (record)
c         call gettext (record,keyword,next)
c         string = record(next:120)
c         if (keyword(1:12) .eq. 'FFT-PACKAGE ') then
c            call getword (record,ffttyp,next)
c         else if (keyword(1:12) .eq. 'EWALD-ALPHA ') then
c            read (string,*,err=20)  aewald
c         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
c            boundary = 'VACUUM'
c         else if (keyword(1:9) .eq. 'PME-GRID ') then
c            ifft1 = 0
c            ifft2 = 0
c            ifft3 = 0
c            read (string,*,err=10,end=10)  ifft1,ifft2,ifft3
c   10       continue
c            if (ifft2 .eq. 0)  ifft2 = ifft1
c            if (ifft3 .eq. 0)  ifft3 = ifft1
c         else if (keyword(1:10) .eq. 'PME-ORDER ') then
c            read (string,*,err=20)  bsorder
c         end if
c   20    continue
c      end do
c
c     grid size must be even, with prime factors of 2, 3 and 5
c
      nfft1 = maxfft
      nfft2 = maxfft
      nfft3 = maxfft
      do i = maxpower, 1, -1
         k = multi(i)
         if (k .le. maxfft) then
            if (k .ge. ifft1)  nfft1 = k
            if (k .ge. ifft2)  nfft2 = k
            if (k .ge. ifft3)  nfft3 = k
         end if
      end do
c
c     set the number of chunks and grid points per chunk
c
      call getchunk
c
c     check the B-spline order and charge grid dimension
c
c      if (bsorder .gt. maxorder) then
c         write (iout,30)
c   30    format (/,' KEWALD  --  B-Spline Order Too Large;',
c     &              ' Increase MAXORDER')
c         call fatal
c      end if
      if (max(nfft1,nfft2,nfft3) .gt. maxfft) then
         write (iout,40)
   40    format (/,' KEWALD  --  FFT Charge Grid Too Large;',
     &              ' Increase MAXFFT')
         call fatal
      else if (nfft1.lt.ifft1 .or. nfft2.lt.ifft2
     &             .or. nfft3.lt.ifft3) then
         write (iout,50)
   50    format (/,' KEWALD  --  Warning, Small Charge Grid',
     &              ' may give Poor Accuracy')
      end if
c
c     perform dynamic allocation of some global arrays
c
c    THETAI1,2,3 BELOW HAVE ALREADY BEEN ALLOCATED AND DO NOT NEED TO BE RESIZED.   
      if (.not. allocated(thetai1))  allocate (thetai1(4,bsorder,n))
      if (.not. allocated(thetai2))  allocate (thetai2(4,bsorder,n))
      if (.not. allocated(thetai3))  allocate (thetai3(4,bsorder,n))  

C   QGRID,QFAC,PMETABLE MUST BE RESIZED BELOW ACCORDING TO THE SUBSYSTEM-SPECIFIC PME GRID.
      if(allocated(qgrid)) deallocate(qgrid)
      if(allocated(qfac)) deallocate(qfac)
      if(allocated(pmetable)) deallocate (pmetable)

      if (.not. allocated(qgrid))  allocate (qgrid(2,nfft1,nfft2,nfft3))
      if (.not. allocated(qfac))  allocate (qfac(nfft1,nfft2,nfft3))
      if (.not. allocated(pmetable))  allocate (pmetable(n,nchunk))
c
c     initialize the PME arrays that can be precomputed
c
      call moduli
      call fftsetup
c
c     print a message listing some of the Ewald parameters
c
      !if (verbose) then
      !   write (iout,60)  aewald,nfft1,nfft2,nfft3,bsorder
   60    format (/,' Smooth Particle Mesh Ewald Parameters :',
     &           //,4x,'Ewald Coefficient',6x,'Charge Grid',
     &              ' Dimensions',6x,'B-Spline Order',
     &           //,8x,f8.4,11x,3i6,12x,i6)
      !end if
      return
      end
c
