c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c
      subroutine udirect1_3b_newbox1bclust(npole3b,pnum,field,cnt)
      use sizes
      use bound
      use boxes
      use boxes1bclust
      use ewald
      use math
      use mpole
      use pme
      implicit none
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 field(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
      integer npole3b,pnum(*),l1,cnt
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole3b))
      allocate (fmp(10,npole3b))
      allocate (cphi(10,npole3b))
      allocate (fphi(20,npole3b))
c
c     copy multipole moments and coordinates to local storage
c
      !do i = 1, npole
      do l1=1,npole3b
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
c     compute B-spline coefficients and spatial decomposition
c
      call bspline_fill_omp1bclust(cnt)
      call table_fill
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp3bclust1b (npole3b,pnum,cmp,fmp,cnt)
      !print*,"Completed cmp_to_fmp3b"
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole3b_new (npole3b,pnum,fmp)
      !print*,"Completed grid_mpole3b_new"
      call fftfront
      !print*,"Completed fftfront"
c
c     make the scalar summation over reciprocal lattice
c
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox1b(cnt)
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      ntot = nfft1 * nfft2 * nfft3
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
         h1 = recip1b(1,1,cnt)*r1 + recip1b(1,2,cnt)*r2 
     &        + recip1b(1,3,cnt)*r3
         h2 = recip1b(2,1,cnt)*r1 + recip1b(2,2,cnt)*r2 
     &        + recip1b(2,3,cnt)*r3
         h3 = recip1b(3,1,cnt)*r1 + recip1b(3,2,cnt)*r2 
     &        + recip1b(3,3,cnt)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox1b(cnt)*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox1b(cnt)
         qfac(1,1,1) = expterm
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
c     perform 3-D FFT backward transform and get field
c
      call fftback
      !print*,"Done w fftback"
      call fphi_mpole3b_new (npole3b,pnum,fphi)
      !print*,"Done w fphi_mpole3b_new"
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi3bclust1b (npole3b,pnum,fphi,cphi,cnt)
      !print*,"Done w fphi_to_cphi3b"

c
c     increment the field at each multipole site
c
      do i = 1, npole3b
         field(1,i) = field(1,i) - cphi(2,i)
         field(2,i) = field(2,i) - cphi(3,i)
         field(3,i) = field(3,i) - cphi(4,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end

