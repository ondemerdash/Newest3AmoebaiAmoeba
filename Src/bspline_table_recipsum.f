c      subroutine udirect1_3b (npole3b,pnum,field,thetai1_3b,
c     &     thetai2_3b,thetai3_3b,qgrid3b,qfac3b,pmetable3b,igrid3b,
c     &     planf3b,planb3b,iprime3b,ffttable3b)
      subroutine bspline_table_recipsum (thetai1_3b,
     &     thetai2_3b,thetai3_3b,qfac3b,pmetable3b,igrid3b
     &     )
      use sizes
      use bound
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use chunks
      use atoms
      use fft
c      use pme, only: nfft1,nfft2,nfft3,maxorder,bsorder,igrid,thetai1,
c     & thetai2,thetai3,bsmod1,bsmod2,bsmod3
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
c      real*8 field(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
c      integer npole3b,pnum(*),l1,l3
c      real*8 field(3,npole3b)
c      real*8 qgrid(2,nfft1,nfft2,nfft3)
c      real*8 qfac(nfft1,nfft2,nfft3)
c      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      real*8 qfac3b(nfft1,nfft2,nfft3)
      real*8 thetai1_3b(4,bsorder,n)
      real*8 thetai2_3b(4,bsorder,n)
      real*8 thetai3_3b(4,bsorder,n)
      integer pmetable3b(n,nchunk)
      integer igrid3b(3,n)
c      integer*8 planf3b
c      integer*8 planb3b
c      integer iprime3b(maxprime,3)
c      real*8 ffttable3b(maxtable,3)

      call bspline_fill3b(igrid3b,thetai1_3b,thetai2_3b,
     &         thetai3_3b)
      !print*,"successful completion of bspline_fill3b" 
      call table_fill3b(igrid3b,pmetable3b)
      !print*,"successful completion of table_fill3b"

      qfac3b(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
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
         end if
         qfac3b(k1,k2,k3) = expterm
      end do

      return
      end

