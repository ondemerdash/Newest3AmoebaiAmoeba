c
      subroutine fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
      use sizes
      use mpole
      use pme
      use mpiparams
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3
      real*8 t0_1,t0_2,t1_1,t1_2
      real*8 t2_1,t2_2,tq_1,tq_2
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu30,tu21,tu12,tu03
      real*8 tu00_1,tu01_1,tu10_1
      real*8 tu00_2,tu01_2,tu10_2
      real*8 tu20_1,tu11_1,tu02_1
      real*8 tu20_2,tu11_2,tu02_2
      real*8 tuv100_1,tuv010_1,tuv001_1
      real*8 tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_1,tuv020_1,tuv002_1
      real*8 tuv110_1,tuv101_1,tuv011_1
      real*8 tuv200_2,tuv020_2,tuv002_2
      real*8 tuv110_2,tuv101_2,tuv011_2
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fdip_phi1(10,npole)
      real*8 fdip_phi2(10,npole)
      real*8 fdip_sum_phi(20,npole)
      integer lstart, lend

      call splitloop(lstart, lend, npole)

      do i = 1, npole
         do k = 1, 10
            fdip_phi1(k,i) = 0.0d0
            fdip_phi2(k,i) = 0.0d0
         end do
      end do

c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,
!$OMP& igrid,bsorder,nfft3,thetai3,nfft2,thetai2,nfft1,
!$OMP& thetai1,qgrid,fdip_phi1,fdip_phi2,fdip_sum_phi)
!$OMP DO
c
c     extract the induced dipole field at each site
c
      do isite = lstart, lend  
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         tuv100_1 = 0.0d0
         tuv010_1 = 0.0d0
         tuv001_1 = 0.0d0
         tuv200_1 = 0.0d0
         tuv020_1 = 0.0d0
         tuv002_1 = 0.0d0
         tuv110_1 = 0.0d0
         tuv101_1 = 0.0d0
         tuv011_1 = 0.0d0
         tuv100_2 = 0.0d0
         tuv010_2 = 0.0d0
         tuv001_2 = 0.0d0
         tuv200_2 = 0.0d0
         tuv020_2 = 0.0d0
         tuv002_2 = 0.0d0
         tuv110_2 = 0.0d0
         tuv101_2 = 0.0d0
         tuv011_2 = 0.0d0
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)
            tu00_1 = 0.0d0
            tu01_1 = 0.0d0
            tu10_1 = 0.0d0
            tu20_1 = 0.0d0
            tu11_1 = 0.0d0
            tu02_1 = 0.0d0
            tu00_2 = 0.0d0
            tu01_2 = 0.0d0
            tu10_2 = 0.0d0
            tu20_2 = 0.0d0
            tu11_2 = 0.0d0
            tu02_2 = 0.0d0
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)
               t0_1 = 0.0d0
               t1_1 = 0.0d0
               t2_1 = 0.0d0
               t0_2 = 0.0d0
               t1_2 = 0.0d0
               t2_2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq_1 = qgrid(1,i,j,k)
                  tq_2 = qgrid(2,i,j,k)
                  t0_1 = t0_1 + tq_1*thetai1(1,it1,iatm)
                  t1_1 = t1_1 + tq_1*thetai1(2,it1,iatm)
                  t2_1 = t2_1 + tq_1*thetai1(3,it1,iatm)
                  t0_2 = t0_2 + tq_2*thetai1(1,it1,iatm)
                  t1_2 = t1_2 + tq_2*thetai1(2,it1,iatm)
                  t2_2 = t2_2 + tq_2*thetai1(3,it1,iatm)
                  t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,iatm)
               end do
               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2
               t0 = t0_1 + t0_2
               t1 = t1_1 + t1_2
               t2 = t2_1 + t2_2
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
            tuv100_1 = tuv100_1 + tu10_1*v0
            tuv010_1 = tuv010_1 + tu01_1*v0
            tuv001_1 = tuv001_1 + tu00_1*v1
            tuv200_1 = tuv200_1 + tu20_1*v0
            tuv020_1 = tuv020_1 + tu02_1*v0
            tuv002_1 = tuv002_1 + tu00_1*v2
            tuv110_1 = tuv110_1 + tu11_1*v0
            tuv101_1 = tuv101_1 + tu10_1*v1
            tuv011_1 = tuv011_1 + tu01_1*v1
            tuv100_2 = tuv100_2 + tu10_2*v0
            tuv010_2 = tuv010_2 + tu01_2*v0
            tuv001_2 = tuv001_2 + tu00_2*v1
            tuv200_2 = tuv200_2 + tu20_2*v0
            tuv020_2 = tuv020_2 + tu02_2*v0
            tuv002_2 = tuv002_2 + tu00_2*v2
            tuv110_2 = tuv110_2 + tu11_2*v0
            tuv101_2 = tuv101_2 + tu10_2*v1
            tuv011_2 = tuv011_2 + tu01_2*v1
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fdip_phi1(2,isite) = tuv100_1
         fdip_phi1(3,isite) = tuv010_1
         fdip_phi1(4,isite) = tuv001_1
         fdip_phi1(5,isite) = tuv200_1
         fdip_phi1(6,isite) = tuv020_1
         fdip_phi1(7,isite) = tuv002_1
         fdip_phi1(8,isite) = tuv110_1
         fdip_phi1(9,isite) = tuv101_1
         fdip_phi1(10,isite) = tuv011_1
         fdip_phi2(2,isite) = tuv100_2
         fdip_phi2(3,isite) = tuv010_2
         fdip_phi2(4,isite) = tuv001_2
         fdip_phi2(5,isite) = tuv200_2
         fdip_phi2(6,isite) = tuv020_2
         fdip_phi2(7,isite) = tuv002_2
         fdip_phi2(8,isite) = tuv110_2
         fdip_phi2(9,isite) = tuv101_2
         fdip_phi2(10,isite) = tuv011_2
         fdip_sum_phi(1,isite) = tuv000
         fdip_sum_phi(2,isite) = tuv100
         fdip_sum_phi(3,isite) = tuv010
         fdip_sum_phi(4,isite) = tuv001
         fdip_sum_phi(5,isite) = tuv200
         fdip_sum_phi(6,isite) = tuv020
         fdip_sum_phi(7,isite) = tuv002
         fdip_sum_phi(8,isite) = tuv110
         fdip_sum_phi(9,isite) = tuv101
         fdip_sum_phi(10,isite) = tuv011
         fdip_sum_phi(11,isite) = tuv300
         fdip_sum_phi(12,isite) = tuv030
         fdip_sum_phi(13,isite) = tuv003
         fdip_sum_phi(14,isite) = tuv210
         fdip_sum_phi(15,isite) = tuv201
         fdip_sum_phi(16,isite) = tuv120
         fdip_sum_phi(17,isite) = tuv021
         fdip_sum_phi(18,isite) = tuv102
         fdip_sum_phi(19,isite) = tuv012
         fdip_sum_phi(20,isite) = tuv111
      end do

c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL

      return
      end
