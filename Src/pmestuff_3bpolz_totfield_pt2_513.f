c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grid_uind  --  put induced dipoles on PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grid_uind" places the fractional induced dipoles onto the
c     particle mesh Ewald grid
c
c
      subroutine grid_uind3b_totfield_513 (npole3b,pnum,fuind,fuinp,
     &         qgrid3b)
      use sizes
      use atoms
      use chunks
      use mpole  
c      use pme, only: nfft1,nfft2,nfft3,maxorder,bsorder,igrid,thetai1,
c     & thetai2,thetai3
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 term01,term11
      real*8 term02,term12
      real*8 fuind(3,*)
      real*8 fuinp(3,*)
      integer npole3b,pnum(*),l1
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      integer k1
      logical doi

c
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid3b(1,i,j,k) = 0.0d0
               qgrid3b(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,u0,u1,term01,term11,term02,term12,t0,t1,l1,
!$OMP& doi,k1)
!$OMP DO
c
c     put the induced dipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, npole
         !do l1=1,npole3b
         !   isite=pnum(l1)
            iatm = ipole(isite)
            doi=.false.
            do k1=1,npole3b
              if(pnum(k1).eq.isite) then
                doi=.true.
                l1=k1
                goto 24
              end if
            end do
   24    continue

            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
c                     term01 = fuind(2,isite)*u1*v0
c     &                           + fuind(3,isite)*u0*v1
c                     term11 = fuind(1,isite)*u0*v0
c                     term02 = fuinp(2,isite)*u1*v0
c     &                           + fuinp(3,isite)*u0*v1
c                     term12 = fuinp(1,isite)*u0*v0
                     if(doi) then
                     term01 = fuind(2,l1)*u1*v0
     &                           + fuind(3,l1)*u0*v1
                     term11 = fuind(1,l1)*u0*v0
                     term02 = fuinp(2,l1)*u1*v0
     &                           + fuinp(3,l1)*u0*v1
                     term12 = fuinp(1,l1)*u0*v0
                     else
                     term01 = 0.0d0 
                     term11 = 0.0d0
                     term02 = 0.0d0
                     term12 = 0.0d0   
                     end if
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
                        qgrid3b(1,i,j,k) = qgrid3b(1,i,j,k) + term01*t0
     &                                      + term11*t1
                        qgrid3b(2,i,j,k) = qgrid3b(2,i,j,k) + term02*t0
     &                                      + term12*t1
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fphi_uind  --  induced potential from grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fphi_uind" extracts the induced dipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uind3b_totfield_513 (npole3b,pnum,fdip_phi1,
     & fdip_phi2,fdip_sum_phi,qgrid3b)
      use sizes
      use mpole
c      use pme, only: nfft1,nfft2,nfft3,maxorder,bsorder,igrid,thetai1,
c     & thetai2,thetai3
      use pme
      use atoms
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
      real*8 fdip_phi1(10,*)
      real*8 fdip_phi2(10,*)
      real*8 fdip_sum_phi(20,*)
      integer npole3b,pnum(*),l1
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      integer k1
      logical doi

c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,
!$OMP& igrid,bsorder,nfft3,thetai3,nfft2,thetai2,nfft1,
!$OMP& thetai1,qgrid3b,fdip_phi1,fdip_phi2,fdip_sum_phi,
!$OMP& npole3b,pnum)
!$OMP DO
c
c     extract the induced dipole field at each site
c
      do isite = 1, npole
      !do l1=1,npole3b
       !  isite=pnum(l1)

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
                  tq_1 = qgrid3b(1,i,j,k)
                  tq_2 = qgrid3b(2,i,j,k)
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

         !fdip_phi1(2,l1) = tuv100_1
         !fdip_phi1(3,l1) = tuv010_1
         !fdip_phi1(4,l1) = tuv001_1
         !fdip_phi1(5,l1) = tuv200_1
         !fdip_phi1(6,l1) = tuv020_1
         !fdip_phi1(7,l1) = tuv002_1
         !fdip_phi1(8,l1) = tuv110_1
         !fdip_phi1(9,l1) = tuv101_1
         !fdip_phi1(10,l1) = tuv011_1
         !fdip_phi2(2,l1) = tuv100_2
         !fdip_phi2(3,l1) = tuv010_2
         !fdip_phi2(4,l1) = tuv001_2
         !fdip_phi2(5,l1) = tuv200_2
         !fdip_phi2(6,l1) = tuv020_2
         !fdip_phi2(7,l1) = tuv002_2
         !fdip_phi2(8,l1) = tuv110_2
         !fdip_phi2(9,l1) = tuv101_2
         !fdip_phi2(10,l1) = tuv011_2
         !fdip_sum_phi(1,l1) = tuv000
         !fdip_sum_phi(2,l1) = tuv100
         !fdip_sum_phi(3,l1) = tuv010
         !fdip_sum_phi(4,l1) = tuv001
         !fdip_sum_phi(5,l1) = tuv200
         !fdip_sum_phi(6,l1) = tuv020
         !fdip_sum_phi(7,l1) = tuv002
         !fdip_sum_phi(8,l1) = tuv110
         !fdip_sum_phi(9,l1) = tuv101
         !fdip_sum_phi(10,l1) = tuv011
         !fdip_sum_phi(11,l1) = tuv300
         !fdip_sum_phi(12,l1) = tuv030
         !fdip_sum_phi(13,l1) = tuv003
         !fdip_sum_phi(14,l1) = tuv210
         !fdip_sum_phi(15,l1) = tuv201
         !fdip_sum_phi(16,l1) = tuv120
         !fdip_sum_phi(17,l1) = tuv021
         !fdip_sum_phi(18,l1) = tuv102
         !fdip_sum_phi(19,l1) = tuv012
         !fdip_sum_phi(20,l1) = tuv111

      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ############################################################
c     ##                                                        ##
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_to_cphi  --  transformation of potential  ##
c     ##                                                            ##
c     ################################################################
c
c
