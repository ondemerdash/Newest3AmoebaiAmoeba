c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine grid_mpole  --  put multipoles on PME grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "grid_mpole" places the fractional atomic multipoles onto
c     the particle mesh Ewald grid
c
c
c      subroutine grid_mpole3b (npole3b,pnum,fmp,qgrid)
      subroutine grid_mpole3b_totfield (fmp,qgrid3b,thetai1_3b,
     &     thetai2_3b,thetai3_3b,pmetable3b,igrid3b)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
c      use pme, only: nfft1,nfft2,nfft3,maxorder,bsorder,igrid,thetai1,
c     & thetai2,thetai3
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
      real*8 v2,u2,t2
      real*8 term0,term1,term2
      real*8 fmp(10,*)
      !integer npole3b,pnum(*),l1
c      real*8 qgrid(2,nfft1,nfft2,nfft3)
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)       
      real*8 thetai1_3b(4,bsorder,n)
      real*8 thetai2_3b(4,bsorder,n)
      real*8 thetai3_3b(4,bsorder,n)
      integer igrid3b(3,n)
      integer pmetable3b(n,nchunk)

c
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
c               qgrid(1,i,j,k) = 0.0d0
c               qgrid(2,i,j,k) = 0.0d0
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
!$OMP& offsetz,v0,v1,v2,u0,u1,u2,term0,term1,term2,t0,t1,t2)
!$OMP DO
c
c     put the permanent multipole moments onto the grid
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
c            if (pmetable(iatm,ichk) .eq. 1) then
            if (pmetable3b(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid3b(1,iatm) + grdoff
               nearpt(2) = igrid3b(2,iatm) + grdoff
               nearpt(3) = igrid3b(3,iatm) + grdoff
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
c                  v0 = thetai3(1,m,iatm)
c                  v1 = thetai3(2,m,iatm)
c                  v2 = thetai3(3,m,iatm)
                  v0 = thetai3_3b(1,m,iatm)
                  v1 = thetai3_3b(2,m,iatm)
                  v2 = thetai3_3b(3,m,iatm)

                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
c                     u0 = thetai2(1,m,iatm)
c                     u1 = thetai2(2,m,iatm)
c                     u2 = thetai2(3,m,iatm)
                     u0 = thetai2_3b(1,m,iatm)
                     u1 = thetai2_3b(2,m,iatm)
                     u2 = thetai2_3b(3,m,iatm)

                     term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
     &                     + fmp(4,isite)*u0*v1 + fmp(6,isite)*u2*v0
     &                     + fmp(7,isite)*u0*v2 + fmp(10,isite)*u1*v1
                     term1 = fmp(2,isite)*u0*v0 + fmp(8,isite)*u1*v0
     &                          + fmp(9,isite)*u0*v1
                     term2 = fmp(5,isite) * u0 * v0

                     !term0 = fmp(1,l1)*u0*v0 + fmp(3,l1)*u1*v0
     &               !      + fmp(4,l1)*u0*v1 + fmp(6,l1)*u2*v0
     &               !      + fmp(7,l1)*u0*v2 + fmp(10,l1)*u1*v1
                     !term1 = fmp(2,l1)*u0*v0 + fmp(8,l1)*u1*v0
     &               !           + fmp(9,l1)*u0*v1
                     !term2 = fmp(5,l1) * u0 * v0


                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
c                        t0 = thetai1(1,m,iatm)
c                        t1 = thetai1(2,m,iatm)
c                        t2 = thetai1(3,m,iatm)
                        t0 = thetai1_3b(1,m,iatm)
                        t1 = thetai1_3b(2,m,iatm)
                        t2 = thetai1_3b(3,m,iatm)

c                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
c     &                                      + term1*t1 + term2*t2
                        qgrid3b(1,i,j,k) = qgrid3b(1,i,j,k) + term0*t0
     &                                      + term1*t1 + term2*t2
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_mpole  --  multipole potential from grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_mpole" extracts the permanent multipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_mpole3b_totfield(fphi,qgrid3b,
     & igrid3b,thetai1_3b,thetai2_3b,thetai3_3b)
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
      real*8 t0,t1,t2,t3,tq
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fphi(20,*)
c      integer npole3b,pnum(*),l1
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      integer igrid3b(3,n)
      real*8 thetai1_3b(4,bsorder,n)
      real*8 thetai2_3b(4,bsorder,n)
      real*8 thetai3_3b(4,bsorder,n)
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,igrid3b,bsorder,
!$OMP& nfft3,thetai3_3b,nfft2,thetai2_3b,nfft1,thetai1_3b,qgrid3b,fphi)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
      do isite = 1, npole
c      do l1=1,npole3b
c         isite=pnum(l1)
         iatm = ipole(isite)
         igrd0 = igrid3b(1,iatm)
         jgrd0 = igrid3b(2,iatm)
         kgrd0 = igrid3b(3,iatm)
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
            v0 = thetai3_3b(1,it3,iatm)
            v1 = thetai3_3b(2,it3,iatm)
            v2 = thetai3_3b(3,it3,iatm)
            v3 = thetai3_3b(4,it3,iatm)
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
               u0 = thetai2_3b(1,it2,iatm)
               u1 = thetai2_3b(2,it2,iatm)
               u2 = thetai2_3b(3,it2,iatm)
               u3 = thetai2_3b(4,it2,iatm)
               t0 = 0.0d0
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq = qgrid3b(1,i,j,k)
                  t0 = t0 + tq*thetai1_3b(1,it1,iatm)
                  t1 = t1 + tq*thetai1_3b(2,it1,iatm)
                  t2 = t2 + tq*thetai1_3b(3,it1,iatm)
                  t3 = t3 + tq*thetai1_3b(4,it1,iatm)
               end do
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
         fphi(1,isite) = tuv000
         fphi(2,isite) = tuv100
         fphi(3,isite) = tuv010
         fphi(4,isite) = tuv001
         fphi(5,isite) = tuv200
         fphi(6,isite) = tuv020
         fphi(7,isite) = tuv002
         fphi(8,isite) = tuv110
         fphi(9,isite) = tuv101
         fphi(10,isite) = tuv011
         fphi(11,isite) = tuv300
         fphi(12,isite) = tuv030
         fphi(13,isite) = tuv003
         fphi(14,isite) = tuv210
         fphi(15,isite) = tuv201
         fphi(16,isite) = tuv120
         fphi(17,isite) = tuv021
         fphi(18,isite) = tuv102
         fphi(19,isite) = tuv012
         fphi(20,isite) = tuv111

      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
