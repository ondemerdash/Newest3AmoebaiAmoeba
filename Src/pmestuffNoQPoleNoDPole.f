      subroutine cmp_to_fmp_noqpole_nodpole (cmp,fmp)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ctf(10,10)
c      real*8 cmp(10,*)
c      real*8 fmp(10,*)
c      real*8 cmp(4,*)
c      real*8 fmp(4,*)
      real*8 cmp(*)
      real*8 fmp(*)

c
c
c     find the matrix to convert Cartesian to fractional
c
c      call cart_to_frac (ctf)
c
c     apply the transformation to get the fractional multipoles
c
      do i = 1, npole
        ! ctf(1,1) is set to 1.0d0 in cart_to_frac above
         fmp(i) = cmp(i)
        ! fmp(1,i) = ctf(1,1) * cmp(1,i)
        ! do j = 2, 4
        !    fmp(j,i) = 0.0d0
        !    do k = 2, 4
        !       fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
        !    end do
        ! end do
      !   do j = 5, 10
      !      fmp(j,i) = 0.0d0
      !      do k = 5, 10
      !         fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
      !      end do
      !   end do
      end do
      return
      end

      subroutine grid_mpole_noqpole_nodpole (cmp)
      use sizes
      use atoms
      use chunks
      use mpole
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
      real*8 v2,u2,t2
      real*8 term0,term1,term2
c      real*8 fmp(10,*)
c      real*8 fmp(4,*)
      real*8 cmp(*) 
c
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
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
            iatm = ipole(isite)
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
                 ! v1 = thetai3(2,m,iatm)
                 ! v2 = thetai3(3,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                   !  u1 = thetai2(2,m,iatm)
                   !  u2 = thetai2(3,m,iatm)
c                     term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
c     &                     + fmp(4,isite)*u0*v1 + fmp(6,isite)*u2*v0
c     &                     + fmp(7,isite)*u0*v2 + fmp(10,isite)*u1*v1
c                     term1 = fmp(2,isite)*u0*v0 + fmp(8,isite)*u1*v0
c     &                          + fmp(9,isite)*u0*v1
c                     term2 = fmp(5,isite) * u0 * v0
                    ! term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
     &              !       + fmp(4,isite)*u0*v1 
                    ! term1 = fmp(2,isite)*u0*v0 
                    term0 = cmp(isite)*u0*v0

                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                       !   t1 = thetai1(2,m,iatm)
                       !   t2 = thetai1(3,m,iatm)
c                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
c     &                                      + term1*t1 + term2*t2
                      !  qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
     &                !                      + term1*t1
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
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

      subroutine fphi_mpole_noqpole_nodpole (fphi)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
c      real*8 v0,v1,v2,v3
c      real*8 u0,u1,u2,u3
      real*8 v0,v1
      real*8 u0,u1
c      real*8 t0,t1,t2,t3,tq
      real*8 t0,t1,tq
c      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu00,tu10,tu01
c      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
c      real*8 tuv200,tuv020,tuv002,tuv110
c      real*8 tuv101,tuv011,tuv300,tuv030
c      real*8 tuv003,tuv210,tuv201,tuv120
c      real*8 tuv021,tuv102,tuv012,tuv111
c      real*8 fphi(20,*)
c      real*8 fphi(10,*)
      real*8 fphi(4,*)
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,igrid,bsorder,
!$OMP& nfft3,thetai3,nfft2,thetai2,nfft1,thetai1,qgrid,fphi)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
        ! tuv200 = 0.0d0
        ! tuv020 = 0.0d0
        ! tuv002 = 0.0d0
        ! tuv110 = 0.0d0
        ! tuv101 = 0.0d0
        ! tuv011 = 0.0d0
c         tuv300 = 0.0d0
c         tuv030 = 0.0d0
c         tuv003 = 0.0d0
c         tuv210 = 0.0d0
c         tuv201 = 0.0d0
c         tuv120 = 0.0d0
c         tuv021 = 0.0d0
c         tuv102 = 0.0d0
c         tuv012 = 0.0d0
c         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            !v2 = thetai3(3,it3,iatm)
c            v3 = thetai3(4,it3,iatm)
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
           ! tu20 = 0.0d0
           ! tu11 = 0.0d0
           ! tu02 = 0.0d0
c            tu30 = 0.0d0
c            tu21 = 0.0d0
c            tu12 = 0.0d0
c            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
             !  u2 = thetai2(3,it2,iatm)
c               u3 = thetai2(4,it2,iatm)
               t0 = 0.0d0
               t1 = 0.0d0
              ! t2 = 0.0d0
              ! t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq = qgrid(1,i,j,k)
                  t0 = t0 + tq*thetai1(1,it1,iatm)
                  t1 = t1 + tq*thetai1(2,it1,iatm)
               !   t2 = t2 + tq*thetai1(3,it1,iatm)
c                  t3 = t3 + tq*thetai1(4,it1,iatm)
               end do
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
              ! tu20 = tu20 + t2*u0
              ! tu11 = tu11 + t1*u1
              ! tu02 = tu02 + t0*u2
c               tu30 = tu30 + t3*u0
c               tu21 = tu21 + t2*u1
c               tu12 = tu12 + t1*u2
c               tu03 = tu03 + t0*u3
            end do
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
           ! tuv200 = tuv200 + tu20*v0
           ! tuv020 = tuv020 + tu02*v0
           ! tuv002 = tuv002 + tu00*v2
           ! tuv110 = tuv110 + tu11*v0
           ! tuv101 = tuv101 + tu10*v1
           ! tuv011 = tuv011 + tu01*v1
c            tuv300 = tuv300 + tu30*v0
c            tuv030 = tuv030 + tu03*v0
c            tuv003 = tuv003 + tu00*v3
c            tuv210 = tuv210 + tu21*v0
c            tuv201 = tuv201 + tu20*v1
c            tuv120 = tuv120 + tu12*v0
c            tuv021 = tuv021 + tu02*v1
c            tuv102 = tuv102 + tu10*v2
c            tuv012 = tuv012 + tu01*v2
c            tuv111 = tuv111 + tu11*v1
         end do
         fphi(1,isite) = tuv000
         fphi(2,isite) = tuv100
         fphi(3,isite) = tuv010
         fphi(4,isite) = tuv001
        ! fphi(5,isite) = tuv200
        ! fphi(6,isite) = tuv020
        ! fphi(7,isite) = tuv002
        ! fphi(8,isite) = tuv110
        ! fphi(9,isite) = tuv101
        ! fphi(10,isite) = tuv011
c         fphi(11,isite) = tuv300
c         fphi(12,isite) = tuv030
c         fphi(13,isite) = tuv003
c         fphi(14,isite) = tuv210
c         fphi(15,isite) = tuv201
c         fphi(16,isite) = tuv120
c         fphi(17,isite) = tuv021
c         fphi(18,isite) = tuv102
c         fphi(19,isite) = tuv012
c         fphi(20,isite) = tuv111
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end

      subroutine fphi_to_cphi_noqpole_nodpole (fphi,cphi)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ftc(10,10)
c      real*8 cphi(10,*)
      real*8 cphi(*)
c      real*8 fphi(20,*)
c      real*8 fphi(10,*)
      real*8 fphi(4,*)      
c
c
c     find the matrix to convert fractional to Cartesian
c
      call frac_to_cart (ftc)
c
c     apply the transformation to get the Cartesian potential
c
      do i = 1, npole
         cphi(i) = ftc(1,1) * fphi(1,i)
c         do j = 2, 4
c            cphi(j,i) = 0.0d0
c            do k = 2, 4
c               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
c            end do
c         end do
c         do j = 5, 10
c            cphi(j,i) = 0.0d0
c            do k = 5, 10
c               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
c            end do
c         end do
      end do
      return
      end

      subroutine bspline_fill_noqpole_nodpole
      use sizes
      use atoms
      use boxes
      use pme
      implicit none
      integer i,ifr
      real*8 xi,yi,zi
      real*8 w,fr,eps
      logical first
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(igrid))  allocate (igrid(3,n))
      end if
c
c     offset used to shift sites off exact lattice bounds
c
      eps = 1.0d-8
c
c     get the B-spline coefficients for each atomic site
c
      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr = dble(nfft1) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(1,i) = ifr - bsorder
         call bsplgen_noqpole_nodpole (w,thetai1(1,1,i))
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(2,i) = ifr - bsorder
         call bsplgen_noqpole_nodpole (w,thetai2(1,1,i))
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr = dble(nfft3) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(3,i) = ifr - bsorder
         call bsplgen_noqpole_nodpole (w,thetai3(1,1,i))
      end do
      return
      end

      subroutine bsplgen_noqpole_nodpole (w,thetai)
      use sizes
      use pme
      use potent
      implicit none
      integer i,j,k
      integer level
      real*8 w,denom
      real*8 thetai(4,*)
      real*8 temp(maxorder,maxorder)
c
c
c     set B-spline depth for partial charges or multipoles
c
      level = 2
c      if (use_mpole .or. use_polar)  level = 4
c
c     initialization to get to 2nd order recursion
c
      temp(2,2) = w
      temp(2,1) = 1.0d0 - w
c
c     perform one pass to get to 3rd order recursion
c
      temp(3,3) = 0.5d0 * w * temp(2,2)
      temp(3,2) = 0.5d0 * ((1.0d0+w)*temp(2,1)+(2.0d0-w)*temp(2,2))
      temp(3,1) = 0.5d0 * (1.0d0-w) * temp(2,1)
c
c     compute standard B-spline recursion to desired order
c
      do i = 4, bsorder
         k = i - 1
         denom = 1.0d0 / dble(k)
         temp(i,i) = denom * w * temp(k,k)
         do j = 1, i-2
            temp(i,i-j) = denom * ((w+dble(j))*temp(k,i-j-1)
     &                             +(dble(i-j)-w)*temp(k,i-j))
         end do
         temp(i,1) = denom * (1.0d0-w) * temp(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      temp(k,bsorder) = temp(k,bsorder-1)
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
c
c     get coefficients for the B-spline second derivative
c

      ! if (level .eq. 4) then
      ! end if
      do i = 1, bsorder
         do j = 1, level
            thetai(j,i) = temp(bsorder-j+1,i)
         end do
      end do
      return
      end

