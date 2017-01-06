c
c
      subroutine empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,eptemp)
      use bound
      use sizes
      use atoms
      use chgpot
      use couple
      use group
      use mpole
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use usage
      use cell
      use boxes
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,fgrp,gfd
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale3i
      real*8 scale5,scale5i
      real*8 scale7,scale7i
      real*8 temp3,temp5,temp7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz,qiyx,qiyy,qiyz,qizx,qizy,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz,qkyx,qkyy,qkyz,qkzx,qkzy,qkzz
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 fdir(3)
      real*8 dixuk(3),dkxui(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qirx,qiry,qirz,qkrx,qkry,qkrz
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqkxz,rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 gli(7),glip(7)
      real*8 sc(10),sci(8),scip(8)
      real*8 gf(7),gfi(6),gti(6)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3
      logical proceed,usei,usek
      character*6 mode

      eptemp = 0.0d0
      eptemp_userep =0.0d0
c
c   Zero out temporary gradient of polarization energy
c
      !mode = 'MPOLE'
      !call switch (mode)
      off2=1.0d12
      !print*,"mpole off2",off2

      call induce0a_3b_PolelecOnly_orig_nopriordir(npole3b,pnum,
     & uind,uinp)

c      do l1=1,npole3b
c         i=pnum(l1)
c         print*,"npole3b i uindx",npole3b,i,uind(1,l1)
c         print*,"npole3b i uindy",npole3b,i,uind(2,l1)
c         print*,"npole3b i uindz",npole3b,i,uind(3,l1)
c      end do

c      call induce0a_3b_PolelecOnly_subtrdir(npole3b,pnum,uind,uinp)

c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(npole3b))
c      allocate (dscale(npole3b))
c      allocate (uscale(npole3b))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, npole3b
         pscale(i) = 1.0d0
c         dscale(i) = 1.0d0
c         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c
c     set scale factors for permanent multipole and induced terms
c

      do l1 = 1, npole3b-1
         i=pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyx = rpole(8,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizx = rpole(11,i)
         qizy = rpole(12,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i12(j,ii)) then
                  pscale(kk)=p2scale
                  goto 31
               end if
            end do 
   31             continue

         end do
         do j = 1, n13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i13(j,ii)) then
                  pscale(kk)=p3scale
                  goto 32
               end if
            end do
   32             continue            
         end do
         do j = 1, n14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i14(j,ii)) then
                  pscale(kk)=p4scale
                  goto 33
               end if
            end do
   33             continue
            do k = 1, np11(ii)
               do kk=1,npole3b             
                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
     &              (pnum(kk).eq.ip11(k,ii)) ) then
                   pscale(kk) = p4scale * p41scale
                   goto 34
                 end if
               end do  
   34             continue
            end do
         end do
         do j = 1, n15(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i15(j,ii)) then
                  pscale(kk)=p5scale
                  goto 35
               end if
            end do
   35             continue
         end do
c         do j = 1, np11(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip11(j,ii)) then
c                  dscale(kk)=d1scale
c                  uscale(kk)=u1scale
c                  goto 36
c               end if
c            end do
c   36             continue
c         end do
c         do j = 1, np12(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip12(j,ii)) then
c                  dscale(kk)=d2scale
c                  uscale(kk)=u2scale
c                  goto 37
c               end if
c            end do
c   37             continue
c         end do
c         do j = 1, np13(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip13(j,ii)) then
c                  dscale(kk)=d3scale
c                  uscale(kk)=u3scale
c                  goto 38
c               end if
c            end do
c   38             continue
c         end do
c         do j = 1, np14(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip14(j,ii)) then
c                  dscale(kk)=d4scale
c                  uscale(kk)=u4scale
c                  goto 39
c               end if
c            end do
c   39             continue
c         end do
         do l3 = l1+1, npole3b
            k=pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr

            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyx = rpole(8,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzx = rpole(11,k)
               qkzy = rpole(12,k)
               qkzz = rpole(13,k)
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
c               scale3i = scale3 * uscale(l3)
c               scale5i = scale5 * uscale(l3)
c               scale7i = scale7 * uscale(l3)
c               dsc3 = scale3 * dscale(l3)
c               dsc5 = scale5 * dscale(l3)
c               dsc7 = scale7 * dscale(l3)
               psc3 = scale3 * pscale(l3)
               psc5 = scale5 * pscale(l3)
               psc7 = scale7 * pscale(l3)
c
c     construct necessary auxiliary vectors
c
c
c     calculate scalar products for permanent components
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qirx*xr + qiry*yr + qirz*zr
               sc(6) = qkrx*xr + qkry*yr + qkrz*zr
               sc(7) = qirx*dkx + qiry*dky + qirz*dkz
               sc(8) = qkrx*dix + qkry*diy + qkrz*diz
               sc(9) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(10) = qixx*qkxx + qixy*qkxy + qixz*qkxz
     &                     + qiyx*qkyx + qiyy*qkyy + qiyz*qkyz
     &                     + qizx*qkzx + qizy*qkzy + qizz*qkzz
c
c     calculate scalar products for induced components
c
               sci(1) = uind(1,l1)*dkx + uind(2,l1)*dky
     &                     + uind(3,l1)*dkz + dix*uind(1,l3)
     &                     + diy*uind(2,l3) + diz*uind(3,l3)
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = qirx*uind(1,l3) + qiry*uind(2,l3)
     &                     + qirz*uind(3,l3)
               sci(8) = qkrx*uind(1,l1) + qkry*uind(2,l1)
     &                     + qkrz*uind(3,l1)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = f * ei
               eptemp = eptemp + ei
c
c     intermediate variables for the induced components
c
c
            end if

   10       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, npole3b
            pscale(j) = 1.0d0
c            dscale(j) = 1.0d0
c            uscale(j) = 1.0d0
         end do

      end do
c
c     perform deallocation of some local arrays
c

c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c       use_replica=.true.
c   NOTE TO SELF:  REMOVED if (use_replica) BECAUSE NCELL WAS
c   ALWAYS EQUAL TO ZERO AFTER PERIODIC BOUNDARY CONDITIONS WERE
c   IMPLEMENTED CORRECTLY.  ALSO, AS AN ADDITIONAL CHECK, THE ENERGY
c   IN THIS SECTION WAS EQUAL TO ZERO.
c
c     perform deallocation of some local arrays
c
c      print*,"eptemp",eptemp
      deallocate (pscale)
c      deallocate (dscale)
c      deallocate (uscale)
      return
      end


