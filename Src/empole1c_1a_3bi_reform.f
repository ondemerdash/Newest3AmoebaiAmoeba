c
c
      subroutine empole1a_3b_Polar(npole3b,pnum,eptemp,
     &   deptemp,virtemp)
ccc      include 'sizes.i'
c      include 'sizes3b.i'
ccc      include 'atoms.i'
c      include 'atoms3b.i'
ccc      include 'chgpot.i'
c      include 'chgpot3b.i'
c      include 'couple.i'
ccc      include 'group.i'
c      include 'group3b.i'
ccc      include 'mpole.i'
c      include 'mpole3b.i'
c      include 'polar2.i'
c      include 'polgrp.i'
c      include 'polpot.i'
c      include 'usage.i'
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
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 dixdk(3),fdir(3)
      real*8 dixuk(3),dkxui(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
c      real*8 rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
c      real*8 rxqkr(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
c      real*8 gl(0:8),gli(7),glip(7)
      real*8 gli(7),glip(7)
      real*8 sc(10),sci(8),scip(8)
      real*8 gf(7),gfi(6),gti(6)
c      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
ccc      real*8 eptemp,deptemp(3,npole),virtemp(3,3)
      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3
      logical proceed,usei,usek
      character*6 mode

      eptemp = 0.0d0
      eptemp_userep =0.0d0

c   Zero out temporary gradient of polarization energy
      do l1 = 1, npole3b
c        i = pnum(l1)
        do j = 1, 3
c          deptemp(j,i) = 0.0d0
          deptemp(j,l1) = 0.0d0
        end do
      end do

c   Zero out temporary virial
      do i=1,3
         do j=1,3
           virtemp(i,j)=0.0d0
         end do
      end do


      call induce0a_3b_PolelecOnly(npole3b,pnum,uind,uinp)

c      do l1 = 1, npole3b
c        i = pnum(l1)
c         print*,"uindx",uind(1,l1)
c         print*,"uindy",uind(2,l1)
c         print*,"uindz",uind(3,l1)
c      end do 
          

c     perform dynamic allocation of some local arrays
c
c      allocate (mscale(n))
c      allocate (pscale(n))
c      allocate (dscale(n))
c      allocate (uscale(n))

c      allocate (mscale(npole3b))
      allocate (pscale(npole3b))
      allocate (dscale(npole3b))
      allocate (uscale(npole3b))

c
c     set arrays needed to scale connected atom interactions
c
c      if (npole .eq. 0)  return
c      do i = 1, n
      do i = 1, npole3b
c         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      call replica(20.0d0)
c      print*,"Polar off2",off2
c      print*,"Polar xcell ycell zcell",xcell,ycell,zcell
c      print*,"Polar xcell2 ycell2 zcell2",xcell2,ycell2,zcell2
c      print*,"Polar xbox2 ybox2 zbox2",xbox2,ybox2,zbox2
c      print*,"use_replica",use_replica

c      mode = 'MPOLE'
c      call switch (mode)

c
c     set scale factors for permanent multipole and induced terms
c
c      if(use_standard) then
c      if(use_replica.eq..false.) then

      do l1 = 1, npole3b-1
         i=pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
c         print*,"ipole(i)",ipole(i)
c         print*,"zaxis(i)",zaxis(i)
c         print*,"xaxis(i)",xaxis(i)
c         print*,"yaxis(i)",yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
c            mscale(i12(j,ii)) = m2scale
            do kk=1,npole3b
               if (pnum(kk).eq.i12(j,ii)) then
                  pscale(kk)=p2scale
                  goto 31
               end if
            end do 
   31             continue

c            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
c            mscale(i13(j,ii)) = m3scale
            do kk=1,npole3b
               if (pnum(kk).eq.i13(j,ii)) then
                  pscale(kk)=p3scale
                  goto 32
               end if
            end do
   32             continue            
c            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
c            mscale(i14(j,ii)) = m4scale
            do kk=1,npole3b
               if (pnum(kk).eq.i14(j,ii)) then
                  pscale(kk)=p4scale
                  goto 33
               end if
            end do
   33             continue
c            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               do kk=1,npole3b             
                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
     &              (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(i14(j,ii)) = p4scale * p41scale
                   pscale(kk) = p4scale * p41scale
                   goto 34
                 end if
               end do  
   34             continue
            end do
         end do
         do j = 1, n15(ii)
c            mscale(i15(j,ii)) = m5scale
            do kk=1,npole3b
               if (pnum(kk).eq.i15(j,ii)) then
                  pscale(kk)=p5scale
                  goto 35
               end if
            end do
   35             continue
c            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
                  dscale(kk)=d1scale
                  uscale(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
c            dscale(ip11(j,ii)) = d1scale
c            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale(kk)=d2scale
                  uscale(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue
c            dscale(ip12(j,ii)) = d2scale
c            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale(kk)=d3scale
                  uscale(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
c            dscale(ip13(j,ii)) = d3scale
c            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale(kk)=d4scale
                  uscale(kk)=u4scale
                  goto 39
               end if
            end do
   39             continue
c            dscale(ip14(j,ii)) = d4scale
c            uscale(ip14(j,ii)) = u4scale
         end do
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
c            print*,"r2 in stdloop empole1a_3b_Polar",r2
c            print*,"xr yr zr in stdloop empole1a_3b_Polar",xr,yr,zr

c            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
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
c               scale3i = scale3 * uscale(kk)
c               scale5i = scale5 * uscale(kk)
c               scale7i = scale7 * uscale(kk)
c               dsc3 = scale3 * dscale(kk)
c               dsc5 = scale5 * dscale(kk)
c               dsc7 = scale7 * dscale(kk)
c               psc3 = scale3 * pscale(kk)
c               psc5 = scale5 * pscale(kk)
c               psc7 = scale7 * pscale(kk)

               scale3i = scale3 * uscale(l3)
               scale5i = scale5 * uscale(l3)
               scale7i = scale7 * uscale(l3)
               dsc3 = scale3 * dscale(l3)
               dsc5 = scale5 * dscale(l3)
               dsc7 = scale7 * dscale(l3)
               psc3 = scale3 * pscale(l3)
               psc5 = scale5 * pscale(l3)
               psc7 = scale7 * pscale(l3)

c
c     construct necessary auxiliary vectors
c

               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)

               dixuk(1) = di(2)*uind(3,l3) - di(3)*uind(2,l3)
               dixuk(2) = di(3)*uind(1,l3) - di(1)*uind(3,l3)
               dixuk(3) = di(1)*uind(2,l3) - di(2)*uind(1,l3)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l3) - di(3)*uinp(2,l3)
               dixukp(2) = di(3)*uinp(1,l3) - di(1)*uinp(3,l3)
               dixukp(3) = di(1)*uinp(2,l3) - di(2)*uinp(1,l3)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l3) + qi(4)*uind(2,l3)
     &                      + qi(7)*uind(3,l3)
               qiuk(2) = qi(2)*uind(1,l3) + qi(5)*uind(2,l3)
     &                      + qi(8)*uind(3,l3)
               qiuk(3) = qi(3)*uind(1,l3) + qi(6)*uind(2,l3)
     &                      + qi(9)*uind(3,l3)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l3) + qi(4)*uinp(2,l3)
     &                       + qi(7)*uinp(3,l3)
               qiukp(2) = qi(2)*uinp(1,l3) + qi(5)*uinp(2,l3)
     &                       + qi(8)*uinp(3,l3)
               qiukp(3) = qi(3)*uinp(1,l3) + qi(6)*uinp(2,l3)
     &                       + qi(9)*uinp(3,l3)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l3)*qir(3) - uind(3,l3)*qir(2)
               ukxqir(2) = uind(3,l3)*qir(1) - uind(1,l3)*qir(3)
               ukxqir(3) = uind(1,l3)*qir(2) - uind(2,l3)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l3)*qir(3) - uinp(3,l3)*qir(2)
               ukxqirp(2) = uinp(3,l3)*qir(1) - uinp(1,l3)*qir(3)
               ukxqirp(3) = uinp(1,l3)*qir(2) - uinp(2,l3)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l3)
     &                     + di(2)*uind(2,l3) + di(3)*uind(3,l3)
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = qir(1)*uind(1,l3) + qir(2)*uind(2,l3)
     &                     + qir(3)*uind(3,l3)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l3)
     &                      + di(2)*uinp(2,l3) + di(3)*uinp(3,l3)
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qir(1)*uinp(1,l3) + qir(2)*uinp(2,l3)
     &                      + qir(3)*uinp(3,l3)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for permanent components
c

c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
c               gl(4) = sc(5)*sc(6)
c               gl(5) = -4.0d0 * sc(9)
c               gl(6) = sc(2)
c               gl(7) = 2.0d0 * (sc(7)-sc(8))
c               gl(8) = 2.0d0 * sc(10)

c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = f * ei
               eptemp = eptemp + ei
c             print*,"eptemp in empole3a=",eptemp
c
c     increment the total intermolecular energy
c

c               if (molcule(ii) .ne. molcule(kk)) then
c                  einter = einter + e + ei
c               end if

c
c     intermediate variables for the permanent components
c

c               gf(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
c     &                    + rr7*(gl(2)+gl(7)+gl(8))
c     &                    + rr9*(gl(3)+gl(5)) + rr11*gl(4)
c               gf(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
c               gf(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
c               gf(4) = 2.0d0 * rr5
c               gf(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
c               gf(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
c               gf(7) = 4.0d0 * rr7

c
c     intermediate variables for the induced components
c
               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
     &                                   + (glip(1)+glip(6))*dsc3
     &                                   + scip(2)*scale3i)
     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0 * rr9 * (gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the permanent force components
c

c               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
c     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
c     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
c               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
c     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
c     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
c               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
c     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
c     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))

c
c     get the induced force components
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     account for partially excluded induced interactions
c

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))

               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
     &                                  +(glip(1)+glip(6))*dscale(l3))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
     &                                  +(glip(2)+glip(7))*dscale(l3))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
     &                                  +glip(3)*dscale(l3))



               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c

c               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(kk)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               temp3 = 0.5d0 * rr3 * uscale(l3) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(l3)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify induced force for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c

c               if (poltyp .eq. 'DIRECT') then
c                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
c     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
c                  temp5 = 0.5d0 * rr5 * scale5i
c                  fdir(1) = gfd*xr + temp5
c     &                         * (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                           +sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
c                  fdir(2) = gfd*yr + temp5
c     &                         * (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                           +sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
c                  fdir(3) = gfd*zr + temp5
c     &                         * (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                           +sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
c                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
c                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
c                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
c               end if


c
c     intermediate terms for induced torque on multipoles
c
               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the permanent torque components
c

c               ttm2(1) = -rr3*dixdk(1) + gf(2)*dixr(1) - gf(5)*rxqir(1)
c     &           + gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
c     &           - gf(7)*(rxqikr(1)+qkrxqir(1))
c               ttm2(2) = -rr3*dixdk(2) + gf(2)*dixr(2) - gf(5)*rxqir(2)
c     &           + gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
c     &           - gf(7)*(rxqikr(2)+qkrxqir(2))
c               ttm2(3) = -rr3*dixdk(3) + gf(2)*dixr(3) - gf(5)*rxqir(3)
c     &           + gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
c     &           - gf(7)*(rxqikr(3)+qkrxqir(3))
c               ttm3(1) = rr3*dixdk(1) + gf(3)*dkxr(1) - gf(6)*rxqkr(1)
c     &           - gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
c     &           - gf(7)*(rxqkir(1)-qkrxqir(1))
c               ttm3(2) = rr3*dixdk(2) + gf(3)*dkxr(2) - gf(6)*rxqkr(2)
c     &           - gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
c     &           - gf(7)*(rxqkir(2)-qkrxqir(2))
c               ttm3(3) = rr3*dixdk(3) + gf(3)*dkxr(3) - gf(6)*rxqkr(3)
c     &           - gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
c     &           - gf(7)*(rxqkir(3)-qkrxqir(3))

c
c     get the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
c                  ftm2(j) = f * ftm2(j) * mscale(kk)
                  ftm2i(j) = f * ftm2i(j)
c                  ttm2(j) = f * ttm2(j) * mscale(kk)
                  ttm2i(j) = f * ttm2i(j)
c                  ttm3(j) = f * ttm3(j) * mscale(kk)
                  ttm3i(j) = f * ttm3i(j)
               end do
c
c     increment gradient due to force and torque on first site
c

c               dem(1,ii) = dem(1,ii) + ftm2(1)
c               dem(2,ii) = dem(2,ii) + ftm2(2)
c               dem(3,ii) = dem(3,ii) + ftm2(3)
c               dep(1,ii) = dep(1,ii) + ftm2i(1)
c               dep(2,ii) = dep(2,ii) + ftm2i(2)
c               dep(3,ii) = dep(3,ii) + ftm2i(3)
c               call torque (i,ttm2,ttm2i,frcxi,frcyi,frczi)

c               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
c               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
c               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
c               call torque_3b (deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               call torque_3b_new(npole3b,pnum,deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)
c
c     increment gradient due to force and torque on second site
c

c               dem(1,kk) = dem(1,kk) - ftm2(1)
c               dem(2,kk) = dem(2,kk) - ftm2(2)
c               dem(3,kk) = dem(3,kk) - ftm2(3)
c               dep(1,kk) = dep(1,kk) - ftm2i(1)
c               dep(2,kk) = dep(2,kk) - ftm2i(2)
c               dep(3,kk) = dep(3,kk) - ftm2i(3)
c               call torque (k,ttm3,ttm3i,frcxk,frcyk,frczk)

c               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
c               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
c               deptemp(3,k) = deptemp(3,k) - ftm2i(3)

               deptemp(1,l3) = deptemp(1,l3) - ftm2i(1)
               deptemp(2,l3) = deptemp(2,l3) - ftm2i(2)
               deptemp(3,l3) = deptemp(3,l3) - ftm2i(3)
c               call torque_3b (deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
               call torque_3b_new(npole3b,pnum,deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)

c
c     increment the internal virial tensor components
c
               iaz = iz
               iax = ix
               iay = iy
               kaz = kz
               kax = kx
               kay = ky
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)

c               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
c     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
c     &                  + xky*frcyk(1) + xkz*frczk(1)
c               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
c     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
c     &                  + yky*frcyk(1) + ykz*frczk(1)
c               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
c     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
c     &                  + zky*frcyk(1) + zkz*frczk(1)
c               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
c     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
c     &                  + yky*frcyk(2) + ykz*frczk(2)
c               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
c     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
c     &                  + zky*frcyk(2) + zkz*frczk(2)
c               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
c     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
c     &                  + zky*frcyk(3) + zkz*frczk(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
c            end if
               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

   10       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c

c         do j = 1, n12(ii)
c            mscale(i12(j,ii)) = 1.0d0
c            pscale(i12(j,ii)) = 1.0d0
c         end do
c         do j = 1, n13(ii)
c            mscale(i13(j,ii)) = 1.0d0
c            pscale(i13(j,ii)) = 1.0d0
c         end do
c         do j = 1, n14(ii)
c            mscale(i14(j,ii)) = 1.0d0
c            pscale(i14(j,ii)) = 1.0d0
c         end do
c         do j = 1, n15(ii)
c            mscale(i15(j,ii)) = 1.0d0
c            pscale(i15(j,ii)) = 1.0d0
c         end do
c         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = 1.0d0
c            uscale(ip11(j,ii)) = 1.0d0
c         end do
c         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = 1.0d0
c            uscale(ip12(j,ii)) = 1.0d0
c         end do
c         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = 1.0d0
c            uscale(ip13(j,ii)) = 1.0d0
c         end do
c         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = 1.0d0
c            uscale(ip14(j,ii)) = 1.0d0
c         end do

         do j = 1, npole3b
c            mscale(i12(j,ii)) = 1.0d0
            pscale(j) = 1.0d0
            dscale(j) = 1.0d0
            uscale(j) = 1.0d0
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

c      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end


      subroutine induce0a_3b_PolelecOnly(npole3b,pnum,uind,uinp)
c      include 'sizes3b.i'
c      include 'atoms3b.i'
c      include 'polar2.i'
      use sizes
      use atoms
      use polar, only: polarity, thole, pdamp
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 field(3,npole3b),off3b
      real*8 fieldp(3,npole3b)
      character*6 mode
      integer l1,l2,l3,k1,k2,i1,i2
      real*8 M_tot(3*npole3b,3*npole3b)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      integer npole3b,pnum(*)


c      if (.not. use_polar)  return
        do l1=1,npole3b
           do j=1,3
              field(j,l1)=0.0d0
              fieldp(j,l1)=0.0d0
           end do
        end do

         do l1 = 1, 3*npole3b
            do l3 = 1, 3*npole3b
               M_tot(l1,l3) = 0.0d0
            end do
         end do

       call field_noewald_umutual_rl_3b (field,fieldp,M_tot,
     &   npole3b,pnum)

         do l1 = 1, npole3b
c            i = pnum(l1)
            do j = 1, 3
               uind(j,l1) = 0.0d0
               uinp(j,l1) = 0.0d0
            end do
         end do

c         call umutual_rl (M_tot)

         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               i1 = 3*(l1-1)+j
c               M_tot(i1,i1) = 1.0d0/polarity(i)
               M_tot(i1,i1) = M_tot(i1,i1)+1.0d0/polarity(i)
            end do
         end do
c         do l1 = 1, 3*npole3b
c            do l3 = 1, 3*npole3b
c               print*,"M_tot",M_tot(l1,l3)
c            end do
c         end do
c         do l1=1,npole3b
c           do j=1,3
c              print*,"Direct field l1 j",l1,j,field(j,l1)
c           end do
c        end do

         call invert(3*npole3b,M_tot)

         do l1 = 1, npole3b
            i = pnum(l1)
            do i1 = 1, 3
               i2 = 3*(l1-1) + i1
               do l3 = 1, npole3b
                  k = pnum(l3)
                  k2 = 3*(l3-1)
c                  uind(i1,i) = uind(i1,i) + M_tot(i2,k2+1)*fdir(1,k) +
c     &                                      M_tot(i2,k2+2)*fdir(2,k) +
c     &                                      M_tot(i2,k2+3)*fdir(3,k)
                  uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*field(1,l3)+
     &                                      M_tot(i2,k2+2)*field(2,l3)+
     &                                      M_tot(i2,k2+3)*field(3,l3)
                  uinp(i1,l1)=uinp(i1,l1)+ M_tot(i2,k2+1)*fieldp(1,l3)+
     &                                     M_tot(i2,k2+2)*fieldp(2,l3)+
     &                                     M_tot(i2,k2+3)*fieldp(3,l3)



               end do
c               print*,i,i1,uind(i1,i)
            end do
         end do

ccc      do l1 = 1, npole3b
ccc         i = pnum(l1)
ccc         do j = 1, 3
ccc            uinp(j,l1)=uind(j,l1)
ccc         end do
ccc      end do


      return
      end
c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
      subroutine field_noewald_umutual_rl_3b(field,fieldp,M_tot,
     &   npole3b,pnum)
c      implicit none
c      include 'sizes3b.i'
c      include 'atoms3b.i'
c      include 'mpole3b.i'
c      include 'couple.i'
c      include 'group3b.i'
c      include 'polar2.i'
c      include 'polgrp.i'
c      include 'polpot.i'
      use sizes
      use atoms
      use mpole
      use couple
      use group
      use polgrp
      use polpot
      use polar, only: pdamp, thole
      implicit none
      real*8 field(3,*),off3b
      real*8 fieldp(3,*)
      real*8, allocatable :: dscale_dir(:)
      real*8, allocatable :: dscale_mut(:)
      real*8, allocatable :: pscale(:)
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3_dir,rr5_dir,rr7_dir
      real*8 rr3,rr5
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3_dir,scale5_dir
      real*8 scale7_dir
      real*8 scale3_mut,scale5_mut
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3),M_tot(3*npole3b,3*npole3b)
      real*8 fip(3),fkp(3)
      logical proceed
      character*6 mode
      integer i,ii,j,l1,k,m,kk
      integer l2,l3,k1,k2,i1,i2
      integer pnum(*),npole3b
      real*8 Txx,Txy,Txz,Tyx
      real*8 Tyy,Tyz,Tzx,Tzy,Tzz
      

c      allocate (dscale_dir(n))
c      allocate (dscale_mut(n))
c      allocate (pscale(n))

      allocate (dscale_dir(npole3b))
      allocate (dscale_mut(npole3b))
      allocate (pscale(npole3b))

c
c     set the switching function coefficients
c

c      mode = 'MPOLE'
c      call switch (mode)

c      mode = 'MPOLE3'
c      call switch (mode,off3b)
c       call replica (off3b)
c       use_replica = .false.

c
c     compute the direct induced dipole moment at each atom
c
c      if(use_replica.eq..false.) then
      do l1 = 1, npole3b-1
         i = pnum(l1)
         i2 = 3*(l1-1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
c         do j = i+1, npole
c            dscale_dir(ipole(j)) = 1.0d0
c            pscale(ipole(j)) = 1.0d0
c            dscale_mut(ipole(j)) = 1.0d0
c         end do
         do j = 1,npole3b
            dscale_dir(j) = 1.0d0
            pscale(j) = 1.0d0
            dscale_mut(j) = 1.0d0
         end do
         do j = 1, n12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i12(j,ii)) then
                  pscale(kk)=p2scale
                  goto 31
               end if
            end do
   31             continue
c            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i13(j,ii)) then
                  pscale(kk)=p3scale
                  goto 32
               end if
            end do
   32             continue
c            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
c            pscale(i14(j,ii)) = p4scale
c            do k = 1, np11(ii)
c               if (i14(j,ii) .eq. ip11(k,ii))
c     &            pscale(i14(j,ii)) = p4scale * p41scale
c            end do
            do kk=1,npole3b
               if (pnum(kk).eq.i14(j,ii)) then
                  pscale(kk)=p4scale
                  goto 33
               end if
            end do
   33             continue
c            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               do kk=1,npole3b
                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
     &              (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(i14(j,ii)) = p4scale * p41scale
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
c            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
                  dscale_dir(kk)=d1scale
                  dscale_mut(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
c            dscale_dir(ip11(j,ii)) = d1scale
c            dscale_mut(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale_dir(kk)=d2scale
                  dscale_mut(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue

c            dscale_dir(ip12(j,ii)) = d2scale
c            dscale_mut(ip12(j,ii)) = u2scale

         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale_dir(kk)=d3scale
                  dscale_mut(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
c            dscale_dir(ip13(j,ii)) = d3scale
c            dscale_mut(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale_dir(kk)=d4scale
                  dscale_mut(kk)=u4scale
                  goto 39
               end if
            end do
   39             continue
c            dscale_dir(ip14(j,ii)) = d4scale
c            dscale_mut(ip14(j,ii)) = u4scale
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
c               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3_dir = 1.0d0
                  scale5_dir = 1.0d0
                  scale7_dir = 1.0d0

c                  scale3_mut = dscale_mut(kk)
c                  scale5_mut = dscale_mut(kk)

                  scale3_mut = dscale_mut(l3)
                  scale5_mut = dscale_mut(l3)

                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3_dir = 1.0d0 - expdamp
                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7_dir = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut = scale3_mut * (1.0d0-expdamp)
                        scale5_mut = scale5_mut * (1.0d0-expdamp
     &                                        *(1.0d0-damp))

                     end if
                  end if
                  rr3_dir = scale3_dir / (r*r2)
                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

                  rr3 = scale3_mut / (r*r2)
                  rr5 = 3.0d0 * scale5_mut / (r*r2*r2)


                  Txx = -(-rr3 + xr*xr*rr5)
                  Txy = -(xr*yr*rr5)
                  Txz = -(xr*zr*rr5)
                  Tyx = Txy
                  Tyy = -(-rr3 + yr*yr*rr5)
                  Tyz = -(yr*zr*rr5)
                  Tzx = Txz
                  Tzy = Tyz
                  Tzz = -(-rr3 + zr*zr*rr5)

                  M_tot(i2+1,k2+1) = Txx
                  M_tot(i2+1,k2+2) = Txy
                  M_tot(i2+1,k2+3) = Txz
                  M_tot(i2+2,k2+1) = Tyx
                  M_tot(i2+2,k2+2) = Tyy
                  M_tot(i2+2,k2+3) = Tyz
                  M_tot(i2+3,k2+1) = Tzx
                  M_tot(i2+3,k2+2) = Tzy
                  M_tot(i2+3,k2+3) = Tzz

                  M_tot(k2+1,i2+1) = Txx
                  M_tot(k2+1,i2+2) = Txy
                  M_tot(k2+1,i2+3) = Txz
                  M_tot(k2+2,i2+1) = Tyx
                  M_tot(k2+2,i2+2) = Tyy
                  M_tot(k2+2,i2+3) = Tyz
                  M_tot(k2+3,i2+1) = Tzx
                  M_tot(k2+3,i2+2) = Tzy
                  M_tot(k2+3,i2+3) = Tzz

                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
                  do j = 1, 3
c                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
c                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
c                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
c                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)

c                     field(j,l1) = field(j,l1) + fid(j)*dscale_dir(kk)
c                     field(j,l3) = field(j,l3) + fkd(j)*dscale_dir(kk)
c                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(kk)
c                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(kk)

                     field(j,l1) = field(j,l1) + fid(j)*dscale_dir(l3)
                     field(j,l3) = field(j,l3) + fkd(j)*dscale_dir(l3)
                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(l3)
                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(l3)


                  end do
c               end if
            end if
         end do
      end do
c      else
c
c     periodic boundary for large cutoffs via replicates method
c
c       use_replica=.true.
c   NOTE TO SELF:  REMOVED if (use_replica) BECAUSE NCELL WAS
c   ALWAYS EQUAL TO ZERO AFTER PERIODIC BOUNDARY CONDITIONS WERE
c   IMPLEMENTED CORRECTLY. 

      deallocate (dscale_dir)
      deallocate (dscale_mut)
      deallocate (pscale)
      return
      end


      subroutine just_chk_rotpole
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
      return
      end

c     ##  COPYRIGHT (C) 2007 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque_3b  --  convert single site torque to force  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque" takes the torque values on a single site defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
c               call torque_3b (pnum2,deptemp,i,
c     &           ttm2,ttm2i,frcxi,frcyi,frczi)

c      subroutine torque_3b(pnum2,deptemp,
c     &    i,trq1,trq2,frcx,frcy,frcz)
      subroutine torque_3b_new(npole3b,pnum,deptemp,
     &    i,trq1,trq2,frcx,frcy,frcz)
c      implicit none
c      include 'sizes3b.i'
c      include 'atoms3b.i'
c      include 'mpole3b.i'
c      include 'mpole3b_2.i'
      use sizes
      use atoms
      use mpole
      implicit none
      integer i,j,l1
      integer ia,ib,ic,id
      real*8 du,dv,dw,random
      real*8 usiz,vsiz,wsiz
      real*8 rsiz,ssiz
      real*8 t1siz,t2siz
      real*8 uvsiz,uwsiz,vwsiz
      real*8 ursiz,ussiz
      real*8 vssiz,wssiz
      real*8 uvcos,uwcos,vwcos
      real*8 urcos,uscos
      real*8 vscos,wscos
      real*8 ut1cos,ut2cos
      real*8 uvsin,uwsin,vwsin
      real*8 ursin,ussin
      real*8 vssin,wssin
      real*8 ut1sin,ut2sin
      real*8 dphidu,dphidv,dphidw
      real*8 dphidr,dphids
      real*8 trq1(3),trq2(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
      character*8 axetyp
      real*8 deptemp(3,*)
c      integer pnum2(*)
      integer ia_l1,ib_l1,ic_l1
      integer id_l1,pnum(*),npole3b
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0d0
         frcx(j) = 0.0d0
         frcy(j) = 0.0d0
      end do
c
c     get the local frame type and the frame-defining atoms
c
c      i=pnum(l1)
      
      ia = zaxis(i)
      ib = ipole(i)
      ic = xaxis(i)
      id = yaxis(i)

         do l1=1,npole3b
           if(pnum(l1).eq.zaxis(i)) then
             ia_l1=l1
             goto 31      
           end if
         end do
   31             continue

         do l1=1,npole3b
           if(pnum(l1).eq.ipole(i)) then
             ib_l1=l1
             goto 32
           end if
         end do 
   32             continue

         do l1=1,npole3b
           if(pnum(l1).eq.xaxis(i)) then
             ic_l1=l1
             goto 33
           end if
         end do 
   33             continue

         do l1=1,npole3b
           if(pnum(l1).eq.yaxis(i)) then
             id_l1=l1
             goto 34
           end if
         end do
   34             continue

         
c      pnum2(ia)=ia_l1
c      pnum2(ib)=ib_l1
c      pnum2(ic)=ic_l1
c      pnum2(id)=id_l1
c      if(i.eq.1) then
c        print*,"zaxis(i)",zaxis(i)
c        print*,"ipole(i)",ipole(i)
c        print*,"xaxis(i)",xaxis(i)
c        print*,"yaxis(i)",yaxis(i)
c      end if
      axetyp = polaxe(i)
      !print*,"i axetyp",i,axetyp
      if (axetyp .eq. 'None')  return
c
c     construct the three rotation axes for the local frame
c
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
      else
         v(1) = random ()
         v(2) = random ()
         v(3) = random ()
      end if
      if (axetyp.eq.'Z-Bisect' .or. axetyp.eq.'3-Fold') then
         w(1) = x(id) - x(ib)
         w(2) = y(id) - y(ib)
         w(3) = z(id) - z(ib)
      else
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
      end if
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      if (axetyp .eq. 'Z-Bisect') then
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         ssiz = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
            s(j) = s(j) / ssiz
         end do
      end if
c
c     find the perpendicular and angle for each pair of axes
c
      uv(1) = v(2)*u(3) - v(3)*u(2)
      uv(2) = v(3)*u(1) - v(1)*u(3)
      uv(3) = v(1)*u(2) - v(2)*u(1)
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
      vwsiz = sqrt(vw(1)*vw(1) + vw(2)*vw(2) + vw(3)*vw(3))
      do j = 1, 3
         uv(j) = uv(j) / uvsiz
         uw(j) = uw(j) / uwsiz
         vw(j) = vw(j) / vwsiz
      end do
      if (axetyp .eq. 'Z-Bisect') then
         ur(1) = r(2)*u(3) - r(3)*u(2)
         ur(2) = r(3)*u(1) - r(1)*u(3)
         ur(3) = r(1)*u(2) - r(2)*u(1)
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
         wssiz = sqrt(ws(1)*ws(1) + ws(2)*ws(2) + ws(3)*ws(3))
         do j = 1, 3
            ur(j) = ur(j) / ursiz
            us(j) = us(j) / ussiz
            vs(j) = vs(j) / vssiz
            ws(j) = ws(j) / wssiz
         end do
      end if
c
c     get sine and cosine of angles between the rotation axes
c
      uvcos = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      uvsin = sqrt(1.0d0 - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0d0 - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0d0 - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0d0 - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0d0 - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0d0 - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0d0 - wscos*wscos)
      end if
c
c     compute the projection of v and w onto the ru-plane
c
      if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
           t1(j) = v(j) - s(j)*vscos
           t2(j) = w(j) - s(j)*wscos
         end do
         t1siz = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         t2siz = sqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
         do j = 1, 3
           t1(j) = t1(j) / t1siz
           t2(j) = t2(j) / t2siz
         end do
         ut1cos = u(1)*t1(1) + u(2)*t1(2) + u(3)*t1(3)
         ut1sin = sqrt(1.0d0 - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0d0 - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq2(1)*u(1) - trq2(2)*u(2) - trq2(3)*u(3)
      dphidv = -trq2(1)*v(1) - trq2(2)*v(2) - trq2(3)*v(3)
      dphidw = -trq2(1)*w(1) - trq2(2)*w(2) - trq2(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq2(1)*r(1) - trq2(2)*r(2) - trq2(3)*r(3)
         dphids = -trq2(1)*s(1) - trq2(2)*s(2) - trq2(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
c            deptemp(j,ia) = deptemp(j,ia) + du
c            deptemp(j,ib) = deptemp(j,ib) - du
            deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
            deptemp(j,ib_l1) = deptemp(j,ib_l1) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
c            deptemp(j,ia) = deptemp(j,ia) + du
c            deptemp(j,ic) = deptemp(j,ic) + dv
c            deptemp(j,ib) = deptemp(j,ib) - du - dv
            deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
            deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
            deptemp(j,ib_l1) = deptemp(j,ib_l1) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
c            deptemp(j,ia) = deptemp(j,ia) + du
c            deptemp(j,ic) = deptemp(j,ic) + dv
c            deptemp(j,ib) = deptemp(j,ib) - du - dv
            deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
            deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
            deptemp(j,ib_l1) = deptemp(j,ib_l1) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Z-Bisect local coordinate method
c
      else if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
            dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &              / (vsiz*(ut1sin+ut2sin))
            dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &              / (wsiz*(ut1sin+ut2sin))
c            deptemp(j,ia) = deptemp(j,ia) + du
c            deptemp(j,ic) = deptemp(j,ic) + dv
c            deptemp(j,id) = deptemp(j,id) + dw
c            deptemp(j,ib) = deptemp(j,ib) - du - dv - dw
            deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
            deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
            deptemp(j,id_l1) = deptemp(j,id_l1) + dw
            deptemp(j,ib_l1) = deptemp(j,ib_l1) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
      else if (axetyp .eq. '3-Fold') then
         do j = 1, 3
            du = uw(j)*dphidw/(usiz*uwsin)
     &              + uv(j)*dphidv/(usiz*uvsin)
     &              - uw(j)*dphidu/(usiz*uwsin)
     &              - uv(j)*dphidu/(usiz*uvsin)
            dv = vw(j)*dphidw/(vsiz*vwsin)
     &              - uv(j)*dphidu/(vsiz*uvsin)
     &              - vw(j)*dphidv/(vsiz*vwsin)
     &              + uv(j)*dphidv/(vsiz*uvsin)
            dw = -uw(j)*dphidu/(wsiz*uwsin)
     &              - vw(j)*dphidv/(wsiz*vwsin)
     &              + uw(j)*dphidw/(wsiz*uwsin)
     &              + vw(j)*dphidw/(wsiz*vwsin)
            du = du / 3.0d0
            dv = dv / 3.0d0
            dw = dw / 3.0d0
c            deptemp(j,ia) = deptemp(j,ia) + du
c            deptemp(j,ic) = deptemp(j,ic) + dv
c            deptemp(j,id) = deptemp(j,id) + dw
c            deptemp(j,ib) = deptemp(j,ib) - du - dv - dw
            deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
            deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
            deptemp(j,id_l1) = deptemp(j,id_l1) + dw
            deptemp(j,ib_l1) = deptemp(j,ib_l1) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
      end if
      return
      end

