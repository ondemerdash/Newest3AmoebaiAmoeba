
c      subroutine ereal1c_3b(aewald3b3b,ewaldcut3b,npole3b,pnum,uind,uinp,
c     &   eptemp,deptemp2,virtemp)
      subroutine ereal1c1_3b_totfieldnpole(npole3b,pnum,uind,uinp,
     &   deptemp2,virtemp)
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use limits
      use ewald
      use math
      use mpole
      use polgrp
      use polpot
      use polar, only: polarity, thole, pdamp      
      use shunt
      implicit none
      integer moli1,moli2,moli3
      integer i,j,k,l1,l3,l2,l,ll
      integer ii,kk,jcell,jj
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
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
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 uind(3,*),eptemp,deptemp2(3,*)
      real*8 virtemp(3,3),uinp(3,*)
      real*8 ewaldcut3b_2
      integer npole3b,pnum(*)
      character*7 mode
      logical flag
      external erfc

      eintra = 0.0d0
c      if (npole3b .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c

c      allocate (mscale(n))
c      allocate (pscale(n))
c      allocate (dscale(n))
c      allocate (uscale(n))

c      allocate (pscale(npole3b))
c      allocate (dscale(npole3b))
      allocate (uscale(npole3b))
c
c     set arrays needed to scale connected atom interactions
c
c      do i = 1, n
c         mscale(i) = 1.0d0
      do i = 1, npole3b
c         pscale(i) = 1.0d0
c         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do

c      do l1=1,npole3b
c         deptemp2(1,l1)=0.0d0
c         deptemp2(2,l1)=0.0d0
c         deptemp2(3,l1)=0.0d0
c      end do
c      do i=1,3
c        do j=1,3
c          virtemp(j,i) = 0.0d0
c        end do
c      end do

c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      mode = 'EWALD'
c      call switch (mode)
c       print*,"ereal1c1_3b_totfieldnpole off2=",off2

c      ewaldcut3b_2=ewaldcut3b*ewaldcut3b
      ewaldcut3b_2= off2     
      aewald3b=aewald
      do l1 = 1, npole3b-1
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
c         ci = rpole(1,i)
c         di(1) = rpole(2,i)
c         di(2) = rpole(3,i)
c         di(3) = rpole(4,i)
c         qi(1) = rpole(5,i)
c         qi(2) = rpole(6,i)
c         qi(3) = rpole(7,i)
c         qi(4) = rpole(8,i)
c         qi(5) = rpole(9,i)
c         qi(6) = rpole(10,i)
c         qi(7) = rpole(11,i)
c         qi(8) = rpole(12,i)
c         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c

c         do j = 1, n12(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i12(j,ii)) then
c                 pscale(kk) = p2scale
c                 goto 31
c               end if
c            end do
c   31   continue 
c         end do
c
c         do j = 1, n13(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i13(j,ii)) then
c                 pscale(kk) = p3scale
c                 goto 32
c               end if
c            end do
c   32   continue
c         end do
c
c         do j = 1, n14(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i14(j,ii)) then
c                 pscale(kk) = p4scale
c                 goto 33
c               end if
c            end do
c   33   continue
c            do k = 1, np11(ii)
c               do kk=1,npole3b 
c                  if ((i14(j,ii) .eq. ip11(k,ii)).and.
c     &               (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(kk) = p4scale * p41scale
c                   goto 34
c                  end if 
c               end do
c   34   continue
c            end do
c         end do
c         do j = 1, n15(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i15(j,ii)) then
c                 pscale(kk) = p5scale
c                 goto 35
c               end if
c            end do      
c   35   continue      
c         end do

         do j = 1, np11(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip11(j,ii)) then
                 uscale(kk) = u1scale                          
                 goto 36
               end if
            end do
   36   continue
         end do

         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = d2scale
c            uscale(ip12(j,ii)) = u2scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip12(j,ii)) then
c                 dscale(kk) = d2scale
                 uscale(kk) = u2scale
                 goto 37
               end if
            end do
   37   continue
         end do
         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = d3scale
c            uscale(ip13(j,ii)) = u3scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip13(j,ii)) then
c                 dscale(kk) = d3scale
                 uscale(kk) = u3scale
                 goto 38
               end if
            end do
   38   continue
         end do
         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = d4scale
c            uscale(ip14(j,ii)) = u4scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip14(j,ii)) then
c                 dscale(kk) = d4scale
                 uscale(kk) = u4scale
                 goto 39
               end if
            end do
   39   continue             
         end do
         do l2 = l1+1, npole3b
            k=pnum(l2)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
c            if (r2 .le. off2) then
            if (r2 .le. ewaldcut3b_2) then
               r = sqrt(r2)
c               ck = rpole(1,k)
c               dk(1) = rpole(2,k)
c               dk(2) = rpole(3,k)
c               dk(3) = rpole(4,k)
c               qk(1) = rpole(5,k)
c               qk(2) = rpole(6,k)
c               qk(3) = rpole(7,k)
c               qk(4) = rpole(8,k)
c               qk(5) = rpole(9,k)
c               qk(6) = rpole(10,k)
c               qk(7) = rpole(11,k)
c               qk(8) = rpole(12,k)
c               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald3b * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald3b**2
               alsq2n = 0.0d0
               if (aewald3b .gt. 0.0d0) 
     &              alsq2n = 1.0d0 / (sqrtpi*aewald3b)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
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
c               dsc3 = 1.0d0 - scale3*dscale(kk)
c               dsc5 = 1.0d0 - scale5*dscale(kk)
c               dsc7 = 1.0d0 - scale7*dscale(kk)
c               psc3 = 1.0d0 - scale3*pscale(kk)
c               psc5 = 1.0d0 - scale5*pscale(kk)
c               psc7 = 1.0d0 - scale7*pscale(kk)
c               usc3 = 1.0d0 - scale3*uscale(kk)
c               usc5 = 1.0d0 - scale5*uscale(kk)

c               dsc3 = 1.0d0 - scale3*dscale(l2)
c               dsc5 = 1.0d0 - scale5*dscale(l2)
c               dsc7 = 1.0d0 - scale7*dscale(l2)
              ! psc3 = 1.0d0 - scale3*pscale(l2)
              ! psc5 = 1.0d0 - scale5*pscale(l2)
              ! psc7 = 1.0d0 - scale7*pscale(l2)
               usc3 = 1.0d0 - scale3*uscale(l2)
               usc5 = 1.0d0 - scale5*uscale(l2)

c
c     construct necessary auxiliary vectors
c

ccc               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
ccc               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
ccc               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
c               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
c               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
c               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
c               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
c               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
c               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
c               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
c               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
c               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
c               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
c               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
c               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
c               dixr(1) = di(2)*zr - di(3)*yr
c               dixr(2) = di(3)*xr - di(1)*zr
c               dixr(3) = di(1)*yr - di(2)*xr
c               dkxr(1) = dk(2)*zr - dk(3)*yr
c               dkxr(2) = dk(3)*xr - dk(1)*zr
c               dkxr(3) = dk(1)*yr - dk(2)*xr
c               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
c               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
c               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
c               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
c               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
c               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
ccc               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
ccc               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
ccc               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
ccc               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
ccc               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
ccc               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
ccc               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
ccc     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
ccc               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
ccc     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
ccc               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
ccc     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
c               rxqir(1) = yr*qir(3) - zr*qir(2)
c               rxqir(2) = zr*qir(1) - xr*qir(3)
c               rxqir(3) = xr*qir(2) - yr*qir(1)
c               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
c               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
c               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
ccc               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
ccc               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
ccc               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
ccc               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
ccc               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
ccc               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
ccc               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
ccc               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
ccc               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
ccc               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
ccc               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
ccc               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
ccc               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
ccc               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
ccc               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
c               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
c     &                      + qi(7)*uind(3,l2)
c               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
c     &                      + qi(8)*uind(3,l2)
c               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
c     &                      + qi(9)*uind(3,l2)
c               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
c     &                      + qk(7)*uind(3,l1)
c               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
c     &                      + qk(8)*uind(3,l1)
c               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
c     &                      + qk(9)*uind(3,l1)
c               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
c     &                       + qi(7)*uinp(3,l2)
c               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
c     &                       + qi(8)*uinp(3,l2)
c               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
c     &                       + qi(9)*uinp(3,l2)
c               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
c     &                       + qk(7)*uinp(3,l1)
c               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
c     &                       + qk(8)*uinp(3,l1)
c               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
c     &                       + qk(9)*uinp(3,l1)
ccc               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
ccc               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
ccc               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
ccc               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
ccc               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
ccc               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
c               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
c               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
c               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
c               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
c               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
c               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
c               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
c               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
c               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
c               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
c               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
c               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
ccc               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
ccc               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
ccc               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
ccc               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
ccc               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
ccc               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
c               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
c               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
c               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
c               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
c               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
c               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
c               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
c               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
c               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
c               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
c               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)

c
c     calculate the scalar products for permanent components
c
c               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
c               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
c               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
c               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
c               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
c               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
c               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
c               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
c               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
c     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
c     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
               !sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &         !            + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
     &         !            + di(2)*uind(2,l2) + di(3)*uind(3,l2)
               sci(1) = 0.0d0
               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               !sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &         !            + qir(3)*uind(3,l2)
               sci(7) = 0.0d0
               !sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &         !            + qkr(3)*uind(3,l1)
               sci(8) = 0.0d0
c               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
c     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
c     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(1) =0.0d0
               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = 0.0d0
               scip(8) = 0.0d0
c               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
c     &                      + qir(3)*uinp(3,l2)
c               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
c     &                      + qkr(3)*uinp(3,l1)

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
c               gli(1) = ck*sci(3) - ci*sci(4)
c               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
c               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c               gli(6) = sci(1)
c               gli(7) = 2.0d0 * (sci(7)-sci(8))
c               glip(1) = ck*scip(3) - ci*scip(4)
c               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
c               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
c               glip(6) = scip(1)
c               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               !ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &         !             + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               !erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &         !          + rr5*(gli(2)+gli(7))*psc5
     &         !          + rr7*gli(3)*psc7)
               !ei = ei - erli
               !ei = f * ei
               !eptemp = eptemp + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c
c     intermediate variables for induced force terms
c

c               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
c     &                     + 0.5d0*bn(2)*scip(2)
c     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
c     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
c     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
c               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
c               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
c               gfi(4) = 2.0d0 * bn(2)
c               gfi(5) = bn(3) * (sci(4)+scip(4))
c               gfi(6) = -bn(3) * (sci(3)+scip(3))
c               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
c     &                            + (glip(1)+glip(6))*dsc3
c     &                            + scip(2)*usc3)
c     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
c     &                            + (glip(7)+glip(2))*dsc5
c     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
c     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
c               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
c               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
c               gfri(4) = 2.0d0 * rr5
c               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
c               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)


               gfi(1) = !0.5d0*bn(2)*(gli(1)+gli(6))
     &                     0.5d0*bn(2)*scip(2)
     &                     !+ 0.5d0*bn(3)*(gli(2)+gli(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     !+ 0.5d0*bn(4)*(gli(3))
               !gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               !gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(2) = 0.0d0
               gfi(3) = 0.0d0
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = 0.0d0 
               gfi(6) = 0.0d0 
               
               gfri(1) = 0.5d0*rr5*(
     &                             scip(2)*usc3)
     &                 + 0.5d0*rr7*(
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.0d0
               gfri(2) = 0.0d0 
               gfri(3) = 0.0d0
               gfri(4) = 2.0d0 * rr5
               gfri(5) = 0.0d0
               gfri(6) = 0.0d0
c
c     get the induced force with screening
c
c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
c     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
c     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
c     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
c     &            + (sci(4)+scip(4))*bn(2)*di(1)
c     &            + (sci(3)+scip(3))*bn(2)*dk(1)
c     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
c     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
c     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
c     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
c     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
c     &            + (sci(4)+scip(4))*bn(2)*di(2)
c     &            + (sci(3)+scip(3))*bn(2)*dk(2)
c     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
c     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
c               ftm2i(3) = gfi(1)*zr + 0.5d0*
c     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
c     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
c     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
c     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
c     &            + (sci(4)+scip(4))*bn(2)*di(3)
c     &            + (sci(3)+scip(3))*bn(2)*dk(3)
c     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
c     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)

               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (!gfi(2)*(uind(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            !+ gfi(3)*(uind(1,l2))
     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            !+ (sci(4))*bn(2)*di(1)
     &            !+ (sci(3))*bn(2)*dk(1)
     &            !+ gfi(4)*(qkui(1)-qiuk(1))
     &              )
     &            !+ gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (!gfi(2)*(uind(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            !+ gfi(3)*(uind(2,l2))
     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            !+ (sci(4))*bn(2)*di(2)
     &            !+ (sci(3))*bn(2)*dk(2)
     &            !+ gfi(4)*(qkui(2)-qiuk(2))
     &              )
     &            !+ gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (!gfi(2)*(uind(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            !+ gfi(3)*(uind(3,l2))
     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            !+ (sci(4))*bn(2)*di(3)
     &            !+ (sci(3))*bn(2)*dk(3)
     &            !+ gfi(4)*(qkui(3)-qiuk(3))
     &              )
     &            !+ gfi(5)*qir(3) + gfi(6)*qkr(3)

c
c     get the induced force without screening
c
c               ftm2ri(1) = gfri(1)*xr + 0.5d0*
c     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
c     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
c     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
c     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
c     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
c     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
c     &            + (qkuip(1)-qiukp(1))*dsc5)
c     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
c               ftm2ri(2) = gfri(1)*yr + 0.5d0*
c     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
c     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
c     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
c     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
c     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
c     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
c     &            + (qkuip(2)-qiukp(2))*dsc5)
c     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
c               ftm2ri(3) = gfri(1)*zr + 0.5d0*
c     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
c     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
c     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
c     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
c     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
c     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
c     &            + (qkuip(3)-qiukp(3))*dsc5)
c     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)


               ftm2ri(1) = gfri(1)*xr + !0.5d0*
     &           !(- rr3*ck*(uind(1,l1)*psc3)
     &           ! + rr5*sc(4)*(uind(1,l1)*psc5)
     &           ! - rr7*sc(6)*(uind(1,l1)*psc7))
     &           ! + (rr3*ci*(uind(1,l2)*psc3)
     &           ! + rr5*sc(3)*(uind(1,l2)*psc5)
     &           ! + rr7*sc(5)*(uind(1,l2)*psc7))*0.5d0
     &             rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &           ! + 0.5d0*(sci(4)*psc5)*rr5*di(1)
     &           ! + 0.5d0*(sci(3)*psc5)*rr5*dk(1)
     &           ! + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &           ! )
     &           ! + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + ! 0.5d0*
     &           !(- rr3*ck*(uind(2,l1)*psc3)
     &           ! + rr5*sc(4)*(uind(2,l1)*psc5)
     &           ! - rr7*sc(6)*(uind(2,l1)*psc7))
     &           ! + (rr3*ci*(uind(2,l2)*psc3)
     &           ! + rr5*sc(3)*(uind(2,l2)*psc5)
     &           ! + rr7*sc(5)*(uind(2,l2)*psc7))*0.5d0
     &             rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            !+ 0.5d0*(sci(4)*psc5)*rr5*di(2)
     &            !+ 0.5d0*(sci(3)*psc5)*rr5*dk(2)
     &           ! + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &           !  )
     &           ! + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + !0.5d0*
     &           !(- rr3*ck*(uind(3,l1)*psc3)
     &           ! + rr5*sc(4)*(uind(3,l1)*psc5)
     &           ! - rr7*sc(6)*(uind(3,l1)*psc7))
     &           ! + (rr3*ci*(uind(3,l2)*psc3)
     &           ! + rr5*sc(3)*(uind(3,l2)*psc5)
     &           ! + rr7*sc(5)*(uind(3,l2)*psc7))*0.5d0
     &             rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &           ! + 0.5d0*(sci(4)*psc5)*rr5*di(3)
     &           ! + 0.5d0*(sci(3)*psc5)*rr5*dk(3)
     &           ! + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &           ! )
     &           ! + gfri(5)*qir(3) + gfri(6)*qkr(3)

c
c     account for partially excluded induced interactions
c

ccc               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
ccc     &                                  +(glip(1)+glip(6))*dscale(kk))
ccc               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
ccc     &                                  +(glip(2)+glip(7))*dscale(kk))
ccc               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
ccc     &                                  +glip(3)*dscale(kk))

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l2)
c     &                                  +(glip(1)+glip(6))*dscale(l2))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l2)
c     &                                  +(glip(2)+glip(7))*dscale(l2))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l2)
c     &                                  +glip(3)*dscale(l2))


               temp3 = 0.0d0
               temp5 = 0.0d0
               temp7 = 0.0d0

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
               temp3 = 0.5d0 * rr3 * uscale(l2) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(l2)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (bn(2)*scip(2)
     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
     &                     - rr7*(scip(3)*sci(4)
     &                           +sci(3)*scip(4))*usc5)
                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &                          +sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &                          +sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &                          +sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &                        + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &                        + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &                        + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if

c
c     intermediate variables for induced torque terms
c
c               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
c               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
c               gti(4) = gfi(4)
c               gti(5) = gfi(5)
c               gti(6) = gfi(6)
c               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
c               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
c               gtri(4) = gfri(4)
c               gtri(5) = gfri(5)
c               gtri(6) = gfri(6)

c               gti(2) = 0.0d0 
c               gti(3) = 0.0d0
c               gti(4) = gfi(4)
c               gti(5) = gfi(5)
c               gti(6) = gfi(6)
c               gtri(2) = 0.0d0
c               gtri(3) = 0.0d0
c               gtri(4) = gfri(4)
c               gtri(5) = gfri(5)
c               gtri(6) = gfri(6)

c
c     get the induced torque with screening
c

c               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
c     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
c     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
c               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
c     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
c     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
c               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
c     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
c     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
c               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
c     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
c     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
c               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
c     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
c     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
c               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
c     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
c     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)


c               ttm2i(1) = -bn(1)*(dixuk(1))*0.5d0
c     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
c     &           )*0.5d0 - gti(5)*rxqir(1)
c               ttm2i(2) = -bn(1)*(dixuk(2))*0.5d0
c     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
c     &           )*0.5d0 - gti(5)*rxqir(2)
c               ttm2i(3) = -bn(1)*(dixuk(3))*0.5d0
c     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
c     &           )*0.5d0 - gti(5)*rxqir(3)
c               ttm3i(1) = -bn(1)*(dkxui(1))*0.5d0
c     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
c     &           )*0.5d0 - gti(6)*rxqkr(1)
c               ttm3i(2) = -bn(1)*(dkxui(2))*0.5d0
c     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
c     &           )*0.5d0 - gti(6)*rxqkr(2)
c               ttm3i(3) = -bn(1)*(dkxui(3))*0.5d0
c     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
c     &           )*0.5d0 - gti(6)*rxqkr(3)

c
c     get the induced torque without screening
c
c               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
c     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
c               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
c     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
c               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
c     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
c               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
c     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
c               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
c     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
c               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
c     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)

c               ttm2ri(1) = -rr3*(dixuk(1)*psc3)*0.5d0
c     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &           )*0.5d0 - gtri(5)*rxqir(1)
c               ttm2ri(2) = -rr3*(dixuk(2)*psc3)*0.5d0
c     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &           )*0.5d0 - gtri(5)*rxqir(2)
c               ttm2ri(3) = -rr3*(dixuk(3)*psc3)*0.5d0
c     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &           )*0.5d0 - gtri(5)*rxqir(3)
c               ttm3ri(1) = -rr3*(dkxui(1)*psc3)*0.5d0
c     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &           )*0.5d0 - gtri(6)*rxqkr(1)
c               ttm3ri(2) = -rr3*(dkxui(2)*psc3)*0.5d0
c     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &           )*0.5d0 - gtri(6)*rxqkr(2)
c               ttm3ri(3) = -rr3*(dkxui(3)*psc3)*0.5d0
c     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &           )*0.5d0 - gtri(6)*rxqkr(3)

c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                 ! ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                 ! ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c
c     increment gradient due to force and torque on first site
c

c               deptemp2(1,i) = deptemp2(1,i) + ftm2i(1)
c               deptemp2(2,i) = deptemp2(2,i) + ftm2i(2)
c               deptemp2(3,i) = deptemp2(3,i) + ftm2i(3)
               deptemp2(1,l1) =deptemp2(1,l1) + ftm2i(1)
               deptemp2(2,l1) =deptemp2(2,l1) + ftm2i(2)
               deptemp2(3,l1) =deptemp2(3,l1) + ftm2i(3)

c               call torque_3b (deptemp2,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               !call torque_3b_new(npole3b,pnum,deptemp2,i,
     &         ! ttm2,ttm2i,frcxi,frcyi,frczi)
c               do j=1,3
c                frcxi(j) = 0.0d0
c                frcyi(j) = 0.0d0
c                frczi(j) = 0.0d0
c               end do

c
c     increment gradient due to force and torque on second site
c

c               deptemp2(1,k) = deptemp2(1,k) - ftm2i(1)
c               deptemp2(2,k) = deptemp2(2,k) - ftm2i(2)
c               deptemp2(3,k) = deptemp2(3,k) - ftm2i(3)
               deptemp2(1,l2) = deptemp2(1,l2) - ftm2i(1)
               deptemp2(2,l2) = deptemp2(2,l2) - ftm2i(2)
               deptemp2(3,l2) = deptemp2(3,l2) - ftm2i(3)

c               call torque_3b (deptemp2,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
               !call torque_3b_new(npole3b,pnum,deptemp2,k,
     &         !   ttm3,ttm3i,frcxk,frcyk,frczk)
c               do j=1,3
c                frcxk(j) = 0.0d0
c                frcyk(j) = 0.0d0
c                frczk(j) = 0.0d0
c               end do

c               iaz = zaxis(i)
c               iax = xaxis(i)
c               iay = yaxis(i)
c               kaz = zaxis(k)
c               kax = xaxis(k)
c               kay = yaxis(k)
c               if (iaz .eq. 0)  iaz = ii
c               if (iax .eq. 0)  iax = ii
c               if (iay .eq. 0)  iay = ii
c               if (kaz .eq. 0)  kaz = kk
c               if (kax .eq. 0)  kax = kk
c               if (kay .eq. 0)  kay = kk
c               xiz = x(iaz) - x(ii)
c               yiz = y(iaz) - y(ii)
c               ziz = z(iaz) - z(ii)
c               xix = x(iax) - x(ii)
c               yix = y(iax) - y(ii)
c               zix = z(iax) - z(ii)
c               xiy = x(iay) - x(ii)
c               yiy = y(iay) - y(ii)
c               ziy = z(iay) - z(ii)
c               xkz = x(kaz) - x(kk)
c               ykz = y(kaz) - y(kk)
c               zkz = z(kaz) - z(kk)
c               xkx = x(kax) - x(kk)
c               ykx = y(kax) - y(kk)
c               zkx = z(kax) - z(kk)
c               xky = x(kay) - x(kk)
c               yky = y(kay) - y(kk)
c               zky = z(kay) - z(kk)

               vxx = -xr*(ftm2i(1)) !+ xix*frcxi(1)
     &                  !+ xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  !+ xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) !+ yix*frcxi(1)
     &                  !+ yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  !+ yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) !+ zix*frcxi(1)
     &                !  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                !  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) !+ yix*frcxi(2)
     &                 ! + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                 ! + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) !+ zix*frcxi(2)
     &                 ! + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                 ! + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) !+ zix*frcxi(3)
     &                 ! + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                 ! + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

            end if
c          end if
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
         do j=1,npole3b
c            pscale(j)=1.0d0
c            dscale(j)=1.0d0
            uscale(j)=1.0d0
         end do
      end do

      if (use_replica) then
c
c     calculate interactions with other unit cells
c
c      print*,'After use_replica!',moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After use_replica l1 Pnum(l1)",l1,pnum(l1)
c      end do 
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
c         ci = rpole(1,i)
c         di(1) = rpole(2,i)
c         di(2) = rpole(3,i)
c         di(3) = rpole(4,i)
c         qi(1) = rpole(5,i)
c         qi(2) = rpole(6,i)
c         qi(3) = rpole(7,i)
c         qi(4) = rpole(8,i)
c         qi(5) = rpole(9,i)
c         qi(6) = rpole(10,i)
c         qi(7) = rpole(11,i)
c         qi(8) = rpole(12,i)
c         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
c         do j = 1, n12(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i12(j,ii)) then
c                 pscale(kk) = p2scale
c                 goto 41
c               end if
c            end do
c   41    continue
c         end do
c         do j = 1, n13(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i13(j,ii)) then
c                 pscale(kk) = p3scale
c                 goto 42
c               end if
c            end do
c   42    continue
c         end do
c         do j = 1, n14(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i14(j,ii)) then
c                 pscale(kk) = p4scale
c                 goto 43
c               end if
c            end do
c   43    continue
c         end do
c         do j = 1, n15(ii)
c            do kk=1,npole3b
c               if(pnum(kk).eq.i15(j,ii)) then
c                 pscale(kk) = p5scale
c                 goto 44
c               end if
c            end do
c   44    continue
c         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip11(j,ii)) then
                 uscale(kk) = u1scale
                 goto 45
               end if
            end do
   45    continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip12(j,ii)) then
                 uscale(kk) = u2scale
                 goto 46
               end if
            end do
   46    continue
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip13(j,ii)) then
                 uscale(kk) = u3scale
                 goto 47
               end if
            end do
   47    continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip14(j,ii)) then
                 uscale(kk) = u4scale
                 goto 48
               end if
            end do
   48    continue
         end do
         do l2 = l1, npole3b
            k=pnum(l2)
            kk = ipole(k)
            do jcell=1,ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not.(use_polymer .and. r2.le.polycut2)) then
c              mscale(kk)=1.0d0
c              pscale(kk)=1.0d0
c              dscale(kk)=1.0d0
c              uscale(kk)=1.0d0
c              pscale(l2)=1.0d0
c              dscale(l2)=1.0d0
              uscale(l2)=1.0d0
            end if
c            if (r2 .le. off2) then
            if (r2 .le. ewaldcut3b_2) then
               r = sqrt(r2)
c               ck = rpole(1,k)
c               dk(1) = rpole(2,k)
c               dk(2) = rpole(3,k)
c               dk(3) = rpole(4,k)
c               qk(1) = rpole(5,k)
c               qk(2) = rpole(6,k)
c               qk(3) = rpole(7,k)
c               qk(4) = rpole(8,k)
c               qk(5) = rpole(9,k)
c               qk(6) = rpole(10,k)
c               qk(7) = rpole(11,k)
c               qk(8) = rpole(12,k)
c               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald3b * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald3b**2
               alsq2n = 0.0d0
               if (aewald3b .gt. 0.0d0)  
     &             alsq2n = 1.0d0 / (sqrtpi*aewald3b)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
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

c               dsc3 = 1.0d0 - scale3*dscale(kk)
c               dsc5 = 1.0d0 - scale5*dscale(kk)
c               dsc7 = 1.0d0 - scale7*dscale(kk)
c               psc3 = 1.0d0 - scale3*pscale(kk)
c               psc5 = 1.0d0 - scale5*pscale(kk)
c               psc7 = 1.0d0 - scale7*pscale(kk)
c               usc3 = 1.0d0 - scale3*uscale(kk)
c               usc5 = 1.0d0 - scale5*uscale(kk)

c               dsc3 = 1.0d0 - scale3*dscale(l2)
c               dsc5 = 1.0d0 - scale5*dscale(l2)
c               dsc7 = 1.0d0 - scale7*dscale(l2)
c               psc3 = 1.0d0 - scale3*pscale(l2)
c               psc5 = 1.0d0 - scale5*pscale(l2)
c               psc7 = 1.0d0 - scale7*pscale(l2)
               usc3 = 1.0d0 - scale3*uscale(l2)
               usc5 = 1.0d0 - scale5*uscale(l2)

c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
c               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
c               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
c               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
c               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
c               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
c               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
c               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
c               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
c               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
c               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
c               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
c               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
c               dixr(1) = di(2)*zr - di(3)*yr
c               dixr(2) = di(3)*xr - di(1)*zr
c               dixr(3) = di(1)*yr - di(2)*xr
c               dkxr(1) = dk(2)*zr - dk(3)*yr
c               dkxr(2) = dk(3)*xr - dk(1)*zr
c               dkxr(3) = dk(1)*yr - dk(2)*xr
c               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
c               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
c               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
c               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
c               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
c               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
c               rxqir(1) = yr*qir(3) - zr*qir(2)
c               rxqir(2) = zr*qir(1) - xr*qir(3)
c               rxqir(3) = xr*qir(2) - yr*qir(1)
c               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
c               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
c               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
c               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
c     &                      + qi(7)*uind(3,l2)
c               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
c     &                      + qi(8)*uind(3,l2)
c               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
c     &                      + qi(9)*uind(3,l2)
c               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
c     &                      + qk(7)*uind(3,l1)
c               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
c     &                      + qk(8)*uind(3,l1)
c               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
c     &                      + qk(9)*uind(3,l1)
c               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
c     &                       + qi(7)*uinp(3,l2)
c               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
c     &                       + qi(8)*uinp(3,l2)
c               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
c     &                       + qi(9)*uinp(3,l2)
c               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
c     &                       + qk(7)*uinp(3,l1)
c               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
c     &                       + qk(8)*uinp(3,l1)
c               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
c     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
c               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
c               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
c               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
c               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
c               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
c               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
c               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
c               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
c               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
c               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
c               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
c               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
c               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
c               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
c               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
c               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
c               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
c               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
c               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
c               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
c               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
c               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
c               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
c               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
c               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
c               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
c               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
c               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
c               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
c               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
c               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
c               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
c     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
c     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &                     + qir(3)*uind(3,l2)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for induced components
c
c               gli(1) = ck*sci(3) - ci*sci(4)
c               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
c               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c               gli(6) = sci(1)
c               gli(7) = 2.0d0 * (sci(7)-sci(8))
c               glip(1) = ck*scip(3) - ci*scip(4)
c               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
c               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
c               glip(6) = scip(1)
c               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
c               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
c     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
c               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
c     &                   + rr5*(gli(2)+gli(7))*psc5
c     &                   + rr7*gli(3)*psc7)
c               ei = ei - erli
c               ei = f * ei
c               if (ii .eq. kk) then
c                 ei = 0.5d0 * ei
c               end if
c               eptemp = eptemp + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l2)
c     &                                  +(glip(1)+glip(6))*dscale(l2))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l2)
c     &                                  +(glip(2)+glip(7))*dscale(l2))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l2)
c     &                                  +glip(3)*dscale(l2))

               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(l2) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(l2)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
              if (use_polymer .and. r2.le.polycut2) then
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              else
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              end if
               if (ii .eq. kk) then
                  do j = 1, 3
                     ftm2i(j) = 0.5d0 * ftm2i(j)
                     ttm2i(j) = 0.5d0 * ttm2i(j)
                     ttm3i(j) = 0.5d0 * ttm3i(j)
                  end do
               end if

c
c     increment gradient due to force and torque on first site
c

c               deptemp2(1,i) = deptemp2(1,i) + ftm2i(1)
c               deptemp2(2,i) = deptemp2(2,i) + ftm2i(2)
c               deptemp2(3,i) = deptemp2(3,i) + ftm2i(3)
c               call torque_3b (deptemp2,i,
c     &           ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

c               deptemp2(1,k) = deptemp2(1,k) - ftm2i(1)
c               deptemp2(2,k) = deptemp2(2,k) - ftm2i(2)
c               deptemp2(3,k) = deptemp2(3,k) - ftm2i(3)
c               call torque_3b (deptemp2,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               deptemp2(1,l1) =deptemp2(1,l1) + ftm2i(1)
               deptemp2(2,l1) =deptemp2(2,l1) + ftm2i(2)
               deptemp2(3,l1) =deptemp2(3,l1) + ftm2i(3)

c               call torque_3b (deptemp2,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               call torque_3b_new(npole3b,pnum,deptemp2,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

c               deptemp2(1,k) = deptemp2(1,k) - ftm2i(1)
c               deptemp2(2,k) = deptemp2(2,k) - ftm2i(2)
c               deptemp2(3,k) = deptemp2(3,k) - ftm2i(3)
               deptemp2(1,l2) = deptemp2(1,l2) - ftm2i(1)
               deptemp2(2,l2) = deptemp2(2,l2) - ftm2i(2)
               deptemp2(3,l2) = deptemp2(3,l2) - ftm2i(3)

c               call torque_3b (deptemp2,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
               call torque_3b_new(npole3b,pnum,deptemp2,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)


               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
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
c               vxx = -xr*ftm2i(1)+xix*frcxi(1)+xiy*frcyi(1)+xiz*frczi(1)
c               vyx = -yr*ftm2i(1)+yix*frcxi(1)+yiy*frcyi(1)+yiz*frczi(1)
c               vzx = -zr*ftm2i(1)+zix*frcxi(1)+ziy*frcyi(1)+ziz*frczi(1)
c               vyy = -yr*ftm2i(2)+yix*frcxi(2)+yiy*frcyi(2)+yiz*frczi(2)
c               vzy = -zr*ftm2i(2)+zix*frcxi(2)+ziy*frcyi(2)+ziz*frczi(2)
c               vzz = -zr*ftm2i(3)+zix*frcxi(3)+ziy*frcyi(3)+ziz*frczi(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
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

            end if
         end do
c         end if
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
         do j=1,npole3b
c            dscale(j)=1.0d0
            uscale(j)=1.0d0
c            pscale(j)=1.0d0
         end do
      end do
      end if
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c
c     perform deallocation of some local arrays
c
c      deallocate (mscale)
c      deallocate (pscale)
c      deallocate (dscale)
      deallocate (uscale)
c      print*,'ep after ereal1c_3b=',eptemp,moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After ereal1c_3b l1 Pnum(l1)",l1,pnum(l1)
c      end do

      return
      end

