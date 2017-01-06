c
c
      subroutine empole1d_3b_Polar_orig_nopriordiromp1b(npole3b,pnum,
     &   eptemp,deptemp,virtemp,cnt)
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
      use cho
      use pcg
      use neigh2clust
      use ewald
      use math
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
      integer npole3b,pnum(*),l1,l3
      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      logical proceed,usei,usek
      character*6 mode
      integer cnt,kkk
      real*8 term,fterm,uix,uiy,uiz,uii,trq(3),trqi(3)
      real*8 frcx,frcy,frcz

      eptemp = 0.0d0
c
c   Zero out temporary gradient of polarization energy
c

      do l1 = 1, npole3b
        do j = 1, 3
          deptemp(j,l1) = 0.0d0
        end do
      end do
c
c   Zero out temporary virial
c
      do i=1,3
         do j=1,3
           virtemp(j,i)=0.0d0
         end do
      end do

      if(usecholesky) then
       call induce0a_3b_cholesky_nopriordir(npole3b,pnum,
     & uind,uinp) 
      else if (usepcgscf) then
       call induce0a_3b_scf_vacomp1bclust(npole3b,pnum,uind,uinp,cnt)
       !print*,"Done w induce0a_3b_scf_vacomp1bclust"
      else
       call induce0a_3b_PolelecOnly_orig_nopriordir(npole3b,pnum,
     & uind,uinp)
      end if


      call emrecip1_3b_Polar_new(npole3b,pnum,uind,uinp,
     & eptemp,deptemp,virtemp)
      !print*,"Done w emrecip1_3b_Polar_new in 1bclust"

      call ereal1d_3b_1bclust(npole3b,pnum,uind,uinp,cnt,
     & eptemp,deptemp,virtemp)
      !print*,"Done w ereal1d_3b_1bclust in 1bclust"

      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      !do i = 1, npole
      do l1 =1, npole3b
         i=pnum(l1)
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
         uix = uind(1,l1)
         uiy = uind(2,l1)
         uiz = uind(3,l1)
c         cii = ci*ci
c         dii = dix*dix + diy*diy + diz*diz
c         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
c     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
c         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
c         em = em + e
         eptemp = eptemp + ei
      end do

      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      !do i = 1, npole
      do l1 = 1, npole3b
         i=pnum(l1)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,l1)+uinp(1,l1))
         uiy = 0.5d0 * (uind(2,l1)+uinp(2,l1))
         uiz = 0.5d0 * (uind(3,l1)+uinp(3,l1))
         trqi(1) = term * (diy*uiz-diz*uiy)
         trqi(2) = term * (diz*uix-dix*uiz)
         trqi(3) = term * (dix*uiy-diy*uix)
         !call torque (i,trq,trqi,frcx,frcy,frcz)
        call torque_3b_new(npole3b,pnum,deptemp,i,
     &            trq,trqi,frcx,frcy,frcz)
      end do
      !print*,"Completed all steps of empole1d 1bclust"
      return
      end
