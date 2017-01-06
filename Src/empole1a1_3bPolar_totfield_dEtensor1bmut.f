      subroutine empole1c_3bPolar_totfield_dEtensor1bmut(
     &  dep3bt1,virep3bt1,uindt,uinpt)
      use bound
      use sizes
      use atoms
      use chgpot
      use couple
      use group
      use mpole
      use polar, only: polarity, thole, pdamp
c      use polar
      use polgrp
      use polpot
      use usage
      use cell
      use boxes
      use ewald     
      use math
      use pme
      use chunks
      use fft
      use totfield
      use limits 
      use dEtensor
      use mpidat
      use molcul
      implicit none 
      integer i,j,ii,l1,l3
      integer moli1,moli2,moli3
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 eptemp,virep3bt(3,3)
      real*8 dep3bt(3,npole)
c      real*8 uind(3,npole3b)
c      real*8 uinp(3,npole3b)
c      integer, allocatable :: pnum(:)
c      real*8, allocatable :: uind(:,:)
c      real*8, allocatable :: uinp(:,:)
      real*8, allocatable :: qgrid3b(:,:,:,:)
      real*8, allocatable :: qfac3b(:,:,:)
      real*8, allocatable :: thetai1_3b(:,:,:)
      real*8, allocatable :: thetai2_3b(:,:,:)
      real*8, allocatable :: thetai3_3b(:,:,:)
      integer, allocatable :: pmetable3b(:,:)
      integer, allocatable :: igrid3b(:,:)
      integer*8 planf3b
      integer*8 planb3b,k
      integer, allocatable :: iprime3b(:,:)
      real*8, allocatable :: ffttable3b(:,:)
      integer, allocatable :: embedlst(:,:)
      integer, allocatable :: nembedlst(:)
      real*8 r2,xr,yr,zr,ewaldcut2
      logical doi,dok
      real*8 t1xxd,t1xyd,t1xzd,t1yxd,t1yyd,t1yzd,t1zxd,t1zyd,t1zzd
      real*8 t2xxd,t2xyd,t2xzd,t2yxd,t2yyd,t2yzd,t2zxd,t2zyd,t2zzd
      real*8 t1xxp,t1xyp,t1xzp,t1yxp,t1yyp,t1yzp,t1zxp,t1zyp,t1zzp
      real*8  t2xxp,t2xyp,t2xzp,t2yxp,t2yyp,t2yzp,t2zxp,t2zyp,t2zzp
      real*8 tau1xxd,tau1xyd,tau1xzd,tau1yxd,tau1yyd,tau1yzd,tau1zxd
      real*8 tau1zyd,tau1zzd
      real*8 tau2xxd,tau2xyd,tau2xzd,tau2yxd,tau2yyd,tau2yzd,tau2zxd
      real*8 tau2zyd,tau2zzd
      real*8 tau1xxp,tau1xyp,tau1xzp,tau1yxp,tau1yyp,tau1yzp,tau1zxp
      real*8 tau1zyp,tau1zzp
      real*8 tau2xxp,tau2xyp,tau2xzp,tau2yxp,tau2yyp,tau2yzp,tau2zxp
      real*8 tau2zyp,tau2zzp
      real*8 frcz1xxd,frcz1xyd,frcz1xzd,frcz1yxd,frcz1yyd,frcz1yzd
      real*8 frcz1zxd,frcz1zyd,frcz1zzd
      real*8 frcz2xxd,frcz2xyd,frcz2xzd,frcz2yxd,frcz2yyd,frcz2yzd
      real*8 frcz2zxd,frcz2zyd,frcz2zzd
      real*8 frcz1xxp,frcz1xyp,frcz1xzp,frcz1yxp,frcz1yyp,frcz1yzp
      real*8 frcz1zxp,frcz1zyp,frcz1zzp
      real*8 frcz2xxp,frcz2xyp,frcz2xzp,frcz2yxp,frcz2yyp,frcz2yzp
      real*8 frcz2zxp,frcz2zyp,frcz2zzp
      real*8 frcx1xxd,frcx1xyd,frcx1xzd,frcx1yxd,frcx1yyd,frcx1yzd
      real*8 frcx1zxd,frcx1zyd,frcx1zzd
      real*8 frcx2xxd,frcx2xyd,frcx2xzd,frcx2yxd,frcx2yyd,frcx2yzd
      real*8 frcx2zxd,frcx2zyd,frcx2zzd
      real*8 frcx1xxp,frcx1xyp,frcx1xzp,frcx1yxp,frcx1yyp,frcx1yzp
      real*8 frcx1zxp,frcx1zyp,frcx1zzp
      real*8 frcx2xxp,frcx2xyp,frcx2xzp,frcx2yxp,frcx2yyp,frcx2yzp
      real*8 frcx2zxp,frcx2zyp,frcx2zzp
      real*8 frcy1xxd,frcy1xyd,frcy1xzd,frcy1yxd,frcy1yyd,frcy1yzd
      real*8 frcy1zxd,frcy1zyd,frcy1zzd
      real*8 frcy2xxd,frcy2xyd,frcy2xzd,frcy2yxd,frcy2yyd,frcy2yzd
      real*8 frcy2zxd,frcy2zyd,frcy2zzd
      real*8 frcy1xxp,frcy1xyp,frcy1xzp,frcy1yxp,frcy1yyp,frcy1yzp
      real*8 frcy1zxp,frcy1zyp,frcy1zzp
      real*8 frcy2xxp,frcy2xyp,frcy2xzp,frcy2yxp,frcy2yyp,frcy2yzp
      real*8 frcy2zxp,frcy2zyp,frcy2zzp
      real*8 frcxi(3),frcyi(3),frczi(3),frcxk(3),frcyk(3),frczk(3)
      real*8 vxx,vyx,vzx,vyy,vzy,vzz
      integer m,k1,l2,length,start,last,moli1rmndr
      real*8 termx,termy,termz
      real*8 taux,tauy,tauz
      real*8, allocatable :: depi(:,:)
      real*8, allocatable :: depk(:,:)
      real*8 virt(3,3)
      integer iaz,iax,iay,kaz,kax,kay 
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      integer npole3b,pnum(30),len
      real*8 deptemp2(3,30),virtemp(3,3)
      real*8 uindt(3,npole),uindtemp(3,30)
      real*8 uinpt(3,npole),uinptemp(3,30)
      real*8 dep3bt1(3,npole),virep3bt1(3,3) 
      character*6 mode
      integer taskid_offset 
      ewaldcut2=ewaldcut*ewaldcut


      f = electric / dielec
      len=3
      if(numtasks_emreal.eq.int((numtasks-2)/2)) then
         taskid_offset=taskid-2
      else
         taskid_offset=taskid
      end if
      mode = 'EWALD'
      call switch (mode)

!$OMP PARALLEL default(private) shared(start_emreal2,last_emreal2,
!$OMP& taskid_offset,len,dep3bt1,virep3bt1,uindt,uinpt)
!$OMP DO reduction(+:dep3bt1,virep3bt1) schedule(guided)
      do k=start_emreal2(taskid_offset),last_emreal2(taskid_offset),len
         npole3b=3
         pnum(1)=k
         pnum(2)=k+1
         pnum(3)=k+2
         call induce0a_3b_PolelecOnly_totfieldnpole(npole3b,pnum,
     &           uindtemp,uinptemp)
         call empole1a1_3b_Polar_totfieldnpole_mod(npole3b,pnum,
     &        uindtemp,uinptemp,deptemp2,virtemp)
         !call induce0c_3b_new(npole3b,pnum,
     &   !          uindtemp,uinptemp)
         !call ereal1c1_3b_totfieldnpole(npole3b,pnum,
     &   !     uindtemp,uinptemp,deptemp2,virtemp)
         !call empole1c_3b_Polar_totfield_901(
     &   ! npole3b,pnum,uindtemp,uinptemp,deptemp2,virtemp) 

         do l1=1,npole3b
            i=pnum(l1)
            dep3bt1(1,i)=dep3bt1(1,i)+deptemp2(1,l1)
            dep3bt1(2,i)=dep3bt1(2,i)+deptemp2(2,l1)
            dep3bt1(3,i)=dep3bt1(3,i)+deptemp2(3,l1)
         end do
         do i=1,3
            do j=1,3
              virep3bt1(j,i)=virep3bt1(j,i)+virtemp(j,i)
            end do
         end do
         do l1=1,npole3b
            i=pnum(l1)
            uindt(1,i)=uindtemp(1,l1)
            uindt(2,i)=uindtemp(2,l1)
            uindt(3,i)=uindtemp(3,l1)
            uinpt(1,i)=uinptemp(1,l1)
            uinpt(2,i)=uinptemp(2,l1)
            uinpt(3,i)=uinptemp(3,l1)
         end do
      end do   
!$OMP END DO
!$OMP END PARALLEL


c      do m=ntpair_start_rmndr(taskid),ntpair_last_rmndr(taskid)

      !aewald3b=aewald
c      term = 2.0d0 * aewald3b * aewald3b
c      fterm = -f * aewald3b / sqrtpi
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         dix = rpole(2,i)
c         diy = rpole(3,i)
c         diz = rpole(4,i)
c         uix = uind(1,i)
c         uiy = uind(2,i)
c         uiz = uind(3,i)
c         uii = dix*uix + diy*uiy + diz*uiz
c         ei = fterm * term * uii / 3.0d0
c         eptemp = eptemp + ei
c      end do

c
c     compute the self-energy torque term due to induced dipole
c

c      trq(1) = 0.0d0
c      trq(2) = 0.0d0
c      trq(3) = 0.0d0
c      term = (4.0d0/3.0d0) * f * aewald3b**3 / sqrtpi
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         dix = rpole(2,i)
c         diy = rpole(3,i)
c         diz = rpole(4,i)
c         uix = 0.5d0 * (uind(1,i)+uinp(1,i))
c         uiy = 0.5d0 * (uind(2,i)+uinp(2,i))
c         uiz = 0.5d0 * (uind(3,i)+uinp(3,i))
c         trqi(1) = term * (diy*uiz-diz*uiy)
c         trqi(2) = term * (diz*uix-dix*uiz)
c         trqi(3) = term * (dix*uiy-diy*uix)
c         call torque_3b_new_npole(deptemp,i,
c     &      trq,trqi,frcx,frcy,frcz)
c      end do
c      deallocate(pnum)
c      deallocate(uind)
c      deallocate(uinp)
c      deallocate(depi)
c      deallocate(depk)
      return
      end

