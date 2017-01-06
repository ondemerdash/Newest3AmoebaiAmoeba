      subroutine empole1c_3b_Polar_totfield_dEtensor1b_nouselist(
     &  deptemp,virtemp)
      use bound
      use sizes
      use atoms
      use chgpot
      use couple
      use group
      use mpole
c      use polar, only: polarity, thole, pdamp
      use polar
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
      real*8 eptemp,virtemp(3,3)
      integer npole3b
      real*8 deptemp(3,npole)
c      real*8 uind(3,npole3b)
c      real*8 uinp(3,npole3b)
      integer, allocatable :: pnum(:)
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
      integer m,k1,l2,length,start,last,moli1rmndr
      real*8 termx,termy,termz
      real*8, allocatable :: depi(:,:)
      real*8, allocatable :: depk(:,:)
      ewaldcut2=ewaldcut*ewaldcut

c      length=last-start+1
c      allocate(pnum(3*length))
c      allocate(uind(3,3*length))
c      allocate(uinp(3,3*length))
      allocate(depi(3,npole))
      allocate(depk(3,npole))
      do i=1,npole
         depi(1,i)=0.0d0
         depi(2,i)=0.0d0
         depi(3,i)=0.0d0
         depk(1,i)=0.0d0
         depk(2,i)=0.0d0
         depk(3,i)=0.0d0      
      end do
c      do l1=0,length-1
c         pnum(3*l1+1)=imol(1,start+l1)
c         pnum(3*l1+2)=imol(1,start+l1)+1
c         pnum(3*l1+3)=imol(2,start+l1)
c      end do      
c      npole3b=3*length

c      do l1 = 1, npole3b
c        i=pnum(l1)
c        do j = 1, 3
c          uind(j,l1) = polarity(i) * fieldnpole(j,i) 
c          uinp(j,l1) = polarity(i) * fieldpnpole(j,i) 
c        end do
c      end do

      f = electric / dielec
      
c       do l1=1,npole3b
c          i=pnum(l1)
c          eptemp = eptemp -0.5d0*f*( uind(1,i)*fieldpnpole(1,i) +
c     &    uind(2,i)*fieldpnpole(2,i) + uind(3,i)*fieldpnpole(3,i))
c       end do

!$OMP PARALLEL default(private) shared(ntpair_start,ntpair_last,
!$OMP& taskid,dEindex,dEd1,dEd2,dEp1,dEp2,x,y,z,
!$OMP& ewaldcut2,depi,depk,uind,uinp)
!$OMP DO reduction(+:depi,depk)
!$OMP& schedule(guided)
      do m=ntpair_start(taskid),ntpair_last(taskid)
         i=dEindex(1,m)
         k=dEindex(2,m)
        
          t2xxd=dEd2(1,m)
          t2xyd=dEd2(2,m)
          t2xzd=dEd2(3,m)
          t2yxd=dEd2(4,m)
          t2yyd=dEd2(5,m)
          t2yzd=dEd2(6,m)
          t2zxd=dEd2(7,m)
          t2zyd=dEd2(8,m)
          t2zzd=dEd2(9,m)

          t2xxp=dEp2(1,m)
          t2xyp=dEp2(2,m)
          t2xzp=dEp2(3,m)
          t2yxp=dEp2(4,m)
          t2yyp=dEp2(5,m)
          t2yzp=dEp2(6,m)
          t2zxp=dEp2(7,m)
          t2zyp=dEp2(8,m)
          t2zzp=dEp2(9,m)

          termx=uind(1,i)*t2xxp+ uind(2,i)*t2xyp + uind(3,i)*t2xzp
          termy=uind(1,i)*t2yxp+ uind(2,i)*t2yyp + uind(3,i)*t2yzp
          termz=uind(1,i)*t2zxp+ uind(2,i)*t2zyp + uind(3,i)*t2zzp

c          deptemp(1,i)=deptemp(1,i)+termx
c          deptemp(2,i)=deptemp(2,i)+termy
c          deptemp(3,i)=deptemp(3,i)+termz

c          deptemp(1,k)=deptemp(1,k)-termx
c          deptemp(2,k)=deptemp(2,k)-termy
c          deptemp(3,k)=deptemp(3,k)-termz

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz

          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

          termx = uinp(1,i)*t2xxd+ uinp(2,i)*t2xyd + uinp(3,i)*t2xzd
          termy = uinp(1,i)*t2yxd+ uinp(2,i)*t2yyd + uinp(3,i)*t2yzd
          termz = uinp(1,i)*t2zxd+ uinp(2,i)*t2zyd + uinp(3,i)*t2zzd

c          deptemp(1,i)=deptemp(1,i)+termx
c          deptemp(2,i)=deptemp(2,i)+termy
c          deptemp(3,i)=deptemp(3,i)+termz
c
c          deptemp(1,k)=deptemp(1,k)-termx
c          deptemp(2,k)=deptemp(2,k)-termy
c          deptemp(3,k)=deptemp(3,k)-termz
          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

c        end if

c        if(dok) then
          t1xxd=dEd1(1,m)
          t1xyd=dEd1(2,m)
          t1xzd=dEd1(3,m)
          t1yxd=dEd1(4,m)
          t1yyd=dEd1(5,m)
          t1yzd=dEd1(6,m)
          t1zxd=dEd1(7,m)
          t1zyd=dEd1(8,m)
          t1zzd=dEd1(9,m)

          t1xxp=dEp1(1,m)
          t1xyp=dEp1(2,m)
          t1xzp=dEp1(3,m)
          t1yxp=dEp1(4,m)
          t1yyp=dEp1(5,m)
          t1yzp=dEp1(6,m)
          t1zxp=dEp1(7,m)
          t1zyp=dEp1(8,m)
          t1zzp=dEp1(9,m)

          termx=uind(1,k)*t1xxp+ uind(2,k)*t1xyp + uind(3,k)*t1xzp
          termy=uind(1,k)*t1yxp+ uind(2,k)*t1yyp + uind(3,k)*t1yzp
          termz=uind(1,k)*t1zxp+ uind(2,k)*t1zyp + uind(3,k)*t1zzp 
c          deptemp(1,i)=deptemp(1,i)+termx
c          deptemp(2,i)=deptemp(2,i)+termy
c          deptemp(3,i)=deptemp(3,i)+termz
c
c          deptemp(1,k)=deptemp(1,k)-termx
c          deptemp(2,k)=deptemp(2,k)-termy
c          deptemp(3,k)=deptemp(3,k)-termz
          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

          termx = uinp(1,k)*t1xxd+ uinp(2,k)*t1xyd + uinp(3,k)*t1xzd
          termy = uinp(1,k)*t1yxd+ uinp(2,k)*t1yyd + uinp(3,k)*t1yzd
          termz = uinp(1,k)*t1zxd+ uinp(2,k)*t1zyd + uinp(3,k)*t1zzd
          
c          deptemp(1,i)=deptemp(1,i)+termx
c          deptemp(2,i)=deptemp(2,i)+termy
c          deptemp(3,i)=deptemp(3,i)+termz
c
c          deptemp(1,k)=deptemp(1,k)-termx
c          deptemp(2,k)=deptemp(2,k)-termy
c          deptemp(3,k)=deptemp(3,k)-termz
          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz
          
c        end if
c       end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      do i=1,npole
        deptemp(1,i)=deptemp(1,i)+depi(1,i)+depk(1,i)
        deptemp(2,i)=deptemp(2,i)+depi(2,i)+depk(2,i)
        deptemp(3,i)=deptemp(3,i)+depi(3,i)+depk(3,i)
      end do     

c      do m=ntpair_start_rmndr(taskid),ntpair_last_rmndr(taskid)
       if(ntpair_start_rmndr(taskid).ne.0) then
         m=ntpair_start_rmndr(taskid)
         i=dEindex(1,m)
         k=dEindex(2,m)

          t2xxd=dEd2(1,m)
          t2xyd=dEd2(2,m)
          t2xzd=dEd2(3,m)
          t2yxd=dEd2(4,m)
          t2yyd=dEd2(5,m)
          t2yzd=dEd2(6,m)
          t2zxd=dEd2(7,m)
          t2zyd=dEd2(8,m)
          t2zzd=dEd2(9,m)
          t2xxp=dEp2(1,m)
          t2xyp=dEp2(2,m)
          t2xzp=dEp2(3,m)
          t2yxp=dEp2(4,m)
          t2yyp=dEp2(5,m)
          t2yzp=dEp2(6,m)
          t2zxp=dEp2(7,m)
          t2zyp=dEp2(8,m)
          t2zzp=dEp2(9,m)

          termx=uind(1,i)*t2xxp+ uind(2,i)*t2xyp + uind(3,i)*t2xzp
          termy=uind(1,i)*t2yxp+ uind(2,i)*t2yyp + uind(3,i)*t2yzp
          termz=uind(1,i)*t2zxp+ uind(2,i)*t2zyp + uind(3,i)*t2zzp

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz

          termx = uinp(1,i)*t2xxd+ uinp(2,i)*t2xyd + uinp(3,i)*t2xzd
          termy = uinp(1,i)*t2yxd+ uinp(2,i)*t2yyd + uinp(3,i)*t2yzd
          termz = uinp(1,i)*t2zxd+ uinp(2,i)*t2zyd + uinp(3,i)*t2zzd

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz

          t1xxd=dEd1(1,m)
          t1xyd=dEd1(2,m)
          t1xzd=dEd1(3,m)
          t1yxd=dEd1(4,m)
          t1yyd=dEd1(5,m)
          t1yzd=dEd1(6,m)
          t1zxd=dEd1(7,m)
          t1zyd=dEd1(8,m)
          t1zzd=dEd1(9,m)

          t1xxp=dEp1(1,m)
          t1xyp=dEp1(2,m)
          t1xzp=dEp1(3,m)
          t1yxp=dEp1(4,m)
          t1yyp=dEp1(5,m)
          t1yzp=dEp1(6,m)
          t1zxp=dEp1(7,m)
          t1zyp=dEp1(8,m)
          t1zzp=dEp1(9,m)

          termx=uind(1,k)*t1xxp+ uind(2,k)*t1xyp + uind(3,k)*t1xzp
          termy=uind(1,k)*t1yxp+ uind(2,k)*t1yyp + uind(3,k)*t1yzp
          termz=uind(1,k)*t1zxp+ uind(2,k)*t1zyp + uind(3,k)*t1zzp
          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz

          termx = uinp(1,k)*t1xxd+ uinp(2,k)*t1xyd + uinp(3,k)*t1xzd
          termy = uinp(1,k)*t1yxd+ uinp(2,k)*t1yyd + uinp(3,k)*t1yzd
          termz = uinp(1,k)*t1zxd+ uinp(2,k)*t1zyd + uinp(3,k)*t1zzd

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz
      end if

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
      deallocate(depi)
      deallocate(depk)
      return
      end

