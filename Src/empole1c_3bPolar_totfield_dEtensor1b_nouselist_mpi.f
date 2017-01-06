      subroutine empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpi(
     &  dep3bt,virep3bt)
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
      real*8 eptemp,virep3bt(3,3)
      integer npole3b
      real*8 dep3bt(3,npole)
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
      real*8 tau1xxd,tau1xyd,tau1xzd,tau1yxd,tau1yyd,tau1yzd,tau1zxd
      real*8 tau1zyd,tau1zzd
      real*8 tau2xxd,tau2xyd,tau2xzd,tau2yxd,tau2yyd,tau2yzd,tau2zxd
      real*8 tau2zyd,tau2zzd
      real*8 tau1xxp,tau1xyp,tau1xzp,tau1yxp,tau1yyp,tau1yzp,tau1zxp
      real*8 tau1zyp,tau1zzp
      real*8 tau2xxp,tau2xyp,tau2xzp,tau2yxp,tau2yyp,tau2yzp,tau2zxp
      real*8 tau2zyp,tau2zzp
      integer m,k1,l2,length,start,last,moli1rmndr
      real*8 termx,termy,termz
      real*8 taux,tauy,tauz
      real*8, allocatable :: depi(:,:)
      real*8, allocatable :: depk(:,:)
      ewaldcut2=ewaldcut*ewaldcut

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

      f = electric / dielec
      

c!$OMP PARALLEL default(private) shared(ntpair_dEtmp,
c!$OMP& taskid,dEindextmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
c!$OMP& depi,depk,uind,uinp,taud1tmp,taud2tmp,taup1tmp,taup2tmp)
c!$OMP DO reduction(+:depi,depk)
c!$OMP& schedule(guided)
c      do m=ntpair_start(taskid),ntpair_last(taskid)

!$OMP PARALLEL default(private) shared(ntpair_dEtmp,
!$OMP& dEindextmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
!$OMP& depi,depk,uind,uinp,tau1indextmp,tau2indextmp,
!$OMP& taud1tmp,taud2tmp,taup1tmp,taup2tmp,ntpair_tau1tmp,
!$OMP& dep3bt,npole)

!$OMP DO reduction(+:depi,depk)
!$OMP& schedule(guided)
      do m=1,ntpair_dEtmp
         i=dEindextmp(1,m)
         k=dEindextmp(2,m)
        
          t2xxd=dEd2tmp(1,m)
          t2xyd=dEd2tmp(2,m)
          t2xzd=dEd2tmp(3,m)
          t2yxd=dEd2tmp(4,m)
          t2yyd=dEd2tmp(5,m)
          t2yzd=dEd2tmp(6,m)
          t2zxd=dEd2tmp(7,m)
          t2zyd=dEd2tmp(8,m)
          t2zzd=dEd2tmp(9,m)

          t2xxp=dEp2tmp(1,m)
          t2xyp=dEp2tmp(2,m)
          t2xzp=dEp2tmp(3,m)
          t2yxp=dEp2tmp(4,m)
          t2yyp=dEp2tmp(5,m)
          t2yzp=dEp2tmp(6,m)
          t2zxp=dEp2tmp(7,m)
          t2zyp=dEp2tmp(8,m)
          t2zzp=dEp2tmp(9,m)


          termx=uind(1,i)*t2xxp+ uind(2,i)*t2xyp + uind(3,i)*t2xzp
          termy=uind(1,i)*t2yxp+ uind(2,i)*t2yyp + uind(3,i)*t2yzp
          termz=uind(1,i)*t2zxp+ uind(2,i)*t2zyp + uind(3,i)*t2zzp

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz

          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

          termx = uinp(1,i)*t2xxd+ uinp(2,i)*t2xyd + uinp(3,i)*t2xzd
          termy = uinp(1,i)*t2yxd+ uinp(2,i)*t2yyd + uinp(3,i)*t2yzd
          termz = uinp(1,i)*t2zxd+ uinp(2,i)*t2zyd + uinp(3,i)*t2zzd

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

          t1xxd=dEd1tmp(1,m)
          t1xyd=dEd1tmp(2,m)
          t1xzd=dEd1tmp(3,m)
          t1yxd=dEd1tmp(4,m)
          t1yyd=dEd1tmp(5,m)
          t1yzd=dEd1tmp(6,m)
          t1zxd=dEd1tmp(7,m)
          t1zyd=dEd1tmp(8,m)
          t1zzd=dEd1tmp(9,m)

          t1xxp=dEp1tmp(1,m)
          t1xyp=dEp1tmp(2,m)
          t1xzp=dEp1tmp(3,m)
          t1yxp=dEp1tmp(4,m)
          t1yyp=dEp1tmp(5,m)
          t1yzp=dEp1tmp(6,m)
          t1zxp=dEp1tmp(7,m)
          t1zyp=dEp1tmp(8,m)
          t1zzp=dEp1tmp(9,m)

          termx=uind(1,k)*t1xxp+ uind(2,k)*t1xyp + uind(3,k)*t1xzp
          termy=uind(1,k)*t1yxp+ uind(2,k)*t1yyp + uind(3,k)*t1yzp
          termz=uind(1,k)*t1zxp+ uind(2,k)*t1zyp + uind(3,k)*t1zzp 

         depi(1,i)=depi(1,i)+termx
         depi(2,i)=depi(2,i)+termy
         depi(3,i)=depi(3,i)+termz

          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

          termx = uinp(1,k)*t1xxd+ uinp(2,k)*t1xyd + uinp(3,k)*t1xzd
          termy = uinp(1,k)*t1yxd+ uinp(2,k)*t1yyd + uinp(3,k)*t1yzd
          termz = uinp(1,k)*t1zxd+ uinp(2,k)*t1zyd + uinp(3,k)*t1zzd

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz
          
      end do
!$OMP END DO

!$OMP DO reduction(+:depi,depk)
!$OMP& schedule(guided)
      do m=1,ntpair_tau1tmp
         i=tau1indextmp(1,m)
         k=tau1indextmp(2,m)

          tau1xxd=taud1tmp(1,m)
          tau1xyd=taud1tmp(2,m)
          tau1xzd=taud1tmp(3,m)
          tau1yxd=taud1tmp(4,m)
          tau1yyd=taud1tmp(5,m)
          tau1yzd=taud1tmp(6,m)
          tau1zxd=taud1tmp(7,m)
          tau1zyd=taud1tmp(8,m)
          tau1zzd=taud1tmp(9,m)

          tau1xxp=taup1tmp(1,m)
          tau1xyp=taup1tmp(2,m)
          tau1xzp=taup1tmp(3,m)
          tau1yxp=taup1tmp(4,m)
          tau1yyp=taup1tmp(5,m)
          tau1yzp=taup1tmp(6,m)
          tau1zxp=taup1tmp(7,m)
          tau1zyp=taup1tmp(8,m)
          tau1zzp=taup1tmp(9,m)

          taux=uind(1,k)*tau1xxp+ uind(2,k)*tau1xyp + uind(3,k)*tau1xzp
          tauy=uind(1,k)*tau1yxp+ uind(2,k)*tau1yyp + uind(3,k)*tau1yzp
          tauz=uind(1,k)*tau1zxp+ uind(2,k)*tau1zyp + uind(3,k)*tau1zzp

          depi(1,i)=depi(1,i)+taux
          depi(2,i)=depi(2,i)+tauy
          depi(3,i)=depi(3,i)+tauz

          taux=uinp(1,k)*tau1xxd+ uinp(2,k)*tau1xyd +uinp(3,k)*tau1xzd
          tauy=uinp(1,k)*tau1yxd+ uinp(2,k)*tau1yyd +uinp(3,k)*tau1yzd
          tauz=uinp(1,k)*tau1zxd+ uinp(2,k)*tau1zyd +uinp(3,k)*tau1zzd

          depi(1,i)=depi(1,i)+taux
          depi(2,i)=depi(2,i)+tauy
          depi(3,i)=depi(3,i)+tauz
c      end do

c      do m=1,ntpair_tau2tmp
         i=tau2indextmp(1,m)
         k=tau2indextmp(2,m)

          tau2xxd=taud2tmp(1,m)
          tau2xyd=taud2tmp(2,m)
          tau2xzd=taud2tmp(3,m)
          tau2yxd=taud2tmp(4,m)
          tau2yyd=taud2tmp(5,m)
          tau2yzd=taud2tmp(6,m)
          tau2zxd=taud2tmp(7,m)
          tau2zyd=taud2tmp(8,m)
          tau2zzd=taud2tmp(9,m)

          tau2xxp=taup2tmp(1,m)
          tau2xyp=taup2tmp(2,m)
          tau2xzp=taup2tmp(3,m)
          tau2yxp=taup2tmp(4,m)
          tau2yyp=taup2tmp(5,m)
          tau2yzp=taup2tmp(6,m)
          tau2zxp=taup2tmp(7,m)
          tau2zyp=taup2tmp(8,m)
          tau2zzp=taup2tmp(9,m)

          taux=uind(1,i)*tau2xxp+ uind(2,i)*tau2xyp + uind(3,i)*tau2xzp
          tauy=uind(1,i)*tau2yxp+ uind(2,i)*tau2yyp + uind(3,i)*tau2yzp
          tauz=uind(1,i)*tau2zxp+ uind(2,i)*tau2zyp + uind(3,i)*tau2zzp

          depk(1,k)=depk(1,k)+taux
          depk(2,k)=depk(2,k)+tauy
          depk(3,k)=depk(3,k)+tauz

          taux=uinp(1,i)*tau2xxd+ uinp(2,i)*tau2xyd +uinp(3,i)*tau2xzd
          tauy=uinp(1,i)*tau2yxd+ uinp(2,i)*tau2yyd +uinp(3,i)*tau2yzd
          tauz=uinp(1,i)*tau2zxd+ uinp(2,i)*tau2zyd +uinp(3,i)*tau2zzd

          depk(1,k)=depk(1,k)+taux
          depk(2,k)=depk(2,k)+tauy
          depk(3,k)=depk(3,k)+tauz
      end do
!$OMP END DO

!$OMP DO
      do i=1,npole
        dep3bt(1,i)=dep3bt(1,i)+depi(1,i)+depk(1,i)
        dep3bt(2,i)=dep3bt(2,i)+depi(2,i)+depk(2,i)
        dep3bt(3,i)=dep3bt(3,i)+depi(3,i)+depk(3,i)
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
      deallocate(depi)
      deallocate(depk)
      return
      end

