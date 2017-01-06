      subroutine empole1c_3bPolar_totfield_dEtensor1bmut(
     &  ep3bt,dep3bt,virep3bt)
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
      integer npole3b,pnum(30)
      real*8 deptemp2(3,30),virtemp(3,3) 
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

      do i=1,3
       do j=1,3
         virt(j,i)=0.0d0
       end do 
      end do

      f = electric / dielec
      lenmol=3

!$OMP PARALLEL default(private) shared(start_emreal2,last_emreal2,
!$OMP& lenmol,depi,virt,uind,uinp)
!$OMP DO reduction(+:depi,virt) schedule(guided)
      do k=start_emreal2(taskid),last_emreal2(taskid),lenmol
         npole3b=3
         pnum(1)=k
         pnum(2)=k+1
         pnum(3)=k+2
         call induce0a_3b_PolelecOnly_totfieldnpole(npole3b,pnum,
     &           uindtemp,uinptemp)
         call empole1a1_3b_Polar_totfieldnpole_mod(npole3b,pnum,
     &        uindtemp,uinptemp,deptemp2,virtemp)
         do l1=1,npole3b
            i=pnum(l1)
            depi(1,i)=depi(1,i)+deptemp2(1,l1)
            depi(2,i)=depi(2,i)+deptemp2(2,l1)
            depi(3,i)=depi(3,i)+deptemp2(3,l1)
         end do
         do i=1,3
            do j=1,3
               virt(j,i)=virt(j,i)+virtemp(j,i)
            end do
         end do
         do l1=1,npole3b
            i=pnum(l1)
            uind(1,i)=uindtemp(1,l1)
            uind(2,i)=uindtemp(2,l1)
            uind(3,i)=uindtemp(3,l1)
            uinp(1,i)=uinptemp(1,l1)
            uinp(2,i)=uinptemp(2,l1)
            uinp(3,i)=uinptemp(3,l1)
         end do
      end do   
!$OMP END DO
!$OMP END PARLLEL

!$OMP PARALLEL default(private) shared(ntpair_dEtmp,
!$OMP& dEindextmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
!$OMP& depi,depk,uind,uinp,tau1indextmp,tau2indextmp,
!$OMP& taud1tmp,taud2tmp,taup1tmp,taup2tmp,ntpair_tau1tmp,
!$OMP& dep3bt,npole,x,y,z,xaxis,yaxis,zaxis,frcztau1d,
!$OMP& frcztau1p,frcxtau1d,frcxtau1p,frcytau1d,frcytau1p,
!$OMP& frcztau2d,frcztau2p,frcxtau2d,frcxtau2p,
!$OMP& frcytau2d,frcytau2p,virt)

!$OMP DO reduction(+:depi,depk,virt)
!$OMP& schedule(guided)
      do m=1,ntpair_dEtmp
         i=dEindextmp(1,m)
         k=dEindextmp(2,m)

         vxx=0.0d0
         vyx=0.0d0
         vzx=0.0d0
         vyy=0.0d0
         vzy=0.0d0
         vzz=0.0d0

           xr = x(k) - x(i)
           yr = y(k) - y(i)
           zr = z(k) - z(i)
           call image (xr,yr,zr)

               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = i
               if (iax .eq. 0)  iax = i
               if (iay .eq. 0)  iay = i
               if (kaz .eq. 0)  kaz = k
               if (kax .eq. 0)  kax = k
               if (kay .eq. 0)  kay = k
               xiz = x(iaz) - x(i)
               yiz = y(iaz) - y(i)
               ziz = z(iaz) - z(i)
               xix = x(iax) - x(i)
               yix = y(iax) - y(i)
               zix = z(iax) - z(i)
               xiy = x(iay) - x(i)
               yiy = y(iay) - y(i)
               ziy = z(iay) - z(i)
               xkz = x(kaz) - x(k)
               ykz = y(kaz) - y(k)
               zkz = z(kaz) - z(k)
               xkx = x(kax) - x(k)
               ykx = y(kax) - y(k)
               zkx = z(kax) - z(k)
               xky = x(kay) - x(k)
               yky = y(kay) - y(k)
               zky = z(kay) - z(k)
        
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

          frcz2xxd=frcztau2d(1,m)
          frcz2xyd=frcztau2d(2,m)
          frcz2xzd=frcztau2d(3,m)
          frcz2yxd=frcztau2d(4,m)
          frcz2yyd=frcztau2d(5,m)
          frcz2yzd=frcztau2d(6,m)
          frcz2zxd=frcztau2d(7,m)
          frcz2zyd=frcztau2d(8,m)
          frcz2zzd=frcztau2d(9,m)

          frcz2xxp=frcztau2p(1,m)
          frcz2xyp=frcztau2p(2,m)
          frcz2xzp=frcztau2p(3,m)
          frcz2yxp=frcztau2p(4,m)
          frcz2yyp=frcztau2p(5,m)
          frcz2yzp=frcztau2p(6,m)
          frcz2zxp=frcztau2p(7,m)
          frcz2zyp=frcztau2p(8,m)
          frcz2zzp=frcztau2p(9,m)

          frcy2xxd=frcytau2d(1,m)
          frcy2xyd=frcytau2d(2,m)
          frcy2xzd=frcytau2d(3,m)
          frcy2yxd=frcytau2d(4,m)
          frcy2yyd=frcytau2d(5,m)
          frcy2yzd=frcytau2d(6,m)
          frcy2zxd=frcytau2d(7,m)
          frcy2zyd=frcytau2d(8,m)
          frcy2zzd=frcytau2d(9,m)

          frcy2xxp=frcytau2p(1,m)
          frcy2xyp=frcytau2p(2,m)
          frcy2xzp=frcytau2p(3,m)
          frcy2yxp=frcytau2p(4,m)
          frcy2yyp=frcytau2p(5,m)
          frcy2yzp=frcytau2p(6,m)
          frcy2zxp=frcytau2p(7,m)
          frcy2zyp=frcytau2p(8,m)
          frcy2zzp=frcytau2p(9,m)

          frcx2xxd=frcxtau2d(1,m)
          frcx2xyd=frcxtau2d(2,m)
          frcx2xzd=frcxtau2d(3,m)
          frcx2yxd=frcxtau2d(4,m)
          frcx2yyd=frcxtau2d(5,m)
          frcx2yzd=frcxtau2d(6,m)
          frcx2zxd=frcxtau2d(7,m)
          frcx2zyd=frcxtau2d(8,m)
          frcx2zzd=frcxtau2d(9,m)

          frcx2xxp=frcxtau2p(1,m)
          frcx2xyp=frcxtau2p(2,m)
          frcx2xzp=frcxtau2p(3,m)
          frcx2yxp=frcxtau2p(4,m)
          frcx2yyp=frcxtau2p(5,m)
          frcx2yzp=frcxtau2p(6,m)
          frcx2zxp=frcxtau2p(7,m)
          frcx2zyp=frcxtau2p(8,m)
          frcx2zzp=frcxtau2p(9,m)

          termx=uind(1,i)*t2xxp+ uind(2,i)*t2xyp + uind(3,i)*t2xzp
          termy=uind(1,i)*t2yxp+ uind(2,i)*t2yyp + uind(3,i)*t2yzp
          termz=uind(1,i)*t2zxp+ uind(2,i)*t2zyp + uind(3,i)*t2zzp

       frczk(1)=uind(1,i)*frcz2xxp+uind(2,i)*frcz2xyp+uind(3,i)*frcz2xzp
       frczk(2)=uind(1,i)*frcz2yxp+uind(2,i)*frcz2yyp+uind(3,i)*frcz2yzp
       frczk(3)=uind(1,i)*frcz2zxp+uind(2,i)*frcz2zyp+uind(3,i)*frcz2zzp

       frcxk(1)=uind(1,i)*frcx2xxp+uind(2,i)*frcx2xyp+uind(3,i)*frcx2xzp
       frcxk(2)=uind(1,i)*frcx2yxp+uind(2,i)*frcx2yyp+uind(3,i)*frcx2yzp
       frcxk(3)=uind(1,i)*frcx2zxp+uind(2,i)*frcx2zyp+uind(3,i)*frcx2zzp

       frcyk(1)=uind(1,i)*frcy2xxp+uind(2,i)*frcy2xyp+uind(3,i)*frcy2xzp
       frcyk(2)=uind(1,i)*frcy2yxp+uind(2,i)*frcy2yyp+uind(3,i)*frcy2yzp
       frcyk(3)=uind(1,i)*frcy2zxp+uind(2,i)*frcy2zyp+uind(3,i)*frcy2zzp

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz

          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz

          vxx =vxx -xr*termx +xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
          vyx =vyx -yr*termx + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1) 
          vzx =vzx -zr*termx + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
          vyy =vyy -yr*termy + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
          vzy =vzy -zr*termy + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
          vzz =vzz -zr*termz + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

          termx = uinp(1,i)*t2xxd+ uinp(2,i)*t2xyd + uinp(3,i)*t2xzd
          termy = uinp(1,i)*t2yxd+ uinp(2,i)*t2yyd + uinp(3,i)*t2yzd
          termz = uinp(1,i)*t2zxd+ uinp(2,i)*t2zyd + uinp(3,i)*t2zzd

       frczk(1)=uinp(1,i)*frcz2xxd+uinp(2,i)*frcz2xyd+uinp(3,i)*frcz2xzd
       frczk(2)=uinp(1,i)*frcz2yxd+uinp(2,i)*frcz2yyd+uinp(3,i)*frcz2yzd
       frczk(3)=uinp(1,i)*frcz2zxd+uinp(2,i)*frcz2zyd+uinp(3,i)*frcz2zzd
          
       frcxk(1)=uinp(1,i)*frcx2xxd+uinp(2,i)*frcx2xyd+uinp(3,i)*frcx2xzd
       frcxk(2)=uinp(1,i)*frcx2yxd+uinp(2,i)*frcx2yyd+uinp(3,i)*frcx2yzd
       frcxk(3)=uinp(1,i)*frcx2zxd+uinp(2,i)*frcx2zyd+uinp(3,i)*frcx2zzd

       frcyk(1)=uinp(1,i)*frcy2xxd+uinp(2,i)*frcy2xyd+uinp(3,i)*frcy2xzd
       frcyk(2)=uinp(1,i)*frcy2yxd+uinp(2,i)*frcy2yyd+uinp(3,i)*frcy2yzd
       frcyk(3)=uinp(1,i)*frcy2zxd+uinp(2,i)*frcy2zyd+uinp(3,i)*frcy2zzd

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz
          vxx =vxx -xr*termx +xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
          vyx =vyx -yr*termx + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
          vzx =vzx -zr*termx + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
          vyy =vyy -yr*termy + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
          vzy =vzy -zr*termy + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
          vzz =vzz -zr*termz + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

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

          frcz1xxd=frcztau1d(1,m)
          frcz1xyd=frcztau1d(2,m)
          frcz1xzd=frcztau1d(3,m)
          frcz1yxd=frcztau1d(4,m)
          frcz1yyd=frcztau1d(5,m)
          frcz1yzd=frcztau1d(6,m)
          frcz1zxd=frcztau1d(7,m)
          frcz1zyd=frcztau1d(8,m)
          frcz1zzd=frcztau1d(9,m)

          frcz1xxp=frcztau1p(1,m)
          frcz1xyp=frcztau1p(2,m)
          frcz1xzp=frcztau1p(3,m)
          frcz1yxp=frcztau1p(4,m)
          frcz1yyp=frcztau1p(5,m)
          frcz1yzp=frcztau1p(6,m)
          frcz1zxp=frcztau1p(7,m)
          frcz1zyp=frcztau1p(8,m)
          frcz1zzp=frcztau1p(9,m)

          frcx1xxd=frcxtau1d(1,m)
          frcx1xyd=frcxtau1d(2,m)
          frcx1xzd=frcxtau1d(3,m)
          frcx1yxd=frcxtau1d(4,m)
          frcx1yyd=frcxtau1d(5,m)
          frcx1yzd=frcxtau1d(6,m)
          frcx1zxd=frcxtau1d(7,m)
          frcx1zyd=frcxtau1d(8,m)
          frcx1zzd=frcxtau1d(9,m)
          
          frcx1xxp=frcxtau1p(1,m)
          frcx1xyp=frcxtau1p(2,m)
          frcx1xzp=frcxtau1p(3,m)
          frcx1yxp=frcxtau1p(4,m)
          frcx1yyp=frcxtau1p(5,m)
          frcx1yzp=frcxtau1p(6,m)
          frcx1zxp=frcxtau1p(7,m)
          frcx1zyp=frcxtau1p(8,m)
          frcx1zzp=frcxtau1p(9,m)

          frcy1xxd=frcytau1d(1,m)
          frcy1xyd=frcytau1d(2,m)
          frcy1xzd=frcytau1d(3,m)
          frcy1yxd=frcytau1d(4,m)
          frcy1yyd=frcytau1d(5,m)
          frcy1yzd=frcytau1d(6,m)
          frcy1zxd=frcytau1d(7,m)
          frcy1zyd=frcytau1d(8,m)
          frcy1zzd=frcytau1d(9,m)
          
          frcy1xxp=frcytau1p(1,m)
          frcy1xyp=frcytau1p(2,m)
          frcy1xzp=frcytau1p(3,m)
          frcy1yxp=frcytau1p(4,m)
          frcy1yyp=frcytau1p(5,m)
          frcy1yzp=frcytau1p(6,m)
          frcy1zxp=frcytau1p(7,m)
          frcy1zyp=frcytau1p(8,m)
          frcy1zzp=frcytau1p(9,m)

          termx=uind(1,k)*t1xxp+ uind(2,k)*t1xyp + uind(3,k)*t1xzp
          termy=uind(1,k)*t1yxp+ uind(2,k)*t1yyp + uind(3,k)*t1yzp
          termz=uind(1,k)*t1zxp+ uind(2,k)*t1zyp + uind(3,k)*t1zzp 

       frczi(1)=uind(1,k)*frcz1xxp+uind(2,k)*frcz1xyp+uind(3,k)*frcz1xzp
       frczi(2)=uind(1,k)*frcz1yxp+uind(2,k)*frcz1yyp+uind(3,k)*frcz1yzp
       frczi(3)=uind(1,k)*frcz1zxp+uind(2,k)*frcz1zyp+uind(3,k)*frcz1zzp

       frcxi(1)=uind(1,k)*frcx1xxp+uind(2,k)*frcx1xyp+uind(3,k)*frcx1xzp
       frcxi(2)=uind(1,k)*frcx1yxp+uind(2,k)*frcx1yyp+uind(3,k)*frcx1yzp
       frcxi(3)=uind(1,k)*frcx1zxp+uind(2,k)*frcx1zyp+uind(3,k)*frcx1zzp

       frcyi(1)=uind(1,k)*frcy1xxp+uind(2,k)*frcy1xyp+uind(3,k)*frcy1xzp
       frcyi(2)=uind(1,k)*frcy1yxp+uind(2,k)*frcy1yyp+uind(3,k)*frcy1yzp
       frcyi(3)=uind(1,k)*frcy1zxp+uind(2,k)*frcy1zyp+uind(3,k)*frcy1zzp

         depi(1,i)=depi(1,i)+termx
         depi(2,i)=depi(2,i)+termy
         depi(3,i)=depi(3,i)+termz

          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz
       vxx=vxx-xr*termx+ xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) 
       vyx=vyx-yr*termx+ yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1)
       vzx=vzx-zr*termx + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1)
       vyy=vyy-yr*termy + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2)
       vzy=vzy-zr*termy + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2)
       vzz=vzz-zr*termz + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3)
          termx = uinp(1,k)*t1xxd+ uinp(2,k)*t1xyd + uinp(3,k)*t1xzd
          termy = uinp(1,k)*t1yxd+ uinp(2,k)*t1yyd + uinp(3,k)*t1yzd
          termz = uinp(1,k)*t1zxd+ uinp(2,k)*t1zyd + uinp(3,k)*t1zzd

       frczi(1)=uinp(1,k)*frcz1xxd+uinp(2,k)*frcz1xyd+uinp(3,k)*frcz1xzd
       frczi(2)=uinp(1,k)*frcz1yxd+uinp(2,k)*frcz1yyd+uinp(3,k)*frcz1yzd
       frczi(3)=uinp(1,k)*frcz1zxd+uinp(2,k)*frcz1zyd+uinp(3,k)*frcz1zzd
         
       frcxi(1)=uinp(1,k)*frcx1xxd+uinp(2,k)*frcx1xyd+uinp(3,k)*frcx1xzd
       frcxi(2)=uinp(1,k)*frcx1yxd+uinp(2,k)*frcx1yyd+uinp(3,k)*frcx1yzd
       frcxi(3)=uinp(1,k)*frcx1zxd+uinp(2,k)*frcx1zyd+uinp(3,k)*frcx1zzd
          
       frcyi(1)=uinp(1,k)*frcy1xxd+uinp(2,k)*frcy1xyd+uinp(3,k)*frcy1xzd
       frcyi(2)=uinp(1,k)*frcy1yxd+uinp(2,k)*frcy1yyd+uinp(3,k)*frcy1yzd
       frcyi(3)=uinp(1,k)*frcy1zxd+uinp(2,k)*frcy1zyd+uinp(3,k)*frcy1zzd

          depi(1,i)=depi(1,i)+termx
          depi(2,i)=depi(2,i)+termy
          depi(3,i)=depi(3,i)+termz
          
          depk(1,k)=depk(1,k)-termx
          depk(2,k)=depk(2,k)-termy
          depk(3,k)=depk(3,k)-termz
       vxx=vxx-xr*termx+ xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1)
       vyx=vyx-yr*termx+ yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1)
       vzx=vzx-zr*termx + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1)
       vyy=vyy-yr*termy + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2)
       vzy=vzy-zr*termy + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2)
       vzz=vzz-zr*termz + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3)
               virt(1,1) = virt(1,1) + vxx
               virt(2,1) = virt(2,1) + vyx
               virt(3,1) = virt(3,1) + vzx
               virt(1,2) = virt(1,2) + vyx
               virt(2,2) = virt(2,2) + vyy
               virt(3,2) = virt(3,2) + vzy
               virt(1,3) = virt(1,3) + vzx
               virt(2,3) = virt(2,3) + vzy
               virt(3,3) = virt(3,3) + vzz
                 
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

      do i=1,3
        do j=1,3
          virep3bt(j,i)=virep3bt(j,i)+virt(j,i)
        end do
      end do

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

