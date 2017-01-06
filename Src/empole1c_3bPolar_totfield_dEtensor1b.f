      subroutine empole1c_3b_Polar_totfield_dEtensor1b(start,last,
     &  moli1rmndr,eptemp,deptemp,virtemp)
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
      real*8, allocatable :: uind(:,:)
      real*8, allocatable :: uinp(:,:)
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

      length=last-start+1
      allocate(pnum(3*length))
      allocate(uind(3,3*length))
      allocate(uinp(3,3*length))
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
      do l1=0,length-1
         pnum(3*l1+1)=imol(1,start+l1)
         pnum(3*l1+2)=imol(1,start+l1)+1
         pnum(3*l1+3)=imol(2,start+l1)
      end do      
      npole3b=3*length

      do l1 = 1, npole3b
        i=pnum(l1)
        do j = 1, 3
          uind(j,l1) = polarity(i) * fieldnpole(j,i) 
          uinp(j,l1) = polarity(i) * fieldpnpole(j,i) 
        end do
      end do

      f = electric / dielec
      
       do l1=1,npole3b
          i=pnum(l1)
          eptemp = eptemp -0.5d0*f*( uind(1,l1)*fieldpnpole(1,i) +
     &    uind(2,l1)*fieldpnpole(2,i) + uind(3,l1)*fieldpnpole(3,i))
       end do

!$OMP PARALLEL default(private) shared(ntpair_start,ntpair_last,
!$OMP& taskid,dEindex,dEd1,dEd2,dEp1,dEp2,pnum,npole3b,x,y,z,
!$OMP& ewaldcut2,depi,depk,uind,uinp)
!$OMP DO reduction(+:depi,depk)
!$OMP& schedule(guided)
      do m=ntpair_start(taskid),ntpair_last(taskid)
         i=dEindex(1,m)
         k=dEindex(2,m)
         doi=.false.
         dok=.false.
            do k1=1,npole3b
               if(pnum(k1).eq.k) then
                 dok=.true.
                 l2=k1
               else if(pnum(k1).eq.i) then
                 doi=.true.
                 l1=k1
               end if
            end do
          xr = x(k) - x(i)
          yr = y(k) - y(i)
          zr = z(k) - z(i)
          call image (xr,yr,zr)
          r2 = xr*xr + yr*yr + zr*zr

       if (r2 .le. ewaldcut2) then
        
        if(doi) then
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

          termx=uind(1,l1)*t2xxp+ uind(2,l1)*t2xyp + uind(3,l1)*t2xzp
          termy=uind(1,l1)*t2yxp+ uind(2,l1)*t2yyp + uind(3,l1)*t2yzp
          termz=uind(1,l1)*t2zxp+ uind(2,l1)*t2zyp + uind(3,l1)*t2zzp

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

          termx = uinp(1,l1)*t2xxd+ uinp(2,l1)*t2xyd + uinp(3,l1)*t2xzd
          termy = uinp(1,l1)*t2yxd+ uinp(2,l1)*t2yyd + uinp(3,l1)*t2yzd
          termz = uinp(1,l1)*t2zxd+ uinp(2,l1)*t2zyd + uinp(3,l1)*t2zzd

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

        end if

        if(dok) then
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

          termx=uind(1,l2)*t1xxp+ uind(2,l2)*t1xyp + uind(3,l2)*t1xzp
          termy=uind(1,l2)*t1yxp+ uind(2,l2)*t1yyp + uind(3,l2)*t1yzp
          termz=uind(1,l2)*t1zxp+ uind(2,l2)*t1zyp + uind(3,l2)*t1zzp 
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

          termx = uinp(1,l2)*t1xxd+ uinp(2,l2)*t1xyd + uinp(3,l2)*t1xzd
          termy = uinp(1,l2)*t1yxd+ uinp(2,l2)*t1yyd + uinp(3,l2)*t1yzd
          termz = uinp(1,l2)*t1zxd+ uinp(2,l2)*t1zyd + uinp(3,l2)*t1zzd
          
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
          
        end if
       end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      do i=1,npole
        deptemp(1,i)=deptemp(1,i)+depi(1,i)+depk(1,i)
        deptemp(2,i)=deptemp(2,i)+depi(2,i)+depk(2,i)
        deptemp(3,i)=deptemp(3,i)+depi(3,i)+depk(3,i)
      end do     
      if(moli1rmndr.ne.0) then
        pnum(1)=imol(1,moli1rmndr)
        pnum(2)=imol(1,moli1rmndr)+1
        pnum(3)=imol(2,moli1rmndr)
       do l1=1,npole3b
         do j=1,3
            uind(j,l1)= 0.0d0
            uinp(j,l1)= 0.0d0
         end do
       end do

      npole3b=3

      do l1 = 1, npole3b
        i=pnum(l1)
        do j = 1, 3
          uind(j,l1) = polarity(i) * fieldnpole(j,i)
          uinp(j,l1) = polarity(i) * fieldpnpole(j,i)
        end do
      end do

       do l1=1,npole3b
          i=pnum(l1)
          eptemp = eptemp -0.5d0*f*( uind(1,l1)*fieldpnpole(1,i) +
     &    uind(2,l1)*fieldpnpole(2,i) + uind(3,l1)*fieldpnpole(3,i))
       end do

      do m=ntpair_start_rmndr(taskid),ntpair_last_rmndr(taskid)
         i=dEindex(1,m)
         k=dEindex(2,m)
         doi=.false.
         dok=.false.
            do k1=1,npole3b
               if(pnum(k1).eq.k) then
                 dok=.true.
                 l2=k1
               else if(pnum(k1).eq.i) then
                 doi=.true.
                 l1=k1
               end if
            end do

          xr = x(k) - x(i)
          yr = y(k) - y(i)
          zr = z(k) - z(i)
          call image (xr,yr,zr)
          r2 = xr*xr + yr*yr + zr*zr

       if (r2 .le. ewaldcut2) then

        if(doi) then
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

          termx=uind(1,l1)*t2xxp+ uind(2,l1)*t2xyp + uind(3,l1)*t2xzp
          termy=uind(1,l1)*t2yxp+ uind(2,l1)*t2yyp + uind(3,l1)*t2yzp
          termz=uind(1,l1)*t2zxp+ uind(2,l1)*t2zyp + uind(3,l1)*t2zzp

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz

          termx = uinp(1,l1)*t2xxd+ uinp(2,l1)*t2xyd + uinp(3,l1)*t2xzd
          termy = uinp(1,l1)*t2yxd+ uinp(2,l1)*t2yyd + uinp(3,l1)*t2yzd
          termz = uinp(1,l1)*t2zxd+ uinp(2,l1)*t2zyd + uinp(3,l1)*t2zzd

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz
        end if

        if(dok) then
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

          termx=uind(1,l2)*t1xxp+ uind(2,l2)*t1xyp + uind(3,l2)*t1xzp
          termy=uind(1,l2)*t1yxp+ uind(2,l2)*t1yyp + uind(3,l2)*t1yzp
          termz=uind(1,l2)*t1zxp+ uind(2,l2)*t1zyp + uind(3,l2)*t1zzp
          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz

          termx = uinp(1,l2)*t1xxd+ uinp(2,l2)*t1xyd + uinp(3,l2)*t1xzd
          termy = uinp(1,l2)*t1yxd+ uinp(2,l2)*t1yyd + uinp(3,l2)*t1yzd
          termz = uinp(1,l2)*t1zxd+ uinp(2,l2)*t1zyd + uinp(3,l2)*t1zzd

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz
        end if
       end if
      end do

      end if

      !aewald3b=aewald
c      term = 2.0d0 * aewald3b * aewald3b
c      fterm = -f * aewald3b / sqrtpi
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         dix = rpole(2,i)
c         diy = rpole(3,i)
c         diz = rpole(4,i)
c         uix = uind(1,l1)
c         uiy = uind(2,l1)
c         uiz = uind(3,l1)
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
c         uix = 0.5d0 * (uind(1,l1)+uinp(1,l1))
c         uiy = 0.5d0 * (uind(2,l1)+uinp(2,l1))
c         uiz = 0.5d0 * (uind(3,l1)+uinp(3,l1))
c         trqi(1) = term * (diy*uiz-diz*uiy)
c         trqi(2) = term * (diz*uix-dix*uiz)
c         trqi(3) = term * (dix*uiy-diy*uix)
c         call torque_3b_new_npole(deptemp,i,
c     &      trq,trqi,frcx,frcy,frcz)
c      end do
      deallocate(pnum)
      deallocate(uind)
      deallocate(uinp)
      deallocate(depi)
      deallocate(depk)
      return
      end

