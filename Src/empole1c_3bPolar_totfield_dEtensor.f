      subroutine empole1c_3b_Polar_totfield_dEtensor(npole3b,pnum,
     &  eptemp,deptemp,virtemp)
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'boxes.i'
c      include 'chgpot.i'
c      include 'ewald.i'
c      include 'inter.i'
c      include 'math.i'
c      include 'mpole.i'
c      include 'polar2.i'
c      include 'polpot.i'
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
c      use pme, only: nfft1, nfft2, nfft3
      use pme
      use chunks
      use fft
      use totfield
      use limits 
      use dEtensor
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
      real*8 eptemp,deptemp2(3,npole3b),virtemp(3,3)
      integer npole3b,pnum(*)
      real*8 deptemp(3,npole)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
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
      integer m,k1,l2
      real*8 termx,termy,termz
      !if (.not. allocated(qgrid3b)) 
     & !          allocate (qgrid3b(2,nfft1,nfft2,nfft3))
      
      !if (.not. allocated(iprime3b)) allocate (iprime3b(maxprime,3))
      !if (.not. allocated(ffttable3b)) allocate (ffttable3b(maxtable,3))
      
      !call fftsetup3b(qgrid3b,planf3b,planb3b,iprime3b,ffttable3b)

c   Zero out temporary gradient of polarization energy
      ewaldcut2=ewaldcut*ewaldcut
      eptemp = 0.0d0

      do l1 = 1, npole3b
        i=pnum(l1)
        do j = 1, 3
          deptemp2(j,l1) = 0.0d0
          uind(j,l1) = 0.0d0
          uinp(j,l1) = 0.0d0
        end do
      end do

      do i=1,npole
         do j=1,3
           deptemp(j,i) = 0.0d0           
         end do 
      end do

      do i=1,3
         do j=1,3
           virtemp(j,i)=0.0d0
         end do
      end do

      f = electric / dielec
            
      !print*,"DETensor Empole1c_3b_totfieldnnpole"

      !!!!!call induce0c_3b_new(npole3b,pnum,uind,uinp)
       !call induce0a_3b(npole3b,pnum,uind,uinp,thetai1_3b,thetai2_3b,
     &  !              thetai3_3b,qgrid3b,qfac3b,pmetable3b,igrid3b,
     &   !              planf3b,planb3b,iprime3b,ffttable3b)   
       call induce0a_3b_totfield(npole3b,pnum,uind,uinp,qgrid3b,
     &                 planf3b,planb3b,iprime3b,ffttable3b)
    
      if(embedtyp.eq.'I') then

      else if (embedtyp.eq.'II') then

      !call emrecip1_1and2_3b_Polar_totfield_parl_513(npole3b,pnum,
     & !  uind,uinp,eptemp,deptemp,virtemp,qgrid3b,
     & !  planf3b,planb3b,
     & !  iprime3b,ffttable3b)

      else
      end if
ccc      call induce0a_3b(npole3b,pnum,uind,uinp,qgrid,qfac)


      !deallocate(qgrid3b)
      !deallocate(iprime3b)
      !deallocate(ffttable3b)

      
      if(embedtyp.eq.'I') then
      else if(embedtyp.eq.'II') then
       if(poltyp.eq.'MUTUAL') then
       call ereal1c1_3b_totfieldnpole(npole3b,
     &  pnum,uind,uinp,deptemp2,virtemp)
       end if
       !call ereal1d2_3b_totfieldnpole(npole3b,
     & !  pnum,uind,uinp,eptemp,deptemp,virtemp,nembedlst,embedlst)
       !call ereal1c2_3b_totfieldnpole_omp(npole3b,pnum,uind,uinp,
     & ! eptemp,deptemp,virtemp)
       do l1=1,npole3b
          i=pnum(l1)
          eptemp = eptemp -0.5d0*f*( uind(1,l1)*fieldpnpole(1,i) +
     &    uind(2,l1)*fieldpnpole(2,i) + uind(3,l1)*fieldpnpole(3,i))
       end do

       do m=1,ntpair_dE
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

        !if(doi.and.(.not.dok)) then
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

c          deptemp(1,k)=deptemp(1,k)+termx
c          deptemp(2,k)=deptemp(2,k)+termy
c          deptemp(3,k)=deptemp(3,k)+termz

          termx = uinp(1,l1)*t2xxd+ uinp(2,l1)*t2xyd + uinp(3,l1)*t2xzd
          termy = uinp(1,l1)*t2yxd+ uinp(2,l1)*t2yyd + uinp(3,l1)*t2yzd
          termz = uinp(1,l1)*t2zxd+ uinp(2,l1)*t2zyd + uinp(3,l1)*t2zzd

          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz

c          deptemp(1,k)=deptemp(1,k)+termx
c          deptemp(2,k)=deptemp(2,k)+termy
c          deptemp(3,k)=deptemp(3,k)+termz

        end if
c        else if(dok.and.(.not.doi)) then
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

c          deptemp(1,i)=deptemp(1,i)+termx
c          deptemp(2,i)=deptemp(2,i)+termy
c          deptemp(3,i)=deptemp(3,i)+termz

          termx = uinp(1,l2)*t1xxd+ uinp(2,l2)*t1xyd + uinp(3,l2)*t1xzd
          termy = uinp(1,l2)*t1yxd+ uinp(2,l2)*t1yyd + uinp(3,l2)*t1yzd
          termz = uinp(1,l2)*t1zxd+ uinp(2,l2)*t1zyd + uinp(3,l2)*t1zzd
          
          deptemp(1,i)=deptemp(1,i)+termx
          deptemp(2,i)=deptemp(2,i)+termy
          deptemp(3,i)=deptemp(3,i)+termz

          deptemp(1,k)=deptemp(1,k)-termx
          deptemp(2,k)=deptemp(2,k)-termy
          deptemp(3,k)=deptemp(3,k)-termz
          
c          deptemp(1,i)=deptemp(1,i)+termx
c          deptemp(2,i)=deptemp(2,i)+termy
c          deptemp(3,i)=deptemp(3,i)+termz
        end if

c        else if(dok.and.doi) then
c          t1xxd=dEd1(1,m)
c          t1xyd=dEd1(2,m)
c          t1xzd=dEd1(3,m)
c          t1yxd=dEd1(4,m)
c          t1yyd=dEd1(5,m)
c          t1yzd=dEd1(6,m)
c          t1zxd=dEd1(7,m)
c          t1zyd=dEd1(8,m)
c          t1zzd=dEd1(9,m)
c      
c          t1xxp=dEp1(1,m)
c          t1xyp=dEp1(2,m)
c          t1xzp=dEp1(3,m)
c          t1yxp=dEp1(4,m)
c          t1yyp=dEp1(5,m)
c          t1yzp=dEp1(6,m)
c          t1zxp=dEp1(7,m)
c          t1zyp=dEp1(8,m)
c          t1zzp=dEp1(9,m)
c
c          t2xxd=dEd2(1,m)
c          t2xyd=dEd2(2,m)
c          t2xzd=dEd2(3,m)
c          t2yxd=dEd2(4,m)
c          t2yyd=dEd2(5,m)
c          t2yzd=dEd2(6,m)
c          t2zxd=dEd2(7,m)
c          t2zyd=dEd2(8,m)
c          t2zzd=dEd2(9,m)

c          t2xxp=dEp2(1,m)
c          t2xyp=dEp2(2,m)
c          t2xzp=dEp2(3,m)
c          t2yxp=dEp2(4,m)
c          t2yyp=dEp2(5,m)
c          t2yzp=dEp2(6,m)
c          t2zxp=dEp2(7,m)
c          t2zyp=dEp2(8,m)
c          t2zzp=dEp2(9,m)

c          termx=uind(1,l1)*t2xxp+ uind(2,l1)*t2xyp + uind(3,l1)*t2xzp
c          termy=uind(1,l1)*t2yxp+ uind(2,l1)*t2yyp + uind(3,l1)*t2yzp
c          termz=uind(1,l1)*t2zxp+ uind(2,l1)*t2zyp + uind(3,l1)*t2zzp

c          deptemp(1,i)=deptemp(1,i)+0.5d0*termx
c          deptemp(2,i)=deptemp(2,i)+0.5d0*termy
c          deptemp(3,i)=deptemp(3,i)+0.5d0*termz
c
c          deptemp(1,k)=deptemp(1,k)-0.5d0*termx
c          deptemp(2,k)=deptemp(2,k)-0.5d0*termy
c          deptemp(3,k)=deptemp(3,k)-0.5d0*termz
c          deptemp(1,k)=deptemp(1,k)+0.5*termx
c          deptemp(2,k)=deptemp(2,k)+0.5*termy
c          deptemp(3,k)=deptemp(3,k)+0.5*termz

c          termx=uinp(1,l1)*t2xxd+ uinp(2,l1)*t2xyd + uinp(3,l1)*t2xzd
c          termy=uinp(1,l1)*t2yxd+ uinp(2,l1)*t2yyd + uinp(3,l1)*t2yzd
c          termz=uinp(1,l1)*t2zxd+ uinp(2,l1)*t2zyd + uinp(3,l1)*t2zzd
         
c          deptemp(1,i)=deptemp(1,i)+0.5d0*termx
c          deptemp(2,i)=deptemp(2,i)+0.5d0*termy
c          deptemp(3,i)=deptemp(3,i)+0.5d0*termz
c
c          deptemp(1,k)=deptemp(1,k)-0.5d0*termx
c          deptemp(2,k)=deptemp(2,k)-0.5d0*termy
c          deptemp(3,k)=deptemp(3,k)-0.5d0*termz
c          deptemp(1,k)=deptemp(1,k)+0.5*termx
c          deptemp(2,k)=deptemp(2,k)+0.5*termy
c          deptemp(3,k)=deptemp(3,k)+0.5*termz

c          termx=uind(1,l2)*t1xxp+ uind(2,l2)*t1xyp + uind(3,l2)*t1xzp
c          termy=uind(1,l2)*t1yxp+ uind(2,l2)*t1yyp + uind(3,l2)*t1yzp
c          termz=uind(1,l2)*t1zxp+ uind(2,l2)*t1zyp + uind(3,l2)*t1zzp
c          deptemp(1,i)=deptemp(1,i)+0.5d0*termx
c          deptemp(2,i)=deptemp(2,i)+0.5d0*termy
c          deptemp(3,i)=deptemp(3,i)+0.5d0*termz
c
c          deptemp(1,k)=deptemp(1,k)-0.5d0*termx
c          deptemp(2,k)=deptemp(2,k)-0.5d0*termy
c          deptemp(3,k)=deptemp(3,k)-0.5d0*termz

c          deptemp(1,i)=deptemp(1,i)+0.5*termx
c          deptemp(2,i)=deptemp(2,i)+0.5*termy
c          deptemp(3,i)=deptemp(3,i)+0.5*termz
c          termx=uinp(1,l2)*t1xxd + uinp(2,l2)*t1xyd + uinp(3,l2)*t1xzd
c          termy=uinp(1,l2)*t1yxd + uinp(2,l2)*t1yyd + uinp(3,l2)*t1yzd
c          termz=uinp(1,l2)*t1zxd + uinp(2,l2)*t1zyd + uinp(3,l2)*t1zzd
c          deptemp(1,i)=deptemp(1,i)+0.5d0*termx
c          deptemp(2,i)=deptemp(2,i)+0.5d0*termy
c          deptemp(3,i)=deptemp(3,i)+0.5d0*termz
c
c          deptemp(1,k)=deptemp(1,k)-0.5d0*termx
c          deptemp(2,k)=deptemp(2,k)-0.5d0*termy
c          deptemp(3,k)=deptemp(3,k)-0.5d0*termz


c          deptemp(1,i)=deptemp(1,i)+0.5*termx
c          deptemp(2,i)=deptemp(2,i)+0.5*termy
c          deptemp(3,i)=deptemp(3,i)+0.5*termz

c        end if


       end do

      else

      end if

      do l1 =1,npole3b
          i =pnum(l1)
          deptemp(1,i)=deptemp2(1,l1)+deptemp(1,i)
          deptemp(2,i)=deptemp2(2,l1)+deptemp(2,i)
          deptemp(3,i)=deptemp2(3,l1)+deptemp(3,i)
      end do

      
      !aewald3b=aewald
      term = 2.0d0 * aewald3b * aewald3b
      fterm = -f * aewald3b / sqrtpi
      do l1 = 1, npole3b
         i=pnum(l1)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,l1)
         uiy = uind(2,l1)
         uiz = uind(3,l1)
         uii = dix*uix + diy*uiy + diz*uiz
         ei = fterm * term * uii / 3.0d0
         eptemp = eptemp + ei
      end do

c
c     compute the self-energy torque term due to induced dipole
c

      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald3b**3 / sqrtpi
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
cccc         call torque_3b (deptemp,i,trq,trqi,frcx,frcy,frcz)
c         call torque_3b_new(npole3b,pnum,deptemp,i,
c     &    trq,trqi,frcx,frcy,frcz)
         call torque_3b_new_npole(deptemp,i,
     &      trq,trqi,frcx,frcy,frcz)
      end do


      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         xup = 0.0d0
         yup = 0.0d0
         zup = 0.0d0
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
c            xu = xu + uind(1,i)
c            yu = yu + uind(2,i)
c            zu = zu + uind(3,i)
c            xup = xup + uinp(1,i)
c            yup = yup + uinp(2,i)
c            zup = zup + uinp(3,i)
            xu = xu + uind(1,l1)
            yu = yu + uind(2,l1)
            zu = zu + uind(3,l1)
            xup = xup + uinp(1,l1)
            yup = yup + uinp(2,l1)
            zup = zup + uinp(3,l1)

         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         eptemp=eptemp + term*(xd*xu+yd*yu+zd*zu)
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
c            dep(1,ii) = dep(1,ii) + term*rpole(1,i)*(xu+xup)
c            dep(2,ii) = dep(2,ii) + term*rpole(1,i)*(yu+yup)
c            dep(3,ii) = dep(3,ii) + term*rpole(1,i)*(zu+zup)
            deptemp(1,l1) = deptemp(1,l1) + term*rpole(1,i)*(xu+xup)
            deptemp(2,l1) = deptemp(2,l1) + term*rpole(1,i)*(yu+yup)
            deptemp(3,l1) = deptemp(3,l1) + term*rpole(1,i)*(zu+zup)

         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do l1 = 1, npole3b
            i=pnum(l1)
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            trqi(1) = rpole(3,i)*zufield - rpole(4,i)*yufield
            trqi(2) = rpole(4,i)*xufield - rpole(2,i)*zufield
            trqi(3) = rpole(2,i)*yufield - rpole(3,i)*xufield
          !  call torque_3b (pnum2,deptemp,i,trq,trqi,frcx,frcy,frcz)
            call torque_3b_new(npole3b,pnum,deptemp,i,
     &    trq,trqi,frcx,frcy,frcz) 
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xq * (xd+0.5d0*(xu+xup))
         yv = yq * (yd+0.5d0*(yu+yup))
         zv = zq * (zd+0.5d0*(zu+zup))
         vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
         virtemp(1,1)=virtemp(1,1)+2.0d0*term*(xq*xq+xv) + vterm
         virtemp(2,1)=virtemp(2,1)+2.0d0*term*(xq*yq+xv)
         virtemp(3,1)=virtemp(3,1)+2.0d0*term*(xq*zq+xv)
         virtemp(1,2)=virtemp(1,2) + 2.0d0*term*(yq*xq+yv)
         virtemp(2,2)=virtemp(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         virtemp(3,2)=virtemp(3,2) + 2.0d0*term*(yq*zq+yv)
         virtemp(1,3)=virtemp(1,3) + 2.0d0*term*(zq*xq+zv)
         virtemp(2,3)=virtemp(2,3) + 2.0d0*term*(zq*yq+zv)
         virtemp(3,3)=virtemp(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
c         if (poltyp .eq. 'DIRECT') then
c            vterm = term * (xu*xup+yu*yup+zu*zup)
c            virtemp(1,1) = virtemp(1,1) + vterm
c            virtemp(2,2) = virtemp(2,2) + vterm
c            virtemp(3,3) = virtemp(3,3) + vterm
c         end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
c      einter = einter + em + eptemp - eintra

c      print*,"End of empole1c_3b_Polar ep=",eptemp,moli1,moli2,moli3
c      do l1=1,npole3b
c         print*,"End of empole1c_3b_Polar l1 Pnum(l1)",l1,pnum(l1)
c      end do
      return
      end

