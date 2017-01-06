c      subroutine empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
c     &  dep3bt,virep3bt)
      subroutine Newempole1c_3b_Polar_totfield_dEtensor2(npole3b,pnum,
     &  eptemp,deptemp,virtemp)
      use bound
      use sizes
      use sizes2
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
      use limits2 
      use dEtensor
      use mpidat
      use molcul
      use dEtensor2
      use neigh2
      use dEtensor3
      implicit none 
      integer i,j,ii,l1,l3,k2,j1,kk,m1,m2,k3
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
      integer npole3b,pnum(*)
      real*8 eptemp,deptemp2(3,npole3b),virtemp(3,3)
      real*8 deptemp(3,npole)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
c      real*8 eptemp,virep3bt(3,3)
c      integer npole3b
c      real*8 dep3bt(3,npole)
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
      integer, allocatable :: oldlistinds(:,:)
      real*8 r2,xr,yr,zr,ewaldcut2
      logical doi,dok,in_npole3b
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
      integer k1,l2,length,start,last,moli1rmndr
      real*8 termx,termy,termz
      real*8 taux,tauy,tauz
      real*8 m
c      real*8, allocatable :: depi(:,:)
c      real*8, allocatable :: depk(:,:)
c      real*8 virt(3,3)
      integer iaz,iax,iay,kaz,kax,kay,kkk 
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      
      ewaldcut2=ewaldcutshort*ewaldcutshort
      

c      allocate(depi(3,npole))
c      allocate(depk(3,npole))
c
c      do i=1,npole
c         depi(1,i)=0.0d0
c         depi(2,i)=0.0d0
c         depi(3,i)=0.0d0
c         depk(1,i)=0.0d0
c         depk(2,i)=0.0d0
c         depk(3,i)=0.0d0      
c      end do

c      do i=1,3
c       do j=1,3
c         virt(j,i)=0.0d0
c       end do 
c      end do

      allocate (nembedlst(npole3b))
      allocate (embedlst(maxelst2,npole3b))
      allocate (oldlistinds(maxelst2,npole3b))

      eptemp = 0.0d0

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

      do l1 = 1, npole3b
        i=pnum(l1)
        nembedlst(l1) = 0
        do j = 1, 3
          deptemp2(j,l1) = 0.0d0
          uind(j,l1) = 0.0d0
          uinp(j,l1) = 0.0d0
        end do
        do j=1,nelst2(i)
          nembedlst(l1)=nembedlst(l1)+1
          embedlst(nembedlst(l1),l1)=elst2(j,i)
          oldlistinds(nembedlst(l1),l1)=j
        end do
        do j=1,i-1
           in_npole3b=.false.
           do k1=1,npole3b
              k2=pnum(k1)
              if(k2.eq.j) then
                in_npole3b=.true.
              goto 21
              end if
           end do
   21    continue
           if(in_npole3b.eqv..false.) then
             do k=1,nelst2(j)
              if(elst2(k,j).eq.i) then
                 nembedlst(l1)=nembedlst(l1)+1
                 embedlst(nembedlst(l1),l1)=j
                 oldlistinds(nembedlst(l1),l1)=k
              end if
             end do
           end if
        end do
      end do

      f = electric / dielec
      !print*,"Newempole1c*fewertens!!"
 
c       call induce0a_3b_totfield(npole3b,pnum,uind,uinp,qgrid3b,
c     &                 planf3b,planb3b,iprime3b,ffttable3b)

c       call ereal1c1_3b_totfieldnpole(npole3b,
c     &  pnum,uind,uinp,deptemp2,virtemp)
      call induce0a_3b_PolelecOnly_totfieldnpole(npole3b,pnum,uind,
     &     uinp)

      call empole1a1_3b_Polar_totfieldnpole_mod(npole3b,pnum,uind,
     &     uinp,deptemp2,virtemp)

       do l1=1,npole3b
          i=pnum(l1)
          eptemp = eptemp -0.5d0*f*( uind(1,l1)*fieldpnpole(1,i) +
     &    uind(2,l1)*fieldpnpole(2,i) + uind(3,l1)*fieldpnpole(3,i))
       end do

c!$OMP PARALLEL default(private) shared(ntpair_dEtmp,
c!$OMP& dEindextmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
c!$OMP& depi,depk,uind,uinp,tau1indextmp,tau2indextmp,
c!$OMP& taud1tmp,taud2tmp,taup1tmp,taup2tmp,ntpair_tau1tmp,
c!$OMP& dep3bt,npole,x,y,z,xaxis,yaxis,zaxis,frcztau1d,
c!$OMP& frcztau1p,frcxtau1d,frcxtau1p,frcytau1d,frcytau1p,
c!$OMP& frcztau2d,frcztau2p,frcxtau2d,frcxtau2p,
c!$OMP& frcytau2d,frcytau2p,virt)
c
c!$OMP DO reduction(+:depi,depk,virt)
c!$OMP& schedule(guided)
c      do m=1,ntpair_dEtmp
c         i=dEindextmp(1,m)
c         k=dEindextmp(2,m)
      do l1=1,npole3b     
         i=pnum(l1)
        do kk=1,nembedlst(l1)
           k=embedlst(kk,l1)       
           kkk=oldlistinds(kk,l1)

          if(i.lt.k) then
          xr = x(k) - x(i)
          yr = y(k) - y(i)
          zr = z(k) - z(i)
          else
          xr = x(i) - x(k)
          yr = y(i) - y(k)
          zr = z(i) - z(k)
          end if
          call image (xr,yr,zr)
          r2=xr*xr + yr*yr +zr*zr
c             NEED TO CORRECT STUFF BELOW
            if(i.lt.k) then
              m=dEd1_3(1,kkk,i)
            else if (i.gt.k) then
              m=dEd1_3(1,kkk,k)
            end if

            dok=.false.
            do k3=1,npole3b
               if(pnum(k3).eq.k) then
                 dok=.true.
                 l3=k3
                 goto 31
               end if
            end do
   31    continue

          if ((r2 .le. ewaldcut2).and.(m.ne.0.0d0)) then

         vxx=0.0d0
         vyx=0.0d0
         vzx=0.0d0
         vyy=0.0d0
         vzy=0.0d0
         vzz=0.0d0

c               if(i.lt.k) then
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
c               else
c               iaz = zaxis(k)
c               iax = xaxis(k)
c               iay = yaxis(k)
c               kaz = zaxis(i)
c               kax = xaxis(i)
c               kay = yaxis(i)
c               if (iaz .eq. 0)  iaz = k
c               if (iax .eq. 0)  iax = k
c               if (iay .eq. 0)  iay = k
c               if (kaz .eq. 0)  kaz = i
c               if (kax .eq. 0)  kax = i
c               if (kay .eq. 0)  kay = i
c               xiz = x(iaz) - x(k)
c               yiz = y(iaz) - y(k)
c               ziz = z(iaz) - z(k)
c               xix = x(iax) - x(k)
c               yix = y(iax) - y(k)
c               zix = z(iax) - z(k)
c               xiy = x(iay) - x(k)
c               yiy = y(iay) - y(k)
c               ziy = z(iay) - z(k)
c               xkz = x(kaz) - x(i)
c               ykz = y(kaz) - y(i)
c               zkz = z(kaz) - z(i)
c               xkx = x(kax) - x(i)
c               ykx = y(kax) - y(i)
c               zkx = z(kax) - z(i)
c               xky = x(kay) - x(i)
c               yky = y(kay) - y(i)
c               zky = z(kay) - z(i)
c               end if

          if(i.lt.k) then        
          t2xxd=dEd2_3(1,kkk,i)
          t2xyd=dEd2_3(2,kkk,i)
          t2xzd=dEd2_3(3,kkk,i)
          t2yxd=dEd2_3(4,kkk,i)
          t2yyd=dEd2_3(5,kkk,i)
          t2yzd=dEd2_3(6,kkk,i)
          t2zxd=dEd2_3(7,kkk,i)
          t2zyd=dEd2_3(8,kkk,i)
          t2zzd=dEd2_3(9,kkk,i)

          t2xxp=dEp2_3(1,kkk,i)
          t2xyp=dEp2_3(2,kkk,i)
          t2xzp=dEp2_3(3,kkk,i)
          t2yxp=dEp2_3(4,kkk,i)
          t2yyp=dEp2_3(5,kkk,i)
          t2yzp=dEp2_3(6,kkk,i)
          t2zxp=dEp2_3(7,kkk,i)
          t2zyp=dEp2_3(8,kkk,i)
          t2zzp=dEp2_3(9,kkk,i)

          frcz2xxd=frcztau2dtot_3(1,kkk,i)
          frcz2xyd=frcztau2dtot_3(2,kkk,i)
          frcz2xzd=frcztau2dtot_3(3,kkk,i)
          frcz2yxd=frcztau2dtot_3(4,kkk,i)
          frcz2yyd=frcztau2dtot_3(5,kkk,i)
          frcz2yzd=frcztau2dtot_3(6,kkk,i)
          frcz2zxd=frcztau2dtot_3(7,kkk,i)
          frcz2zyd=frcztau2dtot_3(8,kkk,i)
          frcz2zzd=frcztau2dtot_3(9,kkk,i)

          frcz2xxp=frcztau2ptot_3(1,kkk,i)
          frcz2xyp=frcztau2ptot_3(2,kkk,i)
          frcz2xzp=frcztau2ptot_3(3,kkk,i)
          frcz2yxp=frcztau2ptot_3(4,kkk,i)
          frcz2yyp=frcztau2ptot_3(5,kkk,i)
          frcz2yzp=frcztau2ptot_3(6,kkk,i)
          frcz2zxp=frcztau2ptot_3(7,kkk,i)
          frcz2zyp=frcztau2ptot_3(8,kkk,i)
          frcz2zzp=frcztau2ptot_3(9,kkk,i)

          frcy2xxd=frcytau2dtot_3(1,kkk,i)
          frcy2xyd=frcytau2dtot_3(2,kkk,i)
          frcy2xzd=frcytau2dtot_3(3,kkk,i)
          frcy2yxd=frcytau2dtot_3(4,kkk,i)
          frcy2yyd=frcytau2dtot_3(5,kkk,i)
          frcy2yzd=frcytau2dtot_3(6,kkk,i)
          frcy2zxd=frcytau2dtot_3(7,kkk,i)
          frcy2zyd=frcytau2dtot_3(8,kkk,i)
          frcy2zzd=frcytau2dtot_3(9,kkk,i)

          frcy2xxp=frcytau2ptot_3(1,kkk,i)
          frcy2xyp=frcytau2ptot_3(2,kkk,i)
          frcy2xzp=frcytau2ptot_3(3,kkk,i)
          frcy2yxp=frcytau2ptot_3(4,kkk,i)
          frcy2yyp=frcytau2ptot_3(5,kkk,i)
          frcy2yzp=frcytau2ptot_3(6,kkk,i)
          frcy2zxp=frcytau2ptot_3(7,kkk,i)
          frcy2zyp=frcytau2ptot_3(8,kkk,i)
          frcy2zzp=frcytau2ptot_3(9,kkk,i)

          frcx2xxd=frcxtau2dtot_3(1,kkk,i)
          frcx2xyd=frcxtau2dtot_3(2,kkk,i)
          frcx2xzd=frcxtau2dtot_3(3,kkk,i)
          frcx2yxd=frcxtau2dtot_3(4,kkk,i)
          frcx2yyd=frcxtau2dtot_3(5,kkk,i)
          frcx2yzd=frcxtau2dtot_3(6,kkk,i)
          frcx2zxd=frcxtau2dtot_3(7,kkk,i)
          frcx2zyd=frcxtau2dtot_3(8,kkk,i)
          frcx2zzd=frcxtau2dtot_3(9,kkk,i)

          frcx2xxp=frcxtau2ptot_3(1,kkk,i)
          frcx2xyp=frcxtau2ptot_3(2,kkk,i)
          frcx2xzp=frcxtau2ptot_3(3,kkk,i)
          frcx2yxp=frcxtau2ptot_3(4,kkk,i)
          frcx2yyp=frcxtau2ptot_3(5,kkk,i)
          frcx2yzp=frcxtau2ptot_3(6,kkk,i)
          frcx2zxp=frcxtau2ptot_3(7,kkk,i)
          frcx2zyp=frcxtau2ptot_3(8,kkk,i)
          frcx2zzp=frcxtau2ptot_3(9,kkk,i)

          t1xxd=dEd1_3(1,kkk,i)
          t1xyd=dEd1_3(2,kkk,i)
          t1xzd=dEd1_3(3,kkk,i)
          t1yxd=dEd1_3(4,kkk,i)
          t1yyd=dEd1_3(5,kkk,i)
          t1yzd=dEd1_3(6,kkk,i)
          t1zxd=dEd1_3(7,kkk,i)
          t1zyd=dEd1_3(8,kkk,i)
          t1zzd=dEd1_3(9,kkk,i)

          t1xxp=dEp1_3(1,kkk,i)
          t1xyp=dEp1_3(2,kkk,i)
          t1xzp=dEp1_3(3,kkk,i)
          t1yxp=dEp1_3(4,kkk,i)
          t1yyp=dEp1_3(5,kkk,i)
          t1yzp=dEp1_3(6,kkk,i)
          t1zxp=dEp1_3(7,kkk,i)
          t1zyp=dEp1_3(8,kkk,i)
          t1zzp=dEp1_3(9,kkk,i)

          frcz1xxd=frcztau1dtot_3(1,kkk,i)
          frcz1xyd=frcztau1dtot_3(2,kkk,i)
          frcz1xzd=frcztau1dtot_3(3,kkk,i)
          frcz1yxd=frcztau1dtot_3(4,kkk,i)
          frcz1yyd=frcztau1dtot_3(5,kkk,i)
          frcz1yzd=frcztau1dtot_3(6,kkk,i)
          frcz1zxd=frcztau1dtot_3(7,kkk,i)
          frcz1zyd=frcztau1dtot_3(8,kkk,i)
          frcz1zzd=frcztau1dtot_3(9,kkk,i)

          frcz1xxp=frcztau1ptot_3(1,kkk,i)
          frcz1xyp=frcztau1ptot_3(2,kkk,i)
          frcz1xzp=frcztau1ptot_3(3,kkk,i)
          frcz1yxp=frcztau1ptot_3(4,kkk,i)
          frcz1yyp=frcztau1ptot_3(5,kkk,i)
          frcz1yzp=frcztau1ptot_3(6,kkk,i)
          frcz1zxp=frcztau1ptot_3(7,kkk,i)
          frcz1zyp=frcztau1ptot_3(8,kkk,i)
          frcz1zzp=frcztau1ptot_3(9,kkk,i)

          frcx1xxd=frcxtau1dtot_3(1,kkk,i)
          frcx1xyd=frcxtau1dtot_3(2,kkk,i)
          frcx1xzd=frcxtau1dtot_3(3,kkk,i)
          frcx1yxd=frcxtau1dtot_3(4,kkk,i)
          frcx1yyd=frcxtau1dtot_3(5,kkk,i)
          frcx1yzd=frcxtau1dtot_3(6,kkk,i)
          frcx1zxd=frcxtau1dtot_3(7,kkk,i)
          frcx1zyd=frcxtau1dtot_3(8,kkk,i)
          frcx1zzd=frcxtau1dtot_3(9,kkk,i)

          frcx1xxp=frcxtau1ptot_3(1,kkk,i)
          frcx1xyp=frcxtau1ptot_3(2,kkk,i)
          frcx1xzp=frcxtau1ptot_3(3,kkk,i)
          frcx1yxp=frcxtau1ptot_3(4,kkk,i)
          frcx1yyp=frcxtau1ptot_3(5,kkk,i)
          frcx1yzp=frcxtau1ptot_3(6,kkk,i)
          frcx1zxp=frcxtau1ptot_3(7,kkk,i)
          frcx1zyp=frcxtau1ptot_3(8,kkk,i)
          frcx1zzp=frcxtau1ptot_3(9,kkk,i)

          frcy1xxd=frcytau1dtot_3(1,kkk,i)
          frcy1xyd=frcytau1dtot_3(2,kkk,i)
          frcy1xzd=frcytau1dtot_3(3,kkk,i)
          frcy1yxd=frcytau1dtot_3(4,kkk,i)
          frcy1yyd=frcytau1dtot_3(5,kkk,i)
          frcy1yzd=frcytau1dtot_3(6,kkk,i)
          frcy1zxd=frcytau1dtot_3(7,kkk,i)
          frcy1zyd=frcytau1dtot_3(8,kkk,i)
          frcy1zzd=frcytau1dtot_3(9,kkk,i)

          frcy1xxp=frcytau1ptot_3(1,kkk,i)
          frcy1xyp=frcytau1ptot_3(2,kkk,i)
          frcy1xzp=frcytau1ptot_3(3,kkk,i)
          frcy1yxp=frcytau1ptot_3(4,kkk,i)
          frcy1yyp=frcytau1ptot_3(5,kkk,i)
          frcy1yzp=frcytau1ptot_3(6,kkk,i)
          frcy1zxp=frcytau1ptot_3(7,kkk,i)
          frcy1zyp=frcytau1ptot_3(8,kkk,i)
          frcy1zzp=frcytau1ptot_3(9,kkk,i)
          else
          t2xxd=dEd2_3(1,kkk,k)
          t2xyd=dEd2_3(2,kkk,k)
          t2xzd=dEd2_3(3,kkk,k)
          t2yxd=dEd2_3(4,kkk,k)
          t2yyd=dEd2_3(5,kkk,k)
          t2yzd=dEd2_3(6,kkk,k)
          t2zxd=dEd2_3(7,kkk,k)
          t2zyd=dEd2_3(8,kkk,k)
          t2zzd=dEd2_3(9,kkk,k)

          t2xxp=dEp2_3(1,kkk,k)
          t2xyp=dEp2_3(2,kkk,k)
          t2xzp=dEp2_3(3,kkk,k)
          t2yxp=dEp2_3(4,kkk,k)
          t2yyp=dEp2_3(5,kkk,k)
          t2yzp=dEp2_3(6,kkk,k)
          t2zxp=dEp2_3(7,kkk,k)
          t2zyp=dEp2_3(8,kkk,k)
          t2zzp=dEp2_3(9,kkk,k)

          frcz2xxd=frcztau2dtot_3(1,kkk,k)
          frcz2xyd=frcztau2dtot_3(2,kkk,k)
          frcz2xzd=frcztau2dtot_3(3,kkk,k)
          frcz2yxd=frcztau2dtot_3(4,kkk,k)
          frcz2yyd=frcztau2dtot_3(5,kkk,k)
          frcz2yzd=frcztau2dtot_3(6,kkk,k)
          frcz2zxd=frcztau2dtot_3(7,kkk,k)
          frcz2zyd=frcztau2dtot_3(8,kkk,k)
          frcz2zzd=frcztau2dtot_3(9,kkk,k)

          frcz2xxp=frcztau2ptot_3(1,kkk,k)
          frcz2xyp=frcztau2ptot_3(2,kkk,k)
          frcz2xzp=frcztau2ptot_3(3,kkk,k)
          frcz2yxp=frcztau2ptot_3(4,kkk,k)
          frcz2yyp=frcztau2ptot_3(5,kkk,k)
          frcz2yzp=frcztau2ptot_3(6,kkk,k)
          frcz2zxp=frcztau2ptot_3(7,kkk,k)
          frcz2zyp=frcztau2ptot_3(8,kkk,k)
          frcz2zzp=frcztau2ptot_3(9,kkk,k)

          frcy2xxd=frcytau2dtot_3(1,kkk,k)
          frcy2xyd=frcytau2dtot_3(2,kkk,k)
          frcy2xzd=frcytau2dtot_3(3,kkk,k)
          frcy2yxd=frcytau2dtot_3(4,kkk,k)
          frcy2yyd=frcytau2dtot_3(5,kkk,k)
          frcy2yzd=frcytau2dtot_3(6,kkk,k)
          frcy2zxd=frcytau2dtot_3(7,kkk,k)
          frcy2zyd=frcytau2dtot_3(8,kkk,k)
          frcy2zzd=frcytau2dtot_3(9,kkk,k)

          frcy2xxp=frcytau2ptot_3(1,kkk,k)
          frcy2xyp=frcytau2ptot_3(2,kkk,k)
          frcy2xzp=frcytau2ptot_3(3,kkk,k)
          frcy2yxp=frcytau2ptot_3(4,kkk,k)
          frcy2yyp=frcytau2ptot_3(5,kkk,k)
          frcy2yzp=frcytau2ptot_3(6,kkk,k)
          frcy2zxp=frcytau2ptot_3(7,kkk,k)
          frcy2zyp=frcytau2ptot_3(8,kkk,k)
          frcy2zzp=frcytau2ptot_3(9,kkk,k)

          frcx2xxd=frcxtau2dtot_3(1,kkk,k)
          frcx2xyd=frcxtau2dtot_3(2,kkk,k)
          frcx2xzd=frcxtau2dtot_3(3,kkk,k)
          frcx2yxd=frcxtau2dtot_3(4,kkk,k)
          frcx2yyd=frcxtau2dtot_3(5,kkk,k)
          frcx2yzd=frcxtau2dtot_3(6,kkk,k)
          frcx2zxd=frcxtau2dtot_3(7,kkk,k)
          frcx2zyd=frcxtau2dtot_3(8,kkk,k)
          frcx2zzd=frcxtau2dtot_3(9,kkk,k)

          frcx2xxp=frcxtau2ptot_3(1,kkk,k)
          frcx2xyp=frcxtau2ptot_3(2,kkk,k)
          frcx2xzp=frcxtau2ptot_3(3,kkk,k)
          frcx2yxp=frcxtau2ptot_3(4,kkk,k)
          frcx2yyp=frcxtau2ptot_3(5,kkk,k)
          frcx2yzp=frcxtau2ptot_3(6,kkk,k)
          frcx2zxp=frcxtau2ptot_3(7,kkk,k)
          frcx2zyp=frcxtau2ptot_3(8,kkk,k)
          frcx2zzp=frcxtau2ptot_3(9,kkk,k)

          t1xxd=dEd1_3(1,kkk,k)
          t1xyd=dEd1_3(2,kkk,k)
          t1xzd=dEd1_3(3,kkk,k)
          t1yxd=dEd1_3(4,kkk,k)
          t1yyd=dEd1_3(5,kkk,k)
          t1yzd=dEd1_3(6,kkk,k)
          t1zxd=dEd1_3(7,kkk,k)
          t1zyd=dEd1_3(8,kkk,k)
          t1zzd=dEd1_3(9,kkk,k)

          t1xxp=dEp1_3(1,kkk,k)
          t1xyp=dEp1_3(2,kkk,k)
          t1xzp=dEp1_3(3,kkk,k)
          t1yxp=dEp1_3(4,kkk,k)
          t1yyp=dEp1_3(5,kkk,k)
          t1yzp=dEp1_3(6,kkk,k)
          t1zxp=dEp1_3(7,kkk,k)
          t1zyp=dEp1_3(8,kkk,k)
          t1zzp=dEp1_3(9,kkk,k)

          frcz1xxd=frcztau1dtot_3(1,kkk,k)
          frcz1xyd=frcztau1dtot_3(2,kkk,k)
          frcz1xzd=frcztau1dtot_3(3,kkk,k)
          frcz1yxd=frcztau1dtot_3(4,kkk,k)
          frcz1yyd=frcztau1dtot_3(5,kkk,k)
          frcz1yzd=frcztau1dtot_3(6,kkk,k)
          frcz1zxd=frcztau1dtot_3(7,kkk,k)
          frcz1zyd=frcztau1dtot_3(8,kkk,k)
          frcz1zzd=frcztau1dtot_3(9,kkk,k)

          frcz1xxp=frcztau1ptot_3(1,kkk,k)
          frcz1xyp=frcztau1ptot_3(2,kkk,k)
          frcz1xzp=frcztau1ptot_3(3,kkk,k)
          frcz1yxp=frcztau1ptot_3(4,kkk,k)
          frcz1yyp=frcztau1ptot_3(5,kkk,k)
          frcz1yzp=frcztau1ptot_3(6,kkk,k)
          frcz1zxp=frcztau1ptot_3(7,kkk,k)
          frcz1zyp=frcztau1ptot_3(8,kkk,k)
          frcz1zzp=frcztau1ptot_3(9,kkk,k)

          frcx1xxd=frcxtau1dtot_3(1,kkk,k)
          frcx1xyd=frcxtau1dtot_3(2,kkk,k)
          frcx1xzd=frcxtau1dtot_3(3,kkk,k)
          frcx1yxd=frcxtau1dtot_3(4,kkk,k)
          frcx1yyd=frcxtau1dtot_3(5,kkk,k)
          frcx1yzd=frcxtau1dtot_3(6,kkk,k)
          frcx1zxd=frcxtau1dtot_3(7,kkk,k)
          frcx1zyd=frcxtau1dtot_3(8,kkk,k)
          frcx1zzd=frcxtau1dtot_3(9,kkk,k)

          frcx1xxp=frcxtau1ptot_3(1,kkk,k)
          frcx1xyp=frcxtau1ptot_3(2,kkk,k)
          frcx1xzp=frcxtau1ptot_3(3,kkk,k)
          frcx1yxp=frcxtau1ptot_3(4,kkk,k)
          frcx1yyp=frcxtau1ptot_3(5,kkk,k)
          frcx1yzp=frcxtau1ptot_3(6,kkk,k)
          frcx1zxp=frcxtau1ptot_3(7,kkk,k)
          frcx1zyp=frcxtau1ptot_3(8,kkk,k)
          frcx1zzp=frcxtau1ptot_3(9,kkk,k)

          frcy1xxd=frcytau1dtot_3(1,kkk,k)
          frcy1xyd=frcytau1dtot_3(2,kkk,k)
          frcy1xzd=frcytau1dtot_3(3,kkk,k)
          frcy1yxd=frcytau1dtot_3(4,kkk,k)
          frcy1yyd=frcytau1dtot_3(5,kkk,k)
          frcy1yzd=frcytau1dtot_3(6,kkk,k)
          frcy1zxd=frcytau1dtot_3(7,kkk,k)
          frcy1zyd=frcytau1dtot_3(8,kkk,k)
          frcy1zzd=frcytau1dtot_3(9,kkk,k)

          frcy1xxp=frcytau1ptot_3(1,kkk,k)
          frcy1xyp=frcytau1ptot_3(2,kkk,k)
          frcy1xzp=frcytau1ptot_3(3,kkk,k)
          frcy1yxp=frcytau1ptot_3(4,kkk,k)
          frcy1yyp=frcytau1ptot_3(5,kkk,k)
          frcy1yzp=frcytau1ptot_3(6,kkk,k)
          frcy1zxp=frcytau1ptot_3(7,kkk,k)
          frcy1zyp=frcytau1ptot_3(8,kkk,k)
          frcy1zzp=frcytau1ptot_3(9,kkk,k)
          end if

c   NEW
        if(i.lt.k) then         
             termx=uind(1,l1)*t2xxp+ uind(2,l1)*t2xyp + uind(3,l1)*t2xzp
             termy=uind(1,l1)*t2yxp+ uind(2,l1)*t2yyp + uind(3,l1)*t2yzp
             termz=uind(1,l1)*t2zxp+ uind(2,l1)*t2zyp + uind(3,l1)*t2zzp

             frczk(1)=uind(1,l1)*frcz2xxp+uind(2,l1)*frcz2xyp
     &                +uind(3,l1)*frcz2xzp
             frczk(2)=uind(1,l1)*frcz2yxp+uind(2,l1)*frcz2yyp
     &                +uind(3,l1)*frcz2yzp
             frczk(3)=uind(1,l1)*frcz2zxp+uind(2,l1)*frcz2zyp
     &                +uind(3,l1)*frcz2zzp

             frcxk(1)=uind(1,l1)*frcx2xxp+uind(2,l1)*frcx2xyp
     &                +uind(3,l1)*frcx2xzp
             frcxk(2)=uind(1,l1)*frcx2yxp+uind(2,l1)*frcx2yyp
     &                +uind(3,l1)*frcx2yzp
             frcxk(3)=uind(1,l1)*frcx2zxp+uind(2,l1)*frcx2zyp
     &                +uind(3,l1)*frcx2zzp

             frcyk(1)=uind(1,l1)*frcy2xxp+uind(2,l1)*frcy2xyp
     &                +uind(3,l1)*frcy2xzp
             frcyk(2)=uind(1,l1)*frcy2yxp+uind(2,l1)*frcy2yyp
     &                +uind(3,l1)*frcy2yzp
             frcyk(3)=uind(1,l1)*frcy2zxp+uind(2,l1)*frcy2zyp
     &                +uind(3,l1)*frcy2zzp

             deptemp(1,i)=deptemp(1,i)+termx
             deptemp(2,i)=deptemp(2,i)+termy
             deptemp(3,i)=deptemp(3,i)+termz

             deptemp(1,k)=deptemp(1,k)-termx
             deptemp(2,k)=deptemp(2,k)-termy
             deptemp(3,k)=deptemp(3,k)-termz

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

             termx=uinp(1,l1)*t2xxd+ uinp(2,l1)*t2xyd + uinp(3,l1)*t2xzd
             termy=uinp(1,l1)*t2yxd+ uinp(2,l1)*t2yyd + uinp(3,l1)*t2yzd
             termz=uinp(1,l1)*t2zxd+ uinp(2,l1)*t2zyd + uinp(3,l1)*t2zzd

             frczk(1)=uinp(1,l1)*frcz2xxd+uinp(2,l1)*frcz2xyd
     &                 +uinp(3,l1)*frcz2xzd
             frczk(2)=uinp(1,l1)*frcz2yxd+uinp(2,l1)*frcz2yyd
     &                 +uinp(3,l1)*frcz2yzd
             frczk(3)=uinp(1,l1)*frcz2zxd+uinp(2,l1)*frcz2zyd
     &                 +uinp(3,l1)*frcz2zzd

             frcxk(1)=uinp(1,l1)*frcx2xxd+uinp(2,l1)*frcx2xyd
     &                  +uinp(3,l1)*frcx2xzd
             frcxk(2)=uinp(1,l1)*frcx2yxd+uinp(2,l1)*frcx2yyd
     &                  +uinp(3,l1)*frcx2yzd
             frcxk(3)=uinp(1,l1)*frcx2zxd+uinp(2,l1)*frcx2zyd
     &                  +uinp(3,l1)*frcx2zzd

             frcyk(1)=uinp(1,l1)*frcy2xxd+uinp(2,l1)*frcy2xyd
     &                  +uinp(3,l1)*frcy2xzd
             frcyk(2)=uinp(1,l1)*frcy2yxd+uinp(2,l1)*frcy2yyd
     &                  +uinp(3,l1)*frcy2yzd
             frcyk(3)=uinp(1,l1)*frcy2zxd+uinp(2,l1)*frcy2zyd
     &                  +uinp(3,l1)*frcy2zzd

             deptemp(1,i)=deptemp(1,i)+termx
             deptemp(2,i)=deptemp(2,i)+termy
             deptemp(3,i)=deptemp(3,i)+termz

             deptemp(1,k)=deptemp(1,k)-termx
             deptemp(2,k)=deptemp(2,k)-termy
             deptemp(3,k)=deptemp(3,k)-termz

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

             if(dok) then
             termx=uind(1,l3)*t1xxp+ uind(2,l3)*t1xyp + uind(3,l3)*t1xzp
             termy=uind(1,l3)*t1yxp+ uind(2,l3)*t1yyp + uind(3,l3)*t1yzp
             termz=uind(1,l3)*t1zxp+ uind(2,l3)*t1zyp + uind(3,l3)*t1zzp

             frczi(1)=uind(1,l3)*frcz1xxp+uind(2,l3)*frcz1xyp
     &                +uind(3,l3)*frcz1xzp
             frczi(2)=uind(1,l3)*frcz1yxp+uind(2,l3)*frcz1yyp
     &                +uind(3,l3)*frcz1yzp
             frczi(3)=uind(1,l3)*frcz1zxp+uind(2,l3)*frcz1zyp
     &                +uind(3,l3)*frcz1zzp

             frcxi(1)=uind(1,l3)*frcx1xxp+uind(2,l3)*frcx1xyp
     &                +uind(3,l3)*frcx1xzp
             frcxi(2)=uind(1,l3)*frcx1yxp+uind(2,l3)*frcx1yyp
     &                +uind(3,l3)*frcx1yzp
             frcxi(3)=uind(1,l3)*frcx1zxp+uind(2,l3)*frcx1zyp
     &                +uind(3,l3)*frcx1zzp

             frcyi(1)=uind(1,l3)*frcy1xxp+uind(2,l3)*frcy1xyp
     &                +uind(3,l3)*frcy1xzp
             frcyi(2)=uind(1,l3)*frcy1yxp+uind(2,l3)*frcy1yyp
     &                +uind(3,l3)*frcy1yzp
             frcyi(3)=uind(1,l3)*frcy1zxp+uind(2,l3)*frcy1zyp
     &                +uind(3,l3)*frcy1zzp

             deptemp(1,i)=deptemp(1,i)+termx
             deptemp(2,i)=deptemp(2,i)+termy
             deptemp(3,i)=deptemp(3,i)+termz

             deptemp(1,k)=deptemp(1,k)-termx
             deptemp(2,k)=deptemp(2,k)-termy
             deptemp(3,k)=deptemp(3,k)-termz
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
             termx=uinp(1,l3)*t1xxd+ uinp(2,l3)*t1xyd + uinp(3,l3)*t1xzd
             termy=uinp(1,l3)*t1yxd+ uinp(2,l3)*t1yyd + uinp(3,l3)*t1yzd
             termz=uinp(1,l3)*t1zxd+ uinp(2,l3)*t1zyd + uinp(3,l3)*t1zzd

             frczi(1)=uinp(1,l3)*frcz1xxd+uinp(2,l3)*frcz1xyd
     &                +uinp(3,l3)*frcz1xzd
             frczi(2)=uinp(1,l3)*frcz1yxd+uinp(2,l3)*frcz1yyd
     &                +uinp(3,l3)*frcz1yzd
             frczi(3)=uinp(1,l3)*frcz1zxd+uinp(2,l3)*frcz1zyd
     &                +uinp(3,l3)*frcz1zzd
         
             frcxi(1)=uinp(1,l3)*frcx1xxd+uinp(2,l3)*frcx1xyd
     &                +uinp(3,l3)*frcx1xzd
             frcxi(2)=uinp(1,l3)*frcx1yxd+uinp(2,l3)*frcx1yyd
     &                +uinp(3,l3)*frcx1yzd
             frcxi(3)=uinp(1,l3)*frcx1zxd+uinp(2,l3)*frcx1zyd
     &                +uinp(3,l3)*frcx1zzd
          
             frcyi(1)=uinp(1,l3)*frcy1xxd+uinp(2,l3)*frcy1xyd
     &                +uinp(3,l3)*frcy1xzd
             frcyi(2)=uinp(1,l3)*frcy1yxd+uinp(2,l3)*frcy1yyd
     &                +uinp(3,l3)*frcy1yzd
             frcyi(3)=uinp(1,l3)*frcy1zxd+uinp(2,l3)*frcy1zyd
     &                +uinp(3,l3)*frcy1zzd

             deptemp(1,i)=deptemp(1,i)+termx
             deptemp(2,i)=deptemp(2,i)+termy
             deptemp(3,i)=deptemp(3,i)+termz
          
             deptemp(1,k)=deptemp(1,k)-termx
             deptemp(2,k)=deptemp(2,k)-termy
             deptemp(3,k)=deptemp(3,k)-termz
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
             end if
        else
             termx=uind(1,l1)*t1xxp+ uind(2,l1)*t1xyp + uind(3,l1)*t1xzp
             termy=uind(1,l1)*t1yxp+ uind(2,l1)*t1yyp + uind(3,l1)*t1yzp
             termz=uind(1,l1)*t1zxp+ uind(2,l1)*t1zyp + uind(3,l1)*t1zzp

             frczk(1)=uind(1,l1)*frcz1xxp+uind(2,l1)*frcz1xyp
     &                +uind(3,l1)*frcz1xzp
             frczk(2)=uind(1,l1)*frcz1yxp+uind(2,l1)*frcz1yyp
     &                +uind(3,l1)*frcz1yzp
             frczk(3)=uind(1,l1)*frcz1zxp+uind(2,l1)*frcz1zyp
     &                +uind(3,l1)*frcz1zzp

             frcxk(1)=uind(1,l1)*frcx1xxp+uind(2,l1)*frcx1xyp
     &                +uind(3,l1)*frcx1xzp
             frcxk(2)=uind(1,l1)*frcx1yxp+uind(2,l1)*frcx1yyp
     &                +uind(3,l1)*frcx1yzp
             frcxk(3)=uind(1,l1)*frcx1zxp+uind(2,l1)*frcx1zyp
     &                +uind(3,l1)*frcx1zzp

             frcyk(1)=uind(1,l1)*frcy1xxp+uind(2,l1)*frcy1xyp
     &                +uind(3,l1)*frcy1xzp
             frcyk(2)=uind(1,l1)*frcy1yxp+uind(2,l1)*frcy1yyp
     &                +uind(3,l1)*frcy1yzp
             frcyk(3)=uind(1,l1)*frcy1zxp+uind(2,l1)*frcy1zyp
     &                +uind(3,l1)*frcy1zzp


c             deptemp(1,i)=deptemp(1,i)+termx
c             deptemp(2,i)=deptemp(2,i)+termy
c             deptemp(3,i)=deptemp(3,i)+termz
c
c             deptemp(1,k)=deptemp(1,k)-termx
c             deptemp(2,k)=deptemp(2,k)-termy
c             deptemp(3,k)=deptemp(3,k)-termz

             deptemp(1,i)=deptemp(1,i)-termx
             deptemp(2,i)=deptemp(2,i)-termy
             deptemp(3,i)=deptemp(3,i)-termz

             deptemp(1,k)=deptemp(1,k)+termx
             deptemp(2,k)=deptemp(2,k)+termy
             deptemp(3,k)=deptemp(3,k)+termz

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

c             vxx =vxx -xr*termx +xix*frcxk(1)
c     &                  + xiy*frcyk(1) + xiz*frczk(1)
c             vyx =vyx -yr*termx + yix*frcxk(1)
c     &                  + yiy*frcyk(1) + yiz*frczk(1)
c             vzx =vzx -zr*termx + zix*frcxk(1)
c     &                  + ziy*frcyk(1) + ziz*frczk(1)
c             vyy =vyy -yr*termy + yix*frcxk(2)
c     &                  + yiy*frcyk(2) + yiz*frczk(2)
c             vzy =vzy -zr*termy + zix*frcxk(2)
c     &                  + ziy*frcyk(2) + ziz*frczk(2)
c             vzz =vzz -zr*termz + zix*frcxk(3)
c     &                  + ziy*frcyk(3) + ziz*frczk(3)


             termx=uinp(1,l1)*t1xxd+ uinp(2,l1)*t1xyd + uinp(3,l1)*t1xzd
             termy=uinp(1,l1)*t1yxd+ uinp(2,l1)*t1yyd + uinp(3,l1)*t1yzd
             termz=uinp(1,l1)*t1zxd+ uinp(2,l1)*t1zyd + uinp(3,l1)*t1zzd

             frczk(1)=uinp(1,l1)*frcz1xxd+uinp(2,l1)*frcz1xyd
     &                +uinp(3,l1)*frcz1xzd
             frczk(2)=uinp(1,l1)*frcz1yxd+uinp(2,l1)*frcz1yyd
     &                +uinp(3,l1)*frcz1yzd
             frczk(3)=uinp(1,l1)*frcz1zxd+uinp(2,l1)*frcz1zyd
     &                +uinp(3,l1)*frcz1zzd

             frcxk(1)=uinp(1,l1)*frcx1xxd+uinp(2,l1)*frcx1xyd
     &                +uinp(3,l1)*frcx1xzd
             frcxk(2)=uinp(1,l1)*frcx1yxd+uinp(2,l1)*frcx1yyd
     &                +uinp(3,l1)*frcx1yzd
             frcxk(3)=uinp(1,l1)*frcx1zxd+uinp(2,l1)*frcx1zyd
     &                +uinp(3,l1)*frcx1zzd

             frcyk(1)=uinp(1,l1)*frcy1xxd+uinp(2,l1)*frcy1xyd
     &                +uinp(3,l1)*frcy1xzd
             frcyk(2)=uinp(1,l1)*frcy1yxd+uinp(2,l1)*frcy1yyd
     &                +uinp(3,l1)*frcy1yzd
             frcyk(3)=uinp(1,l1)*frcy1zxd+uinp(2,l1)*frcy1zyd
     &                +uinp(3,l1)*frcy1zzd

c             deptemp(1,i)=deptemp(1,i)+termx
c             deptemp(2,i)=deptemp(2,i)+termy
c             deptemp(3,i)=deptemp(3,i)+termz
c
c             deptemp(1,k)=deptemp(1,k)-termx
c             deptemp(2,k)=deptemp(2,k)-termy
c             deptemp(3,k)=deptemp(3,k)-termz

             deptemp(1,i)=deptemp(1,i)-termx
             deptemp(2,i)=deptemp(2,i)-termy
             deptemp(3,i)=deptemp(3,i)-termz

             deptemp(1,k)=deptemp(1,k)+termx
             deptemp(2,k)=deptemp(2,k)+termy
             deptemp(3,k)=deptemp(3,k)+termz


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

c             vxx =vxx -xr*termx +xix*frcxk(1)
c     &                  + xiy*frcyk(1) + xiz*frczk(1)
c             vyx =vyx -yr*termx + yix*frcxk(1)
c     &                  + yiy*frcyk(1) + yiz*frczk(1)
c             vzx =vzx -zr*termx + zix*frcxk(1)
c     &                  + ziy*frcyk(1) + ziz*frczk(1)
c             vyy =vyy -yr*termy + yix*frcxk(2)
c     &                  + yiy*frcyk(2) + yiz*frczk(2)
c             vzy =vzy -zr*termy + zix*frcxk(2)
c     &                  + ziy*frcyk(2) + ziz*frczk(2)
c             vzz =vzz -zr*termz + zix*frcxk(3)
c     &                  + ziy*frcyk(3) + ziz*frczk(3)

        end if


               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz
c        end do         
c      end do
c!$OMP END DO
c
c!$OMP DO reduction(+:depi,depk)
c!$OMP& schedule(guided)
c      do m=1,ntpair_tau1tmp
c         i=tau1indextmp(1,m)
c         k=tau1indextmp(2,m)

           !   goto 41

          do j1=1,4
          !print*,"In j1 loop over torque tens"
          if(i.lt.k) then          
            k1=elsttau1_index(j1,kkk,i)
c            if(k1.ne.0) then
c             m1=elsttau1(j1,kkk,i)
c            end if
            k2=elsttau2_index(j1,kkk,i)
c            if(k2.ne.0) then
c             m2=elsttau2(j1,kkk,i)
c            end if
          else
            k1=elsttau1_index(j1,kkk,k)
c            if(k1.ne.0) then
c             m1=elsttau1(j1,kkk,k)
c            end if
            k2=elsttau2_index(j1,kkk,k)
c            if(k2.ne.0) then 
c             m2=elsttau2(j1,kkk,k)
c            end if
          end if         

          if(k1.ne.0) then 
            if(i.lt.k) then
            tau1xxd=taud1_3(1,j1,kkk,i)
            tau1xyd=taud1_3(2,j1,kkk,i)
            tau1xzd=taud1_3(3,j1,kkk,i)
            tau1yxd=taud1_3(4,j1,kkk,i)
            tau1yyd=taud1_3(5,j1,kkk,i)
            tau1yzd=taud1_3(6,j1,kkk,i)
            tau1zxd=taud1_3(7,j1,kkk,i)
            tau1zyd=taud1_3(8,j1,kkk,i)
            tau1zzd=taud1_3(9,j1,kkk,i)

            tau1xxp=taup1_3(1,j1,kkk,i)
            tau1xyp=taup1_3(2,j1,kkk,i)
            tau1xzp=taup1_3(3,j1,kkk,i)
            tau1yxp=taup1_3(4,j1,kkk,i)
            tau1yyp=taup1_3(5,j1,kkk,i)
            tau1yzp=taup1_3(6,j1,kkk,i)
            tau1zxp=taup1_3(7,j1,kkk,i)
            tau1zyp=taup1_3(8,j1,kkk,i)
            tau1zzp=taup1_3(9,j1,kkk,i)
            else
            tau1xxd=taud1_3(1,j1,kkk,k)
            tau1xyd=taud1_3(2,j1,kkk,k)
            tau1xzd=taud1_3(3,j1,kkk,k)
            tau1yxd=taud1_3(4,j1,kkk,k)
            tau1yyd=taud1_3(5,j1,kkk,k)
            tau1yzd=taud1_3(6,j1,kkk,k)
            tau1zxd=taud1_3(7,j1,kkk,k)
            tau1zyd=taud1_3(8,j1,kkk,k)
            tau1zzd=taud1_3(9,j1,kkk,k)

            tau1xxp=taup1_3(1,j1,kkk,k)
            tau1xyp=taup1_3(2,j1,kkk,k)
            tau1xzp=taup1_3(3,j1,kkk,k)
            tau1yxp=taup1_3(4,j1,kkk,k)
            tau1yyp=taup1_3(5,j1,kkk,k)
            tau1yzp=taup1_3(6,j1,kkk,k)
            tau1zxp=taup1_3(7,j1,kkk,k)
            tau1zyp=taup1_3(8,j1,kkk,k)
            tau1zzp=taup1_3(9,j1,kkk,k)
            end if

          end if

          if(k2.ne.0) then

            if(i.lt.k) then
            tau2xxd=taud2_3(1,j1,kkk,i)
            tau2xyd=taud2_3(2,j1,kkk,i)
            tau2xzd=taud2_3(3,j1,kkk,i)
            tau2yxd=taud2_3(4,j1,kkk,i)
            tau2yyd=taud2_3(5,j1,kkk,i)
            tau2yzd=taud2_3(6,j1,kkk,i)
            tau2zxd=taud2_3(7,j1,kkk,i)
            tau2zyd=taud2_3(8,j1,kkk,i)
            tau2zzd=taud2_3(9,j1,kkk,i)

            tau2xxp=taup2_3(1,j1,kkk,i)
            tau2xyp=taup2_3(2,j1,kkk,i)
            tau2xzp=taup2_3(3,j1,kkk,i)
            tau2yxp=taup2_3(4,j1,kkk,i)
            tau2yyp=taup2_3(5,j1,kkk,i)
            tau2yzp=taup2_3(6,j1,kkk,i)
            tau2zxp=taup2_3(7,j1,kkk,i)
            tau2zyp=taup2_3(8,j1,kkk,i)
            tau2zzp=taup2_3(9,j1,kkk,i)
            else
            tau2xxd=taud2_3(1,j1,kkk,k)
            tau2xyd=taud2_3(2,j1,kkk,k)
            tau2xzd=taud2_3(3,j1,kkk,k)
            tau2yxd=taud2_3(4,j1,kkk,k)
            tau2yyd=taud2_3(5,j1,kkk,k)
            tau2yzd=taud2_3(6,j1,kkk,k)
            tau2zxd=taud2_3(7,j1,kkk,k)
            tau2zyd=taud2_3(8,j1,kkk,k)
            tau2zzd=taud2_3(9,j1,kkk,k)

            tau2xxp=taup2_3(1,j1,kkk,k)
            tau2xyp=taup2_3(2,j1,kkk,k)
            tau2xzp=taup2_3(3,j1,kkk,k)
            tau2yxp=taup2_3(4,j1,kkk,k)
            tau2yyp=taup2_3(5,j1,kkk,k)
            tau2yzp=taup2_3(6,j1,kkk,k)
            tau2zxp=taup2_3(7,j1,kkk,k)
            tau2zyp=taup2_3(8,j1,kkk,k)
            tau2zzp=taup2_3(9,j1,kkk,k)
            end if

          end if

c          dok1=.false.
c          do k3=1,npole3b
c            if(pnum(k3).eq.k1) then
c              dok1=.true.
c              l3=k3
c              goto 31
c            end if
c          end do
c   31    continue
c
c          dok2=.false.
c          do k3=1,npole3b
c            if(pnum(k3).eq.k2) then
c              dok2=.true.
c              l3_2=k3
c              goto 31
c            end if
c          end do
c   31    continue


          if((i.lt.k)) then

            if(dok.and.(k1.ne.0)) then                    
c          taux=uind(1,k)*tau1xxp+ uind(2,k)*tau1xyp + uind(3,k)*tau1xzp
c          tauy=uind(1,k)*tau1yxp+ uind(2,k)*tau1yyp + uind(3,k)*tau1yzp
c          tauz=uind(1,k)*tau1zxp+ uind(2,k)*tau1zyp + uind(3,k)*tau1zzp

            taux=uind(1,l3)*tau1xxp+ uind(2,l3)*tau1xyp 
     &           + uind(3,l3)*tau1xzp
            tauy=uind(1,l3)*tau1yxp+ uind(2,l3)*tau1yyp 
     &           + uind(3,l3)*tau1yzp
            tauz=uind(1,l3)*tau1zxp+ uind(2,l3)*tau1zyp 
     &           + uind(3,l3)*tau1zzp

c            depi(1,i)=depi(1,i)+taux
c            depi(2,i)=depi(2,i)+tauy
c            depi(3,i)=depi(3,i)+tauz

            deptemp(1,k1)=deptemp(1,k1)+taux
            deptemp(2,k1)=deptemp(2,k1)+tauy
            deptemp(3,k1)=deptemp(3,k1)+tauz

c          taux=uinp(1,k)*tau1xxd+ uinp(2,k)*tau1xyd +uinp(3,k)*tau1xzd
c          tauy=uinp(1,k)*tau1yxd+ uinp(2,k)*tau1yyd +uinp(3,k)*tau1yzd
c          tauz=uinp(1,k)*tau1zxd+ uinp(2,k)*tau1zyd +uinp(3,k)*tau1zzd

            taux=uinp(1,l3)*tau1xxd+ uinp(2,l3)*tau1xyd 
     &           +uinp(3,l3)*tau1xzd
            tauy=uinp(1,l3)*tau1yxd+ uinp(2,l3)*tau1yyd 
     &           +uinp(3,l3)*tau1yzd
            tauz=uinp(1,l3)*tau1zxd+ uinp(2,l3)*tau1zyd 
     &           +uinp(3,l3)*tau1zzd

c            depi(1,i)=depi(1,i)+taux
c            depi(2,i)=depi(2,i)+tauy
c            depi(3,i)=depi(3,i)+tauz

            deptemp(1,k1)=deptemp(1,k1)+taux
            deptemp(2,k1)=deptemp(2,k1)+tauy
            deptemp(3,k1)=deptemp(3,k1)+tauz
            end if
c      end do
        
c      do m=1,ntpair_tau2tmp
c         i=tau2indextmp(1,m)
c         k=tau2indextmp(2,m)


c          taux=uind(1,i)*tau2xxp+ uind(2,i)*tau2xyp + uind(3,i)*tau2xzp
c          tauy=uind(1,i)*tau2yxp+ uind(2,i)*tau2yyp + uind(3,i)*tau2yzp
c          tauz=uind(1,i)*tau2zxp+ uind(2,i)*tau2zyp + uind(3,i)*tau2zzp
            if(k2.ne.0) then
            taux=uind(1,l1)*tau2xxp+ uind(2,l1)*tau2xyp 
     &         + uind(3,l1)*tau2xzp
            tauy=uind(1,l1)*tau2yxp+ uind(2,l1)*tau2yyp 
     &         + uind(3,l1)*tau2yzp
            tauz=uind(1,l1)*tau2zxp+ uind(2,l1)*tau2zyp 
     &         + uind(3,l1)*tau2zzp

c          depk(1,k)=depk(1,k)+taux
c          depk(2,k)=depk(2,k)+tauy
c          depk(3,k)=depk(3,k)+tauz

            deptemp(1,k2)=deptemp(1,k2)+taux
            deptemp(2,k2)=deptemp(2,k2)+tauy
            deptemp(3,k2)=deptemp(3,k2)+tauz

c          taux=uinp(1,i)*tau2xxd+ uinp(2,i)*tau2xyd +uinp(3,i)*tau2xzd
c          tauy=uinp(1,i)*tau2yxd+ uinp(2,i)*tau2yyd +uinp(3,i)*tau2yzd
c          tauz=uinp(1,i)*tau2zxd+ uinp(2,i)*tau2zyd +uinp(3,i)*tau2zzd

            taux=uinp(1,l1)*tau2xxd+ uinp(2,l1)*tau2xyd
     &           +uinp(3,l1)*tau2xzd
            tauy=uinp(1,l1)*tau2yxd+ uinp(2,l1)*tau2yyd
     &           +uinp(3,l1)*tau2yzd
            tauz=uinp(1,l1)*tau2zxd+ uinp(2,l1)*tau2zyd
     &           +uinp(3,l1)*tau2zzd

c          depk(1,k)=depk(1,k)+taux
c          depk(2,k)=depk(2,k)+tauy
c          depk(3,k)=depk(3,k)+tauz

            deptemp(1,k2)=deptemp(1,k2)+taux
            deptemp(2,k2)=deptemp(2,k2)+tauy
            deptemp(3,k2)=deptemp(3,k2)+tauz
            end if
          else
c            if(dok.and.(k2.ne.0)) then
c            print*,"Tilt! Should not be doing k here!"
c            taux=uind(1,l3)*tau2xxp+ uind(2,l3)*tau2xyp
c     &           +uind(3,l3)*tau2xzp
c            tauy=uind(1,l3)*tau2yxp+ uind(2,l3)*tau2yyp
c     &           +uind(3,l3)*tau2yzp
c            tauz=uind(1,l3)*tau2zxp+ uind(2,l3)*tau2zyp
c     &           +uind(3,l3)*tau2zzp
c
c            deptemp(1,k2)=deptemp(1,k2)+taux
c            deptemp(2,k2)=deptemp(2,k2)+tauy
c            deptemp(3,k2)=deptemp(3,k2)+tauz
c
c
c            taux=uinp(1,l3)*tau2xxd+uinp(2,l3)*tau2xyd
c     &           +uinp(3,l3)*tau2xzd
c            tauy=uinp(1,l3)*tau2yxd+uinp(2,l3)*tau2yyd
c     &           +uinp(3,l3)*tau2yzd
c            tauz=uinp(1,l3)*tau2zxd+uinp(2,l3)*tau2zyd
c     &           +uinp(3,l3)*tau2zzd
c
c
c            deptemp(1,k2)=deptemp(1,k2)+taux
c            deptemp(2,k2)=deptemp(2,k2)+tauy
c            deptemp(3,k2)=deptemp(3,k2)+tauz
c            end if

            if(k1.ne.0) then

          taux=uind(1,l1)*tau1xxp+ uind(2,l1)*tau1xyp+uind(3,l1)*tau1xzp
          tauy=uind(1,l1)*tau1yxp+ uind(2,l1)*tau1yyp+uind(3,l1)*tau1yzp
          tauz=uind(1,l1)*tau1zxp+ uind(2,l1)*tau1zyp+uind(3,l1)*tau1zzp

          deptemp(1,k1)=deptemp(1,k1)+taux
          deptemp(2,k1)=deptemp(2,k1)+tauy
          deptemp(3,k1)=deptemp(3,k1)+tauz

          taux=uinp(1,l1)*tau1xxd+ uinp(2,l1)*tau1xyd+uinp(3,l1)*tau1xzd
          tauy=uinp(1,l1)*tau1yxd+ uinp(2,l1)*tau1yyd+uinp(3,l1)*tau1yzd
          tauz=uinp(1,l1)*tau1zxd+ uinp(2,l1)*tau1zyd+uinp(3,l1)*tau1zzd


          deptemp(1,k1)=deptemp(1,k1)+taux
          deptemp(2,k1)=deptemp(2,k1)+tauy
          deptemp(3,k1)=deptemp(3,k1)+tauz
            end if
          end if
          end do
   41    continue

          end if
        end do
      end do
c!$OMP END DO
c
c!$OMP DO
c      do i=1,npole
c        dep3bt(1,i)=dep3bt(1,i)+depi(1,i)+depk(1,i)
c        dep3bt(2,i)=dep3bt(2,i)+depi(2,i)+depk(2,i)
c        dep3bt(3,i)=dep3bt(3,i)+depi(3,i)+depk(3,i)
c      end do
c!$OMP END DO
c
c!$OMP END PARALLEL
c
c      do i=1,3
c        do j=1,3
c          virep3bt(j,i)=virep3bt(j,i)+virt(j,i)
c        end do
c      end do

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
      do l1 =1,npole3b
          i =pnum(l1)
c      do i =1,npole 
          deptemp(1,i)=deptemp(1,i)+deptemp2(1,l1)
          deptemp(2,i)=deptemp(2,i)+deptemp2(2,l1)
          deptemp(3,i)=deptemp(3,i)+deptemp2(3,l1)
      end do

      deallocate(nembedlst)
      deallocate(embedlst)
      deallocate(oldlistinds)
      return
      end

