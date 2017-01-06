c      subroutine empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
c     &  dep3bt,virep3bt)
      subroutine Newempole1c_3b_Polar_totfield_dEtensor(npole3b,pnum,
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
      integer m,k1,l2,length,start,last,moli1rmndr
      real*8 termx,termy,termz
      real*8 taux,tauy,tauz
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

          xr = x(k) - x(i)
          yr = y(k) - y(i)
          zr = z(k) - z(i)
          call image (xr,yr,zr)
          r2=xr*xr + yr*yr +zr*zr
c             NEED TO CORRECT STUFF BELOW
            if((i.lt.k).and.(elstgrad(kkk,i).ne.0)) then
              m=elstgrad(kkk,i)
            else if ((i.gt.k).and.(elstgrad(kkk,k).ne.0)) then
              m=elstgrad(kkk,k)
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

          if ((r2 .le. ewaldcut2).and.(m.ne.0)) then

         vxx=0.0d0
         vyx=0.0d0
         vzx=0.0d0
         vyy=0.0d0
         vzy=0.0d0
         vzz=0.0d0


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

          frcz2xxd=frcztau2dtot(1,m)
          frcz2xyd=frcztau2dtot(2,m)
          frcz2xzd=frcztau2dtot(3,m)
          frcz2yxd=frcztau2dtot(4,m)
          frcz2yyd=frcztau2dtot(5,m)
          frcz2yzd=frcztau2dtot(6,m)
          frcz2zxd=frcztau2dtot(7,m)
          frcz2zyd=frcztau2dtot(8,m)
          frcz2zzd=frcztau2dtot(9,m)

          frcz2xxp=frcztau2ptot(1,m)
          frcz2xyp=frcztau2ptot(2,m)
          frcz2xzp=frcztau2ptot(3,m)
          frcz2yxp=frcztau2ptot(4,m)
          frcz2yyp=frcztau2ptot(5,m)
          frcz2yzp=frcztau2ptot(6,m)
          frcz2zxp=frcztau2ptot(7,m)
          frcz2zyp=frcztau2ptot(8,m)
          frcz2zzp=frcztau2ptot(9,m)

          frcy2xxd=frcytau2dtot(1,m)
          frcy2xyd=frcytau2dtot(2,m)
          frcy2xzd=frcytau2dtot(3,m)
          frcy2yxd=frcytau2dtot(4,m)
          frcy2yyd=frcytau2dtot(5,m)
          frcy2yzd=frcytau2dtot(6,m)
          frcy2zxd=frcytau2dtot(7,m)
          frcy2zyd=frcytau2dtot(8,m)
          frcy2zzd=frcytau2dtot(9,m)

          frcy2xxp=frcytau2ptot(1,m)
          frcy2xyp=frcytau2ptot(2,m)
          frcy2xzp=frcytau2ptot(3,m)
          frcy2yxp=frcytau2ptot(4,m)
          frcy2yyp=frcytau2ptot(5,m)
          frcy2yzp=frcytau2ptot(6,m)
          frcy2zxp=frcytau2ptot(7,m)
          frcy2zyp=frcytau2ptot(8,m)
          frcy2zzp=frcytau2ptot(9,m)

          frcx2xxd=frcxtau2dtot(1,m)
          frcx2xyd=frcxtau2dtot(2,m)
          frcx2xzd=frcxtau2dtot(3,m)
          frcx2yxd=frcxtau2dtot(4,m)
          frcx2yyd=frcxtau2dtot(5,m)
          frcx2yzd=frcxtau2dtot(6,m)
          frcx2zxd=frcxtau2dtot(7,m)
          frcx2zyd=frcxtau2dtot(8,m)
          frcx2zzd=frcxtau2dtot(9,m)

          frcx2xxp=frcxtau2ptot(1,m)
          frcx2xyp=frcxtau2ptot(2,m)
          frcx2xzp=frcxtau2ptot(3,m)
          frcx2yxp=frcxtau2ptot(4,m)
          frcx2yyp=frcxtau2ptot(5,m)
          frcx2yzp=frcxtau2ptot(6,m)
          frcx2zxp=frcxtau2ptot(7,m)
          frcx2zyp=frcxtau2ptot(8,m)
          frcx2zzp=frcxtau2ptot(9,m)

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

          frcz1xxd=frcztau1dtot(1,m)
          frcz1xyd=frcztau1dtot(2,m)
          frcz1xzd=frcztau1dtot(3,m)
          frcz1yxd=frcztau1dtot(4,m)
          frcz1yyd=frcztau1dtot(5,m)
          frcz1yzd=frcztau1dtot(6,m)
          frcz1zxd=frcztau1dtot(7,m)
          frcz1zyd=frcztau1dtot(8,m)
          frcz1zzd=frcztau1dtot(9,m)

          frcz1xxp=frcztau1ptot(1,m)
          frcz1xyp=frcztau1ptot(2,m)
          frcz1xzp=frcztau1ptot(3,m)
          frcz1yxp=frcztau1ptot(4,m)
          frcz1yyp=frcztau1ptot(5,m)
          frcz1yzp=frcztau1ptot(6,m)
          frcz1zxp=frcztau1ptot(7,m)
          frcz1zyp=frcztau1ptot(8,m)
          frcz1zzp=frcztau1ptot(9,m)

          frcx1xxd=frcxtau1dtot(1,m)
          frcx1xyd=frcxtau1dtot(2,m)
          frcx1xzd=frcxtau1dtot(3,m)
          frcx1yxd=frcxtau1dtot(4,m)
          frcx1yyd=frcxtau1dtot(5,m)
          frcx1yzd=frcxtau1dtot(6,m)
          frcx1zxd=frcxtau1dtot(7,m)
          frcx1zyd=frcxtau1dtot(8,m)
          frcx1zzd=frcxtau1dtot(9,m)

          frcx1xxp=frcxtau1ptot(1,m)
          frcx1xyp=frcxtau1ptot(2,m)
          frcx1xzp=frcxtau1ptot(3,m)
          frcx1yxp=frcxtau1ptot(4,m)
          frcx1yyp=frcxtau1ptot(5,m)
          frcx1yzp=frcxtau1ptot(6,m)
          frcx1zxp=frcxtau1ptot(7,m)
          frcx1zyp=frcxtau1ptot(8,m)
          frcx1zzp=frcxtau1ptot(9,m)

          frcy1xxd=frcytau1dtot(1,m)
          frcy1xyd=frcytau1dtot(2,m)
          frcy1xzd=frcytau1dtot(3,m)
          frcy1yxd=frcytau1dtot(4,m)
          frcy1yyd=frcytau1dtot(5,m)
          frcy1yzd=frcytau1dtot(6,m)
          frcy1zxd=frcytau1dtot(7,m)
          frcy1zyd=frcytau1dtot(8,m)
          frcy1zzd=frcytau1dtot(9,m)

          frcy1xxp=frcytau1ptot(1,m)
          frcy1xyp=frcytau1ptot(2,m)
          frcy1xzp=frcytau1ptot(3,m)
          frcy1yxp=frcytau1ptot(4,m)
          frcy1yyp=frcytau1ptot(5,m)
          frcy1yzp=frcytau1ptot(6,m)
          frcy1zxp=frcytau1ptot(7,m)
          frcy1zyp=frcytau1ptot(8,m)
          frcy1zzp=frcytau1ptot(9,m)

c          termx=uind(1,i)*t2xxp+ uind(2,i)*t2xyp + uind(3,i)*t2xzp
c          termy=uind(1,i)*t2yxp+ uind(2,i)*t2yyp + uind(3,i)*t2yzp
c          termz=uind(1,i)*t2zxp+ uind(2,i)*t2zyp + uind(3,i)*t2zzp
c
c       frczk(1)=uind(1,i)*frcz2xxp+uind(2,i)*frcz2xyp+uind(3,i)*frcz2xzp
c       frczk(2)=uind(1,i)*frcz2yxp+uind(2,i)*frcz2yyp+uind(3,i)*frcz2yzp
c       frczk(3)=uind(1,i)*frcz2zxp+uind(2,i)*frcz2zyp+uind(3,i)*frcz2zzp
c
c       frcxk(1)=uind(1,i)*frcx2xxp+uind(2,i)*frcx2xyp+uind(3,i)*frcx2xzp
c       frcxk(2)=uind(1,i)*frcx2yxp+uind(2,i)*frcx2yyp+uind(3,i)*frcx2yzp
c       frcxk(3)=uind(1,i)*frcx2zxp+uind(2,i)*frcx2zyp+uind(3,i)*frcx2zzp
c
c       frcyk(1)=uind(1,i)*frcy2xxp+uind(2,i)*frcy2xyp+uind(3,i)*frcy2xzp
c       frcyk(2)=uind(1,i)*frcy2yxp+uind(2,i)*frcy2yyp+uind(3,i)*frcy2yzp
c       frcyk(3)=uind(1,i)*frcy2zxp+uind(2,i)*frcy2zyp+uind(3,i)*frcy2zzp
c
c          depi(1,i)=depi(1,i)+termx
c          depi(2,i)=depi(2,i)+termy
c          depi(3,i)=depi(3,i)+termz
c
c          depk(1,k)=depk(1,k)-termx
c          depk(2,k)=depk(2,k)-termy
c          depk(3,k)=depk(3,k)-termz
c
c          vxx =vxx -xr*termx +xkx*frcxk(1)
c     &                  + xky*frcyk(1) + xkz*frczk(1)
c          vyx =vyx -yr*termx + ykx*frcxk(1)
c     &                  + yky*frcyk(1) + ykz*frczk(1) 
c          vzx =vzx -zr*termx + zkx*frcxk(1)
c     &                  + zky*frcyk(1) + zkz*frczk(1)
c          vyy =vyy -yr*termy + ykx*frcxk(2)
c     &                  + yky*frcyk(2) + ykz*frczk(2)
c          vzy =vzy -zr*termy + zkx*frcxk(2)
c     &                  + zky*frcyk(2) + zkz*frczk(2)
c          vzz =vzz -zr*termz + zkx*frcxk(3)
c     &                  + zky*frcyk(3) + zkz*frczk(3)
c
c          termx = uinp(1,i)*t2xxd+ uinp(2,i)*t2xyd + uinp(3,i)*t2xzd
c          termy = uinp(1,i)*t2yxd+ uinp(2,i)*t2yyd + uinp(3,i)*t2yzd
c          termz = uinp(1,i)*t2zxd+ uinp(2,i)*t2zyd + uinp(3,i)*t2zzd
c
c       frczk(1)=uinp(1,i)*frcz2xxd+uinp(2,i)*frcz2xyd+uinp(3,i)*frcz2xzd
c       frczk(2)=uinp(1,i)*frcz2yxd+uinp(2,i)*frcz2yyd+uinp(3,i)*frcz2yzd
c       frczk(3)=uinp(1,i)*frcz2zxd+uinp(2,i)*frcz2zyd+uinp(3,i)*frcz2zzd
c          
c       frcxk(1)=uinp(1,i)*frcx2xxd+uinp(2,i)*frcx2xyd+uinp(3,i)*frcx2xzd
c       frcxk(2)=uinp(1,i)*frcx2yxd+uinp(2,i)*frcx2yyd+uinp(3,i)*frcx2yzd
c       frcxk(3)=uinp(1,i)*frcx2zxd+uinp(2,i)*frcx2zyd+uinp(3,i)*frcx2zzd
c
c       frcyk(1)=uinp(1,i)*frcy2xxd+uinp(2,i)*frcy2xyd+uinp(3,i)*frcy2xzd
c       frcyk(2)=uinp(1,i)*frcy2yxd+uinp(2,i)*frcy2yyd+uinp(3,i)*frcy2yzd
c       frcyk(3)=uinp(1,i)*frcy2zxd+uinp(2,i)*frcy2zyd+uinp(3,i)*frcy2zzd
c
c          depi(1,i)=depi(1,i)+termx
c          depi(2,i)=depi(2,i)+termy
c          depi(3,i)=depi(3,i)+termz
c          
c          depk(1,k)=depk(1,k)-termx
c          depk(2,k)=depk(2,k)-termy
c          depk(3,k)=depk(3,k)-termz
c          vxx =vxx -xr*termx +xkx*frcxk(1)
c     &                  + xky*frcyk(1) + xkz*frczk(1)
c          vyx =vyx -yr*termx + ykx*frcxk(1)
c     &                  + yky*frcyk(1) + ykz*frczk(1)
c          vzx =vzx -zr*termx + zkx*frcxk(1)
c     &                  + zky*frcyk(1) + zkz*frczk(1)
c          vyy =vyy -yr*termy + ykx*frcxk(2)
c     &                  + yky*frcyk(2) + ykz*frczk(2)
c          vzy =vzy -zr*termy + zkx*frcxk(2)
c     &                  + zky*frcyk(2) + zkz*frczk(2)
c          vzz =vzz -zr*termz + zkx*frcxk(3)
c     &                  + zky*frcyk(3) + zkz*frczk(3)

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
c             if(dok) then
c             print*,"Doing k!  Should not be happening!"
c             termx=uind(1,l3)*t2xxp+ uind(2,l3)*t2xyp + uind(3,l3)*t2xzp
c             termy=uind(1,l3)*t2yxp+ uind(2,l3)*t2yyp + uind(3,l3)*t2yzp
c             termz=uind(1,l3)*t2zxp+ uind(2,l3)*t2zyp + uind(3,l3)*t2zzp
c
c             frczi(1)=uind(1,l3)*frcz2xxp+uind(2,l3)*frcz2xyp
c     &                +uind(3,l3)*frcz2xzp
c             frczi(2)=uind(1,l3)*frcz2yxp+uind(2,l3)*frcz2yyp
c     &                +uind(3,l3)*frcz2yzp
c             frczi(3)=uind(1,l3)*frcz2zxp+uind(2,l3)*frcz2zyp
c     &                +uind(3,l3)*frcz2zzp
c
c             frcxi(1)=uind(1,l3)*frcx2xxp+uind(2,l3)*frcx2xyp
c     &                +uind(3,l3)*frcx2xzp
c             frcxi(2)=uind(1,l3)*frcx2yxp+uind(2,l3)*frcx2yyp
c     &                +uind(3,l3)*frcx2yzp
c             frcxi(3)=uind(1,l3)*frcx2zxp+uind(2,l3)*frcx2zyp
c     &                +uind(3,l3)*frcx2zzp
c
c             frcyi(1)=uind(1,l3)*frcy2xxp+uind(2,l3)*frcy2xyp
c     &                +uind(3,l3)*frcy2xzp
c             frcyi(2)=uind(1,l3)*frcy2yxp+uind(2,l3)*frcy2yyp
c     &                +uind(3,l3)*frcy2yzp
c             frcyi(3)=uind(1,l3)*frcy2zxp+uind(2,l3)*frcy2zyp
c     &                +uind(3,l3)*frcy2zzp
c
cc             deptemp(1,i)=deptemp(1,i)+termx
cc             deptemp(2,i)=deptemp(2,i)+termy
cc             deptemp(3,i)=deptemp(3,i)+termz
cc
cc             deptemp(1,k)=deptemp(1,k)-termx
cc             deptemp(2,k)=deptemp(2,k)-termy
cc             deptemp(3,k)=deptemp(3,k)-termz
c
c             deptemp(1,i)=deptemp(1,i)-termx
c             deptemp(2,i)=deptemp(2,i)-termy
c             deptemp(3,i)=deptemp(3,i)-termz
c
c             deptemp(1,k)=deptemp(1,k)+termx
c             deptemp(2,k)=deptemp(2,k)+termy
c             deptemp(3,k)=deptemp(3,k)+termz
c
c
c             vxx=vxx-xr*termx+ xix*frcxi(1)
c     &                  + xiy*frcyi(1) + xiz*frczi(1)
c             vyx=vyx-yr*termx+ yix*frcxi(1)
c     &                  + yiy*frcyi(1) + yiz*frczi(1)
c             vzx=vzx-zr*termx + zix*frcxi(1)
c     &                  + ziy*frcyi(1) + ziz*frczi(1)
c             vyy=vyy-yr*termy + yix*frcxi(2)
c     &                  + yiy*frcyi(2) + yiz*frczi(2)
c             vzy=vzy-zr*termy + zix*frcxi(2)
c     &                  + ziy*frcyi(2) + ziz*frczi(2)
c             vzz=vzz-zr*termz + zix*frcxi(3)
c     &                  + ziy*frcyi(3) + ziz*frczi(3)
c             termx=uinp(1,l3)*t2xxd+ uinp(2,l3)*t2xyd + uinp(3,l3)*t2xzd
c             termy=uinp(1,l3)*t2yxd+ uinp(2,l3)*t2yyd + uinp(3,l3)*t2yzd
c             termz=uinp(1,l3)*t2zxd+ uinp(2,l3)*t2zyd + uinp(3,l3)*t2zzd
c
c             frczi(1)=uinp(1,l3)*frcz2xxd+uinp(2,l3)*frcz2xyd
c     &                +uinp(3,l3)*frcz2xzd
c             frczi(2)=uinp(1,l3)*frcz2yxd+uinp(2,l3)*frcz2yyd
c     &                +uinp(3,l3)*frcz2yzd
c             frczi(3)=uinp(1,l3)*frcz2zxd+uinp(2,l3)*frcz2zyd
c     &                +uinp(3,l3)*frcz2zzd
c
c             frcxi(1)=uinp(1,l3)*frcx2xxd+uinp(2,l3)*frcx2xyd
c     &                +uinp(3,l3)*frcx2xzd
c             frcxi(2)=uinp(1,l3)*frcx2yxd+uinp(2,l3)*frcx2yyd
c     &                +uinp(3,l3)*frcx2yzd
c             frcxi(3)=uinp(1,l3)*frcx2zxd+uinp(2,l3)*frcx2zyd
c     &                +uinp(3,l3)*frcx2zzd
c
c             frcyi(1)=uinp(1,l3)*frcy2xxd+uinp(2,l3)*frcy2xyd
c     &                +uinp(3,l3)*frcy2xzd
c             frcyi(2)=uinp(1,l3)*frcy2yxd+uinp(2,l3)*frcy2yyd
c     &                +uinp(3,l3)*frcy2yzd
c             frcyi(3)=uinp(1,l3)*frcy2zxd+uinp(2,l3)*frcy2zyd
c     &                +uinp(3,l3)*frcy2zzd
c
cc             deptemp(1,i)=deptemp(1,i)+termx
cc             deptemp(2,i)=deptemp(2,i)+termy
cc             deptemp(3,i)=deptemp(3,i)+termz
cc
cc             deptemp(1,k)=deptemp(1,k)-termx
cc             deptemp(2,k)=deptemp(2,k)-termy
cc             deptemp(3,k)=deptemp(3,k)-termz
c
c             deptemp(1,i)=deptemp(1,i)-termx
c             deptemp(2,i)=deptemp(2,i)-termy
c             deptemp(3,i)=deptemp(3,i)-termz
c
c             deptemp(1,k)=deptemp(1,k)+termx
c             deptemp(2,k)=deptemp(2,k)+termy
c             deptemp(3,k)=deptemp(3,k)+termz
c
c             vxx=vxx-xr*termx+ xix*frcxi(1)
c     &                  + xiy*frcyi(1) + xiz*frczi(1)
c             vyx=vyx-yr*termx+ yix*frcxi(1)
c     &                  + yiy*frcyi(1) + yiz*frczi(1)
c             vzx=vzx-zr*termx + zix*frcxi(1)
c     &                  + ziy*frcyi(1) + ziz*frczi(1)
c             vyy=vyy-yr*termy + yix*frcxi(2)
c     &                  + yiy*frcyi(2) + yiz*frczi(2)
c             vzy=vzy-zr*termy + zix*frcxi(2)
c     &                  + ziy*frcyi(2) + ziz*frczi(2)
c             vzz=vzz-zr*termz + zix*frcxi(3)
c     &                  + ziy*frcyi(3) + ziz*frczi(3)
c             end if
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
            if(k1.ne.0) then
             m1=elsttau1(j1,kkk,i)
            end if
            k2=elsttau2_index(j1,kkk,i)
            if(k2.ne.0) then
             m2=elsttau2(j1,kkk,i)
            end if
          else
            k1=elsttau1_index(j1,kkk,k)
            if(k1.ne.0) then
             m1=elsttau1(j1,kkk,k)
            end if
            k2=elsttau2_index(j1,kkk,k)
            if(k2.ne.0) then 
             m2=elsttau2(j1,kkk,k)
            end if
          end if         

          if(k1.ne.0) then 
          tau1xxd=taud1(1,m1)
          tau1xyd=taud1(2,m1)
          tau1xzd=taud1(3,m1)
          tau1yxd=taud1(4,m1)
          tau1yyd=taud1(5,m1)
          tau1yzd=taud1(6,m1)
          tau1zxd=taud1(7,m1)
          tau1zyd=taud1(8,m1)
          tau1zzd=taud1(9,m1)

          tau1xxp=taup1(1,m1)
          tau1xyp=taup1(2,m1)
          tau1xzp=taup1(3,m1)
          tau1yxp=taup1(4,m1)
          tau1yyp=taup1(5,m1)
          tau1yzp=taup1(6,m1)
          tau1zxp=taup1(7,m1)
          tau1zyp=taup1(8,m1)
          tau1zzp=taup1(9,m1)
          end if

          if(k2.ne.0) then
          tau2xxd=taud2(1,m2)
          tau2xyd=taud2(2,m2)
          tau2xzd=taud2(3,m2)
          tau2yxd=taud2(4,m2)
          tau2yyd=taud2(5,m2)
          tau2yzd=taud2(6,m2)
          tau2zxd=taud2(7,m2)
          tau2zyd=taud2(8,m2)
          tau2zzd=taud2(9,m2)

          tau2xxp=taup2(1,m2)
          tau2xyp=taup2(2,m2)
          tau2xzp=taup2(3,m2)
          tau2yxp=taup2(4,m2)
          tau2yyp=taup2(5,m2)
          tau2yzp=taup2(6,m2)
          tau2zxp=taup2(7,m2)
          tau2zyp=taup2(8,m2)
          tau2zzp=taup2(9,m2)
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

