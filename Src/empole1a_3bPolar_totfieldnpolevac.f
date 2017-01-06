c
c
      subroutine empole1a_3b_Polar_totfieldnpolevac(npole3b,pnum,
     &  npole3b1,pnum1,pnum1ref,eptemp1,
     &   deptemp1,virtemp1,field,fieldp)
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
      use totfield
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
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      integer npole3b,npole3b1,pnum(*),pnum1(*),pnum1ref(*)   
      real*8 eptemp1,deptemp1(3,npole3b),virtemp1(3,3)
c      real*8 uind(3,npole3b)
c      real*8 uinp(3,npole3b)
      real*8, allocatable :: uind(:,:)
      real*8, allocatable :: uinp(:,:)
      real*8 off3b,eptemp_userep
      integer l1,l3,l2
      logical proceed,usei,usek
      character*6 mode
      real*8, allocatable :: deptemp2(:,:)
      real*8 field(3,npole3b)
      real*8 fieldp(3,npole3b)

c
c   Zero out temporary gradient of polarization energy
c

      allocate (deptemp2(3,npole3b1))
c      allocate (uind(3,npole))
c      allocate (uinp(3,npole))
      allocate (uind(3,npole3b1))
      allocate (uinp(3,npole3b1))

      eptemp1 = 0.0d0

      do l1=1,npole3b
         do j = 1,3
            deptemp1(j,l1)=0.0d0
         end do
      end do

      !do l1=1,npole3b1
      !   do j = 1,3
      !      deptemp2(j,l1)=0.0d0
      !   end do
      !end do

c
c   Zero out temporary virial
c
      do i=1,3
         do j=1,3
           virtemp1(j,i)=0.0d0
         end do
      end do

c      mode = 'MPOLE'
c      call switch (mode)

c      print*,"off2,use_replica = ",off2,use_replica
      

c      call induce0a_3b_scf_totfieldnpole(npole3b,pnum,uind,uinp) 

      call induce0a_3b_PolelecOnly_totfieldnpolevac(npole3b,pnum,
     &     npole3b1,pnum1,pnum1ref,uind,uinp,field,fieldp)

      call empole1a1_3b_Polar_totfieldnpole_mod(npole3b1,pnum1,uind,
     &     uinp,deptemp2,virtemp1)

      call empole1a2_3b_Polar_totfieldnpolevac(npole3b,pnum,npole3b1,
     &     pnum1,pnum1ref,uind,uinp,eptemp1,deptemp1,virtemp1)

c      call induce0a_3b_PolelecOnly_totfieldnpole3(npole3b,pnum,uind,
c     &     uinp)
c
c      call empole1a1_3b_Polar_totfieldnpole_mod3(npole3b,pnum,uind,
c     &     uinp,deptemp2,virtemp)
c
c      call empole1a2_3b_Polar_totfieldnpole_mod3(npole3b,pnum,uind,uinp,
c     &     eptemp,deptemp,virtemp)
      

      do l1 =1,npole3b1
c          i =pnum(l1)
c      do i =1,npole 
          l2=pnum1ref(l1)
          deptemp1(1,l2)=deptemp1(1,l2)+deptemp2(1,l1)
          deptemp1(2,l2)=deptemp1(2,l2)+deptemp2(2,l1)
          deptemp1(3,l2)=deptemp1(3,l2)+deptemp2(3,l1)
      end do

      deallocate(deptemp2)
      deallocate(uind)
      deallocate(uinp)
      return
      end

