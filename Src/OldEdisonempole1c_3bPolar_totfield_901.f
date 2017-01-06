      subroutine empole1c_3b_Polar_totfield_901(
     &   npole3b,pnum,uind,uinp,deptemp2,virtemp)
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
c      real*8 deptemp(3,npole)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
c      real*8, allocatable :: uind(:,:)
c      real*8, allocatable :: uinp(:,:)
      real*8, allocatable :: qgrid3b(:,:,:,:)
c      real*8, allocatable :: qfac3b(:,:,:)
c      real*8, allocatable :: thetai1_3b(:,:,:)
c      real*8, allocatable :: thetai2_3b(:,:,:)
c      real*8, allocatable :: thetai3_3b(:,:,:)
c      integer, allocatable :: pmetable3b(:,:)
c      integer, allocatable :: igrid3b(:,:)
      integer*8 planf3b
      integer*8 planb3b,k
      integer, allocatable :: iprime3b(:,:)
      real*8, allocatable :: ffttable3b(:,:)
c      integer, allocatable :: embedlst(:,:)
c      integer, allocatable :: nembedlst(:)
      real*8 r2,xr,yr,zr,ewaldcut2
c      integer embedlst(maxelst,9)
c      integer nembedlst(9),kk

c      allocate (nembedlst(npole3b))
c      allocate (embedlst(maxelst,npole3b))
c      if (.not. allocated(thetai1_3b)) 
c     &            allocate (thetai1_3b(4,bsorder,n))
c      if (.not. allocated(thetai2_3b))
c     &            allocate (thetai2_3b(4,bsorder,n))
c      if (.not. allocated(thetai3_3b))
c     &            allocate (thetai3_3b(4,bsorder,n))
      if (.not. allocated(qgrid3b)) 
     &           allocate (qgrid3b(2,nfft1,nfft2,nfft3))
c      if (.not. allocated(qfac3b))  allocate (qfac3b(nfft1,nfft2,nfft3))
c      if (.not. allocated(pmetable3b))  allocate (pmetable3b(n,nchunk))
c         if (.not. allocated(igrid3b))  allocate (igrid3b(3,n))
      
      if (.not. allocated(iprime3b)) allocate (iprime3b(maxprime,3))
      if (.not. allocated(ffttable3b)) allocate (ffttable3b(maxtable,3))
      
      call fftsetup3b(qgrid3b,planf3b,planb3b,iprime3b,ffttable3b)
c      call fftsetup3b(qgrid3b,planf3b,planb3b)

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


      !do i=1,npole
      !   do j=1,3
      !      deptemp(j,i) = 0.0d0
      !   end do
      !end do
      do i=1,3
         do j=1,3
           virtemp(j,i)=0.0d0
         end do
      end do

      f = electric / dielec
            

      !call induce0c_3b_new(npole3b,pnum,uind,uinp)
      ! call induce0a_3b(npole3b,pnum,uind,uinp,thetai1_3b,thetai2_3b,
     &!                 thetai3_3b,qgrid3b,qfac3b,pmetable3b,igrid3b,
     &!                 planf3b,planb3b,iprime3b,ffttable3b)   

       call induce0a_3b_totfield(npole3b,pnum,uind,uinp,qgrid3b,
     &                 planf3b,planb3b,iprime3b,ffttable3b)

c      do l1 = 1, npole3b
c        i=pnum(l1)
c         print*,"i=",i,"uindx,y,z=",uind(1,l1),uind(2,l1),uind(3,l1)
c      end do

       !print*,"Successful completion of TOTFIELDEwald Induce!"
      if(embedtyp.eq.'I') then
c      call erecip1_1_3b_totfieldnpole3b(npole3b,
c     &  pnum,uind,uinp,eptemp,deptemp2,virtemp)

c      call erecip1_2_3b_totfieldnpole3b(npole3b,
c     &  pnum,uind,uinp,deptemp,virtemp)


      else if (embedtyp.eq.'II') then
       !print*,"Embedtyp=",embedtyp
      !call erecip1_1_3b_totfieldnpole(npole3b,
     & ! pnum,uind,uinp,deptemp2,virtemp)
      !call erecip1_2_3b_totfieldnpole(npole3b,
     & ! pnum,uind,uinp,eptemp,deptemp,virtemp)
c      call emrecip1_1_3b_Polar_totfield(npole3b,pnum,uind,uinp,
c     &   deptemp2,virtemp,qgrid3b,qfac3b,thetai1_3b,thetai2_3b,
c     &   thetai3_3b,pmetable3b,igrid3b,planf3b,planb3b,iprime3b,
c     &   ffttable3b)

       call emrecip1_1_3b_Polar_totfield_parl_513(npole3b,pnum,
     &   uind,uinp,eptemp,deptemp2,virtemp,qgrid3b,
     &   planf3b,planb3b,
     &   iprime3b,ffttable3b)

c      call emrecip1_2_3b_Polar_totfield(npole3b,pnum,uind,uinp,
c     &   eptemp,deptemp,virtemp,qgrid3b,qfac3b,thetai1_3b,thetai2_3b,
c     &   thetai3_3b,pmetable3b,igrid3b,planf3b,planb3b,iprime3b,
c     &   ffttable3b)
      !call emrecip1_1and2_3b_Polar_totfield(npole3b,pnum,uind,
     &! uinp,eptemp,deptemp,virtemp,qgrid3b,qfac3b,thetai1_3b,thetai2_3b,
     &!   thetai3_3b,pmetable3b,igrid3b,planf3b,planb3b,iprime3b,
     &!   ffttable3b)
      !call emrecip1_1and2_3b_Polar_totfield_parl(npole3b,pnum,
     & !uind,uinp,eptemp,deptemp,virtemp,qgrid3b,qfac3b,thetai1_3b,
     & !thetai2_3b,thetai3_3b,pmetable3b,igrid3b,planf3b,planb3b,
     & !  iprime3b,ffttable3b)

      !call emrecip1_1and2_3b_Polar_totfield_parl_513(npole3b,pnum,
     & !  uind,uinp,eptemp,deptemp,virtemp,qgrid3b,
     & !  planf3b,planb3b,
     & !  iprime3b,ffttable3b)

      else
         !call erecip1_3b(npole3b,
     &   ! pnum,uind,uinp,eptemp,deptemp,virtemp)
         print*,"IIIo:  Just finished erecip1_3b"
      end if
ccc      call induce0a_3b(npole3b,pnum,uind,uinp,qgrid,qfac)

      !call emrecip1_3b_Polar(npole3b,pnum,uind,uinp,
     &!   eptemp,deptemp,virtemp,qgrid3b,qfac3b,thetai1_3b,thetai2_3b,
     &!   thetai3_3b,pmetable3b,igrid3b,planf3b,planb3b,iprime3b,
     &!   ffttable3b)


      deallocate(qgrid3b)
c      deallocate(qfac3b)
c      deallocate(pmetable3b)
c      deallocate(igrid3b)
c      deallocate(thetai1_3b)
c      deallocate(thetai2_3b)
c      deallocate(thetai3_3b)
      deallocate(iprime3b)
      deallocate(ffttable3b)

      
      if(embedtyp.eq.'I') then
c       print*,"Before ereal1c1_3b_totfieldnpole3b"
c       call ereal1c1_3b_totfieldnpole3b(npole3b,pnum,uind,uinp,
c     &   eptemp,deptemp2,virtemp)
c       print*,"Successful completion of ereal1c1_3b_totfieldnpole3b"
c       print*,"Before ereal1c2_3b_totfieldnpole3b"
c       call ereal1c2_3b_totfieldnpole3b(npole3b,pnum,uind,uinp,
c     &  deptemp,virtemp)
       !print*,"Successful completion of ereal1c2_3b_totfieldnpole3b"
      else if(embedtyp.eq.'II') then
       call ereal1c1_3b_totfieldnpole(npole3b,
     &  pnum,uind,uinp,deptemp2,virtemp)
       !call ereal1c2_3b_totfieldnpole(npole3b,
     & ! pnum,uind,uinp,eptemp,deptemp,virtemp)
        
       
       !call ereal1d2_3b_totfieldnpole(npole3b,
     & ! pnum,uind,uinp,eptemp,deptemp,virtemp,nembedlst,embedlst)

      else

      end if

      !do l1 =1,npole3b
      !    i =pnum(l1)
      !    deptemp(1,i)=deptemp2(1,l1)+deptemp(1,i)
      !    deptemp(2,i)=deptemp2(2,l1)+deptemp(2,i)
      !    deptemp(3,i)=deptemp2(3,l1)+deptemp(3,i)
      !end do

      
      !aewald3b=aewald
      !term = 2.0d0 * aewald3b * aewald3b
      !fterm = -f * aewald3b / sqrtpi
      !do l1 = 1, npole3b
      !   i=pnum(l1)
      !   dix = rpole(2,i)
      !   diy = rpole(3,i)
      !   diz = rpole(4,i)
      !   uix = uind(1,l1)
      !   uiy = uind(2,l1)
      !   uiz = uind(3,l1)
      !   uii = dix*uix + diy*uiy + diz*uiz
      !   ei = fterm * term * uii / 3.0d0
      !   eptemp = eptemp + ei
      !end do

c
c     compute the self-energy torque term due to induced dipole
c

      !trq(1) = 0.0d0
      !trq(2) = 0.0d0
      !trq(3) = 0.0d0
      !term = (4.0d0/3.0d0) * f * aewald3b**3 / sqrtpi
      !do l1 = 1, npole3b
      !   i=pnum(l1)
      !   dix = rpole(2,i)
      !   diy = rpole(3,i)
      !   diz = rpole(4,i)
      !   uix = 0.5d0 * (uind(1,l1)+uinp(1,l1))
      !   uiy = 0.5d0 * (uind(2,l1)+uinp(2,l1))
      !   uiz = 0.5d0 * (uind(3,l1)+uinp(3,l1))
      !   trqi(1) = term * (diy*uiz-diz*uiy)
      !   trqi(2) = term * (diz*uix-dix*uiz)
      !   trqi(3) = term * (dix*uiy-diy*uix)
cccc         call torque_3b (deptemp,i,trq,trqi,frcx,frcy,frcz)
c         call torque_3b_new(npole3b,pnum,deptemp,i,
c     &    trq,trqi,frcx,frcy,frcz)
      !   call torque_3b_new_npole(deptemp,i,
     &!      trq,trqi,frcx,frcy,frcz)
      !end do


      if (boundary .eq. 'VACUUM') then
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

