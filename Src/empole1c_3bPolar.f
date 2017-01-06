      subroutine empole1c_3b_Polar(npole3b,pnum,eptemp,
     &   deptemp,virtemp)
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
      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      integer npole3b,pnum(*)
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
      integer*8 planb3b
      integer, allocatable :: iprime3b(:,:)
      real*8, allocatable :: ffttable3b(:,:)

      if (.not. allocated(thetai1_3b)) 
     &            allocate (thetai1_3b(4,bsorder,n))
      if (.not. allocated(thetai2_3b))
     &            allocate (thetai2_3b(4,bsorder,n))
      if (.not. allocated(thetai3_3b))
     &            allocate (thetai3_3b(4,bsorder,n))
      if (.not. allocated(qgrid3b)) 
     &           allocate (qgrid3b(2,nfft1,nfft2,nfft3))
      if (.not. allocated(qfac3b))  allocate (qfac3b(nfft1,nfft2,nfft3))
      if (.not. allocated(pmetable3b))  allocate (pmetable3b(n,nchunk))
         if (.not. allocated(igrid3b))  allocate (igrid3b(3,n))
      
      if (.not. allocated(iprime3b)) allocate (iprime3b(maxprime,3))
      if (.not. allocated(ffttable3b)) allocate (ffttable3b(maxtable,3))
      
      call fftsetup3b(qgrid3b,planf3b,planb3b,iprime3b,ffttable3b)
c      call fftsetup3b(qgrid3b,planf3b,planb3b)
      eptemp = 0.0d0

c   Zero out temporary gradient of polarization energy
      do l1 = 1, npole3b
c        i = pnum(l1)
        do j = 1, 3
          deptemp(j,l1) = 0.0d0
          uind(j,l1) = 0.0d0
          uinp(j,l1) = 0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
           virtemp(j,i)=0.0d0
         end do
      end do

      f = electric / dielec
            


       ! call induce0c_3b_new(npole3b,pnum,uind,uinp)
       call induce0a_3b(npole3b,pnum,uind,uinp,thetai1_3b,thetai2_3b,
     &                  thetai3_3b,qgrid3b,qfac3b,pmetable3b,igrid3b,
     &                  planf3b,planb3b,iprime3b,ffttable3b)   
c       print*,"Successful completion of Ewald Induce!"
       !call erecip1_3b(npole3b,
     & !  pnum,uind,uinp,eptemp,deptemp,virtemp)

c      call induce0a_3b(npole3b,pnum,uind,uinp,qgrid,qfac)

       call emrecip1_3b_Polar(npole3b,pnum,uind,uinp,
     &    eptemp,deptemp,virtemp,qgrid3b,qfac3b,thetai1_3b,thetai2_3b,
     &    thetai3_3b,pmetable3b,igrid3b,planf3b,planb3b,iprime3b,
     &    ffttable3b)


      deallocate(qgrid3b)
      deallocate(qfac3b)
      deallocate(pmetable3b)
      deallocate(igrid3b)
      deallocate(thetai1_3b)
      deallocate(thetai2_3b)
      deallocate(thetai3_3b)
      deallocate(iprime3b)
      deallocate(ffttable3b)

      call ereal1c_3b(npole3b,pnum,uind,uinp,
     &            eptemp,deptemp,virtemp)
      
      aewald3b=aewald
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
c         call torque_3b (deptemp,i,trq,trqi,frcx,frcy,frcz)
          call torque_3b_new(npole3b,pnum,deptemp,i,
     &    trq,trqi,frcx,frcy,frcz)
       
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

c      subroutine erecip1_3b(aewald3b3b,ewaldcut3b,jmax3b3b,kmax3b3b,lmax3b3b,
c     &   npole3b,pnum,uind,uinp,eptemp,deptemp,virtemp)
      subroutine erecip1_3b(npole3b,pnum,uind,uinp,eptemp,deptemp,
     & virtemp)
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'boxes.i'
c      include 'chgpot.i'
c      include 'ewald.i'
c      include 'ewreg2.i'
c      include 'math.i'
c      include 'mpole.i'
c      include 'polar2.i'
c      include 'polpot.i'
c      include 'units.i'
c      include 'virial.i'
      use sizes
      use atoms
      use chgpot
      use mpole
      use polar, only: polarity, thole, pdamp
      use boxes
      use ewald
      use ewreg3bpolz
      use math 
      implicit none
      integer i,j,k,l,l1
      integer ii,m1,m2
      integer jmin
      integer kmin
      integer lmin
      real*8 e,ei,etot,f,cut
      real*8 expterm,term,eterm
      real*8 uterm,vterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 t1,t2,t3,t4
      real*8 ukp,t3p,t4p
      real*8 de,det1,det2
      real*8 dei,det1i,det2i
      real*8 wterm(3,3)
      real*8 t5(3,3),t6(3,3)
      real*8 t5u(3,3),t5p(3,3)
      real*8 t6u(3,3),t6p(3,3)
      real*8 qt(3,3),dt(3,3)
      real*8 dtu(3,3),dtp(3,3)
c      real*8 ckr(maxatm),skr(maxatm)
c      real*8 cjk(maxatm),sjk(maxatm)
c      real*8 s1(maxatm),s2(maxatm)
c      real*8 s3(maxatm),s4(maxatm)
c      real*8 s3p(maxatm),s4p(maxatm)
c      real*8 cm(maxatm),dm(3,maxatm)
c      real*8 qm(9,maxatm),um(3,maxatm)
c      real*8 trq(3,maxatm),trqi(3,maxatm)
c      real*8 dkx(maxatm),qkx(maxatm)
c      real*8 dky(maxatm),qky(maxatm)
c      real*8 dkz(maxatm),qkz(maxatm)
c      integer npole3b,pnum(*),pnum2(*)
      integer npole3b,pnum(*)
      real*8 ckr(npole3b),skr(npole3b)
      real*8 cjk(npole3b),sjk(npole3b)
      real*8 s1(npole3b),s2(npole3b)
      real*8 s3(npole3b),s4(npole3b)
      real*8 s3p(npole3b),s4p(npole3b)
      real*8 cm(npole3b),dm(3,npole3b)
      real*8 qm(9,npole3b),um(3,npole3b)
c      real*8 trq(3,npole3b),trqi(3,npole3b)
c      real*8 trq(3,npole),trqi(3,npole)
      real*8 trqi(3,npole3b)
      real*8 dkx(npole3b),qkx(npole3b)
      real*8 dky(npole3b),qky(npole3b)
      real*8 dkz(npole3b),qkz(npole3b)
      real*8 ejc(npole3b,0:maxvec)
      real*8 ejs(npole3b,0:maxvec)
      real*8 ekc(npole3b,-maxvec:maxvec)
      real*8 eks(npole3b,-maxvec:maxvec)
      real*8 elc(npole3b,-maxvec:maxvec)
      real*8 els(npole3b,-maxvec:maxvec)
      real*8 uind(3,*),eptemp,deptemp(3,*),virtemp(3,3)
      real*8 frecip,uinp(3,*)
c      real*8 aewald3b,aewald3b3b,ewaldcut3b
c      integer jmax3b3b,kmax3b3b,lmax3b3b
c
c     return if the Ewald coefficient is zero
c
c      aewald3b=aewald3b3b

      if (aewald3b .lt. 1.0d-6)  return

      frecip = 0.5d0
      f = electric / dielec
      term = -0.25d0 / aewald3b**2
      eterm = 4.0d0 * pi * f / volbox
      !print*,"aewald3b=",aewald3b
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1


c      print*,"jmax3b,kmax3b,lmax3b in erecip1_3b",jmax3b,kmax3b,lmax3b

      do l1 = 1, npole3b
         i = pnum(l1)
         cm(l1) = rpole(1,i)
         dm(1,l1) = rpole(2,i)
         dm(2,l1) = rpole(3,i)
         dm(3,l1) = rpole(4,i)
         qm(1,l1) = rpole(5,i)
         qm(2,l1) = rpole(6,i)
         qm(3,l1) = rpole(7,i)
         qm(4,l1) = rpole(8,i)
         qm(5,l1) = rpole(9,i)
         qm(6,l1) = rpole(10,i)
         qm(7,l1) = rpole(11,i)
         qm(8,l1) = rpole(12,i)
         qm(9,l1) = rpole(13,i)
         um(1,l1) = uind(1,l1)
         um(2,l1) = uind(2,l1)
         um(3,l1) = uind(3,l1)
      end do

      do l1 = 1, npole3b
         i = pnum(l1)
c         trqi(1,i) = 0.0d0
c         trqi(2,i) = 0.0d0
c         trqi(3,i) = 0.0d0
         trqi(1,l1) = 0.0d0
         trqi(2,l1) = 0.0d0
         trqi(3,l1) = 0.0d0
      end do

      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(l1,0) = 1.0d0
         ejs(l1,0) = 0.0d0
         ekc(l1,0) = 1.0d0
         eks(l1,0) = 0.0d0
         elc(l1,0) = 1.0d0
         els(l1,0) = 0.0d0
         ejc(l1,1) = cos(xfr)
         ejs(l1,1) = sin(xfr)
         ekc(l1,1) = cos(yfr)
         eks(l1,1) = sin(yfr)
         elc(l1,1) = cos(zfr)
         els(l1,1) = sin(zfr)
         ekc(l1,-1) = ekc(l1,1)
         eks(l1,-1) = -eks(l1,1)
         elc(l1,-1) = elc(l1,1)
         els(l1,-1) = -els(l1,1)
         do j = 2, jmax3b
            ejc(l1,j) = ejc(l1,j-1)*ejc(l1,1) - ejs(l1,j-1)*ejs(l1,1)
            ejs(l1,j) = ejs(l1,j-1)*ejc(l1,1) + ejc(l1,j-1)*ejs(l1,1)
         end do
         do j = 2, kmax3b
            ekc(l1,j) = ekc(l1,j-1)*ekc(l1,1) - eks(l1,j-1)*eks(l1,1)
            eks(l1,j) = eks(l1,j-1)*ekc(l1,1) + ekc(l1,j-1)*eks(l1,1)
            ekc(l1,-j) = ekc(l1,j)
            eks(l1,-j) = -eks(l1,j)
         end do
         do j = 2, lmax3b
            elc(l1,j) = elc(l1,j-1)*elc(l1,1) - els(l1,j-1)*els(l1,1)
            els(l1,j) = els(l1,j-1)*elc(l1,1) + elc(l1,j-1)*els(l1,1)
            elc(l1,-j) = elc(l1,j)
            els(l1,-j) = -els(l1,j)
         end do
      end do

      do j = jmin, jmax3b
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax3b
            rk = 2.0d0 * pi * dble(k)
            do l1 = 1, npole3b
               i = pnum(l1)
c               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
c               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
               cjk(l1) = ejc(l1,j)*ekc(l1,k) - ejs(l1,j)*eks(l1,k)
               sjk(l1) = ejs(l1,j)*ekc(l1,k) + ejc(l1,j)*eks(l1,k)
            end do
            do l = lmin, lmax3b
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
c               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  t3p = 0.0d0
                  t4p = 0.0d0
                  do m2 = 1, 3
                     do m1 = 1, 3
                        t5(m1,m2) = 0.0d0
                        t5u(m1,m2) = 0.0d0
                        t5p(m1,m2) = 0.0d0
                        t6(m1,m2) = 0.0d0
                        t6u(m1,m2) = 0.0d0
                        t6p(m1,m2) = 0.0d0
                     end do
                  end do
                  do l1 = 1, npole3b
                     i = pnum(l1)
                     ckr(l1) = cjk(l1)*elc(l1,l) - sjk(l1)*els(l1,l)
                     skr(l1) = sjk(l1)*elc(l1,l) + cjk(l1)*els(l1,l)
                     ck = cm(l1)
                     dk = h1*dm(1,l1) + h2*dm(2,l1) + h3*dm(3,l1)
                     dkx(l1) = h3*dm(2,l1) - h2*dm(3,l1)
                     dky(l1) = h1*dm(3,l1) - h3*dm(1,l1)
                     dkz(l1) = h2*dm(1,l1) - h1*dm(2,l1)
                     q1 = h1*qm(1,l1) + h2*qm(4,l1) + h3*qm(7,l1)
                     q2 = h1*qm(2,l1) + h2*qm(5,l1) + h3*qm(8,l1)
                     q3 = h1*qm(3,l1) + h2*qm(6,l1) + h3*qm(9,l1)
                     qk = h1*q1 + h2*q2 + h3*q3
                     qkx(l1) = h3*q2 - h2*q3
                     qky(l1) = h1*q3 - h3*q1
                     qkz(l1) = h2*q1 - h1*q2
                     uk = h1*uind(1,l1) + h2*uind(2,l1) + h3*uind(3,l1)
                     ukp = h1*uinp(1,l1) + h2*uinp(2,l1) + h3*uinp(3,l1)
                     s1(l1) = (ck-qk)*skr(l1) + dk*ckr(l1)
                     s2(l1) = (ck-qk)*ckr(l1) - dk*skr(l1)
                     s3(l1) = uk * ckr(l1)
                     s4(l1) = -uk * skr(l1)
                     s3p(l1) = ukp * ckr(l1)
                     s4p(l1) = -ukp * skr(l1)
                     t1 = t1 + s1(l1)
                     t2 = t2 + s2(l1)
                     t3 = t3 + s3(l1)
                     t4 = t4 + s4(l1)
                     t3p = t3p + s3p(l1)
                     t4p = t4p + s4p(l1)
c
c     terms needed for subsequent virial tensor calculation
c
                  qt(1,1)=h1*(h1*qm(1,l1)+h2*qm(4,l1)+h3*qm(7,l1))
                  qt(2,1)=h1*(h1*qm(2,l1)+h2*qm(5,l1)+h3*qm(8,l1))
                  qt(3,1) = h1*(h1*qm(3,l1)+h2*qm(6,l1)+h3*qm(9,l1))
                  qt(1,2) = h2*(h1*qm(1,l1)+h2*qm(4,l1)+h3*qm(7,l1))
                  qt(2,2) = h2*(h1*qm(2,l1)+h2*qm(5,l1)+h3*qm(8,l1))
                  qt(3,2) = h2*(h1*qm(3,l1)+h2*qm(6,l1)+h3*qm(9,l1))
                  qt(1,3) = h3*(h1*qm(1,l1)+h2*qm(4,l1)+h3*qm(7,l1))
                  qt(2,3) = h3*(h1*qm(2,l1)+h2*qm(5,l1)+h3*qm(8,l1))
                  qt(3,3) = h3*(h1*qm(3,l1)+h2*qm(6,l1)+h3*qm(9,l1))
                     dt(1,1) = h1 * dm(1,l1)
                     dt(2,1) = h1 * dm(2,l1)
                     dt(3,1) = h1 * dm(3,l1)
                     dt(1,2) = h2 * dm(1,l1)
                     dt(2,2) = h2 * dm(2,l1)
                     dt(3,2) = h2 * dm(3,l1)
                     dt(1,3) = h3 * dm(1,l1)
                     dt(2,3) = h3 * dm(2,l1)
                     dt(3,3) = h3 * dm(3,l1)
                     dtu(1,1) = h1 * uind(1,l1)
                     dtu(2,1) = h1 * uind(2,l1)
                     dtu(3,1) = h1 * uind(3,l1)
                     dtu(1,2) = h2 * uind(1,l1)
                     dtu(2,2) = h2 * uind(2,l1)
                     dtu(3,2) = h2 * uind(3,l1)
                     dtu(1,3) = h3 * uind(1,l1)
                     dtu(2,3) = h3 * uind(2,l1)
                     dtu(3,3) = h3 * uind(3,l1)
                     dtp(1,1) = h1 * uinp(1,l1)
                     dtp(2,1) = h1 * uinp(2,l1)
                     dtp(3,1) = h1 * uinp(3,l1)
                     dtp(1,2) = h2 * uinp(1,l1)
                     dtp(2,2) = h2 * uinp(2,l1)
                     dtp(3,2) = h2 * uinp(3,l1)
                     dtp(1,3) = h3 * uinp(1,l1)
                     dtp(2,3) = h3 * uinp(2,l1)
                     dtp(3,3) = h3 * uinp(3,l1)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           t5(m1,m2) = t5(m1,m2) - dt(m1,m2)*ckr(l1)
     &                                    + 2.0d0*qt(m1,m2)*skr(l1)
                           t5u(m1,m2) = t5u(m1,m2) - dtu(m1,m2)*ckr(l1)
                           t5p(m1,m2) = t5p(m1,m2) - dtp(m1,m2)*ckr(l1)
                           t6(m1,m2) = t6(m1,m2) + dt(m1,m2)*skr(l1)
     &                                    + 2.0d0*qt(m1,m2)*ckr(l1)
                           t6u(m1,m2) = t6u(m1,m2) + dtu(m1,m2)*skr(l1)
                           t6p(m1,m2) = t6p(m1,m2) + dtp(m1,m2)*skr(l1)
                        end do
                     end do
                  end do
c
c     get the energy contributions for current reciprocal vector
c
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  ei = expterm * (t1*t3+t2*t4)
                  etot = e + ei
                  eptemp = eptemp + ei

                  uterm = expterm * (t1*(t3+t3p) + t3*t3p
     &                                 + t2*(t4+t4p) + t4*t4p)
                  do m2 = 1, 3
                     do m1 = 1, 3
c                        wterm(m1,m2) = 2.0d0 * expterm
c     &                     * (t1*t5(m1,m2) + t2*t6(m1,m2)
c     &                        + 0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
c     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
c     &                        + (t3+t3p)*t5(m1,m2)
c     &                        + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
c     &                        + (t4+t4p)*t6(m1,m2)
c     &                        + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))

c     NEW wterm w/ just the permanent elec part removed (already done)
                        wterm(m1,m2) = 2.0d0 * expterm
     &                     * (
     &                        0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
     &                        + (t3+t3p)*t5(m1,m2)
     &                        + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
     &                        + (t4+t4p)*t6(m1,m2)
     &                        + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))

                     end do
                  end do

                  wterm(2,1) = 0.5d0 * (wterm(2,1)+wterm(1,2))
                  wterm(3,1) = 0.5d0 * (wterm(3,1)+wterm(1,3))
                  wterm(3,2) = 0.5d0 * (wterm(3,2)+wterm(2,3))
                  wterm(1,2) = wterm(2,1)
                  wterm(1,3) = wterm(3,1)
                  wterm(2,3) = wterm(3,2)
                  vterm = 2.0d0 * uterm * (1.0d0-term*hsq) / hsq

                  virtemp(1,1) = virtemp(1,1) + h1*h1*vterm + wterm(1,1)
     &                            - uterm
                  virtemp(2,1) = virtemp(2,1) + h2*h1*vterm + wterm(2,1)
                  virtemp(3,1) = virtemp(3,1) + h3*h1*vterm + wterm(3,1)
                  virtemp(1,2) = virtemp(1,2) + h1*h2*vterm + wterm(1,2)
                  virtemp(2,2) = virtemp(2,2) + h2*h2*vterm + wterm(2,2)
     &                            - uterm
                  virtemp(3,2) = virtemp(3,2) + h3*h2*vterm + wterm(3,2)
                  virtemp(1,3) = virtemp(1,3) + h1*h3*vterm + wterm(1,3)
                  virtemp(2,3) = virtemp(2,3) + h2*h3*vterm + wterm(2,3)
                  virtemp(3,3) = virtemp(3,3) + h3*h3*vterm + wterm(3,3)
     &                             - uterm

c
c     get the force contributions for current reciprocal vector
c
                  expterm = 2.0d0 * expterm
                  do l1 = 1, npole3b
                     i = pnum(l1)
                     ii = ipole(i)
                     dei = 0.5d0 * expterm * ((s4(l1)+s4p(l1))*t1
     &                                       -(s3(l1)+s3p(l1))*t2
     &                                +s2(l1)*(t3+t3p)-s1(l1)*(t4+t4p))
c                     if (poltyp .eq. 'MUTUAL') then
                         dei = dei + 0.5d0 * expterm
     &                            * (s4p(l1)*t3+s4(l1)*t3p
     &                              -s3p(l1)*t4-s3(l1)*t4p)
c                     end if
c                     det1 = expterm * (skr(l1)*t2-ckr(l1)*t1)
c                     det2 = 2.0d0 * expterm * (ckr(l1)*t2+skr(l1)*t1)
                     det1i = 0.5d0 * expterm * (skr(l1)*(t4+t4p)
     &                                         -ckr(l1)*(t3+t3p))
                     det2i =expterm*(ckr(l1)*(t4+t4p)+skr(l1)*(t3+t3p))
c                     deptemp(1,i) = deptemp(1,i) + h1*dei
c                     deptemp(2,i) = deptemp(2,i) + h2*dei
c                     deptemp(3,i) = deptemp(3,i) + h3*dei

                     deptemp(1,l1) = deptemp(1,l1) + h1*dei
                     deptemp(2,l1) = deptemp(2,l1) + h2*dei
                     deptemp(3,l1) = deptemp(3,l1) + h3*dei

c                     trq(1,l1) = trq(1,l1) + dkx(l1)*det1 + qkx(l1)*det2
c                     trq(2,l1) = trq(2,l1) + dky(l1)*det1 + qky(l1)*det2
c                     trq(3,l1) = trq(3,l1) + dkz(l1)*det1 + qkz(l1)*det2

c                     trqi(1,i) = trqi(1,i) + dkx(l1)*det1i
c     &                               + qkx(l1)*det2i
c                     trqi(2,i) = trqi(2,i) + dky(l1)*det1i
c     &                               + qky(l1)*det2i
c                     trqi(3,i) = trqi(3,i) + dkz(l1)*det1i
c     &                               + qkz(l1)*det2i

                     trqi(1,l1) = trqi(1,l1) + dkx(l1)*det1i
     &                               + qkx(l1)*det2i
                     trqi(2,l1) = trqi(2,l1) + dky(l1)*det1i
     &                               + qky(l1)*det2i
                     trqi(3,l1) = trqi(3,l1) + dkz(l1)*det1i
     &                               + qkz(l1)*det2i

                  end do
c               end if
            end do
            lmin = -lmax3b
         end do
         kmin = -kmax3b
      end do

c      call torque_ewreg3b (trqi,deptemp,pnum,npole3b)
      call torque2_3b (npole3b,pnum,trqi,deptemp) 
      return
      end

c      subroutine ereal1c_3b(aewald3b3b,ewaldcut3b,npole3b,pnum,uind,uinp,
c     &   eptemp,deptemp,virtemp)
      subroutine ereal1c_3b(npole3b,pnum,uind,uinp,
     &   eptemp,deptemp,virtemp)
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use limits
      use ewald
      use math
      use mpole
      use polgrp
      use polpot
      use polar, only: polarity, thole, pdamp      
      use shunt
      implicit none
      integer moli1,moli2,moli3
      integer i,j,k,l1,l3,l2,l,ll
      integer ii,kk,jcell,jj
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 uind(3,*),eptemp,deptemp(3,*)
      real*8 virtemp(3,3),uinp(3,*)
      real*8 ewaldcut3b_2
      integer npole3b,pnum(*)
      character*6 mode
      logical flag
      external erfc

      eintra = 0.0d0
c      if (npole3b .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c

c      allocate (mscale(n))
c      allocate (pscale(n))
c      allocate (dscale(n))
c      allocate (uscale(n))

      allocate (pscale(npole3b))
      allocate (dscale(npole3b))
      allocate (uscale(npole3b))

c
c     set arrays needed to scale connected atom interactions
c
c      do i = 1, n
c         mscale(i) = 1.0d0
      do i = 1, npole3b
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c      ewaldcut3b_2=ewaldcut3b*ewaldcut3b
      ewaldcut3b_2= off2     
      aewald3b=aewald
      do l1 = 1, npole3b-1
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i12(j,ii)) then
                 pscale(kk) = p2scale
                 goto 31
               end if
            end do
   31   continue 
c            pscale(i12(j,ii)) = p2scale
         end do

         do j = 1, n13(ii)
c            pscale(i13(j,ii)) = p3scale
            do kk=1,npole3b
               if(pnum(kk).eq.i13(j,ii)) then
                 pscale(kk) = p3scale
                 goto 32
               end if
            end do
   32   continue
         end do

         do j = 1, n14(ii)
c            pscale(i14(j,ii)) = p4scale
            do kk=1,npole3b
               if(pnum(kk).eq.i14(j,ii)) then
                 pscale(kk) = p4scale
                 goto 33
               end if
            end do
   33   continue
            do k = 1, np11(ii)
               do kk=1,npole3b 
                  if ((i14(j,ii) .eq. ip11(k,ii)).and.
     &               (pnum(kk).eq.ip11(k,ii)) ) then
c     &            pscale(i14(j,ii)) = p4scale * p41scale
                   pscale(kk) = p4scale * p41scale
                   goto 34
                  end if 
               end do
   34   continue
            end do
         end do
         do j = 1, n15(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i15(j,ii)) then
                 pscale(kk) = p5scale
                 goto 35
               end if
            end do      
   35   continue      
c            pscale(i15(j,ii)) = p5scale
         end do

         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = d1scale
c            uscale(ip11(j,ii)) = u1scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip11(j,ii)) then
                 dscale(kk) = d1scale
                 uscale(kk) = u1scale                          
                 goto 36
               end if
            end do
   36   continue
         end do

         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = d2scale
c            uscale(ip12(j,ii)) = u2scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip12(j,ii)) then
                 dscale(kk) = d2scale
                 uscale(kk) = u2scale
                 goto 37
               end if
            end do
   37   continue
         end do
         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = d3scale
c            uscale(ip13(j,ii)) = u3scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip13(j,ii)) then
                 dscale(kk) = d3scale
                 uscale(kk) = u3scale
                 goto 38
               end if
            end do
   38   continue
         end do
         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = d4scale
c            uscale(ip14(j,ii)) = u4scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip14(j,ii)) then
                 dscale(kk) = d4scale
                 uscale(kk) = u4scale
                 goto 39
               end if
            end do
   39   continue             
         end do
         do l2 = l1+1, npole3b
            k=pnum(l2)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
c            if (r2 .le. off2) then
            if (r2 .le. ewaldcut3b_2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald3b * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald3b**2
               alsq2n = 0.0d0
               if (aewald3b .gt. 0.0d0) 
     &              alsq2n = 1.0d0 / (sqrtpi*aewald3b)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
c               dsc3 = 1.0d0 - scale3*dscale(kk)
c               dsc5 = 1.0d0 - scale5*dscale(kk)
c               dsc7 = 1.0d0 - scale7*dscale(kk)
c               psc3 = 1.0d0 - scale3*pscale(kk)
c               psc5 = 1.0d0 - scale5*pscale(kk)
c               psc7 = 1.0d0 - scale7*pscale(kk)
c               usc3 = 1.0d0 - scale3*uscale(kk)
c               usc5 = 1.0d0 - scale5*uscale(kk)

               dsc3 = 1.0d0 - scale3*dscale(l2)
               dsc5 = 1.0d0 - scale5*dscale(l2)
               dsc7 = 1.0d0 - scale7*dscale(l2)
               psc3 = 1.0d0 - scale3*pscale(l2)
               psc5 = 1.0d0 - scale5*pscale(l2)
               psc7 = 1.0d0 - scale7*pscale(l2)
               usc3 = 1.0d0 - scale3*uscale(l2)
               usc5 = 1.0d0 - scale5*uscale(l2)

c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
     &                      + qi(7)*uind(3,l2)
               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
     &                      + qi(8)*uind(3,l2)
               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
     &                      + qi(9)*uind(3,l2)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
     &                       + qi(7)*uinp(3,l2)
               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
     &                       + qi(8)*uinp(3,l2)
               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
     &                       + qi(9)*uinp(3,l2)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &                     + qir(3)*uind(3,l2)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for permanent components
c

c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
c               gl(4) = sc(5)*sc(6)
c               gl(5) = -4.0d0 * sc(9)
c               gl(6) = sc(2)
c               gl(7) = 2.0d0 * (sc(7)-sc(8))
c               gl(8) = 2.0d0 * sc(10)

c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = ei - erli
               ei = f * ei
               eptemp = eptemp + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     intermediate variables for permanent force terms
c

c               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
c     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
c     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
c               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
c               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
c               gf(4) = 2.0d0 * bn(2)
c               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
c               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
c               gf(7) = 4.0d0 * bn(3)
c               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
c     &                     + rr7*(gl(2)+gl(7)+gl(8))
c     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
c               gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
c               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
c               gfr(4) = 2.0d0 * rr5
c               gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
c               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
c               gfr(7) = 4.0d0 * rr7

c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))

               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l2)
     &                                  +(glip(1)+glip(6))*dscale(l2))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l2)
     &                                  +(glip(2)+glip(7))*dscale(l2))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l2)
     &                                  +glip(3)*dscale(l2))

               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
c               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(kk)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               temp3 = 0.5d0 * rr3 * uscale(l2) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(l2)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c
c     increment gradient due to force and torque on first site
c

c               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
c               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
c               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
               deptemp(1,l1) =deptemp(1,l1) + ftm2i(1)
               deptemp(2,l1) =deptemp(2,l1) + ftm2i(2)
               deptemp(3,l1) =deptemp(3,l1) + ftm2i(3)

c               call torque_3b (deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               call torque_3b_new(npole3b,pnum,deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

c               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
c               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
c               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
               deptemp(1,l2) = deptemp(1,l2) - ftm2i(1)
               deptemp(2,l2) = deptemp(2,l2) - ftm2i(2)
               deptemp(3,l2) = deptemp(3,l2) - ftm2i(3)

c               call torque_3b (deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
               call torque_3b_new(npole3b,pnum,deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)

               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

            end if
c          end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c

c         do j = 1, n12(ii)
c            mscale(i12(j,ii)) = 1.0d0
c            pscale(i12(j,ii)) = 1.0d0
c         end do
c         do j = 1, n13(ii)
c            mscale(i13(j,ii)) = 1.0d0
c            pscale(i13(j,ii)) = 1.0d0
c         end do
c         do j = 1, n14(ii)
c            mscale(i14(j,ii)) = 1.0d0
c            pscale(i14(j,ii)) = 1.0d0
c         end do
c         do j = 1, n15(ii)
c            mscale(i15(j,ii)) = 1.0d0
c            pscale(i15(j,ii)) = 1.0d0
c         end do
c         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = 1.0d0
c            uscale(ip11(j,ii)) = 1.0d0
c         end do
c         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = 1.0d0
c            uscale(ip12(j,ii)) = 1.0d0
c         end do
c         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = 1.0d0
c            uscale(ip13(j,ii)) = 1.0d0
c         end do
c         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = 1.0d0
c            uscale(ip14(j,ii)) = 1.0d0
c         end do
         do j=1,npole3b
            pscale(j)=1.0d0
            dscale(j)=1.0d0
            uscale(j)=1.0d0
         end do
      end do

      if (use_replica) then
c
c     calculate interactions with other unit cells
c
c      print*,'After use_replica!',moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After use_replica l1 Pnum(l1)",l1,pnum(l1)
c      end do 
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i12(j,ii)) then
                 pscale(kk) = p2scale
                 goto 41
               end if
            end do
   41    continue
c            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i13(j,ii)) then
                 pscale(kk) = p3scale
                 goto 42
               end if
            end do
   42    continue
c            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i14(j,ii)) then
                 pscale(kk) = p4scale
                 goto 43
               end if
            end do
   43    continue
c            pscale(i14(j,ii)) = p4scale
         end do
         do j = 1, n15(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i15(j,ii)) then
                 pscale(kk) = p5scale
                 goto 44
               end if
            end do
   44    continue
c            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip11(j,ii)) then
                 dscale(kk) = d1scale
                 uscale(kk) = u1scale
                 goto 45
               end if
            end do
   45    continue
c            dscale(ip11(j,ii)) = d1scale
c            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip12(j,ii)) then
                 dscale(kk) = d2scale
                 uscale(kk) = u2scale
                 goto 46
               end if
            end do
   46    continue
c            dscale(ip12(j,ii)) = d2scale
c            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip13(j,ii)) then
                 dscale(kk) = d3scale
                 uscale(kk) = u3scale
                 goto 47
               end if
            end do
   47    continue
c            dscale(ip13(j,ii)) = d3scale
c            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.ip14(j,ii)) then
                 dscale(kk) = d4scale
                 uscale(kk) = u4scale
                 goto 48
               end if
            end do
   48    continue
c            dscale(ip14(j,ii)) = d4scale
c            uscale(ip14(j,ii)) = u4scale
         end do
         do l2 = l1, npole3b
            k=pnum(l2)
            kk = ipole(k)
            do jcell=1,ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not.(use_polymer .and. r2.le.polycut2)) then
c              mscale(kk)=1.0d0
c              pscale(kk)=1.0d0
c              dscale(kk)=1.0d0
c              uscale(kk)=1.0d0
              pscale(l2)=1.0d0
              dscale(l2)=1.0d0
              uscale(l2)=1.0d0
            end if
c            if (r2 .le. off2) then
            if (r2 .le. ewaldcut3b_2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald3b * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald3b**2
               alsq2n = 0.0d0
               if (aewald3b .gt. 0.0d0)  
     &             alsq2n = 1.0d0 / (sqrtpi*aewald3b)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if

c               dsc3 = 1.0d0 - scale3*dscale(kk)
c               dsc5 = 1.0d0 - scale5*dscale(kk)
c               dsc7 = 1.0d0 - scale7*dscale(kk)
c               psc3 = 1.0d0 - scale3*pscale(kk)
c               psc5 = 1.0d0 - scale5*pscale(kk)
c               psc7 = 1.0d0 - scale7*pscale(kk)
c               usc3 = 1.0d0 - scale3*uscale(kk)
c               usc5 = 1.0d0 - scale5*uscale(kk)

               dsc3 = 1.0d0 - scale3*dscale(l2)
               dsc5 = 1.0d0 - scale5*dscale(l2)
               dsc7 = 1.0d0 - scale7*dscale(l2)
               psc3 = 1.0d0 - scale3*pscale(l2)
               psc5 = 1.0d0 - scale5*pscale(l2)
               psc7 = 1.0d0 - scale7*pscale(l2)
               usc3 = 1.0d0 - scale3*uscale(l2)
               usc5 = 1.0d0 - scale5*uscale(l2)

c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
     &                      + qi(7)*uind(3,l2)
               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
     &                      + qi(8)*uind(3,l2)
               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
     &                      + qi(9)*uind(3,l2)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
     &                       + qi(7)*uinp(3,l2)
               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
     &                       + qi(8)*uinp(3,l2)
               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
     &                       + qi(9)*uinp(3,l2)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &                     + qir(3)*uind(3,l2)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = ei - erli
               ei = f * ei
               if (ii .eq. kk) then
                 ei = 0.5d0 * ei
               end if
               eptemp = eptemp + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))

               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l2)
     &                                  +(glip(1)+glip(6))*dscale(l2))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l2)
     &                                  +(glip(2)+glip(7))*dscale(l2))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l2)
     &                                  +glip(3)*dscale(l2))

               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(l2) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(l2)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
              if (use_polymer .and. r2.le.polycut2) then
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              else
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              end if
               if (ii .eq. kk) then
                  do j = 1, 3
                     ftm2i(j) = 0.5d0 * ftm2i(j)
                     ttm2i(j) = 0.5d0 * ttm2i(j)
                     ttm3i(j) = 0.5d0 * ttm3i(j)
                  end do
               end if

c
c     increment gradient due to force and torque on first site
c

c               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
c               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
c               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
c               call torque_3b (deptemp,i,
c     &           ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

c               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
c               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
c               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
c               call torque_3b (deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               deptemp(1,l1) =deptemp(1,l1) + ftm2i(1)
               deptemp(2,l1) =deptemp(2,l1) + ftm2i(2)
               deptemp(3,l1) =deptemp(3,l1) + ftm2i(3)

c               call torque_3b (deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               call torque_3b_new(npole3b,pnum,deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

c               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
c               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
c               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
               deptemp(1,l2) = deptemp(1,l2) - ftm2i(1)
               deptemp(2,l2) = deptemp(2,l2) - ftm2i(2)
               deptemp(3,l2) = deptemp(3,l2) - ftm2i(3)

c               call torque_3b (deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
               call torque_3b_new(npole3b,pnum,deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)


               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c               vxx = -xr*ftm2i(1)+xix*frcxi(1)+xiy*frcyi(1)+xiz*frczi(1)
c               vyx = -yr*ftm2i(1)+yix*frcxi(1)+yiy*frcyi(1)+yiz*frczi(1)
c               vzx = -zr*ftm2i(1)+zix*frcxi(1)+ziy*frcyi(1)+ziz*frczi(1)
c               vyy = -yr*ftm2i(2)+yix*frcxi(2)+yiy*frcyi(2)+yiz*frczi(2)
c               vzy = -zr*ftm2i(2)+zix*frcxi(2)+ziy*frcyi(2)+ziz*frczi(2)
c               vzz = -zr*ftm2i(3)+zix*frcxi(3)+ziy*frcyi(3)+ziz*frczi(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

            end if
         end do
c         end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
c         do j = 1, n12(ii)
c            mscale(i12(j,ii)) = 1.0d0
c            pscale(i12(j,ii)) = 1.0d0
c         end do
c         do j = 1, n13(ii)
c            mscale(i13(j,ii)) = 1.0d0
c            pscale(i13(j,ii)) = 1.0d0
c         end do
c         do j = 1, n14(ii)
c            mscale(i14(j,ii)) = 1.0d0
c            pscale(i14(j,ii)) = 1.0d0
c         end do
c         do j = 1, n15(ii)
c            mscale(i15(j,ii)) = 1.0d0
c            pscale(i15(j,ii)) = 1.0d0
c         end do
c         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = 1.0d0
c            uscale(ip11(j,ii)) = 1.0d0
c         end do
c         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = 1.0d0
c            uscale(ip12(j,ii)) = 1.0d0
c         end do
c         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = 1.0d0
c            uscale(ip13(j,ii)) = 1.0d0
c         end do
c         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = 1.0d0
c            uscale(ip14(j,ii)) = 1.0d0
c         end do
         do j=1,npole3b
            dscale(j)=1.0d0
            uscale(j)=1.0d0
            pscale(j)=1.0d0
         end do
      end do
      end if
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c
c     perform deallocation of some local arrays
c
c      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
c      print*,'ep after ereal1c_3b=',eptemp,moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After ereal1c_3b l1 Pnum(l1)",l1,pnum(l1)
c      end do

      return
      end

      subroutine torque_ewreg3b (trq,deptemp,pnum,npole3b)
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'mpole.i'
      use sizes
      use atoms
      use mpole
      implicit none
      integer i,j,ia,ib,ic,l1,l3
      real*8 usiz,vsiz,wsiz
      real*8 upsiz,vpsiz
      real*8 dotdu,dotdv
      real*8 dphidu,dphidv,dphidw
      real*8 c,s,uvdis,vudis,du,dv
      real*8 u(3),v(3),w(3)
      real*8 up(3),vp(3),diff(3)
      real*8 trq(3,*)
      real*8 deptemp(3,*)
c      integer pnum(*),npole3b,pnum2(*)
      integer pnum(*),npole3b
      integer ia_l1,ib_l1,ic_l1,l2
c
c     coordinate frame motion described by rotation about u, v and w
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ia = zaxis(i)
         ib = ipole(i)
         ic = xaxis(i)

         do l2=1,npole3b
            if(pnum(l2).eq.zaxis(i)) then
              ia_l1=l2
              goto 31
            end if
         end do
   31    continue

         do l2=1,npole3b
            if(pnum(l2).eq.ipole(i)) then
              ib_l1=l2
              goto 32
            end if
         end do
   32    continue

         do l2=1,npole3b
            if(pnum(l2).eq.xaxis(i)) then
              ic_l1=l2
              goto 33
            end if
         end do
   33    continue
            
         u(1) = x(ia) - x(ib)
         u(2) = y(ia) - y(ib)
         u(3) = z(ia) - z(ib)
         usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
         wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
         dotdu = 0.0d0
         dotdv = 0.0d0
         do j = 1, 3
            u(j) = u(j) / usiz
            v(j) = v(j) / vsiz
            w(j) = w(j) / wsiz
            diff(j) = v(j) - u(j)
            dotdu = dotdu + u(j)*diff(j)
            dotdv = dotdv + v(j)*diff(j)
         end do
c
c     get perpendiculars to u,v to get direction of motion of
c     u or v due to rotation about the cross product vector w
c
         upsiz = 0.0d0
         vpsiz = 0.0d0
         do j = 1, 3
            up(j) = diff(j) - dotdu*u(j)
            vp(j) = diff(j) - dotdv*v(j)
            upsiz = upsiz + up(j)*up(j)
            vpsiz = vpsiz + vp(j)*vp(j)
         end do
         upsiz = sqrt(upsiz)
         vpsiz = sqrt(vpsiz)
         do j = 1, 3
            up(j) = up(j) / upsiz
            vp(j) = vp(j) / vpsiz
         end do

c
c     negative of dot product of torque with unit vectors along u, v
c     and w give result of infinitesmal rotation along these vectors
c     i.e. dphi/dtheta = dot product, where dphi is work, dtheta is
c     angle; then dphi/dtheta is torque and the dot product is torque
c     component along unit vector
c
c         dphidu = -trq(1,ib)*u(1) - trq(2,ib)*u(2) - trq(3,ib)*u(3)
c         dphidv = -trq(1,ib)*v(1) - trq(2,ib)*v(2) - trq(3,ib)*v(3)
c         dphidw = -trq(1,ib)*w(1) - trq(2,ib)*w(2) - trq(3,ib)*w(3)

         dphidu = -trq(1,ib_l1)*u(1) - trq(2,ib_l1)*u(2)
     &            - trq(3,ib_l1)*u(3)
         dphidv = -trq(1,ib_l1)*v(1) - trq(2,ib_l1)*v(2)
     &            - trq(3,ib_l1)*v(3)
         dphidw = -trq(1,ib_l1)*w(1) - trq(2,ib_l1)*w(2)
     &            - trq(3,ib_l1)*w(3)

c
c     get the projected distances between the vectors
c
         c = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
         s = sqrt(1.0d0 - c*c)
         uvdis = usiz * s
         vudis = vsiz * s
c
c     force distribution for the bisector local coordinate method
c
         if (polaxe(i) .eq. 'Bisector') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/(2.0d0*usiz)
               dv = w(j)*dphidu/vudis + vp(j)*dphidw/(2.0d0*vsiz)
c               deptemp(j,ia) = deptemp(j,ia) + du
c               deptemp(j,ic) = deptemp(j,ic) + dv
c               deptemp(j,ib) = deptemp(j,ib) - dv - du

               deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
               deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
               deptemp(j,ib_l1) = deptemp(j,ib_l1) - dv - du

            end do
c
c     force distribution for the Z-then-X local coordinate method
c
         else if (polaxe(i) .eq. 'Z-then-X') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/usiz
               dv = w(j)*dphidu/vudis
c               deptemp(j,ia) = deptemp(j,ia) + du
c               deptemp(j,ic) = deptemp(j,ic) + dv
c               deptemp(j,ib) = deptemp(j,ib) - dv - du

               deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
               deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
               deptemp(j,ib_l1) = deptemp(j,ib_l1) - dv - du

            end do
         end if
      end do
      return
      end

