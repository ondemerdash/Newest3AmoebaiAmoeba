c
c
      subroutine empole1a_3b_Polar_rtrnfield(npole3b,pnum,eptemp,
     &   deptemp,virtemp,field,fieldp)
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
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3
      real*8 field(3,npole3b),fieldp(3,npole3b)
      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      logical proceed,usei,usek
      character*6 mode

      eptemp = 0.0d0
      eptemp_userep =0.0d0
c
c   Zero out temporary gradient of polarization energy
c
c      mode = 'MPOLE'
c      call switch (mode)
c      print*,"mpole off2",off2
      off2=1.0d12
      do l1 = 1, npole3b
        do j = 1, 3
          deptemp(j,l1) = 0.0d0
        end do
      end do
c
c   Zero out temporary virial
c
      do i=1,3
         do j=1,3
           virtemp(i,j)=0.0d0
         end do
      end do

      call induce0a_3b_PolelecOnly_rtrnfield(npole3b,pnum,uind,uinp,
     &     field,fieldp)

c      do l1=1,npole3b
c         i=pnum(l1)
c         print*,"npole3b i uindx",npole3b,i,uind(1,l1)
c         print*,"npole3b i uindy",npole3b,i,uind(2,l1)
c         print*,"npole3b i uindz",npole3b,i,uind(3,l1)
c      end do

c      call induce0a_3b_PolelecOnly_subtrdir(npole3b,pnum,uind,uinp)

c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(npole3b))
      allocate (dscale(npole3b))
      allocate (uscale(npole3b))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, npole3b
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c
c     set scale factors for permanent multipole and induced terms
c

      do l1 = 1, npole3b-1
         i=pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyx = rpole(8,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizx = rpole(11,i)
         qizy = rpole(12,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i12(j,ii)) then
                  pscale(kk)=p2scale
                  goto 31
               end if
            end do 
   31             continue

         end do
         do j = 1, n13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i13(j,ii)) then
                  pscale(kk)=p3scale
                  goto 32
               end if
            end do
   32             continue            
         end do
         do j = 1, n14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i14(j,ii)) then
                  pscale(kk)=p4scale
                  goto 33
               end if
            end do
   33             continue
            do k = 1, np11(ii)
               do kk=1,npole3b             
                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
     &              (pnum(kk).eq.ip11(k,ii)) ) then
                   pscale(kk) = p4scale * p41scale
                   goto 34
                 end if
               end do  
   34             continue
            end do
         end do
         do j = 1, n15(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i15(j,ii)) then
                  pscale(kk)=p5scale
                  goto 35
               end if
            end do
   35             continue
         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
                  dscale(kk)=d1scale
                  uscale(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale(kk)=d2scale
                  uscale(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale(kk)=d3scale
                  uscale(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale(kk)=d4scale
                  uscale(kk)=u4scale
                  goto 39
               end if
            end do
   39             continue
         end do
         do l3 = l1+1, npole3b
            k=pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr

            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyx = rpole(8,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzx = rpole(11,k)
               qkzy = rpole(12,k)
               qkzz = rpole(13,k)
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
               scale3i = scale3 * uscale(l3)
               scale5i = scale5 * uscale(l3)
               scale7i = scale7 * uscale(l3)
               dsc3 = scale3 * dscale(l3)
               dsc5 = scale5 * dscale(l3)
               dsc7 = scale7 * dscale(l3)
               psc3 = scale3 * pscale(l3)
               psc5 = scale5 * pscale(l3)
               psc7 = scale7 * pscale(l3)
c
c     construct necessary auxiliary vectors
c
               dixuk(1) = diy*uind(3,l3) - diz*uind(2,l3)
               dixuk(2) = diz*uind(1,l3) - dix*uind(3,l3)
               dixuk(3) = dix*uind(2,l3) - diy*uind(1,l3)
               dkxui(1) = dky*uind(3,l1) - dkz*uind(2,l1)
               dkxui(2) = dkz*uind(1,l1) - dkx*uind(3,l1)
               dkxui(3) = dkx*uind(2,l1) - dky*uind(1,l1)
               dixukp(1) = diy*uinp(3,l3) - diz*uinp(2,l3)
               dixukp(2) = diz*uinp(1,l3) - dix*uinp(3,l3)
               dixukp(3) = dix*uinp(2,l3) - diy*uinp(1,l3)
               dkxuip(1) = dky*uinp(3,l1) - dkz*uinp(2,l1)
               dkxuip(2) = dkz*uinp(1,l1) - dkx*uinp(3,l1)
               dkxuip(3) = dkx*uinp(2,l1) - dky*uinp(1,l1)
               dixr(1) = diy*zr - diz*yr
               dixr(2) = diz*xr - dix*zr
               dixr(3) = dix*yr - diy*xr
               dkxr(1) = dky*zr - dkz*yr
               dkxr(2) = dkz*xr - dkx*zr
               dkxr(3) = dkx*yr - dky*xr
               qirx = qixx*xr + qiyx*yr + qizx*zr
               qiry = qixy*xr + qiyy*yr + qizy*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkyx*yr + qkzx*zr
               qkry = qkxy*xr + qkyy*yr + qkzy*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
c               qiqkr(1) = qixx*qkrx + qiyx*qkry + qizx*qkrz
c               qiqkr(2) = qixy*qkrx + qiyy*qkry + qizy*qkrz
c               qiqkr(3) = qixz*qkrx + qiyz*qkry + qizz*qkrz
c               qkqir(1) = qkxx*qirx + qkyx*qiry + qkzx*qirz
c               qkqir(2) = qkxy*qirx + qkyy*qiry + qkzy*qirz
c               qkqir(3) = qkxz*qirx + qkyz*qiry + qkzz*qirz
c               qixqkxx = qixy*qkxz + qiyy*qkyz + qizy*qkzz
c     &                       - qixz*qkxy - qiyz*qkyy - qizz*qkzy
c               qixqkxy = qixz*qkxx + qiyz*qkyx + qizz*qkzx
c     &                       - qixx*qkxz - qiyx*qkyz - qizx*qkzz
c               qixqkxz = qixx*qkxy + qiyx*qkyy + qizx*qkzy
c     &                       - qixy*qkxx - qiyy*qkyx - qizy*qkzx
               rxqir(1) = yr*qirz - zr*qiry
               rxqir(2) = zr*qirx - xr*qirz
               rxqir(3) = xr*qiry - yr*qirx
               rxqkr(1) = yr*qkrz - zr*qkry
               rxqkr(2) = zr*qkrx - xr*qkrz
               rxqkr(3) = xr*qkry - yr*qkrx
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkry*qirz - qkrz*qiry
c               qkrxqir(2) = qkrz*qirx - qkrx*qirz
c               qkrxqir(3) = qkrx*qiry - qkry*qirx
c               qidk(1) = qixx*dkx + qiyx*dky + qizx*dkz
c               qidk(2) = qixy*dkx + qiyy*dky + qizy*dkz
c               qidk(3) = qixz*dkx + qiyz*dky + qizz*dkz
c               qkdi(1) = qkxx*dix + qkyx*diy + qkzx*diz
c               qkdi(2) = qkxy*dix + qkyy*diy + qkzy*diz
c               qkdi(3) = qkxz*dix + qkyz*diy + qkzz*diz
               qiuk(1) = qixx*uind(1,l3) + qiyx*uind(2,l3)
     &                      + qizx*uind(3,l3)
               qiuk(2) = qixy*uind(1,l3) + qiyy*uind(2,l3)
     &                      + qizy*uind(3,l3)
               qiuk(3) = qixz*uind(1,l3) + qiyz*uind(2,l3)
     &                      + qizz*uind(3,l3)
               qkui(1) = qkxx*uind(1,l1) + qkyx*uind(2,l1)
     &                      + qkzx*uind(3,l1)
               qkui(2) = qkxy*uind(1,l1) + qkyy*uind(2,l1)
     &                      + qkzy*uind(3,l1)
               qkui(3) = qkxz*uind(1,l1) + qkyz*uind(2,l1)
     &                      + qkzz*uind(3,l1)
               qiukp(1) = qixx*uinp(1,l3) + qiyx*uinp(2,l3)
     &                       + qizx*uinp(3,l3)
               qiukp(2) = qixy*uinp(1,l3) + qiyy*uinp(2,l3)
     &                       + qizy*uinp(3,l3)
               qiukp(3) = qixz*uinp(1,l3) + qiyz*uinp(2,l3)
     &                       + qizz*uinp(3,l3)
               qkuip(1) = qkxx*uinp(1,l1) + qkyx*uinp(2,l1)
     &                       + qkzx*uinp(3,l1)
               qkuip(2) = qkxy*uinp(1,l1) + qkyy*uinp(2,l1)
     &                       + qkzy*uinp(3,l1)
               qkuip(3) = qkxz*uinp(1,l1) + qkyz*uinp(2,l1)
     &                       + qkzz*uinp(3,l1)
c               dixqkr(1) = diy*qkrz - diz*qkry
c               dixqkr(2) = diz*qkrx - dix*qkrz
c               dixqkr(3) = dix*qkry - diy*qkrx
c               dkxqir(1) = dky*qirz - dkz*qiry
c               dkxqir(2) = dkz*qirx - dkx*qirz
c               dkxqir(3) = dkx*qiry - dky*qirx
               uixqkr(1) = uind(2,l1)*qkrz - uind(3,l1)*qkry
               uixqkr(2) = uind(3,l1)*qkrx - uind(1,l1)*qkrz
               uixqkr(3) = uind(1,l1)*qkry - uind(2,l1)*qkrx
               ukxqir(1) = uind(2,l3)*qirz - uind(3,l3)*qiry
               ukxqir(2) = uind(3,l3)*qirx - uind(1,l3)*qirz
               ukxqir(3) = uind(1,l3)*qiry - uind(2,l3)*qirx
               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
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
c     calculate scalar products for permanent components
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qirx*xr + qiry*yr + qirz*zr
               sc(6) = qkrx*xr + qkry*yr + qkrz*zr
               sc(7) = qirx*dkx + qiry*dky + qirz*dkz
               sc(8) = qkrx*dix + qkry*diy + qkrz*diz
               sc(9) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(10) = qixx*qkxx + qixy*qkxy + qixz*qkxz
     &                     + qiyx*qkyx + qiyy*qkyy + qiyz*qkyz
     &                     + qizx*qkzx + qizy*qkzy + qizz*qkzz
c
c     calculate scalar products for induced components
c
               sci(1) = uind(1,l1)*dkx + uind(2,l1)*dky
     &                     + uind(3,l1)*dkz + dix*uind(1,l3)
     &                     + diy*uind(2,l3) + diz*uind(3,l3)
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = qirx*uind(1,l3) + qiry*uind(2,l3)
     &                     + qirz*uind(3,l3)
               sci(8) = qkrx*uind(1,l1) + qkry*uind(2,l1)
     &                     + qkrz*uind(3,l1)
               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
     &                      + qirz*uinp(3,l3)
               scip(8) = qkrx*uinp(1,l1) + qkry*uinp(2,l1)
     &                      + qkrz*uinp(3,l1)
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
               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = f * ei
               eptemp = eptemp + ei
c
c     intermediate variables for the induced components
c
               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
     &                                   + (glip(1)+glip(6))*dsc3
     &                                   + scip(2)*scale3i)
     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0 * rr9 * (gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force components
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*dix
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dkx
     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qirx + gfi(6)*qkrx
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
     &                                  +(glip(1)+glip(6))*dscale(l3))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
     &                                  +(glip(2)+glip(7))*dscale(l3))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
     &                                  +glip(3)*dscale(l3))



               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c

               temp3 = 0.5d0 * rr3 * uscale(l3) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(l3)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify induced force for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

c
c     correction to convert mutual to direct polarization force
c

               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
                  temp5 = 0.5d0 * rr5 * scale5i
                  fdir(1) = gfd*xr + temp5
     &                         * (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &                           +sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
                  fdir(2) = gfd*yr + temp5
     &                         * (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &                           +sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
                  fdir(3) = gfd*zr + temp5
     &                         * (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &                           +sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if

c
c     intermediate terms for induced torque on multipoles
c
               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * ftm2i(j)
                  ttm2i(j) = f * ttm2i(j)
                  ttm3i(j) = f * ttm3i(j)
               end do
c
c     increment gradient due to force and torque on first site
c
               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
               call torque_3b_new(npole3b,pnum,deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)
c
c     increment gradient due to force and torque on second site
c
               deptemp(1,l3) = deptemp(1,l3) - ftm2i(1)
               deptemp(2,l3) = deptemp(2,l3) - ftm2i(2)
               deptemp(3,l3) = deptemp(3,l3) - ftm2i(3)
               call torque_3b_new(npole3b,pnum,deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)
c
c     increment the internal virial tensor components
c
               iaz = iz
               iax = ix
               iay = iy
               kaz = kz
               kax = kx
               kay = ky
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

   10       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, npole3b
            pscale(j) = 1.0d0
            dscale(j) = 1.0d0
            uscale(j) = 1.0d0
         end do

      end do
c
c     perform deallocation of some local arrays
c

c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c       use_replica=.true.
c   NOTE TO SELF:  REMOVED if (use_replica) BECAUSE NCELL WAS
c   ALWAYS EQUAL TO ZERO AFTER PERIODIC BOUNDARY CONDITIONS WERE
c   IMPLEMENTED CORRECTLY.  ALSO, AS AN ADDITIONAL CHECK, THE ENERGY
c   IN THIS SECTION WAS EQUAL TO ZERO.
c
c     perform deallocation of some local arrays
c
c      print*,"eptemp",eptemp
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end


      subroutine induce0a_3b_PolelecOnly_rtrnfield(npole3b,pnum,uind,
     & uinp,field,fieldp)
      use sizes
      use atoms
      use polar, only: polarity, thole, pdamp
      use polpot
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 field(3,npole3b),off3b
      real*8 fieldp(3,npole3b)
      character*6 mode
      integer l1,l2,l3,k1,k2,i1,i2
      real*8 M_tot(3*npole3b,3*npole3b)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      integer npole3b,pnum(*)

c        print*,"npole3b=",npole3b
        do l1=1,npole3b
           do j=1,3
              field(j,l1)=0.0d0
              fieldp(j,l1)=0.0d0
           end do
        end do

         do l1 = 1, 3*npole3b
            do l3 = 1, 3*npole3b
               M_tot(l1,l3) = 0.0d0
            end do
         end do
         do l1 = 1, npole3b
            do j = 1, 3
               uind(j,l1) = 0.0d0
               uinp(j,l1) = 0.0d0
            end do
         end do

       call field_noewald_umutual_rl_3b (field,fieldp,M_tot,
     &   npole3b,pnum)

      if (poltyp .eq. 'MUTUAL') then
         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               i1 = 3*(l1-1)+j
               M_tot(i1,i1) = M_tot(i1,i1)+1.0d0/polarity(i)
            end do
c             print*,"polarity(i)",polarity(i)
         end do

         call invert(3*npole3b,M_tot)

         do l1 = 1, npole3b
            i = pnum(l1)
            do i1 = 1, 3
               i2 = 3*(l1-1) + i1
               do l3 = 1, npole3b
                  k = pnum(l3)
                  k2 = 3*(l3-1)
                  uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*field(1,l3)+
     &                                      M_tot(i2,k2+2)*field(2,l3)+
     &                                      M_tot(i2,k2+3)*field(3,l3)
                  uinp(i1,l1)=uinp(i1,l1)+ M_tot(i2,k2+1)*fieldp(1,l3)+
     &                                     M_tot(i2,k2+2)*fieldp(2,l3)+
     &                                     M_tot(i2,k2+3)*fieldp(3,l3)

               end do
            end do
         end do
      else
         do l1 = 1, npole3b
            i = pnum(l1)
            do j=1,3
               uind(j,l1)=polarity(i)*field(j,l1)
               uinp(j,l1)=polarity(i)*fieldp(j,l1)
            end do
         end do
      end if

      return
      end
c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
