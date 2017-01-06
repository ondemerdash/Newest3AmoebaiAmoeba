c
c
      subroutine empole1a1_3b_Polar_totfieldnpole_modfij(npole3b,pnum,
     &           uind,uinp,deptemp2,virtemp,deptempmat)
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
c      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      integer npole3b,pnum(*),l1,l3
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 eptemp,deptemp2(3,npole3b),virtemp(3,3)!,deptemp(3,npole3b)
c      real*8 uind(3,npole)
c      real*8 uinp(3,npole)
      real*8 off3b,eptemp_userep
      logical proceed,usei,usek
      character*6 mode
      real*8 deptempmat(3,npole,npole)
c      eptemp = 0.0d0
c      eptemp_userep =0.0d0
c
c   Zero out temporary gradient of polarization energy
c
c      off2=1.0d12
      !print*,"off2=",off2
      do l1 = 1, npole3b
        do j = 1, 3
          deptemp2(j,l1) = 0.0d0
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

c      call induce0a_3b_PolelecOnly_totfieldnpole3b(npole3b,pnum,uind,
c     &     uinp)
c
c     perform dynamic allocation of some local arrays
c


c      allocate (pscale(npole3b))
c      allocate (dscale(npole3b))
      allocate (uscale(npole3b))

c      allocate (pscale(npole))
c      allocate (uscale(npole))
c      allocate (dscale(npole))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, npole3b
c      do i = 1,npole
c         pscale(i) = 1.0d0
c         dscale(i) = 1.0d0
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
c      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
c         ci = rpole(1,i)
c         dix = rpole(2,i)
c         diy = rpole(3,i)
c         diz = rpole(4,i)
c         qixx = rpole(5,i)
c         qixy = rpole(6,i)
c         qixz = rpole(7,i)
c         qiyx = rpole(8,i)
c         qiyy = rpole(9,i)
c         qiyz = rpole(10,i)
c         qizx = rpole(11,i)
c         qizy = rpole(12,i)
c         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
c         do j = 1, n12(ii)
c            mscale(i12(j,ii)) = m2scale
c            pscale(i12(j,ii)) = p2scale
c         end do
c         do j = 1, n13(ii)
c            mscale(i13(j,ii)) = m3scale
c            pscale(i13(j,ii)) = p3scale
c         end do
c         do j = 1, n14(ii)
c            mscale(i14(j,ii)) = m4scale
c            pscale(i14(j,ii)) = p4scale
c            do k = 1, np11(ii)
c                if (i14(j,ii) .eq. ip11(k,ii))
c     &            pscale(i14(j,ii)) = p4scale * p41scale
c            end do
c         end do
c         do j = 1, n15(ii)
c            mscale(i15(j,ii)) = m5scale
c            pscale(i15(j,ii)) = p5scale
c         end do
c         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = d1scale
c            uscale(ip11(j,ii)) = u1scale
c         end do
c         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = d2scale
c            uscale(ip12(j,ii)) = u2scale
c         end do
c         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = d3scale
c            uscale(ip13(j,ii)) = u3scale
c         end do
c         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = d4scale
c            uscale(ip14(j,ii)) = u4scale
c         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
c                  dscale(kk)=d1scale
                  uscale(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
c                  dscale(kk)=d2scale
                  uscale(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
c                  dscale(kk)=d3scale
                  uscale(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
c                  dscale(kk)=d4scale
                  uscale(kk)=u4scale
                  goto 39
               end if
            end do
   39             continue
         end do

         do l3 = l1+1, npole3b
            k=pnum(l3)
c         do k=i+1,npole
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
            !call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr

            !if (r2 .le. off2) then
               r = sqrt(r2)
c               ck = rpole(1,k)
c               dkx = rpole(2,k)
c               dky = rpole(3,k)
c               dkz = rpole(4,k)
c               qkxx = rpole(5,k)
c               qkxy = rpole(6,k)
c               qkxz = rpole(7,k)
c               qkyx = rpole(8,k)
c               qkyy = rpole(9,k)
c               qkyz = rpole(10,k)
c               qkzx = rpole(11,k)
c               qkzy = rpole(12,k)
c               qkzz = rpole(13,k)
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
c               scale3i = scale3 * uscale(k)
c               scale5i = scale5 * uscale(k)
c               scale7i = scale7 * uscale(k)
c               dsc3 = scale3 * dscale(l3)
c               dsc5 = scale5 * dscale(l3)
c               dsc7 = scale7 * dscale(l3)
c               psc3 = scale3 * pscale(k)
c               psc5 = scale5 * pscale(k)
c               psc7 = scale7 * pscale(k)
c
c     construct necessary auxiliary vectors
c
c               dixuk(1) = diy*uind(3,l3) - diz*uind(2,l3)
c               dixuk(2) = diz*uind(1,l3) - dix*uind(3,l3)
c               dixuk(3) = dix*uind(2,l3) - diy*uind(1,l3)
c               dkxui(1) = dky*uind(3,l1) - dkz*uind(2,l1)
c               dkxui(2) = dkz*uind(1,l1) - dkx*uind(3,l1)
c               dkxui(3) = dkx*uind(2,l1) - dky*uind(1,l1)
c               dixukp(1) = diy*uinp(3,l3) - diz*uinp(2,l3)
c               dixukp(2) = diz*uinp(1,l3) - dix*uinp(3,l3)
c               dixukp(3) = dix*uinp(2,l3) - diy*uinp(1,l3)
c               dkxuip(1) = dky*uinp(3,l1) - dkz*uinp(2,l1)
c               dkxuip(2) = dkz*uinp(1,l1) - dkx*uinp(3,l1)
c               dkxuip(3) = dkx*uinp(2,l1) - dky*uinp(1,l1)
c               dixr(1) = diy*zr - diz*yr
c               dixr(2) = diz*xr - dix*zr
c               dixr(3) = dix*yr - diy*xr
c               dkxr(1) = dky*zr - dkz*yr
c               dkxr(2) = dkz*xr - dkx*zr
c               dkxr(3) = dkx*yr - dky*xr
c               qirx = qixx*xr + qiyx*yr + qizx*zr
c               qiry = qixy*xr + qiyy*yr + qizy*zr
c               qirz = qixz*xr + qiyz*yr + qizz*zr
c               qkrx = qkxx*xr + qkyx*yr + qkzx*zr
c               qkry = qkxy*xr + qkyy*yr + qkzy*zr
c               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
ccc               qiqkr(1) = qixx*qkrx + qiyx*qkry + qizx*qkrz
ccc               qiqkr(2) = qixy*qkrx + qiyy*qkry + qizy*qkrz
ccc               qiqkr(3) = qixz*qkrx + qiyz*qkry + qizz*qkrz
ccc               qkqir(1) = qkxx*qirx + qkyx*qiry + qkzx*qirz
ccc               qkqir(2) = qkxy*qirx + qkyy*qiry + qkzy*qirz
ccc               qkqir(3) = qkxz*qirx + qkyz*qiry + qkzz*qirz
ccc               qixqkxx = qixy*qkxz + qiyy*qkyz + qizy*qkzz
ccc     &                       - qixz*qkxy - qiyz*qkyy - qizz*qkzy
ccc               qixqkxy = qixz*qkxx + qiyz*qkyx + qizz*qkzx
ccc     &                       - qixx*qkxz - qiyx*qkyz - qizx*qkzz
ccc               qixqkxz = qixx*qkxy + qiyx*qkyy + qizx*qkzy
ccc     &                       - qixy*qkxx - qiyy*qkyx - qizy*qkzx
c               rxqir(1) = yr*qirz - zr*qiry
c               rxqir(2) = zr*qirx - xr*qirz
c               rxqir(3) = xr*qiry - yr*qirx
c               rxqkr(1) = yr*qkrz - zr*qkry
c               rxqkr(2) = zr*qkrx - xr*qkrz
c               rxqkr(3) = xr*qkry - yr*qkrx
ccc               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
ccc               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
ccc               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
ccc               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
ccc               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
ccc               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
ccc               qkrxqir(1) = qkry*qirz - qkrz*qiry
ccc               qkrxqir(2) = qkrz*qirx - qkrx*qirz
ccc               qkrxqir(3) = qkrx*qiry - qkry*qirx
ccc               qidk(1) = qixx*dkx + qiyx*dky + qizx*dkz
ccc               qidk(2) = qixy*dkx + qiyy*dky + qizy*dkz
ccc               qidk(3) = qixz*dkx + qiyz*dky + qizz*dkz
ccc               qkdi(1) = qkxx*dix + qkyx*diy + qkzx*diz
ccc               qkdi(2) = qkxy*dix + qkyy*diy + qkzy*diz
ccc               qkdi(3) = qkxz*dix + qkyz*diy + qkzz*diz
c               qiuk(1) = qixx*uind(1,l3) + qiyx*uind(2,l3)
c     &                      + qizx*uind(3,l3)
c               qiuk(2) = qixy*uind(1,l3) + qiyy*uind(2,l3)
c     &                      + qizy*uind(3,l3)
c               qiuk(3) = qixz*uind(1,l3) + qiyz*uind(2,l3)
c     &                      + qizz*uind(3,l3)
c               qkui(1) = qkxx*uind(1,l1) + qkyx*uind(2,l1)
c     &                      + qkzx*uind(3,l1)
c               qkui(2) = qkxy*uind(1,l1) + qkyy*uind(2,l1)
c     &                      + qkzy*uind(3,l1)
c               qkui(3) = qkxz*uind(1,l1) + qkyz*uind(2,l1)
c     &                      + qkzz*uind(3,l1)
c               qiukp(1) = qixx*uinp(1,l3) + qiyx*uinp(2,l3)
c     &                       + qizx*uinp(3,l3)
c               qiukp(2) = qixy*uinp(1,l3) + qiyy*uinp(2,l3)
c     &                       + qizy*uinp(3,l3)
c               qiukp(3) = qixz*uinp(1,l3) + qiyz*uinp(2,l3)
c     &                       + qizz*uinp(3,l3)
c               qkuip(1) = qkxx*uinp(1,l1) + qkyx*uinp(2,l1)
c     &                       + qkzx*uinp(3,l1)
c               qkuip(2) = qkxy*uinp(1,l1) + qkyy*uinp(2,l1)
c     &                       + qkzy*uinp(3,l1)
c               qkuip(3) = qkxz*uinp(1,l1) + qkyz*uinp(2,l1)
c     &                       + qkzz*uinp(3,l1)
ccc               dixqkr(1) = diy*qkrz - diz*qkry
ccc               dixqkr(2) = diz*qkrx - dix*qkrz
ccc               dixqkr(3) = dix*qkry - diy*qkrx
ccc               dkxqir(1) = dky*qirz - dkz*qiry
ccc               dkxqir(2) = dkz*qirx - dkx*qirz
ccc               dkxqir(3) = dkx*qiry - dky*qirx
c               uixqkr(1) = uind(2,l1)*qkrz - uind(3,l1)*qkry
c               uixqkr(2) = uind(3,l1)*qkrx - uind(1,l1)*qkrz
c               uixqkr(3) = uind(1,l1)*qkry - uind(2,l1)*qkrx
c               ukxqir(1) = uind(2,l3)*qirz - uind(3,l3)*qiry
c               ukxqir(2) = uind(3,l3)*qirx - uind(1,l3)*qirz
c               ukxqir(3) = uind(1,l3)*qiry - uind(2,l3)*qirx
c               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
c               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
c               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
c               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
c               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
c               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
ccc               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
ccc               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
ccc               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
ccc               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
ccc               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
ccc               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
c               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
c               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
c               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
c               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
c               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
c               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
c               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
c               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
c               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
c               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
c               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate scalar products for permanent components
c
c               sc(2) = dix*dkx + diy*dky + diz*dkz
c               sc(3) = dix*xr + diy*yr + diz*zr
c               sc(4) = dkx*xr + dky*yr + dkz*zr
c               sc(5) = qirx*xr + qiry*yr + qirz*zr
c               sc(6) = qkrx*xr + qkry*yr + qkrz*zr
c               sc(7) = qirx*dkx + qiry*dky + qirz*dkz
c               sc(8) = qkrx*dix + qkry*diy + qkrz*diz
c               sc(9) = qirx*qkrx + qiry*qkry + qirz*qkrz
c               sc(10) = qixx*qkxx + qixy*qkxy + qixz*qkxz
c     &                     + qiyx*qkyx + qiyy*qkyy + qiyz*qkyz
c     &                     + qizx*qkzx + qizy*qkzy + qizz*qkzz
c
c     calculate scalar products for induced components
c

c               sci(1) = uind(1,l1)*dkx + uind(2,l1)*dky
c     &                     + uind(3,l1)*dkz + dix*uind(1,l3)
c     &                     + diy*uind(2,l3) + diz*uind(3,l3)
               sci(1) = 0.0d0             
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
c               sci(7) = qirx*uind(1,l3) + qiry*uind(2,l3)
c     &                     + qirz*uind(3,l3)
               sci(7) =0.0d0 
c               sci(8) = qkrx*uind(1,l1) + qkry*uind(2,l1)
c     &                     + qkrz*uind(3,l1)
               sci(8) =0.0d0
c               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
c     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
c     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
               scip(1)=0.0d0
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
c               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
c     &                      + qirz*uinp(3,l3)
c               scip(8) = qkrx*uinp(1,l1) + qkry*uinp(2,l1)
c     &                      + qkrz*uinp(3,l1)
               scip(7) = 0.0d0
               scip(8) =0.0d0
c
c     calculate the gl functions for induced components
c
c               gli(1) = ck*sci(3) - ci*sci(4)
c               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
c               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c               gli(6) = sci(1)
c               gli(7) = 2.0d0 * (sci(7)-sci(8))
c               glip(1) = ck*scip(3) - ci*scip(4)
c               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
c               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
c               glip(6) = scip(1)
c               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
c               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
c     &                   + rr5*(gli(2)+gli(7))*psc5
c     &                   + rr7*gli(3)*psc7)
c               ei = f * ei
c               eptemp = eptemp + ei
c
c     intermediate variables for the induced components
c
c
c               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
c     &                                   + (glip(1)+glip(6))*dsc3
c     &                                   + scip(2)*scale3i)
c     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
c     &                               + (glip(7)+glip(2))*dsc5
c     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
c     &                + 0.5d0 * rr9 * (gli(3)*psc7+glip(3)*dsc7)
c               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
c               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
c               gfi(4) = 2.0d0 * rr5
c               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
c               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)

c               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
c     &                                   + scip(2)*scale3i)
c     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
c     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
c     &                + 0.5d0 * rr9 * (gli(3)*psc7)
c               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
c               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
c               gfi(4) = 2.0d0 * rr5
c               gfi(5) = rr7 * (sci(4)*psc7)
c               gfi(6) = -rr7 * (sci(3)*psc7)

               gfi(1) = 0.5d0 * rr5 * (
     &                                    scip(2)*scale3i)
     &                + 0.5d0 * rr7 * (
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.0d0
               gfi(2) = 0.0d0
               gfi(3) = 0.0d0
               gfi(4) = 2.0d0 * rr5
               gfi(5) = 0.0d0
               gfi(6) = 0.0d0


c
c     get the induced force components
c

c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
c     &            + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
c     &            + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
c     &            + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*dix
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dkx
c     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
c     &            + (qkuip(1)-qiukp(1))*dsc5)
c     &            + gfi(5)*qirx + gfi(6)*qkrx
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
c     &            + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
c     &            + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
c     &            + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*diy
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dky
c     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
c     &            + (qkuip(2)-qiukp(2))*dsc5)
c     &            + gfi(5)*qiry + gfi(6)*qkry
c               ftm2i(3) = gfi(1)*zr  + 0.5d0*
c     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
c     &            + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
c     &            + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
c     &            + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*diz
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dkz
c     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
c     &            + (qkuip(3)-qiukp(3))*dsc5)
c     &            + gfi(5)*qirz + gfi(6)*qkrz

c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &           (- rr3*ck*(uind(1,l1)*psc3)
c     &            + rr5*sc(4)*(uind(1,l1)*psc5)
c     &            - rr7*sc(6)*(uind(1,l1)*psc7))
c     &            + (rr3*ci*(uind(1,l3)*psc3)
c     &            + rr5*sc(3)*(uind(1,l3)*psc5)
c     &            + rr7*sc(5)*(uind(1,l3)*psc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5)*rr5*dix
c     &            + 0.5d0*(sci(3)*psc5)*rr5*dkx
c     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
c     &            )
c     &            + gfi(5)*qirx + gfi(6)*qkrx
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &           (- rr3*ck*(uind(2,l1)*psc3)
c     &            + rr5*sc(4)*(uind(2,l1)*psc5)
c     &            - rr7*sc(6)*(uind(2,l1)*psc7))
c     &            + (rr3*ci*(uind(2,l3)*psc3)
c     &            + rr5*sc(3)*(uind(2,l3)*psc5)
c     &            + rr7*sc(5)*(uind(2,l3)*psc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5)*rr5*diy
c     &            + 0.5d0*(sci(3)*psc5)*rr5*dky
c     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
c     &            )
c     &            + gfi(5)*qiry + gfi(6)*qkry
c               ftm2i(3) = gfi(1)*zr  + 0.5d0*
c     &           (- rr3*ck*(uind(3,l1)*psc3)
c     &            + rr5*sc(4)*(uind(3,l1)*psc5)
c     &            - rr7*sc(6)*(uind(3,l1)*psc7))
c     &            + (rr3*ci*(uind(3,l3)*psc3)
c     &            + rr5*sc(3)*(uind(3,l3)*psc5)
c     &            + rr7*sc(5)*(uind(3,l3)*psc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5)*rr5*diz
c     &            + 0.5d0*(sci(3)*psc5)*rr5*dkz
c     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
c     &            )
c     &            + gfi(5)*qirz + gfi(6)*qkrz


               ftm2i(1) = gfi(1)*xr 
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
               ftm2i(2) = gfi(1)*yr  
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
               ftm2i(3) = gfi(1)*zr  
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0


c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  +(glip(1)+glip(6))*dscale(l3))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  +(glip(2)+glip(7))*dscale(l3))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  +glip(3)*dscale(l3))

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(k)
c     &                                  )
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(k)
c     &                                  )
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(k)
c     &                                  )

               temp3 =0.0d0
               temp5 =0.0d0
               temp7 =0.0d0
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c

c               temp3 = 0.5d0 * rr3 * uscale(k) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(k)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
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
c               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
c               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
c               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5)
c               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5)
               gti(2) = 0.0d0
               gti(3) = 0.0d0
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the induced torque components
c

c               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
c     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
c               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
c     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
c               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
c     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
c               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
c     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
c               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
c     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
c               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
c     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)


c               ttm2i(1) = -rr3*(dixuk(1)*psc3)*0.5d0
c     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &           )*0.5d0 - gti(5)*rxqir(1)
c               ttm2i(2) = -rr3*(dixuk(2)*psc3)*0.5d0
c     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &           )*0.5d0 - gti(5)*rxqir(2)
c               ttm2i(3) = -rr3*(dixuk(3)*psc3)*0.5d0
c     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &           )*0.5d0 - gti(5)*rxqir(3)
c               ttm3i(1) = -rr3*(dkxui(1)*psc3)*0.5d0
c     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &           )*0.5d0 - gti(6)*rxqkr(1)
c               ttm3i(2) = -rr3*(dkxui(2)*psc3)*0.5d0
c     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &           )*0.5d0 - gti(6)*rxqkr(2)
c               ttm3i(3) = -rr3*(dkxui(3)*psc3)*0.5d0
c     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &           )*0.5d0 - gti(6)*rxqkr(3)


c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * ftm2i(j)
c                  ttm2i(j) = f * ttm2i(j)
c                  ttm3i(j) = f * ttm3i(j)
               end do
c
c     increment gradient due to force and torque on first site
c

c               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
c               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
c               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
c               call torque_3b_new(npole3b,pnum,deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)

               deptemp2(1,l1) = deptemp2(1,l1) + ftm2i(1)
               deptemp2(2,l1) = deptemp2(2,l1) + ftm2i(2)
               deptemp2(3,l1) = deptemp2(3,l1) + ftm2i(3)
               deptempmat(1,k,i)=deptempmat(1,k,i)+ftm2i(1)
               deptempmat(2,k,i)=deptempmat(2,k,i)+ftm2i(2)
               deptempmat(3,k,i)=deptempmat(3,k,i)+ftm2i(3)
c               call torque_3b_new_npole(deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
c               do j=1,3
c                frcxi(j) = 0.0d0
c                frcyi(j) = 0.0d0
c                frczi(j) = 0.0d0
c               end do
c
c     increment gradient due to force and torque on second site
c

c               deptemp(1,l3) = deptemp(1,l3) - ftm2i(1)
c               deptemp(2,l3) = deptemp(2,l3) - ftm2i(2)
c               deptemp(3,l3) = deptemp(3,l3) - ftm2i(3)
c               call torque_3b_new(npole3b,pnum,deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               deptemp2(1,l3) = deptemp2(1,l3) - ftm2i(1)
               deptemp2(2,l3) = deptemp2(2,l3) - ftm2i(2)
               deptemp2(3,l3) = deptemp2(3,l3) - ftm2i(3)
               deptempmat(1,i,k)= deptempmat(1,i,k)-ftm2i(1)
               deptempmat(2,i,k)= deptempmat(2,i,k)-ftm2i(2)
               deptempmat(3,i,k)= deptempmat(3,i,k)-ftm2i(3)
c               call torque_3b_new_npole(deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
c               do j=1,3
c                frcxk(j) = 0.0d0
c                frcyk(j) = 0.0d0
c                frczk(j) = 0.0d0
c               end do
c
c     increment the internal virial tensor components
c
c               iaz = iz
c               iax = ix
c               iay = iy
c               kaz = kz
c               kax = kx
c               kay = ky
c               if (iaz .eq. 0)  iaz = ii
c               if (iax .eq. 0)  iax = ii
c               if (iay .eq. 0)  iay = ii
c               if (kaz .eq. 0)  kaz = kk
c               if (kax .eq. 0)  kax = kk
c               if (kay .eq. 0)  kay = kk
c               xiz = x(iaz) - x(ii)
c               yiz = y(iaz) - y(ii)
c               ziz = z(iaz) - z(ii)
c               xix = x(iax) - x(ii)
c               yix = y(iax) - y(ii)
c               zix = z(iax) - z(ii)
c               xiy = x(iay) - x(ii)
c               yiy = y(iay) - y(ii)
c               ziy = z(iay) - z(ii)
c               xkz = x(kaz) - x(kk)
c               ykz = y(kaz) - y(kk)
c               zkz = z(kaz) - z(kk)
c               xkx = x(kax) - x(kk)
c               ykx = y(kax) - y(kk)
c               zkx = z(kax) - z(kk)
c               xky = x(kay) - x(kk)
c               yky = y(kay) - y(kk)
c               zky = z(kay) - z(kk)
               vxx = -xr*(ftm2i(1)) !+ xix*frcxi(1)
     &             !     + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &             !     + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) !+ yix*frcxi(1)
     &             !     + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &             !     + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) !+ zix*frcxi(1)
     &             !     + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &             !     + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) !+ yix*frcxi(2)
     &             !     + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &             !     + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) !+ zix*frcxi(2)
     &             !     + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &             !     + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) !+ zix*frcxi(3)
     &             !     + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &             !     + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz
            !end if

   10       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         !do j = 1, npole
         do j = 1, npole3b
c            pscale(j) = 1.0d0
c            dscale(j) = 1.0d0
            uscale(j) = 1.0d0
         end do

      end do
c
c

c
c
       use_replica=.false.
      if (use_replica) then

      do l1 = 1, npole3b
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
         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do l3 = l1, npole3b
            k=pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 20
            do jcell = 1, ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            !call imager (xr,yr,zr,jcell)
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
c               scale3i = scale3 * uscale(l3)
c               scale5i = scale5 * uscale(l3)
c               scale7i = scale7 * uscale(l3)
c               dsc3 = scale3 * dscale(l3)
c               dsc5 = scale5 * dscale(l3)
c               dsc7 = scale7 * dscale(l3)
c               psc3 = scale3 * pscale(l3)
c               psc5 = scale5 * pscale(l3)
c               psc7 = scale7 * pscale(l3)

               scale3i = scale3
               scale5i = scale5
               scale7i = scale7
c               dsc3 = scale3
c               dsc5 = scale5
c               dsc7 = scale7
               psc3 = scale3
               psc5 = scale5
               psc7 = scale7
               if (use_polymer) then
                  if (r2 .le. polycut2) then
                     scale3i = scale3i * uscale(k)
                     scale5i = scale5i * uscale(k)
                     scale7i = scale7i * uscale(k)
c                     dsc3 = dsc3 * dscale(kk)
c                     dsc5 = dsc5 * dscale(kk)
c                     dsc7 = dsc7 * dscale(kk)
c                     dsc3 = dsc3 * dscale(l3)
c                     dsc5 = dsc5 * dscale(l3)
c                     dsc7 = dsc7 * dscale(l3)
c                     psc3 = psc3 * pscale(l3)
c                     psc5 = psc5 * pscale(l3)
c                     psc7 = psc7 * pscale(l3)
                  end if
               end if



c
c     construct necessary auxiliary vectors
c
c               dixuk(1) = diy*uind(3,l3) - diz*uind(2,l3)
c               dixuk(2) = diz*uind(1,l3) - dix*uind(3,l3)
c               dixuk(3) = dix*uind(2,l3) - diy*uind(1,l3)
c               dkxui(1) = dky*uind(3,l1) - dkz*uind(2,l1)
c               dkxui(2) = dkz*uind(1,l1) - dkx*uind(3,l1)
c               dkxui(3) = dkx*uind(2,l1) - dky*uind(1,l1)
c               dixukp(1) = diy*uinp(3,l3) - diz*uinp(2,l3)
c               dixukp(2) = diz*uinp(1,l3) - dix*uinp(3,l3)
c               dixukp(3) = dix*uinp(2,l3) - diy*uinp(1,l3)
c               dkxuip(1) = dky*uinp(3,l1) - dkz*uinp(2,l1)
c               dkxuip(2) = dkz*uinp(1,l1) - dkx*uinp(3,l1)
c               dkxuip(3) = dkx*uinp(2,l1) - dky*uinp(1,l1)
c               dixr(1) = diy*zr - diz*yr
c               dixr(2) = diz*xr - dix*zr
c               dixr(3) = dix*yr - diy*xr
c               dkxr(1) = dky*zr - dkz*yr
c               dkxr(2) = dkz*xr - dkx*zr
c               dkxr(3) = dkx*yr - dky*xr
c               qirx = qixx*xr + qiyx*yr + qizx*zr
c               qiry = qixy*xr + qiyy*yr + qizy*zr
c               qirz = qixz*xr + qiyz*yr + qizz*zr
c               qkrx = qkxx*xr + qkyx*yr + qkzx*zr
c               qkry = qkxy*xr + qkyy*yr + qkzy*zr
c               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
ccc               qiqkr(1) = qixx*qkrx + qiyx*qkry + qizx*qkrz
ccc               qiqkr(2) = qixy*qkrx + qiyy*qkry + qizy*qkrz
ccc               qiqkr(3) = qixz*qkrx + qiyz*qkry + qizz*qkrz
ccc               qkqir(1) = qkxx*qirx + qkyx*qiry + qkzx*qirz
ccc               qkqir(2) = qkxy*qirx + qkyy*qiry + qkzy*qirz
ccc               qkqir(3) = qkxz*qirx + qkyz*qiry + qkzz*qirz
ccc               qixqkxx = qixy*qkxz + qiyy*qkyz + qizy*qkzz
ccc     &                       - qixz*qkxy - qiyz*qkyy - qizz*qkzy
ccc               qixqkxy = qixz*qkxx + qiyz*qkyx + qizz*qkzx
ccc     &                       - qixx*qkxz - qiyx*qkyz - qizx*qkzz
ccc               qixqkxz = qixx*qkxy + qiyx*qkyy + qizx*qkzy
ccc     &                       - qixy*qkxx - qiyy*qkyx - qizy*qkzx
c               rxqir(1) = yr*qirz - zr*qiry
c               rxqir(2) = zr*qirx - xr*qirz
c               rxqir(3) = xr*qiry - yr*qirx
c               rxqkr(1) = yr*qkrz - zr*qkry
c               rxqkr(2) = zr*qkrx - xr*qkrz
c               rxqkr(3) = xr*qkry - yr*qkrx
ccc               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
ccc               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
ccc               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
ccc               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
ccc               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
ccc               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
ccc               qkrxqir(1) = qkry*qirz - qkrz*qiry
ccc               qkrxqir(2) = qkrz*qirx - qkrx*qirz
ccc               qkrxqir(3) = qkrx*qiry - qkry*qirx
ccc               qidk(1) = qixx*dkx + qiyx*dky + qizx*dkz
ccc               qidk(2) = qixy*dkx + qiyy*dky + qizy*dkz
ccc               qidk(3) = qixz*dkx + qiyz*dky + qizz*dkz
ccc               qkdi(1) = qkxx*dix + qkyx*diy + qkzx*diz
ccc               qkdi(2) = qkxy*dix + qkyy*diy + qkzy*diz
ccc               qkdi(3) = qkxz*dix + qkyz*diy + qkzz*diz
c               qiuk(1) = qixx*uind(1,l3) + qiyx*uind(2,l3)
c     &                      + qizx*uind(3,l3)
c               qiuk(2) = qixy*uind(1,l3) + qiyy*uind(2,l3)
c     &                      + qizy*uind(3,l3)
c               qiuk(3) = qixz*uind(1,l3) + qiyz*uind(2,l3)
c     &                      + qizz*uind(3,l3)
c               qkui(1) = qkxx*uind(1,l1) + qkyx*uind(2,l1)
c     &                      + qkzx*uind(3,l1)
c               qkui(2) = qkxy*uind(1,l1) + qkyy*uind(2,l1)
c     &                      + qkzy*uind(3,l1)
c               qkui(3) = qkxz*uind(1,l1) + qkyz*uind(2,l1)
c     &                      + qkzz*uind(3,l1)
c               qiukp(1) = qixx*uinp(1,l3) + qiyx*uinp(2,l3)
c     &                       + qizx*uinp(3,l3)
c               qiukp(2) = qixy*uinp(1,l3) + qiyy*uinp(2,l3)
c     &                       + qizy*uinp(3,l3)
c               qiukp(3) = qixz*uinp(1,l3) + qiyz*uinp(2,l3)
c     &                       + qizz*uinp(3,l3)
c               qkuip(1) = qkxx*uinp(1,l1) + qkyx*uinp(2,l1)
c     &                       + qkzx*uinp(3,l1)
c               qkuip(2) = qkxy*uinp(1,l1) + qkyy*uinp(2,l1)
c     &                       + qkzy*uinp(3,l1)
c               qkuip(3) = qkxz*uinp(1,l1) + qkyz*uinp(2,l1)
c     &                       + qkzz*uinp(3,l1)
ccc               dixqkr(1) = diy*qkrz - diz*qkry
ccc               dixqkr(2) = diz*qkrx - dix*qkrz
ccc               dixqkr(3) = dix*qkry - diy*qkrx
ccc               dkxqir(1) = dky*qirz - dkz*qiry
ccc               dkxqir(2) = dkz*qirx - dkx*qirz
ccc               dkxqir(3) = dkx*qiry - dky*qirx
c               uixqkr(1) = uind(2,l1)*qkrz - uind(3,l1)*qkry
c               uixqkr(2) = uind(3,l1)*qkrx - uind(1,l1)*qkrz
c               uixqkr(3) = uind(1,l1)*qkry - uind(2,l1)*qkrx
c               ukxqir(1) = uind(2,l3)*qirz - uind(3,l3)*qiry
c               ukxqir(2) = uind(3,l3)*qirx - uind(1,l3)*qirz
c               ukxqir(3) = uind(1,l3)*qiry - uind(2,l3)*qirx
c               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
c               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
c               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
c               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
c               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
c               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
ccc               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
ccc               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
ccc               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
ccc               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
ccc               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
ccc               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
c               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
c               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
c               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
c               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
c               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
c               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
c               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
c               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
c               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
c               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
c               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate scalar products for permanent components
c

c               sc(2) = dix*dkx + diy*dky + diz*dkz
c               sc(3) = dix*xr + diy*yr + diz*zr
c               sc(4) = dkx*xr + dky*yr + dkz*zr
c               sc(5) = qirx*xr + qiry*yr + qirz*zr
c               sc(6) = qkrx*xr + qkry*yr + qkrz*zr
c               sc(7) = qirx*dkx + qiry*dky + qirz*dkz
c               sc(8) = qkrx*dix + qkry*diy + qkrz*diz
c               sc(9) = qirx*qkrx + qiry*qkry + qirz*qkrz
c               sc(10) = qixx*qkxx + qixy*qkxy + qixz*qkxz
c     &                     + qiyx*qkyx + qiyy*qkyy + qiyz*qkyz
c     &                     + qizx*qkzx + qizy*qkzy + qizz*qkzz

c
c     calculate scalar products for induced components
c

c               sci(1) = uind(1,l1)*dkx + uind(2,l1)*dky
c     &                     + uind(3,l1)*dkz + dix*uind(1,l3)
c     &                     + diy*uind(2,l3) + diz*uind(3,l3)
               sci(1) = 0.0d0
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = 0.0d0
               sci(8) = 0.0d0
c               sci(7) = qirx*uind(1,l3) + qiry*uind(2,l3)
c     &                     + qirz*uind(3,l3)
c               sci(8) = qkrx*uind(1,l1) + qkry*uind(2,l1)
c     &                     + qkrz*uind(3,l1)
c               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
c     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
c     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
               scip(1) = 0.0d0
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
c               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
c     &                      + qirz*uinp(3,l3)
c               scip(8) = qkrx*uinp(1,l1) + qkry*uinp(2,l1)
c     &                      + qkrz*uinp(3,l1)
               scip(7) = 0.0d0
               scip(8) = 0.0d0
c
c     calculate the gl functions for induced components
c
c               gli(1) = ck*sci(3) - ci*sci(4)
c               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
c               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c               gli(6) = sci(1)
c               gli(7) = 2.0d0 * (sci(7)-sci(8))
c               glip(1) = ck*scip(3) - ci*scip(4)
c               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
c               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
c               glip(6) = scip(1)
c               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
c               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
c     &                   + rr5*(gli(2)+gli(7))*psc5
c     &                   + rr7*gli(3)*psc7)
c               ei = f * ei
c               eptemp = eptemp + ei
c
c     intermediate variables for the induced components
c

c               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
c     &                                   + (glip(1)+glip(6))*dsc3
c     &                                   + scip(2)*scale3i)
c     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
c     &                               + (glip(7)+glip(2))*dsc5
c     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
c     &                + 0.5d0 * rr9 * (gli(3)*psc7+glip(3)*dsc7)
c               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
c               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
c               gfi(4) = 2.0d0 * rr5
c               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
c               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)

               gfi(1) = 0.5d0 * rr5 * (
     &                                    scip(2)*scale3i)
     &                + 0.5d0 * rr7 * (
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.0d0
               gfi(2) = 0.0d0
               gfi(3) = 0.0d0
               gfi(4) = 2.0d0 * rr5
               gfi(5) = 0.0d0
               gfi(6) = 0.0d0

c
c     get the induced force components
c

c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
c     &            + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
c     &            + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
c     &            + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*dix
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dkx
c     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
c     &            + (qkuip(1)-qiukp(1))*dsc5)
c     &            + gfi(5)*qirx + gfi(6)*qkrx
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
c     &            + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
c     &            + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
c     &            + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*diy
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dky
c     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
c     &            + (qkuip(2)-qiukp(2))*dsc5)
c     &            + gfi(5)*qiry + gfi(6)*qkry
c               ftm2i(3) = gfi(1)*zr  + 0.5d0*
c     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
c     &            + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
c     &            + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
c     &            + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*diz
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dkz
c     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
c     &            + (qkuip(3)-qiukp(3))*dsc5)
c     &            + gfi(5)*qirz + gfi(6)*qkrz

c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &           (- rr3*ck*(uind(1,l1)*psc3)
c     &            + rr5*sc(4)*(uind(1,l1)*psc5)
c     &            - rr7*sc(6)*(uind(1,l1)*psc7))
c     &            + (rr3*ci*(uind(1,l3)*psc3)
c     &            + rr5*sc(3)*(uind(1,l3)*psc5)
c     &            + rr7*sc(5)*(uind(1,l3)*psc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5)*rr5*dix
c     &            + 0.5d0*(sci(3)*psc5)*rr5*dkx
c     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
c     &            )
c     &            + gfi(5)*qirx + gfi(6)*qkrx
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &           (- rr3*ck*(uind(2,l1)*psc3)
c     &            + rr5*sc(4)*(uind(2,l1)*psc5)
c     &            - rr7*sc(6)*(uind(2,l1)*psc7))
c     &            + (rr3*ci*(uind(2,l3)*psc3)
c     &            + rr5*sc(3)*(uind(2,l3)*psc5)
c     &            + rr7*sc(5)*(uind(2,l3)*psc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5)*rr5*diy
c     &            + 0.5d0*(sci(3)*psc5)*rr5*dky
c     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
c     &            )
c     &            + gfi(5)*qiry + gfi(6)*qkry
c               ftm2i(3) = gfi(1)*zr  + 0.5d0*
c     &           (- rr3*ck*(uind(3,l1)*psc3)
c     &            + rr5*sc(4)*(uind(3,l1)*psc5)
c     &            - rr7*sc(6)*(uind(3,l1)*psc7))
c     &            + (rr3*ci*(uind(3,l3)*psc3)
c     &            + rr5*sc(3)*(uind(3,l3)*psc5)
c     &            + rr7*sc(5)*(uind(3,l3)*psc7))*0.5d0
c     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5)*rr5*diz
c     &            + 0.5d0*(sci(3)*psc5)*rr5*dkz
c     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
c     &            )
c     &            + gfi(5)*qirz + gfi(6)*qkrz

               ftm2i(1) = gfi(1)*xr
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
               ftm2i(2) = gfi(1)*yr
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
               ftm2i(3) = gfi(1)*zr
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0


c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  +(glip(1)+glip(6))*dscale(l3))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  +(glip(2)+glip(7))*dscale(l3))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  +glip(3)*dscale(l3))

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  )
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  )
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  )
               temp3 = 0.0d0
               temp5 = 0.0d0
               temp7 = 0.0d0
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c

               temp3 = 0.5d0 * rr3 * uscale(k) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(k)
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
c               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
c               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
c               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5)
c               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5)
c               gti(4) = gfi(4)
c               gti(5) = gfi(5)
c               gti(6) = gfi(6)

               gti(2) = 0.0d0
               gti(3) = 0.0d0
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the induced torque components
c

c               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
c     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
c               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
c     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
c               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
c     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
c               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
c     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
c               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
c     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
c               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
c     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)


c               ttm2i(1) = -rr3*(dixuk(1)*psc3)*0.5d0
c     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &           )*0.5d0 - gti(5)*rxqir(1)
c               ttm2i(2) = -rr3*(dixuk(2)*psc3)*0.5d0
c     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &           )*0.5d0 - gti(5)*rxqir(2)
c               ttm2i(3) = -rr3*(dixuk(3)*psc3)*0.5d0
c     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &           )*0.5d0 - gti(5)*rxqir(3)
c               ttm3i(1) = -rr3*(dkxui(1)*psc3)*0.5d0
c     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &           )*0.5d0 - gti(6)*rxqkr(1)
c               ttm3i(2) = -rr3*(dkxui(2)*psc3)*0.5d0
c     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &           )*0.5d0 - gti(6)*rxqkr(2)
c               ttm3i(3) = -rr3*(dkxui(3)*psc3)*0.5d0
c     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &           )*0.5d0 - gti(6)*rxqkr(3)

c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * ftm2i(j)
c                  ttm2i(j) = f * ttm2i(j)
c                  ttm3i(j) = f * ttm3i(j)
               end do

               if (ii .eq. kk) then
                  do j = 1, 3
                     ftm2i(j) = 0.5d0 * ftm2i(j)
c                     ttm2i(j) = 0.5d0 * ttm2i(j)
c                     ttm3i(j) = 0.5d0 * ttm3i(j)
                  end do
               end if


c
c     increment gradient due to force and torque on first site
c
               deptemp2(1,l1) = deptemp2(1,l1) + ftm2i(1)
               deptemp2(2,l1) = deptemp2(2,l1) + ftm2i(2)
               deptemp2(3,l1) = deptemp2(3,l1) + ftm2i(3)
c               call torque_3b_new_npole(deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               frcxi = 0.0d0
               frcyi = 0.0d0
               frczi = 0.0d0
c
c     increment gradient due to force and torque on second site
c

c               deptemp(1,l3) = deptemp(1,l3) - ftm2i(1)
c               deptemp(2,l3) = deptemp(2,l3) - ftm2i(2)
c               deptemp(3,l3) = deptemp(3,l3) - ftm2i(3)
c               call torque_3b_new(npole3b,pnum,deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               deptemp2(1,l3) = deptemp2(1,l3) - ftm2i(1)
               deptemp2(2,l3) = deptemp2(2,l3) - ftm2i(2)
               deptemp2(3,l3) = deptemp2(3,l3) - ftm2i(3)
c               call torque_3b_new_npole(deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)
               frcxk = 0.0d0
               frcyk = 0.0d0
               frczk = 0.0d0


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
            end do 
   20       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, npole
c            pscale(j) = 1.0d0
c            dscale(j) = 1.0d0
            uscale(j) = 1.0d0
         end do

      end do
      end if
c
c     perform deallocation of some local arrays
c
c
c
c     perform deallocation of some local arrays
c
c      print*,"eptemp",eptemp
c      deallocate (pscale)
c      deallocate (dscale)
      deallocate (uscale)
      return
      end

