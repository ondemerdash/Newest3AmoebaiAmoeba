c
c
      subroutine empole1a1_3b_Polar_totfieldnpole_mod3(npole3b,pnum,
     &           uind,uinp,deptemp2,virtemp)
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
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 eptemp,deptemp2(3,npole3b),virtemp(3,3),deptemp(3,npole3b)
c      real*8 uind(3,npole)
c      real*8 uinp(3,npole)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3
      logical proceed,usei,usek
      character*6 mode

c      eptemp = 0.0d0
c      eptemp_userep =0.0d0
c
c   Zero out temporary gradient of polarization energy
c
c      do l1 = 1, npole3b
c        do j = 1, 3
c          deptemp(j,l1) = 0.0d0
c        end do
c      end do
c
c   Zero out temporary virial
c
c      do i=1,3
c         do j=1,3
c           virtemp(i,j)=0.0d0
c         end do
c      end do

c      call induce0a_3b_PolelecOnly_totfieldnpole3b(npole3b,pnum,uind,
c     &     uinp)
c
c     perform dynamic allocation of some local arrays
c


c      allocate (pscale(npole3b))
c      allocate (dscale(npole3b))
c      allocate (uscale(npole3b))

      allocate (pscale(npole))
      allocate (uscale(npole))
      allocate (dscale(npole))
c
c     set arrays needed to scale connected atom interactions
c
c      do i = 1, npole3b
      do i = 1,npole
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
c      do i = 1, npole-1
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
               scale3i = scale3 * uscale(k)
               scale5i = scale5 * uscale(k)
               scale7i = scale7 * uscale(k)
c               dsc3 = scale3 * dscale(l3)
c               dsc5 = scale5 * dscale(l3)
c               dsc7 = scale7 * dscale(l3)
               dsc3 = scale3 * dscale(k)
               dsc5 = scale5 * dscale(k)
               dsc7 = scale7 * dscale(k)

c               psc3 = scale3 * pscale(k)
c               psc5 = scale5 * pscale(k)
c               psc7 = scale7 * pscale(k)
c
c     construct necessary auxiliary vectors
c
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
               rxqir(1) = yr*qirz - zr*qiry
               rxqir(2) = zr*qirx - xr*qirz
               rxqir(3) = xr*qiry - yr*qirx
               rxqkr(1) = yr*qkrz - zr*qkry
               rxqkr(2) = zr*qkrx - xr*qkrz
               rxqkr(3) = xr*qkry - yr*qkrx
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
               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
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
               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
c               scip(1)=0.0d0
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
c               gli(1) = ck*sci(3) - ci*sci(4)
c               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
c               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c               gli(6) = sci(1)
c               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
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

               gfi(1) = 0.5d0 * rr5 * (
     &                                    (glip(1)+glip(6))*dsc3
     &                                   + scip(2)*scale3i)
     &                + 0.5d0 * rr7 * (
     &                                (glip(7)+glip(2))*dsc5
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0 * rr9 * (glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (scip(4)*dsc7)
               gfi(6) = -rr7 * (scip(3)*dsc7)


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


               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(1,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*dix
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkx
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qirx + gfi(6)*qkrx
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &           (- rr3*ck*(uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(2,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*(
     &              (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(3,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz


c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  +(glip(1)+glip(6))*dscale(l3))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  +(glip(2)+glip(7))*dscale(l3))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  +glip(3)*dscale(l3))


               temp3 = 0.5d0 * rr3 * (
     &                                  (glip(1)+glip(6))*dscale(k))
               temp5 = 0.5d0 * rr5 * (
     &                                  (glip(2)+glip(7))*dscale(k))
               temp7 = 0.5d0 * rr7 * (
     &                                  glip(3)*dscale(k))

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
               gti(2) = 0.5d0 * rr5 * (scip(4)*dsc5)
               gti(3) = 0.5d0 * rr5 * (scip(3)*dsc5)
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

               ttm2i(1) = -rr3*(dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(
     &           (ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(
     &           (ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(
     &           (ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxuip(1)*dsc3)*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(
     &           (uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxuip(2)*dsc3)*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(
     &           (uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxuip(3)*dsc3)*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(
     &           (uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)



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

c               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
c               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
c               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
c               call torque_3b_new(npole3b,pnum,deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)

               deptemp2(1,l1) = deptemp2(1,l1) + ftm2i(1)
               deptemp2(2,l1) = deptemp2(2,l1) + ftm2i(2)
               deptemp2(3,l1) = deptemp2(3,l1) + ftm2i(3)
c               call torque_3b_new_npole(deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
c               frcxi = 0.0d0
c               frcyi = 0.0d0
c               frczi = 0.0d0
               call torque_3b_new(npole3b,pnum,deptemp2,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

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
c               frcxk = 0.0d0
c               frcyk = 0.0d0
c               frczk = 0.0d0
               call torque_3b_new(npole3b,pnum,deptemp2,k,
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
         do j = 1, npole
            pscale(j) = 1.0d0
            dscale(j) = 1.0d0
            uscale(j) = 1.0d0
         end do

      end do
c
c

c
c
c
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
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
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
c               psc7 = scale7 * pscale(l3)

               scale3i = scale3
               scale5i = scale5
               scale7i = scale7
               dsc3 = scale3
               dsc5 = scale5
               dsc7 = scale7
c               psc3 = scale3
c               psc5 = scale5
c               psc7 = scale7
               if (use_polymer) then
                  if (r2 .le. polycut2) then
                     scale3i = scale3i * uscale(k)
                     scale5i = scale5i * uscale(k)
                     scale7i = scale7i * uscale(k)
                     dsc3 = dsc3 * dscale(kk)
                     dsc5 = dsc5 * dscale(kk)
                     dsc7 = dsc7 * dscale(kk)
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
               rxqir(1) = yr*qirz - zr*qiry
               rxqir(2) = zr*qirx - xr*qirz
               rxqir(3) = xr*qiry - yr*qirx
               rxqkr(1) = yr*qkrz - zr*qkry
               rxqkr(2) = zr*qkrx - xr*qkrz
               rxqkr(3) = xr*qkry - yr*qkrx
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
               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c

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
               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
c               scip(1)=0.0d0
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
c               gli(1) = ck*sci(3) - ci*sci(4)
c               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
c               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c               gli(6) = sci(1)
c               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
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
     &                                    (glip(1)+glip(6))*dsc3
     &                                   + scip(2)*scale3i)
     &                + 0.5d0 * rr7 * (
     &                                (glip(7)+glip(2))*dsc5
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0 * rr9 * (glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (scip(4)*dsc7)
               gfi(6) = -rr7 * (scip(3)*dsc7)

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
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(1,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*dix
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkx
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qirx + gfi(6)*qkrx
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &           (- rr3*ck*(uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(2,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*(
     &              (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(3,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz


c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  +(glip(1)+glip(6))*dscale(l3))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  +(glip(2)+glip(7))*dscale(l3))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  +glip(3)*dscale(l3))
               temp3 = 0.5d0 * rr3 * (
     &                                  (glip(1)+glip(6))*dscale(k))
               temp5 = 0.5d0 * rr5 * (
     &                                  (glip(2)+glip(7))*dscale(k))
               temp7 = 0.5d0 * rr7 * (
     &                                  glip(3)*dscale(k))

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
               gti(2) = 0.5d0 * rr5 * (scip(4)*dsc5)
               gti(3) = 0.5d0 * rr5 * (scip(3)*dsc5)
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

               ttm2i(1) = -rr3*(dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(
     &           (ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(
     &           (ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(
     &           (ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxuip(1)*dsc3)*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(
     &           (uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxuip(2)*dsc3)*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(
     &           (uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxuip(3)*dsc3)*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(
     &           (uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)



c
c     handle the case where scaling is used
c
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * ftm2i(j)
                  ttm2i(j) = f * ttm2i(j)
                  ttm3i(j) = f * ttm3i(j)
               end do

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

c               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
c               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
c               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
c               call torque_3b_new(npole3b,pnum,deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)

               deptemp2(1,l1) = deptemp2(1,l1) + ftm2i(1)
               deptemp2(2,l1) = deptemp2(2,l1) + ftm2i(2)
               deptemp2(3,l1) = deptemp2(3,l1) + ftm2i(3)
c               call torque_3b_new_npole(deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)
c               frcxi = 0.0d0
c               frcyi = 0.0d0
c               frczi = 0.0d0
               call torque_3b_new(npole3b,pnum,deptemp2,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

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
c               frcxk = 0.0d0
c               frcyk = 0.0d0
c               frczk = 0.0d0
               call torque_3b_new(npole3b,pnum,deptemp2,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)


c
c     increment gradient due to force and torque on first site
c

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
            pscale(j) = 1.0d0
            dscale(j) = 1.0d0
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end
