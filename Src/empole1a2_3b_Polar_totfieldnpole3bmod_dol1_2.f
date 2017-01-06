c
c

c      subroutine empole1a2_3b_Polar_totfieldnpole3bmod(npole3b,pnum,
c     &   uinp,deptemp,virtemp)
      subroutine empole1a2_3b_Polar_totfieldnpole3bmod_dol1_2(npole3b,
     &   pnum,uinp,deptemp2,virtemp)
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
      real*8 eptemp,deptemp2(3,npole),virtemp(3,3),deptemp(3,npole3b)
c      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3,k1
      logical proceed,usei,usek
      character*6 mode
      logical doi,dok
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
      allocate (dscale(npole3b))
c      allocate (uscale(npole3b))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, npole3b
c         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
c         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c
c     set scale factors for permanent multipole and induced terms
c

c      do l1 = 1, npole3b-1
c         i=pnum(l1)

c      do i=1,npole-1
      do i=1,npole
         doi=.false.
         dok=.true.
         do k1=1,npole3b
            if(pnum(k1).eq.i) then
              doi=.true.
              l1=k1
              goto 21   
            end if
         end do
   21    continue                  
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
c            do kk=1,npole3b
c               if (pnum(kk).eq.i12(j,ii)) then
c                  pscale(kk)=p2scale
c                  goto 31
c               end if
c            end do 
c   31             continue
c         end do
c         do j = 1, n13(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i13(j,ii)) then
c                  pscale(kk)=p3scale
c                  goto 32
c               end if
c            end do
c   32             continue            
c         end do
c         do j = 1, n14(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i14(j,ii)) then
c                  pscale(kk)=p4scale
c                  goto 33
c               end if
c            end do
c   33             continue
c            do k = 1, np11(ii)
c               do kk=1,npole3b             
c                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
c     &              (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(kk) = p4scale * p41scale
c                   goto 34
c                 end if
c               end do  
c   34             continue
c            end do
c         end do
c         do j = 1, n15(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i15(j,ii)) then
c                  pscale(kk)=p5scale
c                  goto 35
c               end if
c            end do
c   35             continue
c         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
                  dscale(kk)=d1scale
c                  uscale(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale(kk)=d2scale
c                  uscale(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale(kk)=d3scale
c                  uscale(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale(kk)=d4scale
c                  uscale(kk)=u4scale
                  goto 39
               end if
            end do
   39             continue
         end do
c         do l3 = l1+1, npole3b
         do l3 = 1, npole3b          
          k=pnum(l3)
          kk = ipole(k)
c          if(kk.gt.ii) then  
          if(kk.ne.ii) then
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
               dsc3 = scale3 * dscale(l3)
               dsc5 = scale5 * dscale(l3)
               dsc7 = scale7 * dscale(l3)
c               psc3 = scale3 * pscale(l3)
c               psc5 = scale5 * pscale(l3)
c               psc7 = scale7 * pscale(l3)
c
c     construct necessary auxiliary vectors
c

c               dixuk(1) = diy*uind(3,l3) - diz*uind(2,l3)
c               dixuk(2) = diz*uind(1,l3) - dix*uind(3,l3)
c               dixuk(3) = dix*uind(2,l3) - diy*uind(1,l3)
c               dkxui(1) = dky*uind(3,l1) - dkz*uind(2,l1)
c               dkxui(2) = dkz*uind(1,l1) - dkx*uind(3,l1)
c               dkxui(3) = dkx*uind(2,l1) - dky*uind(1,l1)
               dixukp(1) = diy*uinp(3,l3) - diz*uinp(2,l3)
               dixukp(2) = diz*uinp(1,l3) - dix*uinp(3,l3)
               dixukp(3) = dix*uinp(2,l3) - diy*uinp(1,l3)
              ! if(doi.and.(k.gt.i)) then
               if(doi) then
               dkxuip(1) = dky*uinp(3,l1) - dkz*uinp(2,l1)
               dkxuip(2) = dkz*uinp(1,l1) - dkx*uinp(3,l1)
               dkxuip(3) = dkx*uinp(2,l1) - dky*uinp(1,l1)
               else 
               dkxuip(1) = 0.0d0
               dkxuip(2) = 0.0d0
               dkxuip(3) = 0.0d0
               end if 
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
               rxqir(1) = yr*qirz - zr*qiry
               rxqir(2) = zr*qirx - xr*qirz
               rxqir(3) = xr*qiry - yr*qirx
               rxqkr(1) = yr*qkrz - zr*qkry
               rxqkr(2) = zr*qkrx - xr*qkrz
               rxqkr(3) = xr*qkry - yr*qkrx
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
               qiukp(1) = qixx*uinp(1,l3) + qiyx*uinp(2,l3)
     &                       + qizx*uinp(3,l3)
               qiukp(2) = qixy*uinp(1,l3) + qiyy*uinp(2,l3)
     &                       + qizy*uinp(3,l3)
               qiukp(3) = qixz*uinp(1,l3) + qiyz*uinp(2,l3)
     &                       + qizz*uinp(3,l3)
               ! if(doi.and.(k.gt.i)) then
               if(doi) then
               qkuip(1) = qkxx*uinp(1,l1) + qkyx*uinp(2,l1)
     &                       + qkzx*uinp(3,l1)
               qkuip(2) = qkxy*uinp(1,l1) + qkyy*uinp(2,l1)
     &                       + qkzy*uinp(3,l1)
               qkuip(3) = qkxz*uinp(1,l1) + qkyz*uinp(2,l1)
     &                       + qkzz*uinp(3,l1)
               else
               qkuip(1) = 0.0d0
               qkuip(2) = 0.0d0
               qkuip(3) = 0.0d0
               end if
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
              ! if(doi.and.(k.gt.i)) then
               if(doi) then
               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
               else
               uixqkrp(1) = 0.0d0 
               uixqkrp(2) = 0.0d0
               uixqkrp(3) = 0.0d0
               end if
               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
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
c               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
c     &                     + uind(3,l1)*uind(3,l3)
c               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
c               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
c               sci(7) = qirx*uind(1,l3) + qiry*uind(2,l3)
c     &                     + qirz*uind(3,l3)
c               sci(8) = qkrx*uind(1,l1) + qkry*uind(2,l1)
c     &                     + qkrz*uind(3,l1)

             ! if(doi.and.(k.gt.i)) then               
              if(doi) then
               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
c               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
c     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
c     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
     &                      + qirz*uinp(3,l3)
               scip(8) = qkrx*uinp(1,l1) + qkry*uinp(2,l1)
     &                      + qkrz*uinp(3,l1)
c              else if(.not.doi) then
              else
               scip(1) = 
     &                      dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
c               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
c     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
c     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = 0.0d0
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
     &                      + qirz*uinp(3,l3)
               scip(8) = 0.0d0 
              end if
              
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

               gfi(1) = 0.5d0 * rr5 * (
     &                                    (glip(1)+glip(6))*dsc3
     &                                   )
     &                + 0.5d0 * rr7 * (
     &                                (glip(7)+glip(2))*dsc5
     &                                 )
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

             ! if(doi.and.(k.gt.i)) then
              if(doi) then  
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(1,l3)*dsc7))*0.5d0
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
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(3,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*(
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz
c              else if(.not.doi) then
              else
               ftm2i(1) = gfi(1)*xr 
     &            + (rr3*ci*(uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(1,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*dix
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkx
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qirx + gfi(6)*qkrx
               ftm2i(2) = gfi(1)*yr 
     &            + (rr3*ci*(uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(2,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr 
     &            + (rr3*ci*(uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(3,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*(
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz
              end if

c
c     account for partially excluded induced interactions
c     LEFT OFF HERE!!!!
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  +(glip(1)+glip(6))*dscale(l3))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  +(glip(2)+glip(7))*dscale(l3))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  +glip(3)*dscale(l3))

               temp3 = 0.5d0 * rr3 * (
     &                                  (glip(1)+glip(6))*dscale(l3))
               temp5 = 0.5d0 * rr5 * (
     &                                 (glip(2)+glip(7))*dscale(l3))
               temp7 = 0.5d0 * rr7 * (
     &                                  glip(3)*dscale(l3))



               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c

c               temp3 = 0.5d0 * rr3 * uscale(l3) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(l3)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

c               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
c               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
c               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)

c
c     modify induced force for partially excluded interactions
c
c               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
c               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
c               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

               ftm2i(1) = ftm2i(1) - fridmp(1) 
               ftm2i(2) = ftm2i(2) - fridmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3)

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

c              if(doi) then
               ttm2i(1) = -rr3*(dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(
     &           (ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(
     &           (ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(
     &           (ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
c              end if
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


              if(doi) then
               deptemp2(1,i) = deptemp2(1,i) + 0.5d0*ftm2i(1)
               deptemp2(2,i) = deptemp2(2,i) + 0.5d0*ftm2i(2)
               deptemp2(3,i) = deptemp2(3,i) + 0.5d0*ftm2i(3)
               call torque_3b_new_npole(deptemp2,i,
     &          ttm2,0.5d0*ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
               deptemp2(1,k) = deptemp2(1,k) - 0.5d0*ftm2i(1)
               deptemp2(2,k) = deptemp2(2,k) - 0.5d0*ftm2i(2)
               deptemp2(3,k) = deptemp2(3,k) - 0.5d0*ftm2i(3)
               call torque_3b_new_npole(deptemp2,k,
     &            ttm3,0.5d0*ttm3i,frcxk,frcyk,frczk)
              else
               deptemp2(1,i) = deptemp2(1,i) + ftm2i(1)
               deptemp2(2,i) = deptemp2(2,i) + ftm2i(2)
               deptemp2(3,i) = deptemp2(3,i) + ftm2i(3)
               call torque_3b_new_npole(deptemp2,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
               deptemp2(1,k) = deptemp2(1,k) - ftm2i(1)
               deptemp2(2,k) = deptemp2(2,k) - ftm2i(2)
               deptemp2(3,k) = deptemp2(3,k) - ftm2i(3)
               call torque_3b_new_npole(deptemp2,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)
              end if
              
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

             if(doi) then
               vxx = -xr*(0.5d0*ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(0.5d0*ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(0.5d0*ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(0.5d0*ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(0.5d0*ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(0.5d0*ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
             else
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

            !end if
   10       continue
          end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, npole3b
c            pscale(j) = 1.0d0
            dscale(j) = 1.0d0
c            uscale(j) = 1.0d0
         end do

      end do
c
c     perform deallocation of some local arrays
c


c      do l1 = 1, npole3b-1
c         i=pnum(l1)
      if (use_replica) then

      do i=1,npole
         doi=.false.
         do k1=1,npole3b
            if(pnum(k1).eq.i) then
              doi=.true.
              l1=k1
              goto 11   
            end if
         end do
   11    continue                  
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
c            do kk=1,npole3b
c               if (pnum(kk).eq.i12(j,ii)) then
c                  pscale(kk)=p2scale
c                  goto 31
c               end if
c            end do 
c   31             continue
c         end do
c         do j = 1, n13(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i13(j,ii)) then
c                  pscale(kk)=p3scale
c                  goto 32
c               end if
c            end do
c   32             continue            
c         end do
c         do j = 1, n14(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i14(j,ii)) then
c                  pscale(kk)=p4scale
c                  goto 33
c               end if
c            end do
c   33             continue
c            do k = 1, np11(ii)
c               do kk=1,npole3b             
c                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
c     &              (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(kk) = p4scale * p41scale
c                   goto 34
c                 end if
c               end do  
c   34             continue
c            end do
c         end do
c         do j = 1, n15(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i15(j,ii)) then
c                  pscale(kk)=p5scale
c                  goto 35
c               end if
c            end do
c   35             continue
c         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
                  dscale(kk)=d1scale
c                  uscale(kk)=u1scale
                  goto 46
               end if
            end do
   46             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale(kk)=d2scale
c                  uscale(kk)=u2scale
                  goto 47
               end if
            end do
   47             continue
         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale(kk)=d3scale
c                  uscale(kk)=u3scale
                  goto 48
               end if
            end do
   48             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale(kk)=d4scale
c                  uscale(kk)=u4scale
                  goto 49
               end if
            end do
   49             continue
         end do
c         do l3 = l1+1, npole3b
         do l3 = 1, npole3b          
          k=pnum(l3)
          kk = ipole(k)
c          if(kk.ge.ii) then  
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
            call imager (xr,yr,zr,jcell)
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

c               scale3i = scale3
c               scale5i = scale5
c               scale7i = scale7
               dsc3 = scale3
               dsc5 = scale5
               dsc7 = scale7
c               psc3 = scale3
c               psc5 = scale5
c               psc7 = scale7
               if (use_polymer) then
                  if (r2 .le. polycut2) then
c                     scale3i = scale3i * uscale(kk)
c                     scale5i = scale5i * uscale(kk)
c                     scale7i = scale7i * uscale(kk)
c                     dsc3 = dsc3 * dscale(kk)
c                     dsc5 = dsc5 * dscale(kk)
c                     dsc7 = dsc7 * dscale(kk)
                     dsc3 = dsc3 * dscale(l3)
                     dsc5 = dsc5 * dscale(l3)
                     dsc7 = dsc7 * dscale(l3)
c                     psc3 = psc3 * pscale(kk)
c                     psc5 = psc5 * pscale(kk)
c                     psc7 = psc7 * pscale(kk)
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
               dixukp(1) = diy*uinp(3,l3) - diz*uinp(2,l3)
               dixukp(2) = diz*uinp(1,l3) - dix*uinp(3,l3)
               dixukp(3) = dix*uinp(2,l3) - diy*uinp(1,l3)
               if(doi) then
               dkxuip(1) = dky*uinp(3,l1) - dkz*uinp(2,l1)
               dkxuip(2) = dkz*uinp(1,l1) - dkx*uinp(3,l1)
               dkxuip(3) = dkx*uinp(2,l1) - dky*uinp(1,l1)
               else 
               dkxuip(1) = 0.0d0
               dkxuip(2) = 0.0d0
               dkxuip(3) = 0.0d0
               end if 
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
               rxqir(1) = yr*qirz - zr*qiry
               rxqir(2) = zr*qirx - xr*qirz
               rxqir(3) = xr*qiry - yr*qirx
               rxqkr(1) = yr*qkrz - zr*qkry
               rxqkr(2) = zr*qkrx - xr*qkrz
               rxqkr(3) = xr*qkry - yr*qkrx
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
               qiukp(1) = qixx*uinp(1,l3) + qiyx*uinp(2,l3)
     &                       + qizx*uinp(3,l3)
               qiukp(2) = qixy*uinp(1,l3) + qiyy*uinp(2,l3)
     &                       + qizy*uinp(3,l3)
               qiukp(3) = qixz*uinp(1,l3) + qiyz*uinp(2,l3)
     &                       + qizz*uinp(3,l3)
               if(doi) then
               qkuip(1) = qkxx*uinp(1,l1) + qkyx*uinp(2,l1)
     &                       + qkzx*uinp(3,l1)
               qkuip(2) = qkxy*uinp(1,l1) + qkyy*uinp(2,l1)
     &                       + qkzy*uinp(3,l1)
               qkuip(3) = qkxz*uinp(1,l1) + qkyz*uinp(2,l1)
     &                       + qkzz*uinp(3,l1)
               else
               qkuip(1) = 0.0d0
               qkuip(2) = 0.0d0
               qkuip(3) = 0.0d0
               end if
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
               if(doi) then
               uixqkrp(1) = uinp(2,l1)*qkrz - uinp(3,l1)*qkry
               uixqkrp(2) = uinp(3,l1)*qkrx - uinp(1,l1)*qkrz
               uixqkrp(3) = uinp(1,l1)*qkry - uinp(2,l1)*qkrx
               else
               uixqkrp(1) = 0.0d0 
               uixqkrp(2) = 0.0d0
               uixqkrp(3) = 0.0d0
               end if
               ukxqirp(1) = uinp(2,l3)*qirz - uinp(3,l3)*qiry
               ukxqirp(2) = uinp(3,l3)*qirx - uinp(1,l3)*qirz
               ukxqirp(3) = uinp(1,l3)*qiry - uinp(2,l3)*qirx
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
c               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
c     &                     + uind(3,l1)*uind(3,l3)
c               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
c               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
c               sci(7) = qirx*uind(1,l3) + qiry*uind(2,l3)
c     &                     + qirz*uind(3,l3)
c               sci(8) = qkrx*uind(1,l1) + qkry*uind(2,l1)
c     &                     + qkrz*uind(3,l1)

              if(doi) then               
               scip(1) = uinp(1,l1)*dkx + uinp(2,l1)*dky
     &                      + uinp(3,l1)*dkz + dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
c               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
c     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
c     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
     &                      + qirz*uinp(3,l3)
               scip(8) = qkrx*uinp(1,l1) + qkry*uinp(2,l1)
     &                      + qkrz*uinp(3,l1)
              else
               scip(1) = 
     &                      dix*uinp(1,l3)
     &                      + diy*uinp(2,l3) + diz*uinp(3,l3)
c               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
c     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
c     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = 0.0d0
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qirx*uinp(1,l3) + qiry*uinp(2,l3)
     &                      + qirz*uinp(3,l3)
               scip(8) = 0.0d0 
              end if
              
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

               gfi(1) = 0.5d0 * rr5 * (
     &                                    (glip(1)+glip(6))*dsc3
     &                                   )
     &                + 0.5d0 * rr7 * (
     &                                (glip(7)+glip(2))*dsc5
     &                                 )
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

              if(doi) then
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(1,l3)*dsc7))*0.5d0
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
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(3,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*(
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz
              else
               ftm2i(1) = gfi(1)*xr 
     &            + (rr3*ci*(uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(1,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*dix
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkx
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qirx + gfi(6)*qkrx
               ftm2i(2) = gfi(1)*yr 
     &            + (rr3*ci*(uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(2,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diy
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dky
     &            + 0.5d0*gfi(4)*(
     &             (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qiry + gfi(6)*qkry
               ftm2i(3) = gfi(1)*zr 
     &            + (rr3*ci*(uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uinp(3,l3)*dsc7))*0.5d0
     &            + 0.5d0*(scip(4)*dsc5)*rr5*diz
     &            + 0.5d0*(scip(3)*dsc5)*rr5*dkz
     &            + 0.5d0*gfi(4)*(
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qirz + gfi(6)*qkrz
              end if

c
c     account for partially excluded induced interactions
c     LEFT OFF HERE!!!!
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(l3)
c     &                                  +(glip(1)+glip(6))*dscale(l3))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(l3)
c     &                                  +(glip(2)+glip(7))*dscale(l3))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(l3)
c     &                                  +glip(3)*dscale(l3))

               temp3 = 0.5d0 * rr3 * (
     &                                  (glip(1)+glip(6))*dscale(l3))
               temp5 = 0.5d0 * rr5 * (
     &                                 (glip(2)+glip(7))*dscale(l3))
               temp7 = 0.5d0 * rr7 * (
     &                                  glip(3)*dscale(l3))



               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c

c               temp3 = 0.5d0 * rr3 * uscale(l3) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(l3)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

c               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
c               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
c               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)

c
c     modify induced force for partially excluded interactions
c
c               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
c               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
c               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

               ftm2i(1) = ftm2i(1) - fridmp(1) 
               ftm2i(2) = ftm2i(2) - fridmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3)

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
c              if(doi) then    
               ttm2i(1) = -rr3*(dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(
     &           (ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(
     &           (ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(
     &           (ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
c              end if
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
              if(doi) then
               deptemp2(1,i) = deptemp2(1,i) + 0.5*ftm2i(1)
               deptemp2(2,i) = deptemp2(2,i) + 0.5*ftm2i(2)
               deptemp2(3,i) = deptemp2(3,i) + 0.5*ftm2i(3)
               call torque_3b_new_npole(deptemp2,i,
     &          ttm2,0.5*ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
               deptemp2(1,k) = deptemp2(1,k) - 0.5*ftm2i(1)
               deptemp2(2,k) = deptemp2(2,k) - 0.5*ftm2i(2)
               deptemp2(3,k) = deptemp2(3,k) - 0.5*ftm2i(3)
               call torque_3b_new_npole(deptemp2,k,
     &            ttm3,0.5*ttm3i,frcxk,frcyk,frczk)
              else
               deptemp2(1,i) = deptemp2(1,i) + ftm2i(1)
               deptemp2(2,i) = deptemp2(2,i) + ftm2i(2)
               deptemp2(3,i) = deptemp2(3,i) + ftm2i(3)
               call torque_3b_new_npole(deptemp2,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
               deptemp2(1,k) = deptemp2(1,k) - ftm2i(1)
               deptemp2(2,k) = deptemp2(2,k) - ftm2i(2)
               deptemp2(3,k) = deptemp2(3,k) - ftm2i(3)
               call torque_3b_new_npole(deptemp2,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)
              end if

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
c          end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, npole3b
c            pscale(j) = 1.0d0
            dscale(j) = 1.0d0
c            uscale(j) = 1.0d0
         end do

      end do
      end if
c
c     perform deallocation of some local arrays
c

c
c
c
c      print*,"eptemp",eptemp
c      deallocate (pscale)
      deallocate (dscale)
c      deallocate (uscale)
      return
      end


