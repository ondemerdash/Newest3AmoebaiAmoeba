c      subroutine ereal1c_3b(aewald3b3b,ewaldcut3b,npole3b,pnum,uind,uinp,
c     &   eptemp,deptemp,virtemp)
      subroutine ereal1c2_3b_totfieldnpole_omp(npole3b,pnum,uind,uinp,
     &  eptemp,deptemp,virtemp)
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
      integer i,j,k,l1,l3,l2,l,ll,k1
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
      integer npole3b,pnum(*),kkk
      character*7 mode
      logical flag,doi,dorli
      real*8, allocatable :: dept1(:,:)
      real*8, allocatable :: dept2(:,:)
      real*8 virt(3,3),eptt

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
c      allocate (uscale(npole3b))
      allocate (dept1(3,n))
      allocate (dept2(3,n))

      do i = 1, n
         dept1(1,i) = 0.0d0
         dept1(2,i) = 0.0d0
         dept1(3,i) = 0.0d0
         dept2(1,i) = 0.0d0
         dept2(2,i) = 0.0d0
         dept2(3,i) = 0.0d0
      end do
      eptt = 0.0d0

      do i = 1, 3
         do j = 1, 3
            virt(j,i) = 0.0d0
         end do
      end do

c
c     set arrays needed to scale connected atom interactions
c
c      do i = 1, n
c         mscale(i) = 1.0d0
      do i = 1, npole3b
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
c         uscale(i) = 1.0d0
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


!$OMP PARALLEL default(shared) firstprivate(f)
!$OMP& private(i,j,k,ii,kkk,kk,e,ei,bfac,damp,expdamp,
!$OMP& pdi,pti,pgamma,scale3,scale5,scale7,temp3,temp5,temp7,
!$OMP& dsc3,dsc5,dsc7,psc3,psc5,psc7,usc3,usc5,alsq2,alsq2n,
!$OMP& exp2a,ralpha,gfd,gfdr,xr,yr,zr,xix,yix,zix,
!$OMP& xiy,yiy,ziy,xiz,yiz,ziz,xkx,ykx,zkx,xky,yky,zky,
!$OMP& xkz,ykz,zkz,r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
!$OMP& erl,erli,iax,iay,iaz,kax,kay,kaz,vxx,vyy,vzz,vyx,vzx,vzy,
!$OMP& frcxi,frcyi,frczi,frcxk,frcyk,frczk,ci,di,qi,ck,dk,qk,
!$OMP& fridmp,findmp,ftm2,ftm2i,ftm2r,ftm2ri,ttm2,ttm3,
!$OMP& ttm2i,ttm3i,ttm2r,ttm3r,ttm2ri,ttm3ri,fdir,dixdk,
!$OMP& dkxui,dixuk,dixukp,dkxuip,uixqkr,ukxqir,uixqkrp,ukxqirp,
!$OMP& qiuk,qkui,qiukp,qkuip,rxqiuk,rxqkui,rxqiukp,rxqkuip,
!$OMP& qidk,qkdi,qir,qkr,qiqkr,qkqir,qixqk,rxqir,dixr,dkxr,
!$OMP& dixqkr,dkxqir,rxqkr,qkrxqir,rxqikr,rxqkir,rxqidk,rxqkdi,
!$OMP& ddsc3,ddsc5,ddsc7,bn,sc,gl,sci,scip,gli,glip,gf,gfi,
!$OMP& gfr,gfri,gti,gtri,dorli,doi,k1,l1,l2)
!$OMP& firstprivate(pscale,dscale)
!$OMP DO reduction(+:eptt,dept1,dept2,virt)
!$OMP& schedule(guided)
      do i=1,npole
         doi=.false.
         do k1=1,npole3b
            if(pnum(k1).eq.i) then
              doi=.true.
              l1=k1
              goto 21
            end if
         end do
   21    continue  
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
         end do

         do j = 1, n13(ii)
            do kk=1,npole3b
               if(pnum(kk).eq.i13(j,ii)) then
                 pscale(kk) = p3scale
                 goto 32
               end if
            end do
   32   continue
         end do

         do j = 1, n14(ii)
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
         end do

         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = d1scale
c            uscale(ip11(j,ii)) = u1scale
            do kk=1,npole3b
               if(pnum(kk).eq.ip11(j,ii)) then
                 dscale(kk) = d1scale
c                 uscale(kk) = u1scale                          
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
c                 uscale(kk) = u2scale
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
c                 uscale(kk) = u3scale
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
c                 uscale(kk) = u4scale
                 goto 39
               end if
            end do
   39   continue             
         end do
c         do l2 = l1+1, npole3b
         do l2 =1,npole3b
            k=pnum(l2)
            kk = ipole(k)
          if(kk.ne.ii) then
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
c               usc3 = 1.0d0 - scale3*uscale(l2)
c               usc5 = 1.0d0 - scale5*uscale(l2)

c
c     construct necessary auxiliary vectors
c

ccc               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
ccc               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
ccc               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
c               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
c               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
c               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
               if(doi) then
                dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
                dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
                dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
                dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
                dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
                dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               else 
                dkxuip(1) = 0.0d0
                dkxuip(2) = 0.0d0
                dkxuip(3) = 0.0d0  
                dkxui(1) = 0.0d0 
                dkxui(2) = 0.0d0
                dkxui(3) = 0.0d0
               end if
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
ccc               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
ccc               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
ccc               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
ccc               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
ccc               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
ccc               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
ccc               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
ccc     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
ccc               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
ccc     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
ccc               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
ccc     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
ccc               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
ccc               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
ccc               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
ccc               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
ccc               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
ccc               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
ccc               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
ccc               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
ccc               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
ccc               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
ccc               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
ccc               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
ccc               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
ccc               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
ccc               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
     &                      + qi(7)*uind(3,l2)
               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
     &                      + qi(8)*uind(3,l2)
               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
     &                      + qi(9)*uind(3,l2)
c               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
c     &                      + qk(7)*uind(3,l1)
c               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
c     &                      + qk(8)*uind(3,l1)
c               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
c     &                      + qk(9)*uind(3,l1)          
               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
     &                       + qi(7)*uinp(3,l2)
               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
     &                       + qi(8)*uinp(3,l2)
               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
     &                       + qi(9)*uinp(3,l2)
               if(doi) then
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)          
               else
               qkuip(1) = 0.0d0
               qkuip(2) = 0.0d0
               qkuip(3) = 0.0d0
               qkui(1) = 0.0d0 
               qkui(2) = 0.0d0
               qkui(3) = 0.0d0
               end if
ccc               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
ccc               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
ccc               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
ccc               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
ccc               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
ccc               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
c               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
c               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
c               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
               if(doi) then
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               else
               uixqkrp(1) = 0.0d0
               uixqkrp(2) = 0.0d0
               uixqkrp(3) = 0.0d0
               uixqkr(1) = 0.0d0 
               uixqkr(2) = 0.0d0
               uixqkr(3) = 0.0d0
               end if
               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
ccc               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
ccc               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
ccc               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
ccc               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
ccc               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
ccc               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
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
c               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
c     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
c     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
c               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
c     &                     + uind(3,l1)*uind(3,l2)
c               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
c               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
c               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
c     &                     + qir(3)*uind(3,l2)
c               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
c     &                     + qkr(3)*uind(3,l1)
              if(doi) then
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
c               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
c     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
c     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)

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
              else
               scip(1) = 
     &                      di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(3) = 0.0d0
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = 0.0d0

               sci(1) = 
     &                      di(1)*uind(1,l2)
     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
c               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
c     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = 0.0d0 
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &                     + qir(3)*uind(3,l2)
               sci(8) = 0.0d0 
              end if
c
c     calculate the gl functions for permanent components
c

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
               if(doi) then
               eptt = eptt + 0.5d0*ei
               else
               eptt = eptt + ei
               end if
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

c
c     intermediate variables for induced force terms
c
c               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
c     &                     + 0.5d0*bn(2)*scip(2)
c     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
c     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
c     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
c               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
c               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
c               gfi(4) = 2.0d0 * bn(2)
c               gfi(5) = bn(3) * (sci(4)+scip(4))
c               gfi(6) = -bn(3) * (sci(3)+scip(3))
c               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
c     &                            + (glip(1)+glip(6))*dsc3
c     &                            + scip(2)*usc3)
c     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
c     &                            + (glip(7)+glip(2))*dsc5
c     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
c     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
c               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
c               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
c               gfri(4) = 2.0d0 * rr5
c               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
c               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)


               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     !+ 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     !- 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                             )
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                             )
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)


c
c     get the induced force with screening
c
c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
c     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
c     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
c     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
c     &            + (sci(4)+scip(4))*bn(2)*di(1)
c     &            + (sci(3)+scip(3))*bn(2)*dk(1)
c     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
c     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
c     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
c     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
c     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
c     &            + (sci(4)+scip(4))*bn(2)*di(2)
c     &            + (sci(3)+scip(3))*bn(2)*dk(2)
c     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
c     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
c               ftm2i(3) = gfi(1)*zr + 0.5d0*
c     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
c     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
c     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
c     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
c     &            + (sci(4)+scip(4))*bn(2)*di(3)
c     &            + (sci(3)+scip(3))*bn(2)*dk(3)
c     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
c     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
               if(doi) then
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            !+ bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            !+ bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            !+ bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            !+ bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            !+ bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            !+ bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)

               else 
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (!gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            !+ bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            !+ bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (!gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            !+ bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            !+ bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (!gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            !+ bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            !+ bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
               end if               
c
c     get the induced force without screening
c

c               ftm2ri(1) = gfri(1)*xr + 0.5d0*
c     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
c     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
c     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
c     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
c     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
c     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
c     &            + (qkuip(1)-qiukp(1))*dsc5)
c     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
c               ftm2ri(2) = gfri(1)*yr + 0.5d0*
c     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
c     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
c     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
c     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
c     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
c     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
c     &            + (qkuip(2)-qiukp(2))*dsc5)
c     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
c               ftm2ri(3) = gfri(1)*zr + 0.5d0*
c     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
c     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
c     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
c     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
c     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
c     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
c     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
c     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
c     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
c     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
c     &            + (qkuip(3)-qiukp(3))*dsc5)
c     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)

               if(doi) then
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            !+ rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            !+ sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
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
     &            !+ rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            !+ sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
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
     &           ! + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &           ! + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)

               else
               ftm2ri(1) = gfri(1)*xr !+ 0.5d0*
     &           !(- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &           ! + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &           ! - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            !+ rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            !+ sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr !+ 0.5d0*
     &           !(- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &           ! + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &           ! - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
     &            !+ rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            !+ sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr ! + 0.5d0*
     &           !(- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &           ! + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &           ! - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
     &           ! + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &           ! + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)

               end if
c
c     account for partially excluded induced interactions
c
C   LEFT OFF HERE!!!
ccc               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
ccc     &                                  +(glip(1)+glip(6))*dscale(kk))
ccc               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
ccc     &                                  +(glip(2)+glip(7))*dscale(kk))
ccc               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
ccc     &                                  +glip(3)*dscale(kk))

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
ccc               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
ccc               temp5 = -0.5d0 * rr5 * uscale(kk)
ccc     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

c               temp3 = 0.5d0 * rr3 * uscale(l2) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(l2)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))

c               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
c               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
c               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
               findmp(1) = 0.0d0
               findmp(2) = 0.0d0
               findmp(3) = 0.0d0
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
              if(doi) then
               dept1(1,i) = dept1(1,i) + 0.5d0*ftm2i(1)
               dept1(2,i) = dept1(2,i) + 0.5d0*ftm2i(2)
               dept1(3,i) = dept1(3,i) + 0.5d0*ftm2i(3)
               call torque_3b_new_npole(dept1,i,
     &          ttm2,0.5d0*ttm2i,frcxi,frcyi,frczi)

               dept2(1,k) = dept2(1,k) - 0.5d0*ftm2i(1)
               dept2(2,k) = dept2(2,k) - 0.5d0*ftm2i(2)
               dept2(3,k) = dept2(3,k) - 0.5d0*ftm2i(3)
               call torque_3b_new_npole(dept2,k,
     &            ttm3,0.5d0*ttm3i,frcxk,frcyk,frczk)
              else
               dept1(1,i) = dept1(1,i) + ftm2i(1)
               dept1(2,i) = dept1(2,i) + ftm2i(2)
               dept1(3,i) = dept1(3,i) + ftm2i(3)
               call torque_3b_new_npole(dept1,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)
               dept2(1,k) = dept2(1,k) - ftm2i(1)
               dept2(2,k) = dept2(2,k) - ftm2i(2)
               dept2(3,k) = dept2(3,k) - ftm2i(3)
               call torque_3b_new_npole(dept2,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)
              end if

              ! deptemp(1,l1) =deptemp(1,l1) + ftm2i(1)
              ! deptemp(2,l1) =deptemp(2,l1) + ftm2i(2)
              ! deptemp(3,l1) =deptemp(3,l1) + ftm2i(3)

              ! call torque_3b_new(npole3b,pnum,deptemp,i,
     &        !  ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

               ! deptemp(1,l2) = deptemp(1,l2) - ftm2i(1)
               ! deptemp(2,l2) = deptemp(2,l2) - ftm2i(2)
               ! deptemp(3,l2) = deptemp(3,l2) - ftm2i(3)

               ! call torque_3b_new(npole3b,pnum,deptemp,k,
     &         !   ttm3,ttm3i,frcxk,frcyk,frczk)

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

               virt(1,1) = virt(1,1) + vxx
               virt(2,1) = virt(2,1) + vyx
               virt(3,1) = virt(3,1) + vzx
               virt(1,2) = virt(1,2) + vyx
               virt(2,2) = virt(2,2) + vyy
               virt(3,2) = virt(3,2) + vzy
               virt(1,3) = virt(1,3) + vzx
               virt(2,3) = virt(2,3) + vzy
               virt(3,3) = virt(3,3) + vzz

            end if
          end if
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
c            uscale(j)=1.0d0
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL

      eptemp = eptemp +eptt
      do i = 1, n
         do j = 1, 3
            deptemp(j,i) = deptemp(j,i) + dept1(j,i) + dept2(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            virtemp(j,i) = virtemp(j,i) + virt(j,i)
         end do
      end do

c
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (dept1)
      deallocate (dept2)
c      deallocate (uscale)
c      print*,'ep after ereal1c_3b=',eptemp,moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After ereal1c_3b l1 Pnum(l1)",l1,pnum(l1)
c      end do

      return
      end

