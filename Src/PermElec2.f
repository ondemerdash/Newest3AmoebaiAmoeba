c
c      subroutine ereal1d_3b_Perm2(emreal_tmp,demreal_tmp,viremreal_tmp)
      subroutine ereal1d_3b_Perm2
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polar
      use polpot
      use virial
      use couple
      use mplpot
      use neigh
      use limits
      use potent
      implicit none
      include 'mpif.h'
      integer i,j,k,start,ending
      integer ii,kk,kkk
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
c      real*8 damp,expdamp
c      real*8 emreal_tmp,demreal_tmp(3,*),viremreal_tmp(3,3)
c      real*8 pdi,pti,pgamma
c      real*8 scale3,scale5
c      real*8 scale7
c      real*8 temp3,temp5,temp7
c      real*8 dsc3,dsc5,dsc7
c      real*8 psc3,psc5,psc7
c      real*8 usc3,usc5
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
c      real*8 emtt,eptt
      real*8 emtt
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
c      real*8 fridmp(3),findmp(3)
c      real*8 ftm2(3),ftm2i(3)
c      real*8 ftm2r(3),ftm2ri(3)
      real*8 ftm2(3)
      real*8 ftm2r(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
c      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
c      real*8 dkxui(3),dixuk(3)
c      real*8 dixukp(3),dkxuip(3)
c      real*8 uixqkr(3),ukxqir(3)
c      real*8 uixqkrp(3),ukxqirp(3)
c      real*8 qiuk(3),qkui(3)
c      real*8 qiukp(3),qkuip(3)
c      real*8 rxqiuk(3),rxqkui(3)
c      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
c      real*8 ddsc3(3),ddsc5(3)
c      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
c      real*8 sci(8),scip(8)
c      real*8 gli(7),glip(7)
c      real*8 gf(7),gfi(6)
      real*8 gf(7)
c      real*8 gfr(7),gfri(6)
      real*8 gfr(7)
c      real*8 gti(6),gtri(6)
      real*8 viri(3,3)
      real*8, allocatable :: mscale(:)
c      real*8, allocatable :: pscale(:)
c      real*8, allocatable :: dscale(:)
c      real*8, allocatable :: uscale(:)
      real*8, allocatable :: demi(:,:)
      real*8, allocatable :: demk(:,:)
c      real*8, allocatable :: depi(:,:)
c      real*8, allocatable :: depk(:,:)
c      integer atomind
c      real*8 demt(3,npole),viremt(3,3),emt,aewald
      logical dorl,dorli
      character*6 mode
      real*8 off,off2,cut,cut2
      real*8 c0,c1,c2,c3,c4,c5
      real*8 f0,f1,f2,f3,f4,f5,f6,f7
      real*8 denom,term
      real*8 off3,off4,off5
      real*8 off6,off7
      real*8 cut3,cut4,cut5
      real*8 cut6,cut7
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
c      eintra = 0.0d0
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c      allocate (pscale(n))
c      allocate (dscale(n))
c      allocate (uscale(n))
      allocate (demi(3,n))
      allocate (demk(3,n))
c      allocate (depi(3,n))
c      allocate (depk(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
c         pscale(i) = 1.0d0
c         dscale(i) = 1.0d0
c         uscale(i) = 1.0d0
      end do

c      emt=0.0d0
c      do i=1,npole
c         do j=1,3
c           demt(j,i)=0.0d0
c         end do
c      end do
c      do i=1,3
c         do j=1,3
c           viremt(i,j)=0.0d0
c         end do
c      end do

c       print*,"In ereal_perm xcell ycell zcell",xcell,ycell,zcell
c       print *,"In ereal_perm xcell2 ycell2 zcell2",xcell2,ycell2,zcell2
c       print*,"In ereal_perm x(3) y(3) z(3)",x(3),y(3),z(3)

c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      mode = 'EWALD'
c      call switch (mode)

         off = ewaldcut
         cut = ewaldcut

      c0 = 0.0d0
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      c4 = 0.0d0
      c5 = 0.0d0
      f0 = 0.0d0
      f1 = 0.0d0
      f2 = 0.0d0
      f3 = 0.0d0
      f4 = 0.0d0
      f5 = 0.0d0
      f6 = 0.0d0
      f7 = 0.0d0
c
c     store the powers of the switching window cutoffs
c
      off2 = off * off
      off3 = off2 * off
      off4 = off2 * off2
      off5 = off2 * off3
      off6 = off3 * off3
      off7 = off3 * off4
      cut2 = cut * cut
      cut3 = cut2 * cut
      cut4 = cut2 * cut2
      cut5 = cut2 * cut3
      cut6 = cut3 * cut3
      cut6 = cut3 * cut3
      cut7 = cut3 * cut4
      if (cut .lt. off) then
         denom = (off-cut)**5
         c0 = off*off2 * (off2-5.0d0*off*cut+10.0d0*cut2) / denom
         c1 = -30.0d0 * off2*cut2 / denom
         c2 = 30.0d0 * (off2*cut+off*cut2) / denom
         c3 = -10.0d0 * (off2+4.0d0*off*cut+cut2) / denom
         c4 = 15.0d0 * (off+cut) / denom
         c5 = -6.0d0 / denom
      end if


c
c     initialize local variables for OpenMP calculation
c
      emtt = 0.0d0
c      eptt = 0.0d0
      do i = 1, n
         do j = 1, 3
            demi(j,i) = 0.0d0
            demk(j,i) = 0.0d0
c            depi(j,i) = 0.0d0
c            depk(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            viri(j,i) = 0.0d0
         end do
      end do

c
c     set OpenMP directives for the major loop structure
c
cc      do i=start,ending

!$OMP PARALLEL default(shared) firstprivate(f) 
!$OMP& private(i,j,k,ii,kk,kkk,e,bfac,
!$OMP& alsq2,alsq2n,
!$OMP& exp2a,ralpha,gfd,gfdr,xr,yr,zr,xix,yix,zix,
!$OMP& xiy,yiy,ziy,xiz,yiz,ziz,xkx,ykx,zkx,xky,yky,zky,
!$OMP& xkz,ykz,zkz,r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
!$OMP& erl,erli,iax,iay,iaz,kax,kay,kaz,vxx,vyy,vzz,vyx,vzx,vzy,
!$OMP& frcxi,frcyi,frczi,frcxk,frcyk,frczk,ci,di,qi,ck,dk,qk,
!$OMP& ftm2,ftm2r,ttm2,ttm3,
!$OMP& ttm2i,ttm3i,ttm2r,ttm3r,dixdk,
!$OMP& qidk,qkdi,qir,qkr,qiqkr,qkqir,qixqk,rxqir,dixr,dkxr,
!$OMP& dixqkr,dkxqir,rxqkr,qkrxqir,rxqikr,rxqkir,rxqidk,rxqkdi,
!$OMP& bn,sc,gl,gf,
!$OMP& gfr,dorl,dorli)
!$OMP& firstprivate(mscale)
!$OMP DO reduction(+:emtt,viri,demi,demk)
!$OMP& schedule(guided)
      do i = 1, npole
c         i=atomind
         ii = ipole(i)
c         pdi = pdamp(i)
c         pti = thole(i)
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
            mscale(i12(j,ii)) = m2scale
c            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
c            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
c            pscale(i14(j,ii)) = p4scale
c            do k = 1, np11(ii)
c                if (i14(j,ii) .eq. ip11(k,ii))
c     &            pscale(i14(j,ii)) = p4scale * p41scale
c            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
c            pscale(i15(j,ii)) = p5scale
         end do
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
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
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
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
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
c               scale3 = 1.0d0
c               scale5 = 1.0d0
c               scale7 = 1.0d0
c               do j = 1, 3
c                  ddsc3(j) = 0.0d0
c                  ddsc5(j) = 0.0d0
c                  ddsc7(j) = 0.0d0
c               end do
c               damp = pdi * pdamp(k)
c               if (damp .ne. 0.0d0) then
c                  pgamma = min(pti,thole(k))
c                  damp = -pgamma * (r/damp)**3
c                  if (damp .gt. -50.0d0) then
c                     expdamp = exp(damp)
c                     scale3 = 1.0d0 - expdamp
c                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
c                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
c     &                                       *expdamp
c                     temp3 = -3.0d0 * damp * expdamp / r2
c                     temp5 = -damp
c                     temp7 = -0.2d0 - 0.6d0*damp
c                     ddsc3(1) = temp3 * xr
c                     ddsc3(2) = temp3 * yr
c                     ddsc3(3) = temp3 * zr
c                     ddsc5(1) = temp5 * ddsc3(1)
c                     ddsc5(2) = temp5 * ddsc3(2)
c                     ddsc5(3) = temp5 * ddsc3(3)
c                     ddsc7(1) = temp7 * ddsc5(1)
c                     ddsc7(2) = temp7 * ddsc5(2)
c                     ddsc7(3) = temp7 * ddsc5(3)
c                  end if
c               end if
c               dsc3 = 1.0d0 - scale3*dscale(kk)
c               dsc5 = 1.0d0 - scale5*dscale(kk)
c               dsc7 = 1.0d0 - scale7*dscale(kk)
c               psc3 = 1.0d0 - scale3*pscale(kk)
c               psc5 = 1.0d0 - scale5*pscale(kk)
c               psc7 = 1.0d0 - scale7*pscale(kk)
c               usc3 = 1.0d0 - scale3*uscale(kk)
c               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
c               dixuk(1) = di(2)*uind(3,k) - di(3)*uind(2,k)
c               dixuk(2) = di(3)*uind(1,k) - di(1)*uind(3,k)
c               dixuk(3) = di(1)*uind(2,k) - di(2)*uind(1,k)
c               dkxui(1) = dk(2)*uind(3,i) - dk(3)*uind(2,i)
c               dkxui(2) = dk(3)*uind(1,i) - dk(1)*uind(3,i)
c               dkxui(3) = dk(1)*uind(2,i) - dk(2)*uind(1,i)
c               dixukp(1) = di(2)*uinp(3,k) - di(3)*uinp(2,k)
c               dixukp(2) = di(3)*uinp(1,k) - di(1)*uinp(3,k)
c               dixukp(3) = di(1)*uinp(2,k) - di(2)*uinp(1,k)
c               dkxuip(1) = dk(2)*uinp(3,i) - dk(3)*uinp(2,i)
c               dkxuip(2) = dk(3)*uinp(1,i) - dk(1)*uinp(3,i)
c               dkxuip(3) = dk(1)*uinp(2,i) - dk(2)*uinp(1,i)
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
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
c               qiuk(1) = qi(1)*uind(1,k) + qi(4)*uind(2,k)
c     &                      + qi(7)*uind(3,k)
c               qiuk(2) = qi(2)*uind(1,k) + qi(5)*uind(2,k)
c     &                      + qi(8)*uind(3,k)
c               qiuk(3) = qi(3)*uind(1,k) + qi(6)*uind(2,k)
c     &                      + qi(9)*uind(3,k)
c               qkui(1) = qk(1)*uind(1,i) + qk(4)*uind(2,i)
c     &                      + qk(7)*uind(3,i)
c               qkui(2) = qk(2)*uind(1,i) + qk(5)*uind(2,i)
c     &                      + qk(8)*uind(3,i)
c               qkui(3) = qk(3)*uind(1,i) + qk(6)*uind(2,i)
c     &                      + qk(9)*uind(3,i)
c               qiukp(1) = qi(1)*uinp(1,k) + qi(4)*uinp(2,k)
c     &                       + qi(7)*uinp(3,k)
c               qiukp(2) = qi(2)*uinp(1,k) + qi(5)*uinp(2,k)
c     &                       + qi(8)*uinp(3,k)
c               qiukp(3) = qi(3)*uinp(1,k) + qi(6)*uinp(2,k)
c     &                       + qi(9)*uinp(3,k)
c               qkuip(1) = qk(1)*uinp(1,i) + qk(4)*uinp(2,i)
c     &                       + qk(7)*uinp(3,i)
c               qkuip(2) = qk(2)*uinp(1,i) + qk(5)*uinp(2,i)
c     &                       + qk(8)*uinp(3,i)
c               qkuip(3) = qk(3)*uinp(1,i) + qk(6)*uinp(2,i)
c     &                       + qk(9)*uinp(3,i)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
c               uixqkr(1) = uind(2,i)*qkr(3) - uind(3,i)*qkr(2)
c               uixqkr(2) = uind(3,i)*qkr(1) - uind(1,i)*qkr(3)
c               uixqkr(3) = uind(1,i)*qkr(2) - uind(2,i)*qkr(1)
c               ukxqir(1) = uind(2,k)*qir(3) - uind(3,k)*qir(2)
c               ukxqir(2) = uind(3,k)*qir(1) - uind(1,k)*qir(3)
c               ukxqir(3) = uind(1,k)*qir(2) - uind(2,k)*qir(1)
c               uixqkrp(1) = uinp(2,i)*qkr(3) - uinp(3,i)*qkr(2)
c               uixqkrp(2) = uinp(3,i)*qkr(1) - uinp(1,i)*qkr(3)
c               uixqkrp(3) = uinp(1,i)*qkr(2) - uinp(2,i)*qkr(1)
c               ukxqirp(1) = uinp(2,k)*qir(3) - uinp(3,k)*qir(2)
c               ukxqirp(2) = uinp(3,k)*qir(1) - uinp(1,k)*qir(3)
c               ukxqirp(3) = uinp(1,k)*qir(2) - uinp(2,k)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
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

c               sci(1) = uind(1,i)*dk(1) + uind(2,i)*dk(2)
c     &                     + uind(3,i)*dk(3) + di(1)*uind(1,k)
c     &                     + di(2)*uind(2,k) + di(3)*uind(3,k)
c               sci(2) = uind(1,i)*uind(1,k) + uind(2,i)*uind(2,k)
c     &                     + uind(3,i)*uind(3,k)
c               sci(3) = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
c               sci(4) = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
c               sci(7) = qir(1)*uind(1,k) + qir(2)*uind(2,k)
c     &                     + qir(3)*uind(3,k)
c               sci(8) = qkr(1)*uind(1,i) + qkr(2)*uind(2,i)
c     &                     + qkr(3)*uind(3,i)
c               scip(1) = uinp(1,i)*dk(1) + uinp(2,i)*dk(2)
c     &                      + uinp(3,i)*dk(3) + di(1)*uinp(1,k)
c     &                      + di(2)*uinp(2,k) + di(3)*uinp(3,k)
c               scip(2) = uind(1,i)*uinp(1,k)+uind(2,i)*uinp(2,k)
c     &                      + uind(3,i)*uinp(3,k)+uinp(1,i)*uind(1,k)
c     &                      + uinp(2,i)*uind(2,k)+uinp(3,i)*uind(3,k)
c               scip(3) = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
c               scip(4) = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
c               scip(7) = qir(1)*uinp(1,k) + qir(2)*uinp(2,k)
c     &                      + qir(3)*uinp(3,k)
c               scip(8) = qkr(1)*uinp(1,i) + qkr(2)*uinp(2,i)
c     &                      + qkr(3)*uinp(3,i)

c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
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
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
c               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
c     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))

c
c     get the real energy without any screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               erl = erl * (1.0d0-mscale(kk))
c               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
c     &                   + rr5*(gli(2)+gli(7))*psc5
c     &                   + rr7*gli(3)*psc7)
               e = e - erl
c               ei = ei - erli
               e = f * e
c               ei = f * ei
               emtt = emtt + e
c               emt = emt + e

c               eptt = eptt + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + mscale(kk)*erl*f
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     set flags to compute components without screening
c
               dorl = .false.
c               dorli = .false.
               if (mscale(kk) .ne. 1.0d0)  dorl = .true.
c               if (psc3 .ne. 0.0d0)  dorli = .true.
c               if (dsc3 .ne. 0.0d0)  dorli = .true.
c               if (usc3 .ne. 0.0d0)  dorli = .true.

c
c     zero out force and torque components without screening
c
               do j = 1, 3
                  ftm2r(j) = 0.0d0
c                  ftm2ri(j) = 0.0d0
                  ttm2r(j) = 0.0d0
c                  ttm2ri(j) = 0.0d0
                  ttm3r(j) = 0.0d0
c                  ttm3ri(j) = 0.0d0
               end do
c
c     get the permanent force with screening
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
               if (dorl) then
                  gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                        + rr7*(gl(2)+gl(7)+gl(8))
     &                        + rr9*(gl(3)+gl(5)) + rr11*gl(4)
                  gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
                  gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
                  gfr(4) = 2.0d0 * rr5
                  gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
                  gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
                  gfr(7) = 4.0d0 * rr7
                  ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
     &                          + gfr(4)*(qkdi(1)-qidk(1))
     &                          + gfr(5)*qir(1) + gfr(6)*qkr(1)
     &                          + gfr(7)*(qiqkr(1)+qkqir(1))
                  ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
     &                          + gfr(4)*(qkdi(2)-qidk(2))
     &                          + gfr(5)*qir(2) + gfr(6)*qkr(2)
     &                          + gfr(7)*(qiqkr(2)+qkqir(2))
                  ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
     &                          + gfr(4)*(qkdi(3)-qidk(3))
     &                          + gfr(5)*qir(3) + gfr(6)*qkr(3)
     &                          + gfr(7)*(qiqkr(3)+qkqir(3))
               end if
c
c     get the induced force with screening
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
c               ftm2i(1) = gfi(1)*xr + 0.5d0*
c     &             (gfi(2)*(uind(1,i)+uinp(1,i))
c     &            + bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
c     &            + gfi(3)*(uind(1,k)+uinp(1,k))
c     &            + bn(2)*(sci(3)*uinp(1,k)+scip(3)*uind(1,k))
c     &            + (sci(4)+scip(4))*bn(2)*di(1)
c     &            + (sci(3)+scip(3))*bn(2)*dk(1)
c     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
c     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
c               ftm2i(2) = gfi(1)*yr + 0.5d0*
c     &             (gfi(2)*(uind(2,i)+uinp(2,i))
c     &            + bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
c     &            + gfi(3)*(uind(2,k)+uinp(2,k))
c     &            + bn(2)*(sci(3)*uinp(2,k)+scip(3)*uind(2,k))
c     &            + (sci(4)+scip(4))*bn(2)*di(2)
c     &            + (sci(3)+scip(3))*bn(2)*dk(2)
c     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
c     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
c               ftm2i(3) = gfi(1)*zr + 0.5d0*
c     &             (gfi(2)*(uind(3,i)+uinp(3,i))
c     &            + bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
c     &            + gfi(3)*(uind(3,k)+uinp(3,k))
c     &            + bn(2)*(sci(3)*uinp(3,k)+scip(3)*uind(3,k))
c     &            + (sci(4)+scip(4))*bn(2)*di(3)
c     &            + (sci(3)+scip(3))*bn(2)*dk(3)
c     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
c     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)

c
c     get the induced force without screening
c
c               if (dorli) then
c                  gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
c     &                               + (glip(1)+glip(6))*dsc3
c     &                               + scip(2)*usc3)
c     &                    + 0.5d0*rr7*((gli(7)+gli(2))*psc5
c     &                               + (glip(7)+glip(2))*dsc5
c     &                        - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
c     &                    + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
c                  gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
c                  gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
c                  gfri(4) = 2.0d0 * rr5
c                  gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
c                  gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c                  ftm2ri(1) = gfri(1)*xr + 0.5d0*
c     &              (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
c     &               + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
c     &               - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
c     &               + (rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
c     &               + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
c     &               + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
c     &               + rr5*usc5*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
c     &               + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
c     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
c     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
c     &               + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
c     &               + (qkuip(1)-qiukp(1))*dsc5)
c     &               + gfri(5)*qir(1) + gfri(6)*qkr(1)
c                  ftm2ri(2) = gfri(1)*yr + 0.5d0*
c     &              (- rr3*ck*(uind(2,i)*psc3+uinp(2,i)*dsc3)
c     &               + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
c     &               - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
c     &               + (rr3*ci*(uind(2,k)*psc3+uinp(2,k)*dsc3)
c     &               + rr5*sc(3)*(uind(2,k)*psc5+uinp(2,k)*dsc5)
c     &               + rr7*sc(5)*(uind(2,k)*psc7+uinp(2,k)*dsc7))*0.5d0
c     &               + rr5*usc5*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
c     &               + sci(3)*uinp(2,k)+scip(3)*uind(2,k))*0.5d0
c     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
c     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
c     &               + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
c     &               + (qkuip(2)-qiukp(2))*dsc5)
c     &               + gfri(5)*qir(2) + gfri(6)*qkr(2)
c                  ftm2ri(3) = gfri(1)*zr + 0.5d0*
c     &              (- rr3*ck*(uind(3,i)*psc3+uinp(3,i)*dsc3)
c     &               + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
c     &               - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
c     &               + (rr3*ci*(uind(3,k)*psc3+uinp(3,k)*dsc3)
c     &               + rr5*sc(3)*(uind(3,k)*psc5+uinp(3,k)*dsc5)
c     &               + rr7*sc(5)*(uind(3,k)*psc7+uinp(3,k)*dsc7))*0.5d0
c     &               + rr5*usc5*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
c     &               + sci(3)*uinp(3,k)+scip(3)*uind(3,k))*0.5d0
c     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
c     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
c     &               + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
c     &               + (qkuip(3)-qiukp(3))*dsc5)
c     &               + gfri(5)*qir(3) + gfri(6)*qkr(3)
c               end if

c
c     account for partially excluded induced interactions
c

c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))
c               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
c     &                        + temp7*ddsc7(1)
c               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
c     &                        + temp7*ddsc7(2)
c               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c     &                        + temp7*ddsc7(3)

c
c     find some scaling terms for induced-induced force
c

c               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
c               temp5 = -0.5d0 * rr5 * uscale(kk)
c     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
c               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
c               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
c               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)

c
c     modify the forces for partially excluded interactions
c

c               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
c               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
c               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)

c
c     correction to convert mutual to direct polarization force
c
c               if (poltyp .eq. 'DIRECT') then
c                  gfd = 0.5d0 * (bn(2)*scip(2)
c     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
c                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
c     &                     - rr7*(scip(3)*sci(4)
c     &                           +sci(3)*scip(4))*usc5)
c                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
c     &                          +sci(3)*uinp(1,k)+scip(3)*uind(1,k))
c                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
c     &                          +sci(3)*uinp(2,k)+scip(3)*uind(2,k))
c                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
c     &                          +sci(3)*uinp(3,k)+scip(3)*uind(3,k))
c                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
c     &                        + sci(3)*uinp(1,k)+scip(3)*uind(1,k))
c                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
c     &                        + sci(3)*uinp(2,k)+scip(3)*uind(2,k))
c                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
c     &                        + sci(3)*uinp(3,k)+scip(3)*uind(3,k))
c                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
c                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
c                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
c               end if
c
c     get the permanent torque with screening
c
               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
     &                      + gf(4)*(dixqkr(1)+dkxqir(1)
     &                              +rxqidk(1)-2.0d0*qixqk(1))
     &                      - gf(5)*rxqir(1)
     &                      - gf(7)*(rxqikr(1)+qkrxqir(1))
               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
     &                      + gf(4)*(dixqkr(2)+dkxqir(2)
     &                              +rxqidk(2)-2.0d0*qixqk(2))
     &                      - gf(5)*rxqir(2)
     &                      - gf(7)*(rxqikr(2)+qkrxqir(2))
               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
     &                      + gf(4)*(dixqkr(3)+dkxqir(3)
     &                              +rxqidk(3)-2.0d0*qixqk(3))
     &                      - gf(5)*rxqir(3)
     &                      - gf(7)*(rxqikr(3)+qkrxqir(3))
               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
     &                      - gf(4)*(dixqkr(1)+dkxqir(1)
     &                              +rxqkdi(1)-2.0d0*qixqk(1))
     &                      - gf(6)*rxqkr(1)
     &                      - gf(7)*(rxqkir(1)-qkrxqir(1))
               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
     &                      - gf(4)*(dixqkr(2)+dkxqir(2)
     &                              +rxqkdi(2)-2.0d0*qixqk(2))
     &                      - gf(6)*rxqkr(2)
     &                      - gf(7)*(rxqkir(2)-qkrxqir(2))
               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
     &                      - gf(4)*(dixqkr(3)+dkxqir(3)
     &                              +rxqkdi(3)-2.0d0*qixqk(3))
     &                      - gf(6)*rxqkr(3)
     &                      - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
               if (dorl) then
                  ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)
     &                          + gfr(4)*(dixqkr(1)+dkxqir(1)
     &                                   +rxqidk(1)-2.0d0*qixqk(1))
     &                          - gfr(5)*rxqir(1)
     &                          - gfr(7)*(rxqikr(1)+qkrxqir(1))
                  ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)
     &                          + gfr(4)*(dixqkr(2)+dkxqir(2)
     &                                   +rxqidk(2)-2.0d0*qixqk(2))
     &                          - gfr(5)*rxqir(2)
     &                          - gfr(7)*(rxqikr(2)+qkrxqir(2))
                  ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)
     &                          + gfr(4)*(dixqkr(3)+dkxqir(3)
     &                                   +rxqidk(3)-2.0d0*qixqk(3))
     &                          - gfr(5)*rxqir(3)
     &                          - gfr(7)*(rxqikr(3)+qkrxqir(3))
                  ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1)
     &                          - gfr(4)*(dixqkr(1)+dkxqir(1)
     &                                   +rxqkdi(1)-2.0d0*qixqk(1))
     &                          - gfr(6)*rxqkr(1)
     &                          - gfr(7)*(rxqkir(1)-qkrxqir(1))
                  ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2)
     &                          - gfr(4)*(dixqkr(2)+dkxqir(2)
     &                                   +rxqkdi(2)-2.0d0*qixqk(2))
     &                          - gfr(6)*rxqkr(2)
     &                          - gfr(7)*(rxqkir(2)-qkrxqir(2))
                  ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3)
     &                          - gfr(4)*(dixqkr(3)+dkxqir(3)
     &                                   +rxqkdi(3)-2.0d0*qixqk(3))
     &                          - gfr(6)*rxqkr(3)
     &                          - gfr(7)*(rxqkir(3)-qkrxqir(3))
               end if
c
c     get the induced torque with screening
c

c               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
c               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
c               gti(4) = gfi(4)
c               gti(5) = gfi(5)
c               gti(6) = gfi(6)
c               ttm2i(1) = -0.5d0*bn(1)*(dixuk(1)+dixukp(1))
c     &                       + gti(2)*dixr(1) - gti(5)*rxqir(1)
c     &                       + 0.5d0*gti(4)*(ukxqir(1)+rxqiuk(1)
c     &                                      +ukxqirp(1)+rxqiukp(1))
c               ttm2i(2) = -0.5d0*bn(1)*(dixuk(2)+dixukp(2))
c     &                       + gti(2)*dixr(2) - gti(5)*rxqir(2)
c     &                       + 0.5d0*gti(4)*(ukxqir(2)+rxqiuk(2)
c     &                                      +ukxqirp(2)+rxqiukp(2))
c               ttm2i(3) = -0.5d0*bn(1)*(dixuk(3)+dixukp(3))
c     &                       + gti(2)*dixr(3) - gti(5)*rxqir(3)
c     &                       + 0.5d0*gti(4)*(ukxqir(3)+rxqiuk(3)
c     &                                      +ukxqirp(3)+rxqiukp(3))
c               ttm3i(1) = -0.5d0*bn(1)*(dkxui(1)+dkxuip(1))
c     &                       + gti(3)*dkxr(1) - gti(6)*rxqkr(1)
c     &                       - 0.5d0*gti(4)*(uixqkr(1)+rxqkui(1)
c     &                                      +uixqkrp(1)+rxqkuip(1))
c               ttm3i(2) = -0.5d0*bn(1)*(dkxui(2)+dkxuip(2))
c     &                       + gti(3)*dkxr(2) - gti(6)*rxqkr(2)
c     &                       - 0.5d0*gti(4)*(uixqkr(2)+rxqkui(2)
c     &                                       +uixqkrp(2)+rxqkuip(2))
c               ttm3i(3) = -0.5d0*bn(1)*(dkxui(3)+dkxuip(3))
c     &                       + gti(3)*dkxr(3) - gti(6)*rxqkr(3)
c     &                       - 0.5d0*gti(4)*(uixqkr(3)+rxqkui(3)
c     &                                      +uixqkrp(3)+rxqkuip(3))

c
c     get the induced torque without screening
c

c               if (dorli) then
c                  gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
c                  gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
c                  gtri(4) = gfri(4)
c                  gtri(5) = gfri(5)
c                  gtri(6) = gfri(6)
c                  ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
c     &                           + gtri(2)*dixr(1) - gtri(5)*rxqir(1)
c     &                           + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
c     &                             +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0
c                  ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
c     &                           + gtri(2)*dixr(2) - gtri(5)*rxqir(2)
c     &                           + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
c     &                             +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0
c                  ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
c     &                           + gtri(2)*dixr(3) - gtri(5)*rxqir(3)
c     &                           + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
c     &                             +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0
c                  ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
c     &                           + gtri(3)*dkxr(1) - gtri(6)*rxqkr(1)
c     &                           - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
c     &                             +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0
c                  ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
c     &                           + gtri(3)*dkxr(2) - gtri(6)*rxqkr(2)
c     &                           - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
c     &                             +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0
c                  ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
c     &                           + gtri(3)*dkxr(3) - gtri(6)*rxqkr(3)
c     &                           - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
c     &                             +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0
c               end if

c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
c                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
c                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
c                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c     increment gradient due to force and torque on first site
c

               demi(1,ii) = demi(1,ii) + ftm2(1)
               demi(2,ii) = demi(2,ii) + ftm2(2)
               demi(3,ii) = demi(3,ii) + ftm2(3)
c               depi(1,ii) = depi(1,ii) + ftm2i(1)
c               depi(2,ii) = depi(2,ii) + ftm2i(2)
c               depi(3,ii) = depi(3,ii) + ftm2i(3)
c               call torque3 (i,ttm2,ttm2i,frcxi,frcyi,frczi,demi,depi)
c               demt(1,ii) = demt(1,ii) + ftm2(1)
c               demt(2,ii) = demt(2,ii) + ftm2(2)
c               demt(3,ii) = demt(3,ii) + ftm2(3)
c               call torque_3b_Perm_new (demt,
c     &         i,ttm2,ttm2i,frcxi,frcyi,frczi)
               call torque_3b_Perm_new (demi,
     &         i,ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

               demk(1,kk) = demk(1,kk) - ftm2(1)
               demk(2,kk) = demk(2,kk) - ftm2(2)
               demk(3,kk) = demk(3,kk) - ftm2(3)
c               depk(1,kk) = depk(1,kk) - ftm2i(1)
c               depk(2,kk) = depk(2,kk) - ftm2i(2)
c               depk(3,kk) = depk(3,kk) - ftm2i(3)
c               call torque3 (k,ttm3,ttm3i,frcxk,frcyk,frczk,demk,depk)
c               demt(1,kk) = demt(1,kk) - ftm2(1)
c               demt(2,kk) = demt(2,kk) - ftm2(2)
c               demt(3,kk) = demt(3,kk) - ftm2(3)
c               call torque_3b_Perm_new (demt,
c     &          k,ttm3,ttm3i,frcxk,frcyk,frczk)
               call torque_3b_Perm_new (demk,
     &          k,ttm3,ttm3i,frcxk,frcyk,frczk)


c
c     increment the internal virial tensor components
c
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
c               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
c     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
c     &                  + xky*frcyk(1) + xkz*frczk(1)
c               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
c     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
c     &                  + yky*frcyk(1) + ykz*frczk(1)
c               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
c     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
c     &                  + zky*frcyk(1) + zkz*frczk(1)
c               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
c     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
c     &                  + yky*frcyk(2) + ykz*frczk(2)
c               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
c     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
c     &                  + zky*frcyk(2) + zkz*frczk(2)
c               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
c     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
c     &                  + zky*frcyk(3) + zkz*frczk(3)

               vxx = -xr*(ftm2(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               viri(1,1) = viri(1,1) + vxx
               viri(2,1) = viri(2,1) + vyx
               viri(3,1) = viri(3,1) + vzx
               viri(1,2) = viri(1,2) + vyx
               viri(2,2) = viri(2,2) + vyy
               viri(3,2) = viri(3,2) + vzy
               viri(1,3) = viri(1,3) + vzx
               viri(2,3) = viri(2,3) + vzy
               viri(3,3) = viri(3,3) + vzz
c               viremt(1,1) = viremt(1,1) + vxx
c               viremt(2,1) = viremt(2,1) + vyx
c               viremt(3,1) = viremt(3,1) + vzx
c               viremt(1,2) = viremt(1,2) + vyx
c               viremt(2,2) = viremt(2,2) + vyy
c               viremt(3,2) = viremt(3,2) + vzy
c               viremt(1,3) = viremt(1,3) + vzx
c               viremt(2,3) = viremt(2,3) + vzy
c               viremt(3,3) = viremt(3,3) + vzz

            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
c            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
c            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
c            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
c            pscale(i15(j,ii)) = 1.0d0
         end do
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
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
c      print*,"Before add emreal_tmp",emreal_tmp
c      print*,"Before add demreal_tmp",demreal_tmp
      emreal_tmp = emtt
c      print*,"emreal_tmp in ereal1d_perm",emreal_tmp
c      emreal=emreal+emtt 
cccc      ep = ep + eptt
      do i = 1, n
         do j = 1, 3
c            demreal(j,i) = demreal(j,i) + demi(j,i) + demk(j,i)
            demreal_tmp(j,i) =  demi(j,i) + demk(j,i)
cccc            dep(j,i) = dep(j,i) + depi(j,i) + depk(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
c            viremreal(j,i)=viremreal(j,i)+viri(j,i)
            viremreal_tmp(j,i) = viri(j,i)
         end do
      end do

c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
c      deallocate (pscale)
c      deallocate (dscale)
c      deallocate (uscale)
      deallocate (demi)
      deallocate (demk)
c      deallocate (depi)
c      deallocate (depk)
      return
      end
c
      subroutine emrecip1_3b_Perm2
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use virial
      use neigh2clust
      implicit none
      integer i,j,k,ii
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 cphim(4),cphid(4)
      real*8 cphip(4)
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time
      
c
c     derivative indices into the fphi and fphidp arrays
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
c      print*,"In emrecip1_3b_Perm npole",npole
c      print*,"In emrecip1_3b_Perm n",n
c      print*,"In emrecip1_3b_Perm aewald",aewald
      allocate (frc(3,n))
      allocate (trq(3,npole))
c      allocate (fuind(3,npole))
c      allocate (fuinp(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (fphi(20,npole))
c      allocate (fphid(10,npole))
c      allocate (fphip(10,npole))
c      allocate (fphidp(20,npole))
      allocate (cphi(10,npole))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     get the fractional to Cartesian transformation matrix
c

c      call system_clock(t1,clock_rate)
      call frac_to_cart (ftc)
c      call system_clock(t2, clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing frac_to_cart", tot_time
c
c     compute the arrays of B-spline coefficients
c
c      if (.not. use_polar) then

c      call system_clock(t1,clock_rate)
      if(.not.use_ewaldclust) then
      !call bspline_fill
      call bspline_fill_omp
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing bspline_fill", tot_time

c      call system_clock(t1,clock_rate)
      !call table_fill
      call table_fill_omp
      end if
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing table fill", tot_time

c      end if
c
c     perform dynamic allocation of some local arrays
c
c      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c

c      if (use_polar) then
c         end do
c      else
c      call system_clock(t1,clock_rate)
      call cmp_to_fmp (cmp,fmp)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing cmp_to_fmp", tot_time

c      call system_clock(t1,clock_rate)
      call grid_mpole (fmp)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing grid_mpole", tot_time

c      call system_clock(t1,clock_rate)
      call fftfront
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timiing fftfront", tot_time

      ! qgrip does not appear to be necessary for just multipoles  
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
c                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
c               end do
c            end do
c         end do
c      end if

c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)*qgrid(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrid(2,k1,k2,k3)
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vyx = vyx + h2*h1*vterm
            vzx = vzx + h3*h1*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vzy = vzy + h3*h2*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
       ! Condense step below w/ current loop
       !  qfac(k1,k2,k3) = expterm
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
      end do
c
c     assign just the induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c

c      if (use_polar .and. poltyp.eq.'DIRECT') then
c      end if
c
c     perform deallocation of some local arrays
c
c      deallocate (qgrip)
c
c     transform permanent multipoles without induced dipoles
c

c      if (use_polar) then
c      end if

c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = 0.5d0 * expterm * struc2
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c

c      do k = 1, nfft3
c         do j = 1, nfft2
c            do i = 1, nfft1
c               term = qfac(i,j,k)
c               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
c               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
c            end do
c         end do
c      end do


c
c     perform 3-D FFT backward transform and get potential
c
c      call system_clock(t1,clock_rate)
      call fftback
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing fftback", tot_time
      
c      call system_clock(t1,clock_rate)
      call fphi_mpole (fphi)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing fphi_mpole", tot_time

      do i = 1, npole
         do j = 1, 20
            fphi(j,i) = electric * fphi(j,i)
         end do
      end do

c      call system_clock(t1,clock_rate)
      call fphi_to_cphi (fphi,cphi)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing fphi_to_cphi", tot_time
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
        ! Condensed from below
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + frc(1,i)
         dem(2,ii) = dem(2,ii) + frc(2,i)
         dem(3,ii) = dem(3,ii) + frc(3,i)

         trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &                 + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &                 + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &                 - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
         trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &                 + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &                 + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &                 - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
         trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &                 + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &                 + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &                 - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)

         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0

         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
     &             - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
         vyx = vyx - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
     &             - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &             - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
     &             - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
         vzx = vzx - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
     &             - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &             - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &             - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
     &             - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
         vzy = vzy - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
     &             - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &             - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &             - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
     &             - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)

      end do
      e = 0.5d0 * e
      em = em + e

      !do i = 1, npole
      !   ii = ipole(i)
      !   dem(1,ii) = dem(1,ii) + frc(1,i)
      !   dem(2,ii) = dem(2,ii) + frc(2,i)
      !   dem(3,ii) = dem(3,ii) + frc(3,i)
      !end do
c
c     distribute torques into the permanent multipole gradient
c
      !do i = 1, npole
      !   trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
c     &                 + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
c     &                 + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
c     &                 - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
c         trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
c     &                 + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
c     &                 + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
c     &                 - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
c         trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
c     &                 + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
c     &                 + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
c     &                 - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
c      end do
      !do i = 1, n
      !   frc(1,i) = 0.0d0
      !   frc(2,i) = 0.0d0
      !   frc(3,i) = 0.0d0
      !end do

c      call system_clock(t1,clock_rate)

      !call torque2 (trq,frc)
      call torque2_omp (trq,frc)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "timing torque2", tot_time

      do i = 1, n
         dem(1,i) = dem(1,i) + frc(1,i)
         dem(2,i) = dem(2,i) + frc(2,i)
         dem(3,i) = dem(3,i) + frc(3,i)
      end do
c
c     permanent multipole contribution to the internal virial
c
c      do i = 1, npole
c         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
c     &             - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
c         vyx = vyx - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
c     &             - (cmp(5,i)+cmp(6,i))*cphi(8,i)
c     &             - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
c     &             - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
c         vzx = vzx - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
c     &             - (cmp(5,i)+cmp(7,i))*cphi(9,i)
c     &             - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
c     &             - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
c         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
c     &             - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
c         vzy = vzy - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
c     &             - (cmp(6,i)+cmp(7,i))*cphi(10,i)
c     &             - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
c     &             - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
c         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
c     &             - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)
c      end do
c
c     convert Cartesian induced dipoles to fractional coordinates
c

c      if (use_polar) then
c      end if

c
c     increment the internal virial tensor components
c
      viremrecip(1,1) = viremrecip(1,1) + vxx
      viremrecip(2,1) = viremrecip(2,1) + vyx
      viremrecip(3,1) = viremrecip(3,1) + vzx
      viremrecip(1,2) = viremrecip(1,2) + vyx
      viremrecip(2,2) = viremrecip(2,2) + vyy
      viremrecip(3,2) = viremrecip(3,2) + vzy
      viremrecip(1,3) = viremrecip(1,3) + vzx
      viremrecip(2,3) = viremrecip(2,3) + vzy
      viremrecip(3,3) = viremrecip(3,3) + vzz
c
c     perform deallocation of some local arrays
c
c      print*,"End of emrecip1_3b_Perm",em
      deallocate (frc)
      deallocate (trq)
c      deallocate (fuind)
c      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
c      deallocate (fphid)
c      deallocate (fphip)
c      deallocate (fphidp)
      deallocate (cphi)
      return
      end
c     ################################################################
c     ##                                                            ##
c     ##              Subroutine empole1c_3b_Perm                   ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c_3b_Perm" calculates the multipole energy and derivatives 
c     with respect to Cartesian coordinates using regular Ewald 
c
c
      subroutine empole1d_3b_Perm_selfeng
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polar
      use polpot
      use virial
      use mplpot
      implicit none
      integer i,j,ii
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
      real*8 em_recip

c     set the energy unit conversion factor
c
      f = electric / dielec

      em=em+emreal

      do i=1,npole
         do j=1,3
            dem(j,i)=dem(j,i)+demreal(j,i)
         end do 
      end do 
c      print*,"After ereal1c_3b_Perm em real",em-em_recip
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c

      if (boundary .eq. 'VACUUM') then
        trqi(1) = 0.0d0
        trqi(2) = 0.0d0
        trqi(3) = 0.0d0
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xu = 0.0d0
        yu = 0.0d0
        zu = 0.0d0
        xup = 0.0d0
        yup = 0.0d0
        zup = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
          yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
          zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
        end do
        term = (2.0d0/3.0d0) * f * (pi/volbox)
        em = em + term*(xd*xd+yd*yd+zd*zd)
        do i = 1, npole
          ii = ipole(i)
          dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
          dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
          dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
        end do
        xdfield = -2.0d0 * term * xd
        ydfield = -2.0d0 * term * yd
        zdfield = -2.0d0 * term * zd
        do i = 1, npole
          trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
          trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
          trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
c          call torque_3b_Perm (i,trq,trqi,frcx,frcy,frcz)
c          call torque_3b_Perm_new (demk,
c     &          k,ttm3,ttm3i,frcxk,frcyk,frczk)
          call torque_3b_Perm_new (dem,
     &          i,trq,trqi,frcx,frcy,frcz)

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
        do i = 1, npole
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
        vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
        vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
        vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
        vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
        vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
        vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
        vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
        vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
        vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
        if (poltyp .eq. 'DIRECT') then
          vterm = term * (xu*xup+yu*yup+zu*zup)
          vir(1,1) = vir(1,1) + vterm
          vir(2,2) = vir(2,2) + vterm
          vir(3,3) = vir(3,3) + vterm
        end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      return
      end
c
c
      subroutine empole1d_3b_Perm_selfeng_bcast
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polar
      use polpot
      use virial
      use mplpot
      implicit none
      integer i,j,ii
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
      real*8 em_recip

c     set the energy unit conversion factor
c
      f = electric / dielec

      em=em+emreal_tmp

      do i=1,npole
         do j=1,3
            dem(j,i)=dem(j,i)+demreal_tmp(j,i)
         end do 
      end do 
c      print*,"After ereal1c_3b_Perm em real",em-em_recip
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c

      if (boundary .eq. 'VACUUM') then
        trqi(1) = 0.0d0
        trqi(2) = 0.0d0
        trqi(3) = 0.0d0
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xu = 0.0d0
        yu = 0.0d0
        zu = 0.0d0
        xup = 0.0d0
        yup = 0.0d0
        zup = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
          yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
          zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
        end do
        term = (2.0d0/3.0d0) * f * (pi/volbox)
        em = em + term*(xd*xd+yd*yd+zd*zd)
        do i = 1, npole
          ii = ipole(i)
          dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
          dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
          dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
        end do
        xdfield = -2.0d0 * term * xd
        ydfield = -2.0d0 * term * yd
        zdfield = -2.0d0 * term * zd
        do i = 1, npole
          trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
          trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
          trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
c          call torque_3b_Perm (i,trq,trqi,frcx,frcy,frcz)
c          call torque_3b_Perm_new (demk,
c     &          k,ttm3,ttm3i,frcxk,frcyk,frczk)
          call torque_3b_Perm_new (dem,
     &          i,trq,trqi,frcx,frcy,frcz)

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
        do i = 1, npole
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
        vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
        vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
        vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
        vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
        vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
        vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
        vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
        vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
        vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
        if (poltyp .eq. 'DIRECT') then
          vterm = term * (xu*xup+yu*yup+zu*zup)
          vir(1,1) = vir(1,1) + vterm
          vir(2,2) = vir(2,2) + vterm
          vir(3,3) = vir(3,3) + vterm
        end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      return
      end
c
