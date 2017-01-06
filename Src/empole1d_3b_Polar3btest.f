      subroutine empole1d_3b_Polar3btest(npole3b,pnum,
     &   eptemp,deptemp,virtemp,cnt,
     &   emtemp,demtemp)
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
      use polpot
      use polar, only: polarity, thole, pdamp
      use cho
      use pcg
      use neigh2clust
      use mpidat
      implicit none
      include 'mpif.h'
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
      integer npole3b
      integer pnum(*),l1,l3
      real*8 eptemp,deptemp(3,*),virtemp(3,3)
      real*8 emtemp,demtemp(3,*),viremtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      integer cnt
      real*8 t1,t2,t3,t4,t5,t6,t7,t8

            t1=mpi_Wtime()

      eptemp = 0.0d0
      emtemp = 0.0d0
c
c   Zero out temporary gradient of polarization energy
c

c      do l1 = 1, npole3b
      do i = 1, npole
        do j = 1, 3
          deptemp(j,i) = 0.0d0
          demtemp(j,i) = 0.0d0
        end do
      end do

c
c     set the energy unit conversion factor
c
      f = electric / dielec

c
c   Zero out temporary virial
c
      do i=1,3
         do j=1,3
           virtemp(j,i)=0.0d0
           viremtemp(j,i)=0.0d0
         end do
      end do

c      if(usecholesky) then
c       call induce0a_3b_cholesky_nopriordir(npole3b,pnum,
c     & uind,uinp) 
c      else if (usepcgscf) then
       call induce0a_3b_scf_vacomp3bclust(npole3b,pnum,uind,uinp,cnt)
c       print*,"Done w induce0a_3b_scf_vacomp1bclust"
c      else
c       call induce0a_3b_PolelecOnly_orig_nopriordir(npole3b,pnum,
c     & uind,uinp)
c      end if


      call emrecip1_3b_Polartest(npole3b,pnum,uind,uinp,
     & eptemp,deptemp,virtemp,emtemp,demtemp)
c       call emrecip1_3b_PolarPermtest(npole3b,pnum,uind,uinp,
c     & eptemp,deptemp,virtemp,emtemp,demtemp) 
      !print*,"Done w emrecip1_3b_Polar_new in 1bclust"

      call ereal1d_3b_3bclustPolartest(npole3b,pnum,uind,uinp,cnt,
     & eptemp,deptemp,virtemp,emtemp,demtemp)

c      print*,"Done w ereal1d_3b_3bclust in 3bclust"

      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      !do i = 1, npole
      do l1 =1, npole3b
         i=pnum(l1)
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
         uix = uind(1,l1)
         uiy = uind(2,l1)
         uiz = uind(3,l1)
c         cii = ci*ci
c         dii = dix*dix + diy*diy + diz*diz
c         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
c     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
c         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
c         emtemp = emtemp + e
         eptemp = eptemp + ei
      end do

      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      !do i = 1, npole
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
         !call torque (i,trq,trqi,frcx,frcy,frcz)
c        if(i.eq.0) then
c       print*,"Tilt!empole1d_PolarPerm3b b4 SelfTorq",taskid,cnt,npole3b
c        else 
c       print*,"In empole1d_PolarPerm3b b4 SelfTorque",taskid,cnt,npole3b
c        end if
c         call torque3_PolarPerm (npole3b,pnum,
c     &    i,trq,trqi,frcx,frcy,frcz,demtemp,deptemp)
         call torque3_Polar (
     &    i,trq,trqi,frcx,frcy,frcz,deptemp)
      end do
            t2=mpi_Wtime()
       !print*,"tsk=",taskid,"TotalCost Within empole1d_3b_3btest",t2-t1
c      print*,"Completed all steps of empole1d 3bclust"
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ereal1d  --  ewald real space derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1d" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c     and dipole polarizability
c
c
      subroutine ereal1d_3b_3bclustPolartest(npole3b,pnum,uind,uinp,
     & cnt,eptemp,deptemp,virtemp,emtemp,demtemp)
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use limits
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use shunt
      use virial
      use neigh2clust 
      implicit none
      integer i,j,k
      integer ii,kk,kkk
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
      real*8 emo,epo,eintrao
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
      real*8 viro(3,3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: demo1(:,:)
      real*8, allocatable :: demo2(:,:)
      real*8, allocatable :: depo1(:,:)
      real*8, allocatable :: depo2(:,:)
      logical dorl,dorli
      character*6 mode
      integer npole3b,pnum(*),cnt,l1,l3
      real*8 uind(3,*),eptemp,deptemp(3,*)
      real*8 virtemp(3,3),uinp(3,*)
      real*8 emtemp,demtemp(3,*),viremtemp(3,3)
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(npole))
      allocate (pscale(npole))
      allocate (dscale(npole))
      allocate (uscale(npole))
c      allocate (demo1(3,npole))
c      allocate (demo2(3,npole))
      allocate (depo1(3,npole))
      allocate (depo2(3,npole))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, npole
         mscale(i) = 1.0d0
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
c
c     initialize local variables for OpenMP calculation
c
c      emo = 0.0d0
      epo = 0.0d0
c      eintrao = eintra
      do i = 1, npole
         do j = 1, 3
c            demo1(j,i) = 0.0d0
c            demo2(j,i) = 0.0d0
            depo1(j,i) = 0.0d0
            depo2(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            viro(j,i) = 0.0d0
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) firstprivate(f)
!$OMP& private(i,j,k,ii,kk,kkk,e,ei,bfac,damp,expdamp,
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
!$OMP& gfr,gfri,gti,gtri,dorl,dorli,l1,l3)
!$OMP& firstprivate(mscale,pscale,dscale,uscale)
!$OMP DO reduction(+:epo,depo1,depo2,
!$OMP& viro)
!$OMP& schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      !do i = 1, npole
      do l1 =1, npole3b
         i=pnum(l1)
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
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
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

         !do kkk = 1, nelst(i)
         !   k = elst(kkk,i)
         do kkk = 1,nelst3b(l1,cnt)
            l3=elst3b(kkk,l1,cnt)
            k=pnum(l3)
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
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l3) - di(3)*uind(2,l3)
               dixuk(2) = di(3)*uind(1,l3) - di(1)*uind(3,l3)
               dixuk(3) = di(1)*uind(2,l3) - di(2)*uind(1,l3)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l3) - di(3)*uinp(2,l3)
               dixukp(2) = di(3)*uinp(1,l3) - di(1)*uinp(3,l3)
               dixukp(3) = di(1)*uinp(2,l3) - di(2)*uinp(1,l3)
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
               qiuk(1) = qi(1)*uind(1,l3) + qi(4)*uind(2,l3)
     &                      + qi(7)*uind(3,l3)
               qiuk(2) = qi(2)*uind(1,l3) + qi(5)*uind(2,l3)
     &                      + qi(8)*uind(3,l3)
               qiuk(3) = qi(3)*uind(1,l3) + qi(6)*uind(2,l3)
     &                      + qi(9)*uind(3,l3)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l3) + qi(4)*uinp(2,l3)
     &                       + qi(7)*uinp(3,l3)
               qiukp(2) = qi(2)*uinp(1,l3) + qi(5)*uinp(2,l3)
     &                       + qi(8)*uinp(3,l3)
               qiukp(3) = qi(3)*uinp(1,l3) + qi(6)*uinp(2,l3)
     &                       + qi(9)*uinp(3,l3)
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
               ukxqir(1) = uind(2,l3)*qir(3) - uind(3,l3)*qir(2)
               ukxqir(2) = uind(3,l3)*qir(1) - uind(1,l3)*qir(3)
               ukxqir(3) = uind(1,l3)*qir(2) - uind(2,l3)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l3)*qir(3) - uinp(3,l3)*qir(2)
               ukxqirp(2) = uinp(3,l3)*qir(1) - uinp(1,l3)*qir(3)
               ukxqirp(3) = uinp(1,l3)*qir(2) - uinp(2,l3)*qir(1)
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
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l3)
     &                     + di(2)*uind(2,l3) + di(3)*uind(3,l3)
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = qir(1)*uind(1,l3) + qir(2)*uind(2,l3)
     &                     + qir(3)*uind(3,l3)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l3)
     &                      + di(2)*uinp(2,l3) + di(3)*uinp(3,l3)
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                     + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                     + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qir(1)*uinp(1,l3) + qir(2)*uinp(2,l3)
     &                      + qir(3)*uinp(3,l3)
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
c               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
c     &                + bn(2)*(gl(2)+gl(7)+gl(8))
c     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
c               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
c     &                  + rr5*(gl(2)+gl(7)+gl(8))
c     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
c               e = e - erl * (1.0d0-mscale(kk))
               ei = ei - erli
c               e = f * e
               ei = f * ei
c               emo = emo + e
               epo = epo + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintrao = eintrao + mscale(kk)*erl*f
c                  eintrao = eintrao + 0.5d0*pscale(kk)
c     &                         * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     set flags to compute components without screening
c
               dorl = .false.
               dorli = .false.
               if (mscale(kk) .ne. 1.0d0)  dorl = .true.
               if (psc3 .ne. 0.0d0)  dorli = .true.
               if (dsc3 .ne. 0.0d0)  dorli = .true.
               if (usc3 .ne. 0.0d0)  dorli = .true.
c
c     zero out force and torque components without screening
c
               do j = 1, 3
                  ftm2r(j) = 0.0d0
                  ftm2ri(j) = 0.0d0
                  ttm2r(j) = 0.0d0
                  ttm2ri(j) = 0.0d0
                  ttm3r(j) = 0.0d0
                  ttm3ri(j) = 0.0d0

                  ftm2(j) = 0.0d0
                  ttm2(j) = 0.0d0
                  ttm3(j) = 0.0d0
               end do
c
c     get the permanent force with screening
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
c               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
c     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
c     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
c               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
c     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
c     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
c               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
c     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
c     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
c               if (dorl) then
c                  gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
c     &                        + rr7*(gl(2)+gl(7)+gl(8))
c     &                        + rr9*(gl(3)+gl(5)) + rr11*gl(4)
c                  gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
c                  gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
c                  gfr(4) = 2.0d0 * rr5
c                  gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
c                  gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
c                  gfr(7) = 4.0d0 * rr7
c                  ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
c     &                          + gfr(4)*(qkdi(1)-qidk(1))
c     &                          + gfr(5)*qir(1) + gfr(6)*qkr(1)
c     &                          + gfr(7)*(qiqkr(1)+qkqir(1))
c                  ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
c     &                          + gfr(4)*(qkdi(2)-qidk(2))
c     &                          + gfr(5)*qir(2) + gfr(6)*qkr(2)
c     &                          + gfr(7)*(qiqkr(2)+qkqir(2))
c                  ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
c     &                          + gfr(4)*(qkdi(3)-qidk(3))
c     &                          + gfr(5)*qir(3) + gfr(6)*qkr(3)
c     &                          + gfr(7)*(qiqkr(3)+qkqir(3))
c               end if
c
c     get the induced force with screening
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
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l3)+uinp(1,l3))
     &            + bn(2)*(sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l3)+uinp(2,l3))
     &            + bn(2)*(sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l3)+uinp(3,l3))
     &            + bn(2)*(sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               if (dorli) then
                  gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                               + (glip(1)+glip(6))*dsc3
     &                               + scip(2)*usc3)
     &                    + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                        - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                    + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
                  gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
                  gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
                  gfri(4) = 2.0d0 * rr5
                  gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
                  gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
                  ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &              (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &               + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &               - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &               + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
     &               + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
     &             + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
     &               + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &               + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &               + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &               + (qkuip(1)-qiukp(1))*dsc5)
     &               + gfri(5)*qir(1) + gfri(6)*qkr(1)
                  ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &              (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &               + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &               - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &               + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
     &               + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
     &             + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
     &               + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &               + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &               + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &               + (qkuip(2)-qiukp(2))*dsc5)
     &               + gfri(5)*qir(2) + gfri(6)*qkr(2)
                  ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &              (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &               + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &               - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &               + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
     &               + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
     &             + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
     &               + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &               + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &               + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &               + (qkuip(3)-qiukp(3))*dsc5)
     &               + gfri(5)*qir(3) + gfri(6)*qkr(3)
               end if
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
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
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (bn(2)*scip(2)
     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
     &                     - rr7*(scip(3)*sci(4)
     &                           +sci(3)*scip(4))*usc5)
                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &                          +sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &                          +sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &                          +sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &                        + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &                        + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &                        + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if
c
c     get the permanent torque with screening
c
c               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
c     &                      + gf(4)*(dixqkr(1)+dkxqir(1)
c     &                              +rxqidk(1)-2.0d0*qixqk(1))
c     &                      - gf(5)*rxqir(1)
c     &                      - gf(7)*(rxqikr(1)+qkrxqir(1))
c               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
c     &                      + gf(4)*(dixqkr(2)+dkxqir(2)
c     &                              +rxqidk(2)-2.0d0*qixqk(2))
c     &                      - gf(5)*rxqir(2)
c     &                      - gf(7)*(rxqikr(2)+qkrxqir(2))
c               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
c     &                      + gf(4)*(dixqkr(3)+dkxqir(3)
c     &                              +rxqidk(3)-2.0d0*qixqk(3))
c     &                      - gf(5)*rxqir(3)
c     &                      - gf(7)*(rxqikr(3)+qkrxqir(3))
c               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
c     &                      - gf(4)*(dixqkr(1)+dkxqir(1)
c     &                              +rxqkdi(1)-2.0d0*qixqk(1))
c     &                      - gf(6)*rxqkr(1)
c     &                      - gf(7)*(rxqkir(1)-qkrxqir(1))
c               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
c     &                      - gf(4)*(dixqkr(2)+dkxqir(2)
c     &                              +rxqkdi(2)-2.0d0*qixqk(2))
c     &                      - gf(6)*rxqkr(2)
c     &                      - gf(7)*(rxqkir(2)-qkrxqir(2))
c               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
c     &                      - gf(4)*(dixqkr(3)+dkxqir(3)
c     &                              +rxqkdi(3)-2.0d0*qixqk(3))
c     &                      - gf(6)*rxqkr(3)
c     &                      - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
c               if (dorl) then
c                  ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)
c     &                          + gfr(4)*(dixqkr(1)+dkxqir(1)
c     &                                   +rxqidk(1)-2.0d0*qixqk(1))
c     &                          - gfr(5)*rxqir(1)
c     &                          - gfr(7)*(rxqikr(1)+qkrxqir(1))
c                  ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)
c     &                          + gfr(4)*(dixqkr(2)+dkxqir(2)
c     &                                   +rxqidk(2)-2.0d0*qixqk(2))
c     &                          - gfr(5)*rxqir(2)
c     &                          - gfr(7)*(rxqikr(2)+qkrxqir(2))
c                  ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)
c     &                          + gfr(4)*(dixqkr(3)+dkxqir(3)
c     &                                   +rxqidk(3)-2.0d0*qixqk(3))
c     &                          - gfr(5)*rxqir(3)
c     &                          - gfr(7)*(rxqikr(3)+qkrxqir(3))
c                  ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1)
c     &                          - gfr(4)*(dixqkr(1)+dkxqir(1)
c     &                                   +rxqkdi(1)-2.0d0*qixqk(1))
c     &                          - gfr(6)*rxqkr(1)
c     &                          - gfr(7)*(rxqkir(1)-qkrxqir(1))
c                  ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2)
c     &                          - gfr(4)*(dixqkr(2)+dkxqir(2)
c     &                                   +rxqkdi(2)-2.0d0*qixqk(2))
c     &                          - gfr(6)*rxqkr(2)
c     &                          - gfr(7)*(rxqkir(2)-qkrxqir(2))
c                  ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3)
c     &                          - gfr(4)*(dixqkr(3)+dkxqir(3)
c     &                                   +rxqkdi(3)-2.0d0*qixqk(3))
c     &                          - gfr(6)*rxqkr(3)
c     &                          - gfr(7)*(rxqkir(3)-qkrxqir(3))
c               end if
c
c     get the induced torque with screening
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               ttm2i(1) = -0.5d0*bn(1)*(dixuk(1)+dixukp(1))
     &                       + gti(2)*dixr(1) - gti(5)*rxqir(1)
     &                       + 0.5d0*gti(4)*(ukxqir(1)+rxqiuk(1)
     &                                      +ukxqirp(1)+rxqiukp(1))
               ttm2i(2) = -0.5d0*bn(1)*(dixuk(2)+dixukp(2))
     &                       + gti(2)*dixr(2) - gti(5)*rxqir(2)
     &                       + 0.5d0*gti(4)*(ukxqir(2)+rxqiuk(2)
     &                                      +ukxqirp(2)+rxqiukp(2))
               ttm2i(3) = -0.5d0*bn(1)*(dixuk(3)+dixukp(3))
     &                       + gti(2)*dixr(3) - gti(5)*rxqir(3)
     &                       + 0.5d0*gti(4)*(ukxqir(3)+rxqiuk(3)
     &                                      +ukxqirp(3)+rxqiukp(3))
               ttm3i(1) = -0.5d0*bn(1)*(dkxui(1)+dkxuip(1))
     &                       + gti(3)*dkxr(1) - gti(6)*rxqkr(1)
     &                       - 0.5d0*gti(4)*(uixqkr(1)+rxqkui(1)
     &                                      +uixqkrp(1)+rxqkuip(1))
               ttm3i(2) = -0.5d0*bn(1)*(dkxui(2)+dkxuip(2))
     &                       + gti(3)*dkxr(2) - gti(6)*rxqkr(2)
     &                       - 0.5d0*gti(4)*(uixqkr(2)+rxqkui(2)
     &                                       +uixqkrp(2)+rxqkuip(2))
               ttm3i(3) = -0.5d0*bn(1)*(dkxui(3)+dkxuip(3))
     &                       + gti(3)*dkxr(3) - gti(6)*rxqkr(3)
     &                       - 0.5d0*gti(4)*(uixqkr(3)+rxqkui(3)
     &                                      +uixqkrp(3)+rxqkuip(3))
c
c     get the induced torque without screening
c
               if (dorli) then
                  gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
                  gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
                  gtri(4) = gfri(4)
                  gtri(5) = gfri(5)
                  gtri(6) = gfri(6)
                  ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(1) - gtri(5)*rxqir(1)
     &                           + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &                             +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0
                  ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(2) - gtri(5)*rxqir(2)
     &                           + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &                             +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0
                  ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(3) - gtri(5)*rxqir(3)
     &                           + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &                             +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0
                  ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(1) - gtri(6)*rxqkr(1)
     &                           - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &                             +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0
                  ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(2) - gtri(6)*rxqkr(2)
     &                           - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &                             +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0
                  ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(3) - gtri(6)*rxqkr(3)
     &                           - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &                             +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0
               end if
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c     increment gradient due to force and torque on first site
c
c               demo1(1,ii) = demo1(1,ii) + ftm2(1)
c               demo1(2,ii) = demo1(2,ii) + ftm2(2)
c               demo1(3,ii) = demo1(3,ii) + ftm2(3)
               depo1(1,ii) = depo1(1,ii) + ftm2i(1)
               depo1(2,ii) = depo1(2,ii) + ftm2i(2)
               depo1(3,ii) = depo1(3,ii) + ftm2i(3)
c               call torque3 (
c     &   i,ttm2,ttm2i,frcxi,frcyi,frczi,demo1,depo1)
               call torque3_Polar (
     &   i,ttm2,ttm2i,frcxi,frcyi,frczi,depo1)

c
c     increment gradient due to force and torque on second site
c
c               demo2(1,kk) = demo2(1,kk) - ftm2(1)
c               demo2(2,kk) = demo2(2,kk) - ftm2(2)
c               demo2(3,kk) = demo2(3,kk) - ftm2(3)
               depo2(1,kk) = depo2(1,kk) - ftm2i(1)
               depo2(2,kk) = depo2(2,kk) - ftm2i(2)
               depo2(3,kk) = depo2(3,kk) - ftm2i(3)
c               call torque3 (
c     &    k,ttm3,ttm3i,frcxk,frcyk,frczk,demo2,depo2)
               call torque3_Polar(
     &    k,ttm3,ttm3i,frcxk,frcyk,frczk,depo2)
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
               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
               viro(1,1) = viro(1,1) + vxx
               viro(2,1) = viro(2,1) + vyx
               viro(3,1) = viro(3,1) + vzx
               viro(1,2) = viro(1,2) + vyx
               viro(2,2) = viro(2,2) + vyy
               viro(3,2) = viro(3,2) + vzy
               viro(1,3) = viro(1,3) + vzx
               viro(2,3) = viro(2,3) + vzy
               viro(3,3) = viro(3,3) + vzz
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
c      emtemp = emtemp + emo
      eptemp = eptemp + epo
c      eintra = eintrao
      do l1 = 1, npole3b
         i=pnum(l1)
         do j = 1, 3
ccc            demtemp(j,l1) = demtemp(j,l1) + demo1(j,i) + demo2(j,i)
ccc            deptemp(j,l1) = deptemp(j,l1) + depo1(j,i) + depo2(j,i)
c            demtemp(j,i) = demtemp(j,i) + demo1(j,i) + demo2(j,i)
            deptemp(j,i) = deptemp(j,i) + depo1(j,i) + depo2(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            virtemp(j,i) = virtemp(j,i) + viro(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
c      deallocate (demo1)
c      deallocate (demo2)
      deallocate (depo1)
      deallocate (depo2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  
c     ##                                                              ##
c     ##################################################################
c
      subroutine emrecip1_3b_Polartest(npole3b,pnum,uind,uinp,
     & eptemp,deptemp,virtemp,emtemp,demtemp)
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
      use polar, only: polarity, thole, pdamp
      use polpot
      use potent
      use virial
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
      integer npole3b,pnum(*),cnt,l1
      real*8 uind(3,*),eptemp,deptemp(3,*)
      real*8 virtemp(3,3),uinp(3,*)
      real*8 emtemp,demtemp(3,*),viremtemp(3,3)
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
      allocate (frc(3,n))
      allocate (trq(3,n))
      allocate (fuind(3,npole3b))
      allocate (fuinp(3,npole3b))
      allocate (cmp(10,npole3b))
      allocate (fmp(10,npole3b))
      allocate (fphi(20,npole3b))
      allocate (fphid(10,npole3b))
      allocate (fphip(10,npole3b))
      allocate (fphidp(20,npole3b))
      allocate (cphi(10,npole3b))
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
      !do i = 1, npole
      do l1 = 1, npole3b
         i=pnum(l1)
         cmp(1,l1) = rpole(1,i)
         cmp(2,l1) = rpole(2,i)
         cmp(3,l1) = rpole(3,i)
         cmp(4,l1) = rpole(4,i)
         cmp(5,l1) = rpole(5,i)
         cmp(6,l1) = rpole(9,i)
         cmp(7,l1) = rpole(13,i)
         cmp(8,l1) = 2.0d0 * rpole(6,i)
         cmp(9,l1) = 2.0d0 * rpole(7,i)
         cmp(10,l1) = 2.0d0 * rpole(10,i)
      end do
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     compute the arrays of B-spline coefficients
c
      if (.not. use_polar) then
         call bspline_fill
         call table_fill
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
      if (use_polar) then
c         do i = 1, npole
         do l1=1,npole3b
            i=pnum(l1)
            do j = 2, 4
               cmp(j,l1) = cmp(j,l1) + uinp(j-1,l1)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
c         do i = 1, npole
         do l1=1,npole3b
            i=pnum(l1)
            do j = 2, 4
               cmp(j,l1) = cmp(j,l1) + uind(j-1,l1) - uinp(j-1,l1)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         !do i = 1, npole
         do l1 = 1, npole3b
            do j = 2, 4
               cmp(j,l1) = cmp(j,l1) - uind(j-1,l1)
            end do
         end do
      else
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
      end if
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
            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vyx = vyx + h2*h1*vterm
            vzx = vzx + h3*h1*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vzy = vzy + h3*h2*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     assign just the induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
      if (use_polar .and. poltyp.eq.'DIRECT') then
         do i = 1, npole3b
            do j = 1, 10
               cmp(j,i) = 0.0d0
            end do
            do j = 2, 4
               cmp(j,i) = uinp(j-1,i)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
         do i = 1, npole3b
            do j = 2, 4
               cmp(j,i) = uind(j-1,i)
            end do
         end do
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
         do l1 = 1, npole3b
            i=pnum(l1)
            cmp(1,l1) = rpole(1,i)
            cmp(2,l1) = rpole(2,i)
            cmp(3,l1) = rpole(3,i)
            cmp(4,l1) = rpole(4,i)
            cmp(5,l1) = rpole(5,i)
            cmp(6,l1) = rpole(9,i)
            cmp(7,l1) = rpole(13,i)
            cmp(8,l1) = 2.0d0 * rpole(6,i)
            cmp(9,l1) = 2.0d0 * rpole(7,i)
            cmp(10,l1) = 2.0d0 * rpole(10,i)
         end do
c
c     make the scalar summation over reciprocal lattice
c
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
               struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                     + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
               eterm = 0.5d0 * electric * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx - h1*h1*vterm + eterm
               vyx = vyx - h2*h1*vterm
               vzx = vzx - h3*h1*vterm
               vyy = vyy - h2*h2*vterm + eterm
               vzy = vzy - h3*h2*vterm
               vzz = vzz - h3*h3*vterm + eterm
            end if
         end do
      end if


c         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
c         call grid_mpole3b_new (npole3b,pnum,fmp)
c         call fftfront
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
c                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
c               end do
c            end do
c         end do

c
c     Need to subtract permanent multipole contribution to reciprocal lattice
c
c      ntot = nfft1 * nfft2 * nfft3
c      pterm = (pi/aewald)**2
c      volterm = pi * volbox
c      nff = nfft1 * nfft2
c      nf1 = (nfft1+1) / 2
c      nf2 = (nfft2+1) / 2
c      nf3 = (nfft3+1) / 2
c      do i = 1, ntot-1
c         k3 = i/nff + 1
c         j = i - (k3-1)*nff
c         k2 = j/nfft1 + 1
c         k1 = j - (k2-1)*nfft1 + 1
c         m1 = k1 - 1
c         m2 = k2 - 1
c         m3 = k3 - 1
c         if (k1 .gt. nf1)  m1 = m1 - nfft1
c         if (k2 .gt. nf2)  m2 = m2 - nfft2
c         if (k3 .gt. nf3)  m3 = m3 - nfft3
c         r1 = dble(m1)
c         r2 = dble(m2)
c         r3 = dble(m3)
c         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c         hsq = h1*h1 + h2*h2 + h3*h3
c         term = -pterm * hsq
c         expterm = 0.0d0
c         if (term .gt. -50.0d0) then
c            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c            expterm = exp(term) / denom
c            if (.not. use_bounds) then
c               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
c            else if (octahedron) then
c               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
c            end if
c            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
c     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
c            eterm = 0.5d0 * electric * expterm * struc2
c            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
c            vxx = vxx - (h1*h1*vterm - eterm)
c            vyx = vyx - (h2*h1*vterm)
c            vzx = vzx - (h3*h1*vterm)
c            vyy = vyy - (h2*h2*vterm - eterm)
c            vzy = vzy - (h3*h2*vterm)
c            vzz = vzz - (h3*h3*vterm - eterm)
c         end if
cc         qfac(k1,k2,k3) = expterm
c      end do

c     p3rform deallocation of some local arrays
c
      deallocate (qgrip)
c
c     transform permanent multipoles without induced dipoles
c
      if (use_polar) then
         call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
         call grid_mpole3b_new (npole3b,pnum,fmp)
         call fftfront
      end if
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
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole3b_new (npole3b,pnum,fphi)
      do i = 1, npole3b
         do j = 1, 20
            fphi(j,i) = electric * fphi(j,i)
         end do
      end do
      call fphi_to_cphi3b (npole3b,pnum,fphi,cphi)
c
c     increment the permanent multipole energy and gradient
c
c      e = 0.0d0
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         f1 = 0.0d0
c         f2 = 0.0d0
c         f3 = 0.0d0
c         do k = 1, 10
c            e = e + fmp(k,l1)*fphi(k,l1)
c            f1 = f1 + fmp(k,l1)*fphi(deriv1(k),l1)
c            f2 = f2 + fmp(k,l1)*fphi(deriv2(k),l1)
c            f3 = f3 + fmp(k,l1)*fphi(deriv3(k),l1)
c         end do
c         f1 = dble(nfft1) * f1
c         f2 = dble(nfft2) * f2
c         f3 = dble(nfft3) * f3
c         frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
c         frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
c         frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
c      end do
c      e = 0.5d0 * e
c      emtemp = emtemp + e
c      do l1 = 1, npole3b
c         i=pnum(l1)
c      !   ii = ipole(i)
cc         demtemp(1,l1) = demtemp(1,l1) + frc(1,i)
cc         demtemp(2,l1) = demtemp(2,l1) + frc(2,i)
cc         demtemp(3,l1) = demtemp(3,l1) + frc(3,i)
c
c         demtemp(1,i) = demtemp(1,i) + frc(1,i)
c         demtemp(2,i) = demtemp(2,i) + frc(2,i)
c         demtemp(3,i) = demtemp(3,i) + frc(3,i)
c      end do
c
c     distribute torques into the permanent multipole gradient
c
c      do l1 = 1, npole3b
c         i=pnum(l1)
c         trq(1,i) = cmp(4,l1)*cphi(3,l1) - cmp(3,l1)*cphi(4,l1)
c     &                 + 2.0d0*(cmp(7,l1)-cmp(6,l1))*cphi(10,l1)
c     &                 + cmp(9,l1)*cphi(8,l1) + cmp(10,l1)*cphi(6,l1)
c     &                 - cmp(8,l1)*cphi(9,l1) - cmp(10,l1)*cphi(7,l1)
c         trq(2,i) = cmp(2,l1)*cphi(4,l1) - cmp(4,l1)*cphi(2,l1)
c     &                 + 2.0d0*(cmp(5,l1)-cmp(7,l1))*cphi(9,l1)
c     &                 + cmp(8,l1)*cphi(10,l1) + cmp(9,l1)*cphi(7,l1)
c     &                 - cmp(9,l1)*cphi(5,l1) - cmp(10,l1)*cphi(8,l1)
c         trq(3,i) = cmp(3,l1)*cphi(2,l1) - cmp(2,l1)*cphi(3,l1)
c     &                 + 2.0d0*(cmp(6,l1)-cmp(5,l1))*cphi(8,l1)
c     &                 + cmp(8,l1)*cphi(5,l1) + cmp(10,l1)*cphi(9,l1)
c     &                 - cmp(8,l1)*cphi(6,l1) - cmp(9,l1)*cphi(10,l1)
c      end do
c      do i = 1, n
c         frc(1,i) = 0.0d0
c         frc(2,i) = 0.0d0
c         frc(3,i) = 0.0d0
c      end do
c      call torque2_3bnew (npole3b,pnum,trq,frc)
c      do l1 = 1, npole3b
c         i=pnum(l1)
cc         demtemp(1,l1) = demtemp(1,l1) + frc(1,i)
cc         demtemp(2,l1) = demtemp(2,l1) + frc(2,i)
cc         demtemp(3,l1) = demtemp(3,l1) + frc(3,i)
c
c         demtemp(1,i) = demtemp(1,i) + frc(1,i)
c         demtemp(2,i) = demtemp(2,i) + frc(2,i)
c         demtemp(3,i) = demtemp(3,i) + frc(3,i)
c      end do
c
c     permanent multipole contribution to the internal virial
c
c      do i = 1, npole3b
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
      if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         do i = 1, npole3b
            do j = 1, 3
               fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
     &                          + a(j,3)*uind(3,i)
               fuinp(j,i) = a(j,1)*uinp(1,i) + a(j,2)*uinp(2,i)
     &                          + a(j,3)*uinp(3,i)
            end do
         end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
         call grid_uind3b_new (npole3b,pnum,fuind,fuinp)
         call fftfront
c
c     account for the zeroth grid point for a finite system
c
         if (.not. use_bounds) then
            expterm = 0.5d0 * pi / xbox
            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
            e = 0.5d0 * expterm * struc2
         end if
c
c     complete the transformation of the PME grid
c
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  term = qfac(i,j,k)
                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
               end do
            end do
         end do
c
c     perform 3-D FFT backward transform and get potential
c
         call fftback
         call fphi_uind3b_new (npole3b,pnum,fphid,fphip,fphidp)
         do i = 1, npole3b
            do j = 1, 10
               fphid(j,i) = electric * fphid(j,i)
               fphip(j,i) = electric * fphip(j,i)
            end do
            do j = 1, 20
               fphidp(j,i) = electric * fphidp(j,i)
            end do
         end do
c
c     increment the induced dipole energy and gradient
c
         e = 0.0d0
         do l1 = 1, npole3b
            i=pnum(l1)
            f1 = 0.0d0
            f2 = 0.0d0
            f3 = 0.0d0
            do k = 1, 3
               j1 = deriv1(k+1)
               j2 = deriv2(k+1)
               j3 = deriv3(k+1)
               e = e + fuind(k,l1)*fphi(k+1,l1)
               f1 = f1 + (fuind(k,l1)+fuinp(k,l1))*fphi(j1,l1)
     &                 + fuind(k,l1)*fphip(j1,l1)
     &                 + fuinp(k,l1)*fphid(j1,l1)
               f2 = f2 + (fuind(k,l1)+fuinp(k,l1))*fphi(j2,l1)
     &                 + fuind(k,l1)*fphip(j2,l1)
     &                 + fuinp(k,l1)*fphid(j2,l1)
               f3 = f3 + (fuind(k,l1)+fuinp(k,l1))*fphi(j3,l1)
     &                 + fuind(k,l1)*fphip(j3,l1)
     &                 + fuinp(k,l1)*fphid(j3,l1)
               if (poltyp .eq. 'DIRECT') then
                  f1 = f1 - fuind(k,l1)*fphip(j1,l1)
     &                    - fuinp(k,l1)*fphid(j1,l1)
                  f2 = f2 - fuind(k,l1)*fphip(j2,l1)
     &                    - fuinp(k,l1)*fphid(j2,l1)
                  f3 = f3 - fuind(k,l1)*fphip(j3,l1)
     &                    - fuinp(k,l1)*fphid(j3,l1)
               end if
            end do
            do k = 1, 10
               f1 = f1 + fmp(k,l1)*fphidp(deriv1(k),l1)
               f2 = f2 + fmp(k,l1)*fphidp(deriv2(k),l1)
               f3 = f3 + fmp(k,l1)*fphidp(deriv3(k),l1)
            end do
            f1 = 0.5d0 * dble(nfft1) * f1
            f2 = 0.5d0 * dble(nfft2) * f2
            f3 = 0.5d0 * dble(nfft3) * f3
            frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
            frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
            frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         end do
         e = 0.5d0 * e
         eptemp = eptemp + e
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
c            deptemp(1,l1) = deptemp(1,l1) + frc(1,i)
c            deptemp(2,l1) = deptemp(2,l1) + frc(2,i)
c            deptemp(3,l1) = deptemp(3,l1) + frc(3,i)

            deptemp(1,i) = deptemp(1,i) + frc(1,i)
            deptemp(2,i) = deptemp(2,i) + frc(2,i)
            deptemp(3,i) = deptemp(3,i) + frc(3,i)
         end do
c
c     set the potential to be the induced dipole average
c
         do i = 1, npole3b
            do k = 1, 10
               fphidp(k,i) = 0.5d0 * fphidp(k,i)
            end do
         end do
         call fphi_to_cphi3b (npole3b,pnum,fphidp,cphi)
c
c     distribute torques into the induced dipole gradient
c
         do l1 = 1, npole3b
            i=pnum(l1)
            trq(1,i) = cmp(4,l1)*cphi(3,l1) - cmp(3,l1)*cphi(4,l1)
     &                    + 2.0d0*(cmp(7,l1)-cmp(6,l1))*cphi(10,l1)
     &                    + cmp(9,l1)*cphi(8,l1) + cmp(10,l1)*cphi(6,l1)
     &                    - cmp(8,l1)*cphi(9,l1) - cmp(10,l1)*cphi(7,l1)
            trq(2,i) = cmp(2,l1)*cphi(4,l1) - cmp(4,l1)*cphi(2,l1)
     &                    + 2.0d0*(cmp(5,l1)-cmp(7,l1))*cphi(9,l1)
     &                    + cmp(8,l1)*cphi(10,l1) + cmp(9,l1)*cphi(7,l1)
     &                    - cmp(9,l1)*cphi(5,l1) - cmp(10,l1)*cphi(8,l1)
            trq(3,i) = cmp(3,l1)*cphi(2,l1) - cmp(2,l1)*cphi(3,l1)
     &                    + 2.0d0*(cmp(6,l1)-cmp(5,l1))*cphi(8,l1)
     &                    + cmp(8,l1)*cphi(5,l1) + cmp(10,l1)*cphi(9,l1)
     &                    - cmp(8,l1)*cphi(6,l1) - cmp(9,l1)*cphi(10,l1)
         end do
         do i = 1, n
            frc(1,i) = 0.0d0
            frc(2,i) = 0.0d0
            frc(3,i) = 0.0d0
         end do
         call torque2_3bnew (npole3b,pnum,trq,frc)
         do l1 = 1, npole3b
            i=pnum(l1)
c            deptemp(1,l1) = deptemp(1,l1) + frc(1,i)
c            deptemp(2,l1) = deptemp(2,l1) + frc(2,i)
c            deptemp(3,l1) = deptemp(3,l1) + frc(3,i)

            deptemp(1,i) = deptemp(1,i) + frc(1,i)
            deptemp(2,i) = deptemp(2,i) + frc(2,i)
            deptemp(3,i) = deptemp(3,i) + frc(3,i)
         end do
c
c     induced dipole contribution to the internal virial
c
         do i = 1, npole3b
            l1=i
            do j = 2, 4
               cphim(j) = 0.0d0
               cphid(j) = 0.0d0
               cphip(j) = 0.0d0
               do k = 2, 4
                  cphim(j) = cphim(j) + ftc(j,k)*fphi(k,i)
                  cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
                  cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
               end do
            end do
            vxx = vxx - cphi(2,i)*cmp(2,i)
     &                - 0.5d0*(cphim(2)*(uind(1,l1)+uinp(1,l1))
     &                        +cphid(2)*uinp(1,l1)+cphip(2)*uind(1,l1))
            vyx = vyx - 0.5d0*(cphi(2,i)*cmp(3,i)+cphi(3,i)*cmp(2,i))
     &                - 0.25d0*(cphim(2)*(uind(2,l1)+uinp(2,l1))
     &                         +cphim(3)*(uind(1,l1)+uinp(1,l1))
     &                         +cphid(2)*uinp(2,l1)+cphip(2)*uind(2,l1)
     &                         +cphid(3)*uinp(1,l1)+cphip(3)*uind(1,l1))
            vzx = vzx - 0.5d0*(cphi(2,i)*cmp(4,i)+cphi(4,i)*cmp(2,i))
     &                - 0.25d0*(cphim(2)*(uind(3,l1)+uinp(3,l1))
     &                         +cphim(4)*(uind(1,l1)+uinp(1,l1))
     &                         +cphid(2)*uinp(3,l1)+cphip(2)*uind(3,l1)
     &                         +cphid(4)*uinp(1,l1)+cphip(4)*uind(1,l1))
            vyy = vyy - cphi(3,i)*cmp(3,i)
     &                - 0.5d0*(cphim(3)*(uind(2,l1)+uinp(2,l1))
     &                        +cphid(3)*uinp(2,l1)+cphip(3)*uind(2,l1))
            vzy = vzy - 0.5d0*(cphi(3,i)*cmp(4,i)+cphi(4,i)*cmp(3,i))
     &                - 0.25d0*(cphim(3)*(uind(3,l1)+uinp(3,l1))
     &                         +cphim(4)*(uind(2,l1)+uinp(2,l1))
     &                         +cphid(3)*uinp(3,l1)+cphip(3)*uind(3,l1)
     &                         +cphid(4)*uinp(2,l1)+cphip(4)*uind(2,l1))
            vzz = vzz - cphi(4,i)*cmp(4,i)
     &                - 0.5d0*(cphim(4)*(uind(3,l1)+uinp(3,l1))
     &                        +cphid(4)*uinp(3,l1)+cphip(4)*uind(3,l1))
            vxx = vxx - 2.0d0*cmp(5,i)*cphi(5,i) - cmp(8,i)*cphi(8,i)
     &                - cmp(9,i)*cphi(9,i)
            vyx = vyx - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &                - 0.5d0*(cmp(8,i)*(cphi(6,i)+cphi(5,i))
     &                     +cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
            vzx = vzx - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &                - 0.5d0*(cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &                     +cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
            vyy = vyy - 2.0d0*cmp(6,i)*cphi(6,i) - cmp(8,i)*cphi(8,i)
     &                - cmp(10,i)*cphi(10,i)
            vzy = vzy - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &                - 0.5d0*(cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &                     +cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
            vzz = vzz - 2.0d0*cmp(7,i)*cphi(7,i) - cmp(9,i)*cphi(9,i)
     &                - cmp(10,i)*cphi(10,i)
            if (poltyp .eq. 'DIRECT') then
             vxx = vxx + 0.5d0*(cphid(2)*uinp(1,l1)+cphip(2)*uind(1,l1))
             vyx = vyx + 0.25d0*(cphid(2)*uinp(2,l1)+cphip(2)*uind(2,l1)
     &                         +cphid(3)*uinp(1,l1)+cphip(3)*uind(1,l1))
             vzx = vzx + 0.25d0*(cphid(2)*uinp(3,l1)+cphip(2)*uind(3,l1)
     &                         +cphid(4)*uinp(1,l1)+cphip(4)*uind(1,l1))
             vyy = vyy + 0.5d0*(cphid(3)*uinp(2,l1)+cphip(3)*uind(2,l1))
             vzy = vzy + 0.25d0*(cphid(3)*uinp(3,l1)+cphip(3)*uind(3,l1)
     &                         +cphid(4)*uinp(2,l1)+cphip(4)*uind(2,l1))
             vzz = vzz + 0.5d0*(cphid(4)*uinp(3,l1)+cphip(4)*uind(3,l1))
            end if
         end do
      end if
c
c     increment the internal virial tensor components
c
      virtemp(1,1) = virtemp(1,1) + vxx
      virtemp(2,1) = virtemp(2,1) + vyx
      virtemp(3,1) = virtemp(3,1) + vzx
      virtemp(1,2) = virtemp(1,2) + vyx
      virtemp(2,2) = virtemp(2,2) + vyy
      virtemp(3,2) = virtemp(3,2) + vzy
      virtemp(1,3) = virtemp(1,3) + vzx
      virtemp(2,3) = virtemp(2,3) + vzy
      virtemp(3,3) = virtemp(3,3) + vzz
c
c     perform deallocation of some local arrays
c
      deallocate (frc)
      deallocate (trq)
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
      deallocate (cphi)
      return
      end
