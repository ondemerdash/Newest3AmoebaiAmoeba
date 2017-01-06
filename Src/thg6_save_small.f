c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ereal1d  --  ewald real space derivs via list  ##
c     ##               Omar Demerdash June 2015                     ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1d" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c     and dipole polarizability
c
c
c      subroutine LoadBal_ereal1d_3b_Perm_sentlist2_totfield_dEtensorOmp2
c     & (emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet,
c     &  ntpair_dEtmp,dEindextmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp)
      subroutine thg6_save_small
     & (emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
      use sizes
      use atoms
      use boxes
      use chgpot
      use ewald
      use inter
      use math
      use mpole
      use polpot
      use couple
      use mplpot
      use neigh
      use limits
      use potent
      use paremneigh 
      use mpidat
      use polgrp
      use polar, only: polarity, thole, pdamp
      use openmp
      use dEtensor
      use shunt
      implicit none
      include 'mpif.h'
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 temp3,temp5,temp7
      real*8 fact1,fact2,fact3,fact4,fact5,fact6,fact7,fact8
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      integer sz,m
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
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
      real*8 ftm2(3)
      real*8 ftm2r(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 fdir(3),dixdk(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 gf(7)
      real*8 gfr(7)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: demi(:,:)
      real*8, allocatable :: demk(:,:)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      integer atomind,iter3
      real*8 demrealt(3,npole),viremrealt(3,3),emrealt
      logical dorl,dorli
      character*6 mode
      real*8 fieldnpolet(3,npole),fieldpnpolet(3,npole),bcn(3)
      real*8 fimd(3),fkmd(3),fimp(3),fkmp(3),rr5_field,rr7_field
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5,scale7    
      integer nlocal,maxlocal,maxlocal_allthread
      integer tid,toffset0
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8, allocatable :: dlocal1d(:,:) 
      real*8, allocatable :: dlocal2d(:,:)
      real*8, allocatable :: dlocal1p(:,:)
      real*8, allocatable :: dlocal2p(:,:)
      real*8 gfi(6),gfri_d(6),gfri_p(6)
      real*8 t1xxd,t1xyd,t1xzd,t1yxd,t1yyd,t1yzd,t1zxd,t1zyd,t1zzd
      real*8 t2xxd,t2xyd,t2xzd,t2yxd,t2yyd,t2yzd,t2zxd,t2zyd,t2zzd
      real*8 t1xxp,t1xyp,t1xzp,t1yxp,t1yyp,t1yzp,t1zxp,t1zyp,t1zzp
      real*8  t2xxp,t2xyp,t2xzp,t2yxp,t2yyp,t2yzp,t2zxp,t2zyp,t2zzp
c      real*8 tau1d(1),tau1d(2),tau1d(3),tau1d(4),tau1d(5),tau1d(6)
c      real*8 tau1d(7),tau1d(8),tau1d(9)
c      real*8 tau2d(1),tau2d(2),tau2d(3),tau2d(4),tau2d(5),tau2d(6)
c      real*8 tau2d(7),tau2d(8),tau2d(9)
c      real*8 tau1p(1),tau1p(2),tau1p(3),tau1p(4),tau1p(5),tau1p(6)
c      real*8 tau1p(7),tau1p(8),tau1p(9)
c      real*8 tau2p(1),tau2p(2),tau2p(3),tau2p(4),tau2p(5),tau2p(6)
c      real*8 tau2p(7),tau2p(8),tau2p(9)
      real*8 tau1d(9),tau1p(9),tau2d(9),tau2p(9)
      real*8, allocatable :: taulocal1d(:,:)
      real*8, allocatable :: taulocal2d(:,:)
      real*8, allocatable :: taulocal1p(:,:)
      real*8, allocatable :: taulocal2p(:,:)
      integer ik,ierr
      logical done_ik
      integer nlocaltau1,nlocaltau2
      integer, allocatable :: ilocaltau1(:,:)
      integer, allocatable :: ilocaltau2(:,:)
c      integer ntpair_dEtmp,dEindextmp(2,*)
c      real*8 dEd1tmp(9,*),dEd2tmp(9,*),dEp1tmp(9,*),dEp2tmp(9,*)
      integer toffset0tau1,toffset0tau2,maxlocaltau
      integer, allocatable :: toffsettau1(:)
      integer, allocatable :: toffsettau2(:)
      external erfc
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (demi(3,n))
      allocate (demk(3,n))
      allocate (fieldt(3,n))
      allocate (fieldtp(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do

      emrealt=0.0d0
      do i=1,npole
         do j=1,3
           demi(j,i) = 0.0d0
           demk(j,i) = 0.0d0
           fieldt(j,i) = 0.0d0
           fieldtp(j,i) = 0.0d0
         end do
      end do
      do i=1,3
         do j=1,3
           viremrealt(j,i)=0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
      nlocal = 0
      nlocaltau1 = 0
      nlocaltau2 = 0

      toffset0 = 0
      toffset0tau1 = 0
      toffset0tau2 = 0
      sz=last_emreal2(taskid)-start_emreal2(taskid)+1
      allocate (toffset(0:nthread-1))
      allocate (toffsettau1(0:nthread-1))
      allocate (toffsettau2(0:nthread-1))
c      print*,"taskid",taskid,"maxsiz_elst(taskid)=",maxsize_elst(taskid)
c      maxlocal = int(dble(sz)*dble(maxsize_elst(taskid))/dble(nthread))
      maxlocal = int(dble(sz) * dble(maxelst)/dble(nthread))
      maxlocaltau=3*maxlocal

       !print*,"THG6_sm"

c
c     set OpenMP directives for the major loop structure
c

c!$OMP PARALLEL default(shared) firstprivate(f) 
c!$OMP& private(i,j,k,ii,kk,kkk,e,bfac,damp,expdamp,
c!$OMP& pdi,pti,pgamma,scale3,scale5,scale7,temp3,temp5,temp7,
c!$OMP& dsc3,dsc5,dsc7,psc3,psc5,psc7,alsq2,alsq2n,
c!$OMP& exp2a,ralpha,gfd,gfdr,xr,yr,zr,xix,yix,zix,
c!$OMP& xiy,yiy,ziy,xiz,yiz,ziz,xkx,ykx,zkx,xky,yky,zky,
c!$OMP& xkz,ykz,zkz,r,r2,rr1,rr3,rr5,rr7,rr9,rr11,rr5_field,rr7_field,
c!$OMP& erl,erli,iax,iay,iaz,kax,kay,kaz,vxx,vyy,vzz,vyx,vzx,vzy,
c!$OMP& frcxi,frcyi,frczi,frcxk,frcyk,frczk,ci,di,qi,ck,dk,qk,
c!$OMP& ftm2,ftm2r,ttm2,ttm3,
c!$OMP& ttm2i,ttm3i,ttm2r,ttm3r,dixdk,
c!$OMP& qidk,qkdi,qir,qkr,qiqkr,qkqir,qixqk,rxqir,dixr,dkxr,
c!$OMP& dixqkr,dkxqir,rxqkr,qkrxqir,rxqikr,rxqkir,rxqidk,rxqkdi,
c!$OMP& ddsc3,ddsc5,ddsc7,bn,sc,gl,gf,
c!$OMP& gfr,dorl,dorli,iter3,bcn,fimd,fkmd,fimp,fkmp,
c!$OMP& ilocal,dlocal1d,dlocal2d,dlocal1p,dlocal2p,gfi,gfri_d,gfri_p,
c!$OMP& t1xxd,t1xyd,t1xzd,t1yxd,t1yyd,t1yzd,t1zxd,t1zyd,t1zzd,
c!$OMP& t2xxd,t2xyd,t2xzd,t2yxd,t2yyd,t2yzd,t2zxd,t2zyd,t2zzd,
c!$OMP& t1xxp,t1xyp,t1xzp,t1yxp,t1yyp,t1yzp,t1zxp,t1zyp,t1zzp,
c!$OMP& t2xxp,t2xyp,t2xzp,t2yxp,t2yyp,t2yzp,t2zxp,t2zyp,t2zzp,m,
c!$OMP& done_ik,ik,tau1p,tau1d,tau2p,tau2d,tid,
c!$OMP& tau1d(1),tau1d(2),tau1d(3),tau1d(4),tau1d(5),tau1d(6),
c!$OMP& tau1d(7),tau1d(8),tau1d(9),
c!$OMP& tau2d(1),tau2d(2),tau2d(3),tau2d(4),tau2d(5),tau2d(6),
c!$OMP& tau2d(7),tau2d(8),tau2d(9),
c!$OMP& tau1p(1),tau1p(2),tau1p(3),tau1p(4),tau1p(5),tau1p(6),
c!$OMP& tau1p(7),tau1p(8),tau1p(9),
c!$OMP& tau2p(1),tau2p(2),tau2p(3),tau2p(4),tau2p(5),tau2p(6),
c!$OMP& tau2p(7),tau2p(8),tau2p(9),nlocal)
c!$OMP& firstprivate(mscale,pscale,dscale)

!$OMP PARALLEL default(private) shared(demi,demk,fieldt,fieldtp,
!$OMP& viremrealt,emrealt,npole,electric,dielec,off2,toffset,
!$OMP& toffset0,dEindextmp,ntpair_dEtmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
!$OMP& nthread,n,x,y,z,rpole,ipole,pdamp,thole,n12,i12,maxlocal,
!$OMP& m2scale,p2scale,n13,m3scale,p3scale,i13,n14,i14,p4scale,
!$OMP& m4scale,np11,ip11,p41scale,n15,m5scale,p5scale,i15,
!$OMP& np12,ip12,d1scale,d2scale,np13,ip13,d3scale,np14,ip14,d4scale,
!$OMP& nelst_recv,elst_recv,xaxis,yaxis,zaxis,f,aewald,
!$OMP& start_emreal2,last_emreal2,taskid,tau1indextmp,tau2indextmp,
!$OMP& taud1tmp,taud2tmp,taup1tmp,taup2tmp,ntpair_tau1tmp,
!$OMP& ntpair_tau2tmp,toffsettau1,toffset0tau1,toffsettau2,
!$OMP& toffset0tau2,maxlocaltau)
!$OMP& firstprivate(pscale,dscale,mscale,nlocal,nlocaltau1,
!$OMP& nlocaltau2)

      allocate (ilocal(2,maxlocal))
      allocate (dlocal1d(9,maxlocal))
      allocate (dlocal2d(9,maxlocal))
      allocate (dlocal1p(9,maxlocal))
      allocate (dlocal2p(9,maxlocal))
      allocate (taulocal1d(9,maxlocaltau))
      allocate (taulocal2d(9,maxlocaltau))
      allocate (taulocal1p(9,maxlocaltau))
      allocate (taulocal2p(9,maxlocaltau))
      allocate (ilocaltau1(2,maxlocaltau))
      allocate (ilocaltau2(2,maxlocaltau))

         !do m=1,maxlocal
         !   do j=1,9
         !      dlocal1d(j,m)=0.0d0
         !      dlocal1p(j,m)=0.0d0
         !      dlocal2d(j,m)=0.0d0
         !      dlocal2p(j,m)=0.0d0
         !   end do
         !end do
!$OMP DO reduction(+:emrealt,viremrealt,demi,demk,fieldt,fieldtp)
!$OMP& schedule(guided)
      !do i=1,npole
      do i=start_emreal2(taskid),last_emreal2(taskid)
         iter3=i-start_emreal2(taskid)+1

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
        end do
        do j = 1, np12(ii)
          dscale(ip12(j,ii)) = d2scale
        end do
        do j = 1, np13(ii)
          dscale(ip13(j,ii)) = d3scale
        end do
        do j = 1, np14(ii)
          dscale(ip14(j,ii)) = d4scale
        end do
        !do kkk = 1, nelst(i)
        !  k = elst(kkk,i)
        do kkk = 1, nelst_recv(iter3)
           k = elst_recv(kkk,iter3)
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
            qk(1) = rpole(5,k) ! qkxx
            qk(2) = rpole(6,k) ! qkxy
            qk(3) = rpole(7,k) ! qkxz
            qk(4) = rpole(8,k) ! qkyx
            qk(5) = rpole(9,k) ! qkyy
            qk(6) = rpole(10,k) ! qkyz
            qk(7) = rpole(11,k) ! qkzx
            qk(8) = rpole(12,k) ! qkzy
            qk(9) = rpole(13,k) ! qkzz
            t1xxd=0.0d0
            t1xyd=0.0d0
            t1xzd=0.0d0
            t1yxd=0.0d0
            t1yyd=0.0d0
            t1yzd=0.0d0
            t1zxd=0.0d0
            t1zyd=0.0d0
            t1zzd=0.0d0
            t2xxd=0.0d0
            t2xyd=0.0d0
            t2xzd=0.0d0
            t2yxd=0.0d0
            t2yyd=0.0d0
            t2yzd=0.0d0
            t2zxd=0.0d0
            t2zyd=0.0d0
            t2zzd=0.0d0
            t1xxp=0.0d0
            t1xyp=0.0d0
            t1xzp=0.0d0
            t1yxp=0.0d0
            t1yyp=0.0d0
            t1yzp=0.0d0
            t1zxp=0.0d0
            t1zyp=0.0d0
            t1zzp=0.0d0
            t2xxp=0.0d0
            t2xyp=0.0d0
            t2xzp=0.0d0
            t2yxp=0.0d0
            t2yyp=0.0d0
            t2yzp=0.0d0
            t2zxp=0.0d0
            t2zyp=0.0d0
            t2zzp=0.0d0
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
            rr5_field = rr3 / r2
            rr7_field = rr5_field / r2
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
                scale7 = 1.0d0-(1.0d0-damp+0.6d0*damp**2)*expdamp
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
c
c     construct necessary auxiliary vectors
c
            dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
            dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
            dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
            dixr(1) = di(2)*zr - di(3)*yr
            dixr(2) = di(3)*xr - di(1)*zr
            dixr(3) = di(1)*yr - di(2)*xr
            dkxr(1) = dk(2)*zr - dk(3)*yr
            dkxr(2) = dk(3)*xr - dk(1)*zr
            dkxr(3) = dk(1)*yr - dk(2)*xr
            qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr ! qix
            qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr ! qiy
            qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr ! qiz
            qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr ! qkx
            qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr ! qky
            qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr ! qkz
            qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
            qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
            qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
            qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
            qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
            qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
            qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)-
     &                 qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
            qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)-
     &                 qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
            qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)-
     &                 qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
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
            dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
            dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
            dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
            dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
            dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
            dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
            rxqidk(1) = yr*qidk(3) - zr*qidk(2)
            rxqidk(2) = zr*qidk(1) - xr*qidk(3)
            rxqidk(3) = xr*qidk(2) - yr*qidk(1)
            rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
            rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
            rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c
c     calculate the scalar products for permanent components
c
            sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
            sc(3) = di(1)*xr + di(2)*yr + di(3)*zr ! dir
            sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr ! dkr
            sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr ! qir
            sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr ! qkr
            sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
            sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
            sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
            sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)+
     &               qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)+
     &               qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)


            bcn(1) = bn(1)-(1.0d0-scale3*dscale(kk))*rr3
            bcn(2) = bn(2)-3.0d0*(1.0d0-scale5*dscale(kk))*rr5_field
            bcn(3) = bn(3)-15.0d0*(1.0d0-scale7*dscale(kk))*rr7_field
            fimd(1) = -xr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &                     - bcn(1)*dk(1) + 2.0d0*bcn(2)*qkr(1)
            fimd(2) = -yr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &                     - bcn(1)*dk(2) + 2.0d0*bcn(2)*qkr(2)
            fimd(3) = -zr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &                     - bcn(1)*dk(3) + 2.0d0*bcn(2)*qkr(3)
            fkmd(1) = xr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &                     - bcn(1)*di(1) - 2.0d0*bcn(2)*qir(1)
            fkmd(2) = yr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &                     - bcn(1)*di(2) - 2.0d0*bcn(2)*qir(2)
            fkmd(3) = zr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &                     - bcn(1)*di(3) - 2.0d0*bcn(2)*qir(3)
            bcn(1) = bn(1)-(1.0d0-scale3*pscale(kk))*rr3
            bcn(2) = bn(2)-3.0d0*(1.0d0-scale5*pscale(kk))*rr5_field
            bcn(3) = bn(3)-15.0d0*(1.0d0-scale7*pscale(kk))*rr7_field
            fimp(1) = -xr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &                     - bcn(1)*dk(1) + 2.0d0*bcn(2)*qkr(1)
            fimp(2) = -yr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &                     - bcn(1)*dk(2) + 2.0d0*bcn(2)*qkr(2)
            fimp(3) = -zr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &                     - bcn(1)*dk(3) + 2.0d0*bcn(2)*qkr(3)
            fkmp(1) = xr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &                     - bcn(1)*di(1) - 2.0d0*bcn(2)*qir(1)
            fkmp(2) = yr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &                     - bcn(1)*di(2) - 2.0d0*bcn(2)*qir(2)
            fkmp(3) = zr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &                     - bcn(1)*di(3) - 2.0d0*bcn(2)*qir(3)
            do j = 1, 3
              fieldt(j,i) = fieldt(j,i) + fimd(j)
              fieldt(j,k) = fieldt(j,k) + fkmd(j)
              fieldtp(j,i) = fieldtp(j,i) + fimp(j)
              fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
            end do
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
c     compute the energy contributions for this interaction
c
            e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)

c
c     get the real energy without any screening function
c
            erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
            erl = erl * (1.0d0-mscale(kk))
            e = e - erl
            e = f * e
            emrealt = emrealt + e
c
c     set flags to compute components without screening
c
            dorl = .false.
            if (mscale(kk) .ne. 1.0d0)  dorl = .true.
c
c     zero out force and torque components without screening
c
            do j = 1, 3
              ftm2r(j) = 0.0d0
              ttm2r(j) = 0.0d0
              ttm3r(j) = 0.0d0
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
            ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)+
     &                gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)+
     &                gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
            ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)+
     &                gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)+
     &                gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
            ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)+
     &                gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)+
     &                gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
            if (dorl) then
              gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))+
     &                 rr7*(gl(2)+gl(7)+gl(8))+
     &                 rr9*(gl(3)+gl(5)) + rr11*gl(4)
              gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
              gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
              gfr(4) = 2.0d0 * rr5
              gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
              gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
              gfr(7) = 4.0d0 * rr7
              ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)+
     &                   gfr(4)*(qkdi(1)-qidk(1))+
     &                   gfr(5)*qir(1) + gfr(6)*qkr(1)+
     &                   gfr(7)*(qiqkr(1)+qkqir(1))
              ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)+
     &                   gfr(4)*(qkdi(2)-qidk(2))+
     &                   gfr(5)*qir(2) + gfr(6)*qkr(2)+
     &                   gfr(7)*(qiqkr(2)+qkqir(2))
              ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)+
     &                   gfr(4)*(qkdi(3)-qidk(3))+
     &                   gfr(5)*qir(3) + gfr(6)*qkr(3)+
     &                   gfr(7)*(qiqkr(3)+qkqir(3))
            end if
c
c     get the induced force with screening
c
C 1ST TERM W/ SCREENING
c
            fact1=0.5d0*bn(2)
            fact2=0.5d0*bn(3)
            fact3=fact1*ck-fact2*sc(4)+0.5d0*bn(4)*sc(6)
            gfi(2) = 0.5d0*(-ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3))
            t2xxd = t2xxd + xr*(xr*fact3+fact1*dk(1)-bn(3)*qkr(1))+
     &              (fact1*dk(1)-bn(3)*qkr(1))*xr+bn(2)*qk(1)+gfi(2)
            t2xyd = t2xyd + xr*(yr*fact3+fact1*dk(2)-bn(3)*qkr(2))+
     &                      (fact1*dk(1)-bn(3)*qkr(1))*yr+bn(2)*qk(4)
            t2xzd = t2xzd + xr*(zr*fact3+fact1*dk(3)-bn(3)*qkr(3))+
     &                      (fact1*dk(1)-bn(3)*qkr(1))*zr+bn(2)*qk(7)
            t2yxd = t2yxd + yr*(xr*fact3+fact1*dk(1)-bn(3)*qkr(1))+
     &                      (fact1*dk(2)-bn(3)*qkr(2))*xr+bn(2)*qk(2)
            t2yyd = t2yyd + yr*(yr*fact3+fact1*dk(2)-bn(3)*qkr(2))+
     &              (fact1*dk(2)-bn(3)*qkr(2))*yr+bn(2)*qk(5)+gfi(2)
            t2yzd = t2yzd + yr*(zr*fact3+fact1*dk(3)-bn(3)*qkr(3))+
     &                      (fact1*dk(2)-bn(3)*qkr(2))*zr+bn(2)*qk(8)
            t2zxd = t2zxd + zr*(xr*fact3+fact1*dk(1)-bn(3)*qkr(1))+
     &                      (fact1*dk(3)-bn(3)*qkr(3))*xr+bn(2)*qk(3)
            t2zyd = t2zyd + zr*(yr*fact3+fact1*dk(2)-bn(3)*qkr(2))+
     &                      (fact1*dk(3)-bn(3)*qkr(3))*yr+bn(2)*qk(6)
            t2zzd = t2zzd + zr*(zr*fact3+fact1*dk(3)-bn(3)*qkr(3))+
     &              (fact1*dk(3)-bn(3)*qkr(3))*zr+bn(2)*qk(9)+gfi(2)

            fact3=fact1*ci+fact2*sc(3)+0.5d0*bn(4)*sc(5)
            gfi(3) = 0.5d0*(ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3))
            t1xxd = t1xxd - xr*(xr*fact3-fact1*di(1)-bn(3)*qir(1))+
     &              (fact1*di(1)+bn(3)*qir(1))*xr-bn(2)*qi(1)+gfi(3)
            t1xyd = t1xyd - xr*(yr*fact3-fact1*di(2)-bn(3)*qir(2))+
     &                      (fact1*di(1)+bn(3)*qir(1))*yr-bn(2)*qi(4)
            t1xzd = t1xzd - xr*(zr*fact3-fact1*di(3)-bn(3)*qir(3))+
     &                      (fact1*di(1)+bn(3)*qir(1))*zr-bn(2)*qi(7)
            t1yxd = t1yxd - yr*(xr*fact3-fact1*di(1)-bn(3)*qir(1))+
     &                      (fact1*di(2)+bn(3)*qir(2))*xr-bn(2)*qi(2)
            t1yyd = t1yyd - yr*(yr*fact3-fact1*di(2)-bn(3)*qir(2))+
     &              (fact1*di(2)+bn(3)*qir(2))*yr-bn(2)*qi(5)+gfi(3)
            t1yzd = t1yzd - yr*(zr*fact3-fact1*di(3)-bn(3)*qir(3))+
     &                      (fact1*di(2)+bn(3)*qir(2))*zr-bn(2)*qi(8)
            t1zxd = t1zxd - zr*(xr*fact3-fact1*di(1)-bn(3)*qir(1))+
     &                      (fact1*di(3)+bn(3)*qir(3))*xr-bn(2)*qi(3)
            t1zyd = t1zyd - zr*(yr*fact3-fact1*di(2)-bn(3)*qir(2))+
     &                      (fact1*di(3)+bn(3)*qir(3))*yr-bn(2)*qi(6)
            t1zzd = t1zzd - zr*(zr*fact3-fact1*di(3)-bn(3)*qir(3))+
     &              (fact1*di(3)+bn(3)*qir(3))*zr-bn(2)*qi(9)+gfi(3)

            t2xxp = t2xxd
            t2xyp = t2xyd
            t2xzp = t2xzd
            t2yxp = t2yxd
            t2yyp = t2yyd
            t2yzp = t2yzd
            t2zxp = t2zxd
            t2zyp = t2zyd
            t2zzp = t2zzd

            t1xxp = t1xxd
            t1xyp = t1xyd
            t1xzp = t1xzd
            t1yxp = t1yxd
            t1yyp = t1yyd
            t1yzp = t1yzd
            t1zxp = t1zxd
            t1zyp = t1zyd
            t1zzp = t1zzd
c
c     get the p/dscaled tensor terms without screening
c
c
C 1ST TERM W/OUT SCREENING
c
            fact1=-0.5d0*rr5*dsc3         
            fact2=-0.5d0*rr5*psc3         
            fact3=0.5d0*rr7*dsc5*sc(4)
            fact4=0.5d0*rr7*psc5*sc(4)
            fact5=-0.5d0*rr9*dsc7*sc(6)
            fact6=-0.5d0*rr9*psc7*sc(6)
            fact7=-rr7*dsc5
            fact8=-rr7*psc5
            t2xxd=t2xxd+xr*((fact1*ck+fact5+fact3)*xr+fact1*dk(1)-
     &                       fact7*qkr(1))
            t2xxp=t2xxp+xr*((fact2*ck+fact6+fact4)*xr+fact2*dk(1)-
     &                       fact8*qkr(1))
            t2xyd=t2xyd+xr*((fact1*ck+fact5+fact3)*yr+fact1*dk(2)-
     &                       fact7*qkr(2))
            t2xyp=t2xyp+xr*((fact2*ck+fact6+fact4)*yr+fact2*dk(2)-
     &                       fact8*qkr(2))
            t2xzd=t2xzd+xr*((fact1*ck+fact5+fact3)*zr+fact1*dk(3)-
     &                       fact7*qkr(3))
            t2xzp=t2xzp+xr*((fact2*ck+fact6+fact4)*zr+fact2*dk(3)-
     &                       fact8*qkr(3))

            t2yxd=t2yxd+yr*((fact1*ck+fact5+fact3)*xr+fact1*dk(1)-
     &                       fact7*qkr(1))
            t2yxp=t2yxp+yr*((fact2*ck+fact6+fact4)*xr+fact2*dk(1)-
     &                       fact8*qkr(1))
            t2yyd=t2yyd+yr*((fact1*ck+fact5+fact3)*yr+fact1*dk(2)-
     &                       fact7*qkr(2))
            t2yyp=t2yyp+yr*((fact2*ck+fact6+fact4)*yr+fact2*dk(2)-
     &                       fact8*qkr(2))
            t2yzd=t2yzd+yr*((fact1*ck+fact5+fact3)*zr+fact1*dk(3)-
     &                       fact7*qkr(3))
            t2yzp=t2yzp+yr*((fact2*ck+fact6+fact4)*zr+fact2*dk(3)-
     &                       fact8*qkr(3))

            t2zxd=t2zxd+zr*((fact1*ck+fact5+fact3)*xr+fact1*dk(1)-
     &                       fact7*qkr(1))
            t2zxp=t2zxp+zr*((fact2*ck+fact6+fact4)*xr+fact2*dk(1)-
     &                       fact8*qkr(1))
            t2zyd=t2zyd+zr*((fact1*ck+fact5+fact3)*yr+fact1*dk(2)-
     &                       fact7*qkr(2))
            t2zyp=t2zyp+zr*((fact2*ck+fact6+fact4)*yr+fact2*dk(2)-
     &                       fact8*qkr(2))
            t2zzd=t2zzd+zr*((fact1*ck+fact5+fact3)*zr+fact1*dk(3)-
     &                       fact7*qkr(3))
            t2zzp=t2zzp+zr*((fact2*ck+fact6+fact4)*zr+fact2*dk(3)-
     &                       fact8*qkr(3))

            fact3=0.5d0*rr7*dsc5*sc(3)
            fact4=0.5d0*rr7*psc5*sc(3)
            fact5=-0.5d0*rr9*dsc7*sc(5)
            fact6=-0.5d0*rr9*psc7*sc(5)
            fact7=-rr7*dsc5
            fact8=-rr7*psc5
            t1xxd=t1xxd-xr*((fact1*ci+fact5-fact3)*xr-fact1*di(1)-
     &                       fact7*qir(1))
            t1xxp=t1xxp-xr*((fact2*ci+fact6-fact4)*xr-fact2*di(1)-
     &                       fact8*qir(1))
            t1xyd=t1xyd-xr*((fact1*ci+fact5-fact3)*yr-fact1*di(2)-
     &                       fact7*qir(2))
            t1xyp=t1xyp-xr*((fact2*ci+fact6-fact4)*yr-fact2*di(2)-
     &                       fact8*qir(2))
            t1xzd=t1xzd-xr*((fact1*ci+fact5-fact3)*zr-fact1*di(3)-
     &                       fact7*qir(3))
            t1xzp=t1xzp-xr*((fact2*ci+fact6-fact4)*zr-fact2*di(3)-
     &                       fact8*qir(3))

            t1yxd=t1yxd-yr*((fact1*ci+fact5-fact3)*xr-fact1*di(1)-
     &                       fact7*qir(1))
            t1yxp=t1yxp-yr*((fact2*ci+fact6-fact4)*xr-fact2*di(1)-
     &                       fact8*qir(1))
            t1yyd=t1yyd-yr*((fact1*ci+fact5-fact3)*yr-fact1*di(2)-
     &                       fact7*qir(2))
            t1yyp=t1yyp-yr*((fact2*ci+fact6-fact4)*yr-fact2*di(2)-
     &                       fact8*qir(2))
            t1yzd=t1yzd-yr*((fact1*ci+fact5-fact3)*zr-fact1*di(3)-
     &                       fact7*qir(3))
            t1yzp=t1yzp-yr*((fact2*ci+fact6-fact4)*zr-fact2*di(3)-
     &                       fact8*qir(3))

            t1zxd=t1zxd-zr*((fact1*ci+fact5-fact3)*xr-fact1*di(1)-
     &                       fact7*qir(1))
            t1zxp=t1zxp-zr*((fact2*ci+fact6-fact4)*xr-fact2*di(1)-
     &                       fact8*qir(1))
            t1zyd=t1zyd-zr*((fact1*ci+fact5-fact3)*yr-fact1*di(2)-
     &                       fact7*qir(2))
            t1zyp=t1zyp-zr*((fact2*ci+fact6-fact4)*yr-fact2*di(2)-
     &                       fact8*qir(2))
            t1zzd=t1zzd-zr*((fact1*ci+fact5-fact3)*zr-fact1*di(3)-
     &                       fact7*qir(3))
            t1zzp=t1zzp-zr*((fact2*ci+fact6-fact4)*zr-fact2*di(3)-
     &                       fact8*qir(3))
c
C 4TH TERM, 7TH W/OUT SCREENING; tensor term for temp5
c
            fact1=-0.5d0*rr5*dsc5
            fact2=-0.5d0*rr5*psc5
            fact3=-rr7*dsc7
            fact4=-rr7*psc7
            fact5=0.5d0*rr7*dscale(kk)*sc(5)*ddsc7(1)
            fact6=0.5d0*rr7*pscale(kk)*sc(5)*ddsc7(1)
            fact7=-rr5*dsc5
            fact8=-rr5*psc5
            t1xxd=t1xxd+xr*(fact1*di(1)+fact3*qir(1)+fact5)-
     &                      fact7*qi(1)
            t1xxp=t1xxp+xr*(fact2*di(1)+fact4*qir(1)+fact6)-
     &                      fact8*qi(1)
            t1xyd=t1xyd+yr*(fact1*di(1)+fact3*qir(1)+fact5)-
     &                      fact7*qi(4)
            t1xyp=t1xyp+yr*(fact2*di(1)+fact4*qir(1)+fact6)-
     &                      fact8*qi(4)
            t1xzd=t1xzd+zr*(fact1*di(1)+fact3*qir(1)+fact5)-
     &                      fact7*qi(7)
            t1xzp=t1xzp+zr*(fact2*di(1)+fact4*qir(1)+fact6)-
     &                      fact8*qi(7)

            fact5=0.5d0*rr7*dscale(kk)*sc(5)*ddsc7(2)
            fact6=0.5d0*rr7*pscale(kk)*sc(5)*ddsc7(2)
            t1yxd=t1yxd+xr*(fact1*di(2)+fact3*qir(2)+fact5)-
     &                      fact7*qi(2)
            t1yxp=t1yxp+xr*(fact2*di(2)+fact4*qir(2)+fact6)-
     &                      fact8*qi(2)
            t1yyd=t1yyd+yr*(fact1*di(2)+fact3*qir(2)+fact5)-
     &                      fact7*qi(5)
            t1yyp=t1yyp+yr*(fact2*di(2)+fact4*qir(2)+fact6)-
     &                      fact8*qi(5)
            t1yzd=t1yzd+zr*(fact1*di(2)+fact3*qir(2)+fact5)-
     &                      fact7*qi(8)
            t1yzp=t1yzp+zr*(fact2*di(2)+fact4*qir(2)+fact6)-
     &                      fact8*qi(8)

            fact5=0.5d0*rr7*dscale(kk)*sc(5)*ddsc7(3)
            fact6=0.5d0*rr7*pscale(kk)*sc(5)*ddsc7(3)
            t1zxd=t1zxd+xr*(fact1*di(3)+fact3*qir(3)+fact5)-
     &                      fact7*qi(3)
            t1zxp=t1zxp+xr*(fact2*di(3)+fact4*qir(3)+fact6)-
     &                      fact8*qi(3)
            t1zyd=t1zyd+yr*(fact1*di(3)+fact3*qir(3)+fact5)-
     &                      fact7*qi(6)
            t1zyp=t1zyp+yr*(fact2*di(3)+fact4*qir(3)+fact6)-
     &                      fact8*qi(6)
            t1zzd=t1zzd+zr*(fact1*di(3)+fact3*qir(3)+fact5)-
     &                      fact7*qi(9)
            t1zzp=t1zzp+zr*(fact2*di(3)+fact4*qir(3)+fact6)-
     &                      fact8*qi(9)
c
C 5TH TERM W/OUT SCREENING
c
            fact5=-0.5d0*rr7*dscale(kk)*sc(6)*ddsc7(1)
            fact6=-0.5d0*rr7*pscale(kk)*sc(6)*ddsc7(1)
            t2xxd=t2xxd+xr*(fact1*dk(1)-fact3*qkr(1)+fact5)+
     &                      fact7*qk(1)
            t2xxp=t2xxp+xr*(fact2*dk(1)-fact4*qkr(1)+fact6)+
     &                      fact8*qk(1)
            t2xyd=t2xyd+yr*(fact1*dk(1)-fact3*qkr(1)+fact5)+
     &                      fact7*qk(4)
            t2xyp=t2xyp+yr*(fact2*dk(1)-fact4*qkr(1)+fact6)+
     &                      fact8*qk(4)
            t2xzd=t2xzd+zr*(fact1*dk(1)-fact3*qkr(1)+fact5)+
     &                      fact7*qk(7)
            t2xzp=t2xzp+zr*(fact2*dk(1)-fact4*qkr(1)+fact6)+
     &                      fact8*qk(7)

            fact5=-0.5d0*rr7*dscale(kk)*sc(6)*ddsc7(2)
            fact6=-0.5d0*rr7*pscale(kk)*sc(6)*ddsc7(2)
            t2yxd=t2yxd+xr*(fact1*dk(2)-fact3*qkr(2)+fact5)+
     &                      fact7*qk(2)
            t2yxp=t2yxp+xr*(fact2*dk(2)-fact4*qkr(2)+fact6)+
     &                      fact8*qk(2)
            t2yyd=t2yyd+yr*(fact1*dk(2)-fact3*qkr(2)+fact5)+
     &                      fact7*qk(5)
            t2yyp=t2yyp+yr*(fact2*dk(2)-fact4*qkr(2)+fact6)+
     &                      fact8*qk(5)
            t2yzd=t2yzd+zr*(fact1*dk(2)-fact3*qkr(2)+fact5)+
     &                      fact7*qk(8)
            t2yzp=t2yzp+zr*(fact2*dk(2)-fact4*qkr(2)+fact6)+
     &                      fact8*qk(8)

            fact5=-0.5d0*rr7*dscale(kk)*sc(6)*ddsc7(3)
            fact6=-0.5d0*rr7*pscale(kk)*sc(6)*ddsc7(3)
            t2zxd=t2zxd+xr*(fact1*dk(3)-fact3*qkr(3)+fact5)+
     &                      fact7*qk(3)
            t2zxp=t2zxp+xr*(fact2*dk(3)-fact4*qkr(3)+fact6)+
     &                      fact8*qk(3)
            t2zyd=t2zyd+yr*(fact1*dk(3)-fact3*qkr(3)+fact5)+
     &                      fact7*qk(6)
            t2zyp=t2zyp+yr*(fact2*dk(3)-fact4*qkr(3)+fact6)+
     &                      fact8*qk(6)
            t2zzd=t2zzd+zr*(fact1*dk(3)-fact3*qkr(3)+fact5)+
     &                      fact7*qk(9)
            t2zzp=t2zzp+zr*(fact2*dk(3)-fact4*qkr(3)+fact6)+
     &                      fact8*qk(9)
c
c
C 6TH TERM W/OUT SCREENING
C  TENSOR TERMS CORRESPONDING TO TEMP3 and TEMP5 FOR 'fridmp(j)' terms
c
            fact1=-0.5d0*rr3*dscale(kk)*ddsc3(1)
            fact2=-0.5d0*rr3*pscale(kk)*ddsc3(1)
            fact3=rr5*dscale(kk)*ddsc5(1)
            fact4=rr5*pscale(kk)*ddsc5(1)
            fact5=0.5d0*sc(3)
            fact6=0.5d0*sc(4)
            t2xxd=t2xxd+(fact1*ck+fact3*fact6)*xr+fact1*dk(1)+
     &                   fact3*qkr(1)
            t2xxp=t2xxp+(fact2*ck+fact4*fact6)*xr+fact2*dk(1)+
     &                   fact4*qkr(1)
            t2xyd=t2xyd+(fact1*ck+fact3*fact6)*yr+fact1*dk(2)+
     &                   fact3*qkr(2)
            t2xyp=t2xyp+(fact2*ck+fact4*fact6)*yr+fact2*dk(2)+
     &                   fact4*qkr(2)
            t2xzd=t2xzd+(fact1*ck+fact3*fact6)*zr+fact1*dk(3)+
     &                   fact3*qkr(3)
            t2xzp=t2xzp+(fact2*ck+fact4*fact6)*zr+fact2*dk(3)+
     &                   fact4*qkr(3)
            t1xxd=t1xxd-(fact1*ci-fact3*fact5)*xr+fact1*di(1)-
     &                   fact3*qir(1)
            t1xxp=t1xxp-(fact2*ci-fact4*fact5)*xr+fact2*di(1)-
     &                   fact4*qir(1)
            t1xyd=t1xyd-(fact1*ci-fact3*fact5)*yr+fact1*di(2)-
     &                   fact3*qir(2)
            t1xyp=t1xyp-(fact2*ci-fact4*fact5)*yr+fact2*di(2)-
     &                   fact4*qir(2)
            t1xzd=t1xzd-(fact1*ci-fact3*fact5)*zr+fact1*di(3)-
     &                   fact3*qir(3)
            t1xzp=t1xzp-(fact2*ci-fact4*fact5)*zr+fact2*di(3)-
     &                   fact4*qir(3)

            fact1=-0.5d0*rr3*dscale(kk)*ddsc3(2)
            fact2=-0.5d0*rr3*pscale(kk)*ddsc3(2)
            fact3=rr5*dscale(kk)*ddsc5(2)
            fact4=rr5*pscale(kk)*ddsc5(2)
            t2yxd=t2yxd+(fact1*ck+fact3*fact6)*xr+fact1*dk(1)+
     &                   fact3*qkr(1)
            t2yxp=t2yxp+(fact2*ck+fact4*fact6)*xr+fact2*dk(1)+
     &                   fact4*qkr(1)
            t2yyd=t2yyd+(fact1*ck+fact3*fact6)*yr+fact1*dk(2)+
     &                   fact3*qkr(2)
            t2yyp=t2yyp+(fact2*ck+fact4*fact6)*yr+fact2*dk(2)+
     &                   fact4*qkr(2)
            t2yzd=t2yzd+(fact1*ck+fact3*fact6)*zr+fact1*dk(3)+
     &                   fact3*qkr(3)
            t2yzp=t2yzp+(fact2*ck+fact4*fact6)*zr+fact2*dk(3)+
     &                   fact4*qkr(3)
            t1yxd=t1yxd-(fact1*ci-fact3*fact5)*xr+fact1*di(1)-
     &                   fact3*qir(1)
            t1yxp=t1yxp-(fact2*ci-fact4*fact5)*xr+fact2*di(1)-
     &                   fact4*qir(1)
            t1yyd=t1yyd-(fact1*ci-fact3*fact5)*yr+fact1*di(2)-
     &                   fact3*qir(2)
            t1yyp=t1yyp-(fact2*ci-fact4*fact5)*yr+fact2*di(2)-
     &                   fact4*qir(2)
            t1yzd=t1yzd-(fact1*ci-fact3*fact5)*zr+fact1*di(3)-
     &                   fact3*qir(3)
            t1yzp=t1yzp-(fact2*ci-fact4*fact5)*zr+fact2*di(3)-
     &                   fact4*qir(3)
            fact1=-0.5d0*rr3*dscale(kk)*ddsc3(3)
            fact2=-0.5d0*rr3*pscale(kk)*ddsc3(3)
            fact3=rr5*dscale(kk)*ddsc5(3)
            fact4=rr5*pscale(kk)*ddsc5(3)
            t2zxd=t2zxd+(fact1*ck+fact3*fact6)*xr+fact1*dk(1)+
     &                   fact3*qkr(1)
            t2zxp=t2zxp+(fact2*ck+fact4*fact6)*xr+fact2*dk(1)+
     &                   fact4*qkr(1)
            t2zyd=t2zyd+(fact1*ck+fact3*fact6)*yr+fact1*dk(2)+
     &                   fact3*qkr(2)
            t2zyp=t2zyp+(fact2*ck+fact4*fact6)*yr+fact2*dk(2)+
     &                   fact4*qkr(2)
            t2zzd=t2zzd+(fact1*ck+fact3*fact6)*zr+fact1*dk(3)+
     &                   fact3*qkr(3)
            t2zzp=t2zzp+(fact2*ck+fact4*fact6)*zr+fact2*dk(3)+
     &                   fact4*qkr(3)
            t1zxd=t1zxd-(fact1*ci-fact3*fact5)*xr+fact1*di(1)-
     &                   fact3*qir(1)
            t1zxp=t1zxp-(fact2*ci-fact4*fact5)*xr+fact2*di(1)-
     &                   fact4*qir(1)
            t1zyd=t1zyd-(fact1*ci-fact3*fact5)*yr+fact1*di(2)-
     &                   fact3*qir(2)
            t1zyp=t1zyp-(fact2*ci-fact4*fact5)*yr+fact2*di(2)-
     &                   fact4*qir(2)
            t1zzd=t1zzd-(fact1*ci-fact3*fact5)*zr+fact1*di(3)-
     &                   fact3*qir(3)
            t1zzp=t1zzp-(fact2*ci-fact4*fact5)*zr+fact2*di(3)-
     &                   fact4*qir(3)
c
C 2ND TERM W/OUT SCREENING
c
            gfri_d(2) = -rr3*ck*dsc3 + rr5*sc(4)*dsc5 - rr7*sc(6)*dsc7
            gfri_p(2) = -rr3*ck*psc3 + rr5*sc(4)*psc5 - rr7*sc(6)*psc7
            t2xxd = t2xxd - 0.5d0*gfri_d(2)
            t2xxp = t2xxp - 0.5d0*gfri_p(2)
            t2yyd = t2yyd - 0.5d0*gfri_d(2)
            t2yyp = t2yyp - 0.5d0*gfri_p(2)
            t2zzd = t2zzd - 0.5d0*gfri_d(2)
            t2zzp = t2zzp - 0.5d0*gfri_p(2)
c
C 3RD TERM W/OUT SCREENING
c
            gfri_d(3) =  rr3*ci*dsc3 + rr5*sc(3)*dsc5 + rr7*sc(5)*dsc7
            gfri_p(3) =  rr3*ci*psc3 + rr5*sc(3)*psc5 + rr7*sc(5)*psc7
            t1xxd = t1xxd - 0.5d0*gfri_d(3)
            t1xxp = t1xxp - 0.5d0*gfri_p(3)
            t1yyd = t1yyd - 0.5d0*gfri_d(3)
            t1yyp = t1yyp - 0.5d0*gfri_p(3)
            t1zzd = t1zzd - 0.5d0*gfri_d(3)
            t1zzp = t1zzp - 0.5d0*gfri_p(3)
c
c    NO: modify the forces for partially excluded interactions
c
C     Apply 'f' factor term
         
             t1xxd=f*t1xxd
             t1xyd=f*t1xyd
             t1xzd=f*t1xzd
             t1yxd=f*t1yxd
             t1yyd=f*t1yyd
             t1yzd=f*t1yzd
             t1zxd=f*t1zxd
             t1zyd=f*t1zyd
             t1zzd=f*t1zzd
             t2xxd=f*t2xxd
             t2xyd=f*t2xyd
             t2xzd=f*t2xzd
             t2yxd=f*t2yxd
             t2yyd=f*t2yyd
             t2yzd=f*t2yzd
             t2zxd=f*t2zxd
             t2zyd=f*t2zyd
             t2zzd=f*t2zzd
             t1xxp=f*t1xxp
             t1xyp=f*t1xyp
             t1xzp=f*t1xzp
             t1yxp=f*t1yxp
             t1yyp=f*t1yyp
             t1yzp=f*t1yzp
             t1zxp=f*t1zxp
             t1zyp=f*t1zyp
             t1zzp=f*t1zzp
             t2xxp=f*t2xxp
             t2xyp=f*t2xyp
             t2xzp=f*t2xzp
             t2yxp=f*t2yxp
             t2yyp=f*t2yyp
             t2yzp=f*t2yzp
             t2zxp=f*t2zxp
             t2zyp=f*t2zyp
             t2zzp=f*t2zzp

             nlocal = nlocal + 1
             ilocal(1,nlocal) = i
             ilocal(2,nlocal) = k
             dlocal1d(1,nlocal)=t1xxd
             dlocal1d(2,nlocal)=t1xyd
             dlocal1d(3,nlocal)=t1xzd
             dlocal1d(4,nlocal)=t1yxd
             dlocal1d(5,nlocal)=t1yyd
             dlocal1d(6,nlocal)=t1yzd
             dlocal1d(7,nlocal)=t1zxd
             dlocal1d(8,nlocal)=t1zyd
             dlocal1d(9,nlocal)=t1zzd
       
             dlocal2d(1,nlocal)=t2xxd
             dlocal2d(2,nlocal)=t2xyd
             dlocal2d(3,nlocal)=t2xzd
             dlocal2d(4,nlocal)=t2yxd
             dlocal2d(5,nlocal)=t2yyd
             dlocal2d(6,nlocal)=t2yzd
             dlocal2d(7,nlocal)=t2zxd
             dlocal2d(8,nlocal)=t2zyd
             dlocal2d(9,nlocal)=t2zzd
       
             dlocal1p(1,nlocal)=t1xxp
             dlocal1p(2,nlocal)=t1xyp
             dlocal1p(3,nlocal)=t1xzp
             dlocal1p(4,nlocal)=t1yxp
             dlocal1p(5,nlocal)=t1yyp
             dlocal1p(6,nlocal)=t1yzp
             dlocal1p(7,nlocal)=t1zxp
             dlocal1p(8,nlocal)=t1zyp
             dlocal1p(9,nlocal)=t1zzp
       
             dlocal2p(1,nlocal)=t2xxp
             dlocal2p(2,nlocal)=t2xyp
             dlocal2p(3,nlocal)=t2xzp
             dlocal2p(4,nlocal)=t2yxp
             dlocal2p(5,nlocal)=t2yyp
             dlocal2p(6,nlocal)=t2yzp
             dlocal2p(7,nlocal)=t2zxp
             dlocal2p(8,nlocal)=t2zyp
             dlocal2p(9,nlocal)=t2zzp

C   NOW, HAVE TO ACCUMULATE ALL THE TERMS FOR THE TORQUE CONTRIBUTIONS TO THE TENSOR!
      tau1d(1)=0.0d0
      tau1d(2)=0.0d0
      tau1d(3)=0.0d0
      tau1d(4)=0.0d0
      tau1d(5)=0.0d0
      tau1d(6)=0.0d0
      tau1d(7)=0.0d0
      tau1d(8)=0.0d0
      tau1d(9)=0.0d0

      tau1p(1)=0.0d0
      tau1p(2)=0.0d0
      tau1p(3)=0.0d0
      tau1p(4)=0.0d0
      tau1p(5)=0.0d0
      tau1p(6)=0.0d0
      tau1p(7)=0.0d0
      tau1p(8)=0.0d0
      tau1p(9)=0.0d0

      tau2d(1)=0.0d0
      tau2d(2)=0.0d0
      tau2d(3)=0.0d0
      tau2d(4)=0.0d0
      tau2d(5)=0.0d0
      tau2d(6)=0.0d0
      tau2d(7)=0.0d0
      tau2d(8)=0.0d0
      tau2d(9)=0.0d0

      tau2p(1)=0.0d0
      tau2p(2)=0.0d0
      tau2p(3)=0.0d0
      tau2p(4)=0.0d0
      tau2p(5)=0.0d0
      tau2p(6)=0.0d0
      tau2p(7)=0.0d0
      tau2p(8)=0.0d0
      tau2p(9)=0.0d0

C   FIRST, TORQUE TERMS W/ EWALD SCREENING

      tau1d(1)= tau1d(1)+ 0.5d0*bn(2)*dixr(1)*xr
     &        +0.5d0*2.0d0*bn(2)*(yr*qi(3)-zr*qi(2))
     &         - bn(3)*rxqir(1)*xr
      tau1p(1)=tau1d(1)
      tau1d(2)= tau1d(2) - 0.5d0*bn(1)*(-di(3)) + 0.5d0*bn(2)*dixr(1)*yr
     &          + 0.5d0*2.0d0*bn(2)*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -bn(3)*rxqir(1)*yr
      tau1p(2)=tau1d(2)
      tau1d(3)=tau1d(3) -0.5d0*bn(1)*di(2) + 0.5d0*bn(2)*dixr(1)*zr
     &          + 0.5d0*2.0d0*bn(2)*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -bn(3)*rxqir(1)*zr
      tau1p(3)=tau1d(3)

      tau1d(4)= tau1d(4) -0.5d0*bn(1)*di(3) + 0.5d0*bn(2)*dixr(2)*xr
     &          +0.5d0*2.0d0*bn(2)*( -qir(3) + zr*qi(1) - xr*qi(3))
     &          -bn(3)*rxqir(2)*xr
      tau1p(4)=tau1d(4)
      tau1d(5)= tau1d(5)+0.5d0*bn(2)*dixr(2)*yr
     &         +0.5d0*2.0d0*bn(2)*( zr*qi(4)-xr*qi(6))
     &         -bn(3)*rxqir(2)*yr
      tau1p(5)=tau1d(5)
      tau1d(6)=tau1d(6) -0.5d0*bn(1)*(-di(1))+0.5d0*bn(2)*dixr(2)*zr
     &         +0.5d0*2.0d0*bn(2)*( qir(1) + zr*qi(7)-xr*qi(9) )
     &         -bn(3)*rxqir(2)*zr
      tau1p(6)=tau1d(6)

      tau1d(7)=tau1d(7) -0.5d0*bn(1)*(-di(2)) + 0.5d0*bn(2)*dixr(3)*xr
     &         +0.5d0*2.0d0*bn(2)*(qir(2) + xr*qi(2)-yr*qi(1) )
     &         -bn(3)*rxqir(3)*xr
      tau1p(7)=tau1d(7)
      tau1d(8)=tau1d(8) -0.5d0*bn(1)*di(1) + 0.5d0*bn(2)*dixr(3)*yr
     &        +0.5d0*2.0d0*bn(2)*(-qir(1) + xr*qi(5)-yr*qi(4) )
     &        -bn(3)*rxqir(3)*yr
      tau1p(8)=tau1d(8)
      tau1d(9)=tau1d(9)+0.5d0*bn(2)*dixr(3)*zr
     &        +0.5d0*2.0d0*bn(2)*(xr*qi(8)-yr*qi(7))
     &        -bn(3)*rxqir(3)*zr
      tau1p(9)=tau1d(9)

      tau2d(1)=tau2d(1) + 0.5d0*bn(2)*dkxr(1)*xr
     &        -0.5d0*2.0d0*bn(2)*(yr*qk(3)-zr*qk(2))
     &        -(-bn(3))*rxqkr(1)*xr
      tau2p(1)=tau2d(1)
      tau2d(2)=tau2d(2) -0.5d0*bn(1)*(-dk(3)) + 0.5d0*bn(2)*dkxr(1)*yr
     &       -0.5d0*2.0d0*bn(2)*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-bn(3))*rxqkr(1)*yr
      tau2p(2)=tau2d(2)
      tau2d(3)=tau2d(3) -0.5d0*bn(1)*dk(2) + 0.5d0*bn(2)*dkxr(1)*zr
     &        -0.5d0*2.0d0*bn(2)*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-bn(3))*rxqkr(1)*zr
      tau2p(3)=tau2d(3)


      tau2d(4)= tau2d(4) -0.5d0*bn(1)*dk(3) + 0.5d0*bn(2)*dkxr(2)*xr
     &          -0.5d0*2.0d0*bn(2)*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-bn(3))*rxqkr(2)*xr
      tau2p(4)=tau2d(4)
      tau2d(5)= tau2d(5)+0.5d0*bn(2)*dkxr(2)*yr
     &         -0.5d0*2.0d0*bn(2)*( zr*qk(4)-xr*qk(6))
     &         -(-bn(3))*rxqkr(2)*yr
      tau2p(5)=tau2d(5)
      tau2d(6)=tau2d(6) -0.5d0*bn(1)*(-dk(1))+0.5d0*bn(2)*dkxr(2)*zr
     &         -0.5d0*2.0d0*bn(2)*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-bn(3))*rxqkr(2)*zr
      tau2p(6)=tau2d(6)

      tau2d(7)=tau2d(7) -0.5d0*bn(1)*(-dk(2)) + 0.5d0*bn(2)*dkxr(3)*xr
     &         -0.5d0*2.0d0*bn(2)*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-bn(3))*rxqkr(3)*xr
      tau2p(7)=tau2d(7)
      tau2d(8)=tau2d(8) -0.5d0*bn(1)*dk(1) + 0.5d0*bn(2)*dkxr(3)*yr
     &        -0.5d0*2.0d0*bn(2)*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-bn(3))*rxqkr(3)*yr
      tau2p(8)=tau2d(8)
      tau2d(9)=tau2d(9)+0.5d0*bn(2)*dkxr(3)*zr
     &        -0.5d0*2.0d0*bn(2)*(xr*qk(8)-yr*qk(7))
     &        -(-bn(3))*rxqkr(3)*zr
      tau2p(9)=tau2d(9)

C   NOW W/OUT SCREENING

      tau1d(1)=tau1d(1)+ (0.5d0*rr5*dsc5*dixr(1)*xr
     &        +0.5d0*2.0d0*rr5*dsc5*(yr*qi(3)-zr*qi(2))
     &         - rr7*dsc7*rxqir(1)*xr)*(-1.0d0)
      tau1d(2)=tau1d(2)+(-0.5d0*rr3*dsc3*(-di(3))
     &             +0.5d0*rr5*dsc5*dixr(1)*yr
     &          + 0.5d0*2.0d0*rr5*dsc5*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -rr7*dsc7*rxqir(1)*yr)*(-1.0d0)
      tau1d(3)=tau1d(3)+ (-0.5d0*rr3*dsc3*di(2)
     &           + 0.5d0*rr5*dsc5*dixr(1)*zr
     &          + 0.5d0*2.0d0*rr5*dsc5*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -rr7*dsc7*rxqir(1)*zr)*(-1.0d0)


      tau1d(4)= tau1d(4)+ (-0.5d0*rr3*dsc3*di(3)
     &             + 0.5d0*rr5*dsc5*dixr(2)*xr
     &          +0.5d0*2.0d0*rr5*dsc5*( -qir(3) + zr*qi(1)-xr*qi(3))
     &          -rr7*dsc7*rxqir(2)*xr)*(-1.0d0)
      tau1d(5)= tau1d(5)+ (0.5d0*rr5*dsc5*dixr(2)*yr
     &         +0.5d0*2.0d0*rr5*dsc5*( zr*qi(4)-xr*qi(6))
     &         -rr7*dsc7*rxqir(2)*yr)*(-1.0d0)
      tau1d(6)=tau1d(6)+ (-0.5d0*rr3*dsc3*(-di(1))
     &               +0.5d0*rr5*dsc5*dixr(2)*zr
     &         +0.5d0*2.0d0*rr5*dsc5*( qir(1) + zr*qi(7)-xr*qi(9))
     &         -rr7*dsc7*rxqir(2)*zr)*(-1.0d0)

      tau1d(7)=tau1d(7)+ (-0.5d0*rr3*dsc3*(-di(2))
     &         + 0.5d0*rr5*dsc5*dixr(3)*xr
     &         +0.5d0*2.0d0*rr5*dsc5*(qir(2) + xr*qi(2)-yr*qi(1))
     &         -rr7*dsc7*rxqir(3)*xr)*(-1.0d0)
      tau1d(8)=tau1d(8)+ (-0.5d0*rr3*dsc3*di(1)
     &         + 0.5d0*rr5*dsc5*dixr(3)*yr
     &        +0.5d0*2.0d0*rr5*dsc5*(-qir(1) + xr*qi(5)-yr*qi(4))
     &        -rr7*dsc7*rxqir(3)*yr)*(-1.0d0)
      tau1d(9)=tau1d(9)+ (0.5d0*rr5*dsc5*dixr(3)*zr
     &        +0.5d0*2.0d0*rr5*dsc5*(xr*qi(8)-yr*qi(7))
     &        -rr7*dsc7*rxqir(3)*zr)*(-1.0d0)


      tau2d(1)=tau2d(1) + (0.5d0*rr5*dsc5*dkxr(1)*xr
     &        -0.5d0*2.0d0*rr5*dsc5*(yr*qk(3)-zr*qk(2))
     &        -(-rr7*dsc7)*rxqkr(1)*xr)*(-1.0d0)
      tau2d(2)=tau2d(2)+ (-0.5d0*rr3*dsc3*(-dk(3))
     &           +0.5d0*rr5*dsc5*dkxr(1)*yr
     &       -0.5d0*2.0d0*rr5*dsc5*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-rr7*dsc7)*rxqkr(1)*yr)*(-1.0d0)
      tau2d(3)=tau2d(3)+ (-0.5d0*rr3*dsc3*dk(2)
     &          + 0.5d0*rr5*dsc5*dkxr(1)*zr
     &        -0.5d0*2.0d0*rr5*dsc5*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-rr7*dsc7)*rxqkr(1)*zr)*(-1.0d0)
c
c
      tau2d(4)= tau2d(4)+ (-0.5d0*rr3*dsc3*dk(3)
     &                 + 0.5d0*rr5*dsc5*dkxr(2)*xr
     &          -0.5d0*2.0d0*rr5*dsc5*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-rr7*dsc7)*rxqkr(2)*xr)*(-1.0d0)
      tau2d(5)= tau2d(5)+ (0.5d0*rr5*dsc5*dkxr(2)*yr
     &         -0.5d0*2.0d0*rr5*dsc5*( zr*qk(4)-xr*qk(6))
     &         -(-rr7*dsc7)*rxqkr(2)*yr)*(-1.0d0)
      tau2d(6)=tau2d(6)+ (-0.5d0*rr3*dsc3*(-dk(1))
     &                  +0.5d0*rr5*dsc5*dkxr(2)*zr
     &         -0.5d0*2.0d0*rr5*dsc5*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-rr7*dsc7)*rxqkr(2)*zr)*(-1.0d0)

      tau2d(7)=tau2d(7)+ (-0.5d0*rr3*dsc3*(-dk(2))
     &                     +0.5d0*rr5*dsc5*dkxr(3)*xr
     &         -0.5d0*2.0d0*rr5*dsc5*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-rr7*dsc7)*rxqkr(3)*xr)*(-1.0d0)
      tau2d(8)=tau2d(8)+ (-0.5d0*rr3*dsc3*dk(1)
     &                     + 0.5d0*rr5*dsc5*dkxr(3)*yr
     &        -0.5d0*2.0d0*rr5*dsc5*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-rr7*dsc7)*rxqkr(3)*yr)*(-1.0d0)
      tau2d(9)=tau2d(9)+(0.5d0*rr5*dsc5*dkxr(3)*zr
     &        -0.5d0*2.0d0*rr5*dsc5*(xr*qk(8)-yr*qk(7))
     &        -(-rr7*dsc7)*rxqkr(3)*zr)*(-1.0d0)

      tau1p(1)=tau1p(1)+ (0.5d0*rr5*psc5*dixr(1)*xr
     &        +0.5d0*2.0d0*rr5*psc5*(yr*qi(3)-zr*qi(2))
     &         - rr7*psc7*rxqir(1)*xr)*(-1.0d0)
      tau1p(2)=tau1p(2)+(-0.5d0*rr3*psc3*(-di(3))
     &             +0.5d0*rr5*psc5*dixr(1)*yr
     &          + 0.5d0*2.0d0*rr5*psc5*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -rr7*psc7*rxqir(1)*yr)*(-1.0d0)
      tau1p(3)=tau1p(3)+ (-0.5d0*rr3*psc3*di(2)
     &           + 0.5d0*rr5*psc5*dixr(1)*zr
     &          + 0.5d0*2.0d0*rr5*psc5*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -rr7*psc7*rxqir(1)*zr)*(-1.0d0)


      tau1p(4)= tau1p(4)+ (-0.5d0*rr3*psc3*di(3)
     &             + 0.5d0*rr5*psc5*dixr(2)*xr
     &          +0.5d0*2.0d0*rr5*psc5*( -qir(3) + zr*qi(1) - xr*qi(3))
     &          -rr7*psc7*rxqir(2)*xr)*(-1.0d0)
      tau1p(5)= tau1p(5)+ (0.5d0*rr5*psc5*dixr(2)*yr
     &         +0.5d0*2.0d0*rr5*psc5*( zr*qi(4)-xr*qi(6))
     &         -rr7*psc7*rxqir(2)*yr)*(-1.0d0)
      tau1p(6)=tau1p(6)+ (-0.5d0*rr3*psc3*(-di(1))
     &               +0.5d0*rr5*psc5*dixr(2)*zr
     &         +0.5d0*2.0d0*rr5*psc5*( qir(1) + zr*qi(7)-xr*qi(9) )
     &         -rr7*psc7*rxqir(2)*zr)*(-1.0d0)

      tau1p(7)=tau1p(7)+ (-0.5d0*rr3*psc3*(-di(2))
     &         + 0.5d0*rr5*psc5*dixr(3)*xr
     &         +0.5d0*2.0d0*rr5*psc5*(qir(2) + xr*qi(2)-yr*qi(1) )
     &         -rr7*psc7*rxqir(3)*xr)*(-1.0d0)
      tau1p(8)=tau1p(8)+ (-0.5d0*rr3*psc3*di(1)
     &         + 0.5d0*rr5*psc5*dixr(3)*yr
     &        +0.5d0*2.0d0*rr5*psc5*(-qir(1) + xr*qi(5)-yr*qi(4) )
     &        -rr7*psc7*rxqir(3)*yr)*(-1.0d0)
      tau1p(9)=tau1p(9)+ (0.5d0*rr5*psc5*dixr(3)*zr
     &        +0.5d0*2.0d0*rr5*psc5*(xr*qi(8)-yr*qi(7))
     &        -rr7*psc7*rxqir(3)*zr)*(-1.0d0)

      tau2p(1)=tau2p(1) + (0.5d0*rr5*psc5*dkxr(1)*xr
     &        -0.5d0*2.0d0*rr5*psc5*(yr*qk(3)-zr*qk(2))
     &        -(-rr7*psc7)*rxqkr(1)*xr)*(-1.0d0)
      tau2p(2)=tau2p(2)+ (-0.5d0*rr3*psc3*(-dk(3))
     &           +0.5d0*rr5*psc5*dkxr(1)*yr
     &       -0.5d0*2.0d0*rr5*psc5*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-rr7*psc7)*rxqkr(1)*yr)*(-1.0d0)
      tau2p(3)=tau2p(3)+ (-0.5d0*rr3*psc3*dk(2)
     &          + 0.5d0*rr5*psc5*dkxr(1)*zr
     &        -0.5d0*2.0d0*rr5*psc5*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-rr7*psc7)*rxqkr(1)*zr)*(-1.0d0)


      tau2p(4)= tau2p(4)+ (-0.5d0*rr3*psc3*dk(3)
     &                 + 0.5d0*rr5*psc5*dkxr(2)*xr
     &          -0.5d0*2.0d0*rr5*psc5*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-rr7*psc7)*rxqkr(2)*xr)*(-1.0d0)
      tau2p(5)= tau2p(5)+ (0.5d0*rr5*psc5*dkxr(2)*yr
     &         -0.5d0*2.0d0*rr5*psc5*( zr*qk(4)-xr*qk(6))
     &         -(-rr7*psc7)*rxqkr(2)*yr)*(-1.0d0)
      tau2p(6)=tau2p(6)+ (-0.5d0*rr3*psc3*(-dk(1))
     &                  +0.5d0*rr5*psc5*dkxr(2)*zr
     &         -0.5d0*2.0d0*rr5*psc5*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-rr7*psc7)*rxqkr(2)*zr)*(-1.0d0)

      tau2p(7)=tau2p(7)+ (-0.5d0*rr3*psc3*(-dk(2))
     &                     +0.5d0*rr5*psc5*dkxr(3)*xr
     &         -0.5d0*2.0d0*rr5*psc5*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-rr7*psc7)*rxqkr(3)*xr)*(-1.0d0)
      tau2p(8)=tau2p(8)+ (-0.5d0*rr3*psc3*dk(1)
     &                     + 0.5d0*rr5*psc5*dkxr(3)*yr
     &        -0.5d0*2.0d0*rr5*psc5*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-rr7*psc7)*rxqkr(3)*yr)*(-1.0d0)
      tau2p(9)=tau2p(9)+(0.5d0*rr5*psc5*dkxr(3)*zr
     &        -0.5d0*2.0d0*rr5*psc5*(xr*qk(8)-yr*qk(7))
     &        -(-rr7*psc7)*rxqkr(3)*zr)*(-1.0d0)

      tau1d(1)=f*tau1d(1)
      tau1d(2)=f*tau1d(2)
      tau1d(3)=f*tau1d(3)

      tau1d(4)=f*tau1d(4)
      tau1d(5)=f*tau1d(5)
      tau1d(6)=f*tau1d(6)

      tau1d(7)=f*tau1d(7)
      tau1d(8)=f*tau1d(8)
      tau1d(9)=f*tau1d(9)

      tau1p(1)=f*tau1p(1)
      tau1p(2)=f*tau1p(2)
      tau1p(3)=f*tau1p(3)

      tau1p(4)=f*tau1p(4)
      tau1p(5)=f*tau1p(5)
      tau1p(6)=f*tau1p(6)

      tau1p(7)=f*tau1p(7)
      tau1p(8)=f*tau1p(8)
      tau1p(9)=f*tau1p(9)

      tau2d(1)=f*tau2d(1)
      tau2d(2)=f*tau2d(2)
      tau2d(3)=f*tau2d(3)

      tau2d(4)=f*tau2d(4)
      tau2d(5)=f*tau2d(5)
      tau2d(6)=f*tau2d(6)

      tau2d(7)=f*tau2d(7)
      tau2d(8)=f*tau2d(8)
      tau2d(9)=f*tau2d(9)

      tau2p(1)=f*tau2p(1)
      tau2p(2)=f*tau2p(2)
      tau2p(3)=f*tau2p(3)

      tau2p(4)=f*tau2p(4)
      tau2p(5)=f*tau2p(5)
      tau2p(6)=f*tau2p(6)

      tau2p(7)=f*tau2p(7)
      tau2p(8)=f*tau2p(8)
      tau2p(9)=f*tau2p(9)

c      tau1d(1)=tau1d(1)
c      tau1d(2)=tau1d(2)
c      tau1d(3)=tau1d(3)
c      tau1d(4)=tau1d(4)
c      tau1d(5)=tau1d(5)
c      tau1d(6)=tau1d(6)
c      tau1d(7)=tau1d(7)
c      tau1d(8)=tau1d(8)
c      tau1d(9)=tau1d(9)
c
c      tau2d(1)=tau2d(1)
c      tau2d(2)=tau2d(2)
c      tau2d(3)=tau2d(3)
c      tau2d(4)=tau2d(4)
c      tau2d(5)=tau2d(5)
c      tau2d(6)=tau2d(6)
c      tau2d(7)=tau2d(7)
c      tau2d(8)=tau2d(8)
c      tau2d(9)=tau2d(9)
c
c      tau1p(1)=tau1p(1)
c      tau1p(2)=tau1p(2)
c      tau1p(3)=tau1p(3)
c      tau1p(4)=tau1p(4)
c      tau1p(5)=tau1p(5)
c      tau1p(6)=tau1p(6)
c      tau1p(7)=tau1p(7)
c      tau1p(8)=tau1p(8)
c      tau1p(9)=tau1p(9)
c
c      tau2p(1)=tau2p(1)
c      tau2p(2)=tau2p(2)
c      tau2p(3)=tau2p(3)
c      tau2p(4)=tau2p(4)
c      tau2p(5)=tau2p(5)
c      tau2p(6)=tau2p(6)
c      tau2p(7)=tau2p(7)
c      tau2p(8)=tau2p(8)
c      tau2p(9)=tau2p(9)


c
c     get the permanent torque with screening
c
             ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)+
     &             gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))-
     &                  gf(5)*rxqir(1) - gf(7)*(rxqikr(1)+qkrxqir(1))
             ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)+
     &             gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))-
     &                  gf(5)*rxqir(2) - gf(7)*(rxqikr(2)+qkrxqir(2))
             ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)+
     &             gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))-
     &                  gf(5)*rxqir(3) - gf(7)*(rxqikr(3)+qkrxqir(3))
             ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)-
     &             gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))-
     &                 gf(6)*rxqkr(1) - gf(7)*(rxqkir(1)-qkrxqir(1))
             ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)-
     &             gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))-
     &                 gf(6)*rxqkr(2) - gf(7)*(rxqkir(2)-qkrxqir(2))
             ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)-
     &             gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))-
     &                 gf(6)*rxqkr(3) - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
             if (dorl) then
               ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)+
     &            gfr(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))-
     &                   gfr(5)*rxqir(1) - gfr(7)*(rxqikr(1)+qkrxqir(1))
               ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)+
     &            gfr(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))-
     &                   gfr(5)*rxqir(2) - gfr(7)*(rxqikr(2)+qkrxqir(2))
               ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)+
     &            gfr(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))-
     &                   gfr(5)*rxqir(3) - gfr(7)*(rxqikr(3)+qkrxqir(3))
               ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1)-
     &            gfr(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))-
     &                  gfr(6)*rxqkr(1) - gfr(7)*(rxqkir(1)-qkrxqir(1))
               ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2)-
     &           gfr(4)*(dixqkr(2)+dkxqir(2) +rxqkdi(2)-2.0d0*qixqk(2))-
     &                  gfr(6)*rxqkr(2) - gfr(7)*(rxqkir(2)-qkrxqir(2))
               ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3)-
     &            gfr(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))-
     &                   gfr(6)*rxqkr(3) - gfr(7)*(rxqkir(3)-qkrxqir(3))
             end if
c
c     handle the case where scaling is used
c
             do j = 1, 3
               ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
               ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
               ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
             end do
c
c     increment gradient due to force and torque on first site
c
             demi(1,ii) = demi(1,ii) + ftm2(1)
             demi(2,ii) = demi(2,ii) + ftm2(2)
             demi(3,ii) = demi(3,ii) + ftm2(3)
           ! call torque_3b_Perm_new(demi,i,ttm2,ttm2i,frcxi,frcyi,frczi)
           call torque3_dEtensor_72015_taunew2_pt1(i,ttm2,ttm2i,frcxi,
     &     frcyi,frczi,demi,tau1d,tau1p,nlocaltau1,ilocaltau1,
     &     taulocal1d,taulocal1p,k)

c
c     increment gradient due to force and torque on second site
c
             demk(1,kk) = demk(1,kk) - ftm2(1)
             demk(2,kk) = demk(2,kk) - ftm2(2)
             demk(3,kk) = demk(3,kk) - ftm2(3)
          ! call torque_3b_Perm_new (demk,k,ttm3,ttm3i,frcxk,frcyk,frczk)
           call torque3_dEtensor_72015_taunew2_pt2(i,ttm3,ttm3i,frcxk,
     &       frcyk,frczk,demk,tau2d,tau2p,nlocaltau2,ilocaltau2,
     &       taulocal2d,taulocal2p,k)

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

             viremrealt(1,1) = viremrealt(1,1) + vxx
             viremrealt(2,1) = viremrealt(2,1) + vyx
             viremrealt(3,1) = viremrealt(3,1) + vzx
             viremrealt(1,2) = viremrealt(1,2) + vyx
             viremrealt(2,2) = viremrealt(2,2) + vyy
             viremrealt(3,2) = viremrealt(3,2) + vzy
             viremrealt(1,3) = viremrealt(1,3) + vzx
             viremrealt(2,3) = viremrealt(2,3) + vzy
             viremrealt(3,3) = viremrealt(3,3) + vzz
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
         end do
         do j = 1, np12(ii)
           dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
           dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
           dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO

!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = toffset0
      toffsettau1(tid)=toffset0tau1
      toffsettau2(tid)=toffset0tau2
c      print*,"tid=",tid,"nlocal=",nlocal
      toffset0 = toffset0 + nlocal
      toffset0tau1 = toffset0tau1 + nlocaltau1
      toffset0tau2 = toffset0tau2 + nlocaltau2
      ntpair_dEtmp = toffset0
      ntpair_tau1tmp = toffset0tau1
      ntpair_tau2tmp = toffset0tau2
      !ntpair_dEtmp = nlocal
      !ntpair_tau1tmp = nlocaltau1
      !ntpair_tau2tmp = nlocaltau2

c      print*,"taskid=",taskid,"ntpair_dEtmp=",ntpair_dEtmp
!$OMP END CRITICAL
         k = toffset(tid)
         do i = 1, nlocal
            m = k + i
           ! m = i
            dEindextmp(1,m) = ilocal(1,i)
            dEindextmp(2,m) = ilocal(2,i)
            do j = 1, 9
               dEd1tmp(j,m) = dlocal1d(j,i)
               dEd2tmp(j,m) = dlocal2d(j,i)
               dEp1tmp(j,m) = dlocal1p(j,i)
               dEp2tmp(j,m) = dlocal2p(j,i)
            end do
         end do

         k = toffsettau1(tid)
         do i = 1, nlocaltau1
            m = k + i
           ! m = i
            tau1indextmp(1,m) = ilocaltau1(1,i)
            tau1indextmp(2,m) = ilocaltau1(2,i)
            do j = 1, 9
               taud1tmp(j,m) = taulocal1d(j,i)
               taup1tmp(j,m) = taulocal1p(j,i)
            end do
         end do

         k = toffsettau2(tid)
         do i = 1, nlocaltau2
            m = k + i
           ! m = i
            tau2indextmp(1,m) = ilocaltau2(1,i)
            tau2indextmp(2,m) = ilocaltau2(2,i)
            do j = 1, 9
               taud2tmp(j,m) = taulocal2d(j,i)
               taup2tmp(j,m) = taulocal2p(j,i)
            end do
         end do

         deallocate (ilocal)
         deallocate (ilocaltau1)
         deallocate (ilocaltau2)
         deallocate (dlocal1d)
         deallocate (dlocal2d)
         deallocate (dlocal1p)
         deallocate (dlocal2p)
         deallocate (taulocal1d)
         deallocate (taulocal2d)
         deallocate (taulocal1p)
         deallocate (taulocal2p)


!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
       
      do i = 1, n
         do j = 1, 3
c            dem(j,i) = dem(j,i) + demi(j,i) + demk(j,i)
            demrealt(j,i) = demrealt(j,i) + demi(j,i) + demk(j,i)
            fieldnpolet(j,i) = fieldnpolet(j,i) + fieldt(j,i)
            fieldpnpolet(j,i) = fieldpnpolet(j,i) + fieldtp(j,i)
         end do
      end do
c      print*,"ntpair_dEtmp=",ntpair_dEtmp
c      print*,"ntpair_tau1tmp=",ntpair_tau1tmp
c      print*,"ntpair_tau2tmp=",ntpair_tau2tmp
c      print*,"End of MPI TensorLoad Bal Innerloop EMreal"
c
c     perform deallocation of some local arrays
c
        deallocate (mscale)
        deallocate (pscale)
        deallocate (dscale)
        deallocate (fieldt)
        deallocate (fieldtp)
        deallocate (demi)
        deallocate (demk)
        deallocate (toffset)
        deallocate (toffsettau1)
        deallocate (toffsettau2)
      return
      end
