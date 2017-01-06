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
      subroutine thg5_small_vir_nolistsend_Nn_noewtens 
     & (fieldnpolet,fieldpnpolet)
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
      use neigh2
      use limits
      use potent
      use neigh
      use paremneigh 
      use mpidat
      use polgrp
      use polar, only: polarity, thole, pdamp
      use openmp
      use dEtensor
      use dEtensor2
      use shunt
      use usage
      implicit none
      include 'mpif.h'
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 temp3,temp5,temp7
      real*8 fact1,fact2,fact3,fact4
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
c      real*8 demrealt(3,npole),viremrealt(3,3),emrealt
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
      real*8 t2xxp,t2xyp,t2xzp,t2yxp,t2yyp,t2yzp,t2zxp,t2zyp,t2zzp
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
      real*8, allocatable :: frcztaulocal1d(:,:)
      real*8, allocatable :: frcytaulocal1d(:,:)
      real*8, allocatable :: frcxtaulocal1d(:,:)
      real*8, allocatable :: frcztaulocal1p(:,:)
      real*8, allocatable :: frcytaulocal1p(:,:)
      real*8, allocatable :: frcxtaulocal1p(:,:)
      real*8, allocatable :: frcztaulocal2d(:,:)
      real*8, allocatable :: frcytaulocal2d(:,:)
      real*8, allocatable :: frcxtaulocal2d(:,:)
      real*8, allocatable :: frcztaulocal2p(:,:)
      real*8, allocatable :: frcytaulocal2p(:,:)
      real*8, allocatable :: frcxtaulocal2p(:,:)
      integer taskid_offset
      real*8 dir,dix,diy,diz,dkr,dkx,dky,dkz,qix,qixx
      real*8 fid(3),fkd(3)
c      real*8 fid1,fid2,fid3,fkd1,fkd2,fkd3
      integer ix,iy,iz,j1
      integer kx,ky,kz
      real*8 fgrp
      logical proceed,use_group,use_intra,usei,usek
      real*8 qixy,qixz,qiy,qiyy,qiyz,qiz,qizz,qkx,qkxx,qkxy,qkxz
      real*8 qky,qkyy,qkyz,qkz,qkzz,rr3_field,qir_scalar,qkr_scalar

c
c     perform dynamic allocation of some local arrays
c
c      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
c      allocate (demi(3,n))
c      allocate (demk(3,n))
      allocate (fieldt(3,n))
      allocate (fieldtp(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
c         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do

c      emrealt=0.0d0
      do i=1,npole
         do j=1,3
c           demi(j,i) = 0.0d0
c           demk(j,i) = 0.0d0
           fieldt(j,i) = 0.0d0
           fieldtp(j,i) = 0.0d0
         end do
      end do
c      do i=1,3
c         do j=1,3
c           viremrealt(j,i)=0.0d0
c         end do
c      end do

c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD2'
      call switch2 (mode)
      nlocal = 0
      nlocaltau1 = 0
      nlocaltau2 = 0

c      toffset0 = 0
c      toffset0tau1 = 0
c      toffset0tau2 = 0

      !if(numtasks_emreal.eq.int((numtasks-2)/2)) then
      if(numtasks_emreal.lt.numtasks) then
         taskid_offset=taskid-2
      else
         taskid_offset=taskid
      end if
      
      sz=last_emreal2(taskid_offset)-start_emreal2(taskid_offset)+1

      do i=last_emreal2(taskid_offset),start_emreal2(taskid_offset)
         iter3=i-start_emreal2(taskid_offset)+1
         do kkk=1,nelst2(i)
            elstgrad_tmp(kkk,iter3)=0
            do j1=1,4
              elsttau1_index_tmp(j1,kkk,iter3)=0
              elsttau2_index_tmp(j1,kkk,iter3)=0
            end do
         end do
      end do

c      allocate (toffset(0:nthread-1))
c      allocate (toffsettau1(0:nthread-1))
c      allocate (toffsettau2(0:nthread-1))
c      print*,"taskid",taskid,"maxsiz_elst(taskid)=",maxsize_elst(taskid)
c      maxlocal=0
c      maxlocal = int(dble(maxlocal)/dble(nthread))
c      maxlocal = int(dble(sz) * dble(maxelst2)/dble(nthread))
c      maxlocal = int(dble(sz) * dble(maxelst2))
c      maxlocaltau=3*maxlocal


       !print*,"VIRTHG5"

c      print*,"taskid,taskoffset",taskid,taskid_offset,
c     & "start_emreal2(taskid_offset)", start_emreal2(taskid_offset)
c      print*,"taskid,taskoffset",taskid,taskid_offset,
c     & "last_emreal2(taskid_offset)", last_emreal2(taskid_offset)

c     c print*,"In Tens taskid=",taskid,"maxlocal=",maxlocal
c      print*,"In Tens taskid=",taskid,"maxlocaltau=",maxlocaltau
c
c     set OpenMP directives for the major loop structure
c
      !allocate (ilocal(2,maxlocal))
      !allocate (dlocal1d(9,maxlocal))
      !allocate (dlocal2d(9,maxlocal))
      !allocate (dlocal1p(9,maxlocal))
      !allocate (dlocal2p(9,maxlocal))
      !allocate (taulocal1d(9,maxlocaltau))
      !allocate (taulocal2d(9,maxlocaltau))
      !allocate (taulocal1p(9,maxlocaltau))
      !allocate (taulocal2p(9,maxlocaltau))
      !allocate (ilocaltau1(2,maxlocaltau))
      !allocate (ilocaltau2(2,maxlocaltau))

      !allocate (frcztaulocal1d(9,maxlocal))
      !allocate (frcytaulocal1d(9,maxlocal))
      !allocate (frcxtaulocal1d(9,maxlocal))
      !allocate (frcztaulocal1p(9,maxlocal))
      !allocate (frcytaulocal1p(9,maxlocal))
      !allocate (frcxtaulocal1p(9,maxlocal))
      !allocate (frcztaulocal2d(9,maxlocal))
      !allocate (frcytaulocal2d(9,maxlocal))
      !allocate (frcxtaulocal2d(9,maxlocal))
      !allocate (frcztaulocal2p(9,maxlocal))
      !allocate (frcytaulocal2p(9,maxlocal))
      !allocate (frcxtaulocal2p(9,maxlocal))



c!$OMP PARALLEL default(private) shared(fieldt,fieldtp,
c!$OMP& npole,electric,dielec,off2,toffset,
c!$OMP& toffset0,dEindextmp,ntpair_dEtmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
c!$OMP& nthread,n,x,y,z,rpole,ipole,pdamp,thole,n12,i12,maxlocal,
c!$OMP& m2scale,p2scale,n13,m3scale,p3scale,i13,n14,i14,p4scale,
c!$OMP& m4scale,np11,ip11,p41scale,n15,m5scale,p5scale,i15,
c!$OMP& np12,ip12,d1scale,d2scale,np13,ip13,d3scale,np14,ip14,d4scale,
c!$OMP& nelst2,elst2,xaxis,yaxis,zaxis,f,aewald,
c!$OMP& start_emreal3,last_emreal3,taskid,tau1indextmp,tau2indextmp,
c!$OMP& taud1tmp,taud2tmp,taup1tmp,taup2tmp,ntpair_tau1tmp,
c!$OMP& ntpair_tau2tmp,toffsettau1,toffset0tau1,toffsettau2,
c!$OMP& toffset0tau2,maxlocaltau,taskid_offset,
c!$OMP& frcztau1d, frcztau1p,frcxtau1d,frcxtau1p,frcytau1d,frcytau1p,
c!$OMP& frcztau2d,frcztau2p,frcxtau2d,frcxtau2p,frcytau2d,frcytau2p,
c!$OMP& sz)
c!$OMP& firstprivate(pscale,dscale,nlocal,nlocaltau1,
c!$OMP& nlocaltau2)

c      allocate (ilocal(2,maxlocal))
c      allocate (dlocal1d(9,maxlocal))
c      allocate (dlocal2d(9,maxlocal))
c      allocate (dlocal1p(9,maxlocal))
c      allocate (dlocal2p(9,maxlocal))
c      allocate (taulocal1d(9,maxlocaltau))
c      allocate (taulocal2d(9,maxlocaltau))
c      allocate (taulocal1p(9,maxlocaltau))
c      allocate (taulocal2p(9,maxlocaltau))
c      allocate (ilocaltau1(2,maxlocaltau))
c      allocate (ilocaltau2(2,maxlocaltau))
c
c      allocate (frcztaulocal1d(9,maxlocal))
c      allocate (frcytaulocal1d(9,maxlocal))
c      allocate (frcxtaulocal1d(9,maxlocal))
c      allocate (frcztaulocal1p(9,maxlocal))
c      allocate (frcytaulocal1p(9,maxlocal))
c      allocate (frcxtaulocal1p(9,maxlocal))
c      allocate (frcztaulocal2d(9,maxlocal))
c      allocate (frcytaulocal2d(9,maxlocal))
c      allocate (frcxtaulocal2d(9,maxlocal))
c      allocate (frcztaulocal2p(9,maxlocal))
c      allocate (frcytaulocal2p(9,maxlocal))
c      allocate (frcxtaulocal2p(9,maxlocal))
c
c      allocate (elstgrad_local(maxelst2,sz))
c      allocate (elsttau1_local(3*maxelst2,sz))
c      allocate (elsttau2_local(3*maxelst2,sz))

!$OMP PARALLEL default(private) shared(fieldt,fieldtp,
!$OMP& npole,electric,dielec,off2,toffset,
!$OMP& toffset0,dEindextmp,ntpair_dEtmp,dEd1tmp,dEd2tmp,dEp1tmp,dEp2tmp,
!$OMP& nthread,n,x,y,z,rpole,ipole,pdamp,thole,n12,i12,maxlocal,
!$OMP& m2scale,p2scale,n13,m3scale,p3scale,i13,n14,i14,p4scale,
!$OMP& m4scale,np11,ip11,p41scale,n15,m5scale,p5scale,i15,
!$OMP& np12,ip12,d1scale,d2scale,np13,ip13,d3scale,np14,ip14,d4scale,
!$OMP& nelst2,elst2,xaxis,yaxis,zaxis,f,aewald,
!$OMP& taskid,tau1indextmp,tau2indextmp,
!$OMP& taud1tmp,taud2tmp,taup1tmp,taup2tmp,ntpair_tau1tmp,
!$OMP& ntpair_tau2tmp,toffsettau1,toffset0tau1,toffsettau2,
!$OMP& toffset0tau2,maxlocaltau,taskid_offset,
!$OMP& frcztau1d, frcztau1p,frcxtau1d,frcxtau1p,frcytau1d,frcytau1p,
!$OMP& frcztau2d,frcztau2p,frcxtau2d,frcxtau2p,frcytau2d,frcytau2p,
!$OMP& sz,elstgrad_tmp,elsttau1_tmp,elsttau2_tmp,nlocal,
!$OMP& start_emreal2,last_emreal2,
!$OMP& nlocaltau1,nlocaltau2)
!$OMP& firstprivate(pscale,dscale)

!$OMP DO reduction(+:fieldt,fieldtp)
!$OMP& schedule(guided)
      do i=start_emreal2(taskid_offset),last_emreal2(taskid_offset)
         iter3=i-start_emreal2(taskid_offset)+1
        ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
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
        dix = di(1)
        diy = di(2)
        diz = di(3)
        qixx = qi(1)
        qixy = qi(2)
        qixz = qi(3)
        qiyy = qi(5)
        qiyz = qi(6)
        qizz = qi(9)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))

c
c     set interaction scaling coefficients for connected atoms
c
        do j = 1, n12(ii)
c          mscale(i12(j,ii)) = m2scale
          pscale(i12(j,ii)) = p2scale
        end do
        do j = 1, n13(ii)
c          mscale(i13(j,ii)) = m3scale
          pscale(i13(j,ii)) = p3scale
        end do
        do j = 1, n14(ii)
c          mscale(i14(j,ii)) = m4scale
          pscale(i14(j,ii)) = p4scale
          do k = 1, np11(ii)
            if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
          end do
        end do
        do j = 1, n15(ii)
c          mscale(i15(j,ii)) = m5scale
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
        do kkk = 1, nelst2(i)
          k = elst2(kkk,i)
        !do kkk = 1, nelst_recv(iter3)
        !   k = elst_recv(kkk,iter3)
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

            dkx = dk(1)
            dky = dk(2)
            dkz = dk(3)
            qkxx = qk(1)
            qkxy = qk(2)
            qkxz = qk(3)
            qkyy = qk(5)
            qkyz = qk(6)
            qkzz = qk(9)

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
            !ralpha = aewald * r
            !bn(0) = erfc(ralpha) / r
            !alsq2 = 2.0d0 * aewald**2
            !alsq2n = 0.0d0
            !if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
            !exp2a = exp(-ralpha**2)
            !do j = 1, 5
            !  bfac = dble(2*j-1)
            !  alsq2n = alsq2 * alsq2n
            !  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
            !end do
c
c     apply Thole polarization damping to scale factors
c
            rr1 = 1.0d0 / r
            rr3 = rr1 / r2
            rr5 = 3.0d0 * rr3 / r2
            rr7 = 5.0d0 * rr5 / r2
            rr9 = 7.0d0 * rr7 / r2
            rr11 = 9.0d0 * rr9 / r2
            !rr5_field = rr3 / r2
            !rr7_field = rr5_field / r2
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
            !dsc3 = 1.0d0 - scale3*dscale(kk)
            !dsc5 = 1.0d0 - scale5*dscale(kk)
            !dsc7 = 1.0d0 - scale7*dscale(kk)
            !psc3 = 1.0d0 - scale3*pscale(kk)
            !psc5 = 1.0d0 - scale5*pscale(kk)
            !psc7 = 1.0d0 - scale7*pscale(kk)

               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)

                  rr3_field = scale3 / (r*r2)
                  rr5_field = 3.0d0 * scale5 / (r*r2*r2)
                  rr7_field = 15.0d0 * scale7 / (r*r2*r2*r2)

c
c     construct necessary auxiliary vectors
c
            !dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
            !dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
            !dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
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
            !qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
            !qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
            !qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
            !qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
            !qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
            !qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
            !qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)-
     &      !           qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
            !qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)-
     &      !           qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
            !qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)-
     &      !           qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
            rxqir(1) = yr*qir(3) - zr*qir(2)
            rxqir(2) = zr*qir(1) - xr*qir(3)
            rxqir(3) = xr*qir(2) - yr*qir(1)
            rxqkr(1) = yr*qkr(3) - zr*qkr(2)
            rxqkr(2) = zr*qkr(1) - xr*qkr(3)
            rxqkr(3) = xr*qkr(2) - yr*qkr(1)
            !rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
            !rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
            !rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
            !rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
            !rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
            !rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
            !qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
            !qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
            !qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
            !qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
            !qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
            !qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
            !qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
            !qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
            !qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
            !dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
            !dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
            !dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
            !dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
            !dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
            !dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
            !rxqidk(1) = yr*qidk(3) - zr*qidk(2)
            !rxqidk(2) = zr*qidk(1) - xr*qidk(3)
            !rxqidk(3) = xr*qidk(2) - yr*qidk(1)
            !rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
            !rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
            !rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
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


            !bcn(1) = bn(1)-(1.0d0-scale3*dscale(kk))*rr3
            !bcn(2) = bn(2)-3.0d0*(1.0d0-scale5*dscale(kk))*rr5_field
            !bcn(3) = bn(3)-15.0d0*(1.0d0-scale7*dscale(kk))*rr7_field
            !fimd(1) = -xr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &      !               - bcn(1)*dk(1) + 2.0d0*bcn(2)*qkr(1)
            !fimd(2) = -yr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &      !               - bcn(1)*dk(2) + 2.0d0*bcn(2)*qkr(2)
            !fimd(3) = -zr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &      !               - bcn(1)*dk(3) + 2.0d0*bcn(2)*qkr(3)
            !fkmd(1) = xr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &      !               - bcn(1)*di(1) - 2.0d0*bcn(2)*qir(1)
            !fkmd(2) = yr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &      !               - bcn(1)*di(2) - 2.0d0*bcn(2)*qir(2)
            !fkmd(3) = zr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &      !               - bcn(1)*di(3) - 2.0d0*bcn(2)*qir(3)
            !bcn(1) = bn(1)-(1.0d0-scale3*pscale(kk))*rr3
            !bcn(2) = bn(2)-3.0d0*(1.0d0-scale5*pscale(kk))*rr5_field
            !bcn(3) = bn(3)-15.0d0*(1.0d0-scale7*pscale(kk))*rr7_field
            !fimp(1) = -xr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &      !               - bcn(1)*dk(1) + 2.0d0*bcn(2)*qkr(1)
            !fimp(2) = -yr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &      !               - bcn(1)*dk(2) + 2.0d0*bcn(2)*qkr(2)
            !fimp(3) = -zr*(bcn(1)*ck-bcn(2)*sc(4)+bcn(3)*sc(6))
     &      !               - bcn(1)*dk(3) + 2.0d0*bcn(2)*qkr(3)
            !fkmp(1) = xr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &      !               - bcn(1)*di(1) - 2.0d0*bcn(2)*qir(1)
            !fkmp(2) = yr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &      !               - bcn(1)*di(2) - 2.0d0*bcn(2)*qir(2)
            !fkmp(3) = zr*(bcn(1)*ci+bcn(2)*sc(3)+bcn(3)*sc(5))
     &      !               - bcn(1)*di(3) - 2.0d0*bcn(2)*qir(3)
            !do j = 1, 3
            !  fieldt(j,i) = fieldt(j,i) + fimd(j)
            !  fieldt(j,k) = fieldt(j,k) + fkmd(j)
            !  fieldtp(j,i) = fieldtp(j,i) + fimp(j)
            !  fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
            !end do

                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir_scalar = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr_scalar = qkx*xr + qky*yr + qkz*zr

            fid(1)=-xr*(rr3_field*ck-rr5_field*dkr+rr7_field*qkr_scalar)
     &                        - rr3_field*dkx + 2.0d0*rr5_field*qkx
            fid(2)=-yr*(rr3_field*ck-rr5_field*dkr+rr7_field*qkr_scalar)
     &                       - rr3_field*dky + 2.0d0*rr5_field*qky
            fid(3)=-zr*(rr3_field*ck-rr5_field*dkr+rr7_field*qkr_scalar)
     &                        - rr3_field*dkz + 2.0d0*rr5_field*qkz
            fkd(1)=xr*(rr3_field*ci+rr5_field*dir+rr7_field*qir_scalar)
     &                        - rr3_field*dix - 2.0d0*rr5_field*qix
            fkd(2)=yr*(rr3_field*ci+rr5_field*dir+rr7_field*qir_scalar)
     &                        - rr3_field*diy - 2.0d0*rr5_field*qiy
            fkd(3)=zr*(rr3_field*ci+rr5_field*dir+rr7_field*qir_scalar)
     &                        - rr3_field*diz - 2.0d0*rr5_field*qiz

                  do j = 1, 3
                   fieldt(j,i) = fieldt(j,i) + fid(j)*dscale(kk)
                   fieldt(j,k) = fieldt(j,k) + fkd(j)*dscale(kk)
                   fieldtp(j,i) = fieldtp(j,i) + fid(j)*pscale(kk)
                   fieldtp(j,k) = fieldtp(j,k) + fkd(j)*pscale(kk)
                  end do

c                   fieldt(1,i) = fieldt(1,i) + fid1*dscale(kk)
c                   fieldt(1,k) = fieldt(1,k) + fkd1*dscale(kk)
c                   fieldtp(1,i) = fieldtp(1,i) + fid1*pscale(kk)
c                   fieldtp(1,k) = fieldtp(1,k) + fkd1*pscale(kk)
c
c                   fieldt(2,i) = fieldt(2,i) + fid2*dscale(kk)
c                   fieldt(2,k) = fieldt(2,k) + fkd2*dscale(kk)
c                   fieldtp(2,i) = fieldtp(2,i) + fid2*pscale(kk)
c                   fieldtp(2,k) = fieldtp(2,k) + fkd2*pscale(kk)
c
c                   fieldt(3,i) = fieldt(3,i) + fid3*dscale(kk)
c                   fieldt(3,k) = fieldt(3,k) + fkd3*dscale(kk)
c                   fieldtp(3,i) = fieldtp(3,i) + fid3*pscale(kk)
c                   fieldtp(3,k) = fieldtp(3,k) + fkd3*pscale(kk)
                 
c
c     calculate the gl functions for permanent components
c
c
c     compute the energy contributions for this interaction
c

c
c     get the real energy without any screening function
c
c
c     set flags to compute components without screening
c
c
c     zero out force and torque components without screening
c
c
c     get the permanent force with screening
c
c
c     get the permanent force without screening
c
c
c     get the induced force with screening
c
C 1ST TERM W/ SCREENING
c

c
C 2ND TERM W SCREENING
c
c
C 3RD TERM W SCREENING
c
c
C 4TH TERM W SCREENING
c

c
C 5TH TERM W SCREENING
c


c
C 6TH TERM W SCREENING
c
c
C 7TH TERM W SCREENING
c

c
C 8TH TERM W SCREENING
c


c
c     get the p/dscaled tensor terms without screening
c

c
C 1ST TERM W/OUT SCREENING
c
C     fact1 and fact2 originally have a minus sign, since ftm2ri is
C     being subtracted.   
            !fact1=-0.5d0*rr5*dsc3         
            !fact2=-0.5d0*rr5*psc3         

            fact1= 0.5d0*rr5*dsc3
            fact2= 0.5d0*rr5*psc3

            t2xxd = t2xxd + fact1*xr*(ck*xr + dk(1)) 
            t2xxp = t2xxp + fact2*xr*(ck*xr + dk(1)) 
            t2xyd = t2xyd + fact1*xr*(ck*yr + dk(2))
            t2xyp = t2xyp + fact2*xr*(ck*yr + dk(2))
            t2xzd = t2xzd + fact1*xr*(ck*zr + dk(3))
            t2xzp = t2xzp + fact2*xr*(ck*zr + dk(3))

            t2yxd = t2yxd + fact1*yr*(ck*xr + dk(1))
            t2yxp = t2yxp + fact2*yr*(ck*xr + dk(1))
            t2yyd = t2yyd + fact1*yr*(ck*yr + dk(2)) 
            t2yyp = t2yyp + fact2*yr*(ck*yr + dk(2)) 
            t2yzd = t2yzd + fact1*yr*(ck*zr + dk(3))
            t2yzp = t2yzp + fact2*yr*(ck*zr + dk(3))

            t2zxd = t2zxd + fact1*zr*(ck*xr + dk(1))
            t2zxp = t2zxp + fact2*zr*(ck*xr + dk(1))
            t2zyd = t2zyd + fact1*zr*(ck*yr + dk(2))
            t2zyp = t2zyp + fact2*zr*(ck*yr + dk(2))
            t2zzd = t2zzd + fact1*zr*(ck*zr + dk(3)) 
            t2zzp = t2zzp + fact2*zr*(ck*zr + dk(3)) 

            t1xxd = t1xxd - fact1*xr*(ci*xr - di(1)) 
            t1xxp = t1xxp - fact2*xr*(ci*xr - di(1)) 
            t1xyd = t1xyd - fact1*xr*(ci*yr - di(2))
            t1xyp = t1xyp - fact2*xr*(ci*yr - di(2))
            t1xzd = t1xzd - fact1*xr*(ci*zr - di(3))
            t1xzp = t1xzp - fact2*xr*(ci*zr - di(3))

            t1yxd = t1yxd - fact1*yr*(ci*xr - di(1))
            t1yxp = t1yxp - fact2*yr*(ci*xr - di(1))
            t1yyd = t1yyd - fact1*yr*(ci*yr - di(2)) 
            t1yyp = t1yyp - fact2*yr*(ci*yr - di(2)) 
            t1yzd = t1yzd - fact1*yr*(ci*zr - di(3))
            t1yzp = t1yzp - fact2*yr*(ci*zr - di(3))

            t1zxd = t1zxd - fact1*zr*(ci*xr - di(1))
            t1zxp = t1zxp - fact2*zr*(ci*xr - di(1))
            t1zyd = t1zyd - fact1*zr*(ci*yr - di(2))
            t1zyp = t1zyp - fact2*zr*(ci*yr - di(2))
            t1zzd = t1zzd - fact1*zr*(ci*zr - di(3)) 
            t1zzp = t1zzp - fact2*zr*(ci*zr - di(3)) 
C  As above, we must flip the signs on fact1 and fact2 in the case of no
C  Ewald
            !fact1=0.5d0*rr7*dsc5*sc(3)
            !fact2=0.5d0*rr7*psc5*sc(3)

            fact1=-0.5d0*rr7*dsc5*sc(3)
            fact2=-0.5d0*rr7*psc5*sc(3)
            
            t1xxd = t1xxd + fact1*xr*xr
            t1xxp = t1xxp + fact2*xr*xr
            t1xyd = t1xyd + fact1*xr*yr
            t1xyp = t1xyp + fact2*xr*yr
            t1xzd = t1xzd + fact1*xr*zr
            t1xzp = t1xzp + fact2*xr*zr

            t1yxd = t1yxd + fact1*yr*xr
            t1yxp = t1yxp + fact2*yr*xr
            t1yyd = t1yyd + fact1*yr*yr
            t1yyp = t1yyp + fact2*yr*yr
            t1yzd = t1yzd + fact1*yr*zr
            t1yzp = t1yzp + fact2*yr*zr

            t1zxd = t1zxd + fact1*zr*xr
            t1zxp = t1zxp + fact2*zr*xr
            t1zyd = t1zyd + fact1*zr*yr
            t1zyp = t1zyp + fact2*zr*yr
            t1zzd = t1zzd + fact1*zr*zr
            t1zzp = t1zzp + fact2*zr*zr

            !fact1=0.5d0*rr7*dsc5*sc(4)
            !fact2=0.5d0*rr7*psc5*sc(4)

            fact1=-0.5d0*rr7*dsc5*sc(4)
            fact2=-0.5d0*rr7*psc5*sc(4)

            t2xxd = t2xxd + fact1*xr*xr
            t2xxp = t2xxp + fact2*xr*xr
            t2xyd = t2xyd + fact1*xr*yr
            t2xyp = t2xyp + fact2*xr*yr
            t2xzd = t2xzd + fact1*xr*zr
            t2xzp = t2xzp + fact2*xr*zr

            t2yxd = t2yxd + fact1*yr*xr
            t2yxp = t2yxp + fact2*yr*xr
            t2yyd = t2yyd + fact1*yr*yr
            t2yyp = t2yyp + fact2*yr*yr
            t2yzd = t2yzd + fact1*yr*zr
            t2yzp = t2yzp + fact2*yr*zr

            t2zxd = t2zxd + fact1*zr*xr
            t2zxp = t2zxp + fact2*zr*xr
            t2zyd = t2zyd + fact1*zr*yr
            t2zyp = t2zyp + fact2*zr*yr
            t2zzd = t2zzd + fact1*zr*zr
            t2zzp = t2zzp + fact2*zr*zr

            !fact1=-rr7*dsc5
            !fact2=-rr7*psc5

            fact1=rr7*dsc5
            fact2=rr7*psc5

            t1xxd = t1xxd + fact1*qir(1)*xr
            t1xxp = t1xxp + fact2*qir(1)*xr
            t1xyd = t1xyd + fact1*qir(2)*xr
            t1xyp = t1xyp + fact2*qir(2)*xr
            t1xzd = t1xzd + fact1*qir(3)*xr
            t1xzp = t1xzp + fact2*qir(3)*xr

            t1yxd = t1yxd + fact1*qir(1)*yr
            t1yxp = t1yxp + fact2*qir(1)*yr
            t1yyd = t1yyd + fact1*qir(2)*yr
            t1yyp = t1yyp + fact2*qir(2)*yr
            t1yzd = t1yzd + fact1*qir(3)*yr
            t1yzp = t1yzp + fact2*qir(3)*yr

            t1zxd = t1zxd + fact1*qir(1)*zr
            t1zxp = t1zxp + fact2*qir(1)*zr
            t1zyd = t1zyd + fact1*qir(2)*zr
            t1zyp = t1zyp + fact2*qir(2)*zr
            t1zzd = t1zzd + fact1*qir(3)*zr
            t1zzp = t1zzp + fact2*qir(3)*zr

            t2xxd = t2xxd - fact1*qkr(1)*xr
            t2xxp = t2xxp - fact2*qkr(1)*xr
            t2xyd = t2xyd - fact1*qkr(2)*xr
            t2xyp = t2xyp - fact2*qkr(2)*xr
            t2xzd = t2xzd - fact1*qkr(3)*xr
            t2xzp = t2xzp - fact2*qkr(3)*xr

            t2yxd = t2yxd - fact1*qkr(1)*yr
            t2yxp = t2yxp - fact2*qkr(1)*yr
            t2yyd = t2yyd - fact1*qkr(2)*yr
            t2yyp = t2yyp - fact2*qkr(2)*yr
            t2yzd = t2yzd - fact1*qkr(3)*yr
            t2yzp = t2yzp - fact2*qkr(3)*yr

            t2zxd = t2zxd - fact1*qkr(1)*zr
            t2zxp = t2zxp - fact2*qkr(1)*zr
            t2zyd = t2zyd - fact1*qkr(2)*zr
            t2zyp = t2zyp - fact2*qkr(2)*zr
            t2zzd = t2zzd - fact1*qkr(3)*zr
            t2zzp = t2zzp - fact2*qkr(3)*zr

            !fact1=-0.5d0*rr9*dsc7*sc(6)
            !fact2=-0.5d0*rr9*psc7*sc(6)
            fact1=0.5d0*rr9*dsc7*sc(6)
            fact2=0.5d0*rr9*psc7*sc(6)
            t2xxd = t2xxd + fact1*xr*xr
            t2xxp = t2xxp + fact2*xr*xr
            t2xyd = t2xyd + fact1*yr*xr
            t2xyp = t2xyp + fact2*yr*xr
            t2xzd = t2xzd + fact1*zr*xr
            t2xzp = t2xzp + fact2*zr*xr

            t2yxd = t2yxd + fact1*xr*yr
            t2yxp = t2yxp + fact2*xr*yr
            t2yyd = t2yyd + fact1*yr*yr
            t2yyp = t2yyp + fact2*yr*yr
            t2yzd = t2yzd + fact1*zr*yr
            t2yzp = t2yzp + fact2*zr*yr

            t2zxd = t2zxd + fact1*xr*zr
            t2zxp = t2zxp + fact2*xr*zr
            t2zyd = t2zyd + fact1*yr*zr
            t2zyp = t2zyp + fact2*yr*zr
            t2zzd = t2zzd + fact1*zr*zr
            t2zzp = t2zzp + fact2*zr*zr
      
            !fact1=-0.5d0*rr9*dsc7*sc(5)
            !fact2=-0.5d0*rr9*psc7*sc(5)
            fact1=0.5d0*rr9*dsc7*sc(5)
            fact2=0.5d0*rr9*psc7*sc(5)
            t1xxd = t1xxd - fact1*xr*xr
            t1xxp = t1xxp - fact2*xr*xr
            t1xyd = t1xyd - fact1*yr*xr
            t1xyp = t1xyp - fact2*yr*xr
            t1xzd = t1xzd - fact1*zr*xr
            t1xzp = t1xzp - fact2*zr*xr

            t1yxd = t1yxd - fact1*xr*yr
            t1yxp = t1yxp - fact2*xr*yr
            t1yyd = t1yyd - fact1*yr*yr
            t1yyp = t1yyp - fact2*yr*yr
            t1yzd = t1yzd - fact1*zr*yr
            t1yzp = t1yzp - fact2*zr*yr

            t1zxd = t1zxd - fact1*xr*zr
            t1zxp = t1zxp - fact2*xr*zr
            t1zyd = t1zyd - fact1*yr*zr
            t1zyp = t1zyp - fact2*yr*zr
            t1zzd = t1zzd - fact1*zr*zr
            t1zzp = t1zzp - fact2*zr*zr
c
C 2ND TERM W/OUT SCREENING
c
            gfri_d(2) = -rr3*ck*dsc3 + rr5*sc(4)*dsc5 - rr7*sc(6)*dsc7
            gfri_p(2) = -rr3*ck*psc3 + rr5*sc(4)*psc5 - rr7*sc(6)*psc7
C 
            !t2xxd = t2xxd - 0.5d0*gfri_d(2)
            !t2xxp = t2xxp - 0.5d0*gfri_p(2)
            !t2yyd = t2yyd - 0.5d0*gfri_d(2)
            !t2yyp = t2yyp - 0.5d0*gfri_p(2)
            !t2zzd = t2zzd - 0.5d0*gfri_d(2)
            !t2zzp = t2zzp - 0.5d0*gfri_p(2)

            t2xxd = t2xxd + 0.5d0*gfri_d(2)
            t2xxp = t2xxp + 0.5d0*gfri_p(2)
            t2yyd = t2yyd + 0.5d0*gfri_d(2)
            t2yyp = t2yyp + 0.5d0*gfri_p(2)
            t2zzd = t2zzd + 0.5d0*gfri_d(2)
            t2zzp = t2zzp + 0.5d0*gfri_p(2)

c
C 3RD TERM W/OUT SCREENING
c
            gfri_d(3) =  rr3*ci*dsc3 + rr5*sc(3)*dsc5 + rr7*sc(5)*dsc7
            gfri_p(3) =  rr3*ci*psc3 + rr5*sc(3)*psc5 + rr7*sc(5)*psc7
            !t1xxd = t1xxd - 0.5d0*gfri_d(3)
            !t1xxp = t1xxp - 0.5d0*gfri_p(3)
            !t1yyd = t1yyd - 0.5d0*gfri_d(3)
            !t1yyp = t1yyp - 0.5d0*gfri_p(3)
            !t1zzd = t1zzd - 0.5d0*gfri_d(3)
            !t1zzp = t1zzp - 0.5d0*gfri_p(3)

            t1xxd = t1xxd + 0.5d0*gfri_d(3)
            t1xxp = t1xxp + 0.5d0*gfri_p(3)
            t1yyd = t1yyd + 0.5d0*gfri_d(3)
            t1yyp = t1yyp + 0.5d0*gfri_p(3)
            t1zzd = t1zzd + 0.5d0*gfri_d(3)
            t1zzp = t1zzp + 0.5d0*gfri_p(3)

c
C 4TH TERM W/OUT SCREENING
c
            !fact1=-0.5d0*rr5*dsc5
            !fact2=-0.5d0*rr5*psc5

            fact1=0.5d0*rr5*dsc5
            fact2=0.5d0*rr5*psc5

            t1xxd = t1xxd + fact1*di(1)*xr
            t1xxp = t1xxp + fact2*di(1)*xr
            t1xyd = t1xyd + fact1*di(1)*yr
            t1xyp = t1xyp + fact2*di(1)*yr
            t1xzd = t1xzd + fact1*di(1)*zr
            t1xzp = t1xzp + fact2*di(1)*zr

            t1yxd = t1yxd + fact1*di(2)*xr
            t1yxp = t1yxp + fact2*di(2)*xr
            t1yyd = t1yyd + fact1*di(2)*yr
            t1yyp = t1yyp + fact2*di(2)*yr
            t1yzd = t1yzd + fact1*di(2)*zr
            t1yzp = t1yzp + fact2*di(2)*zr

            t1zxd = t1zxd + fact1*di(3)*xr
            t1zxp = t1zxp + fact2*di(3)*xr
            t1zyd = t1zyd + fact1*di(3)*yr
            t1zyp = t1zyp + fact2*di(3)*yr
            t1zzd = t1zzd + fact1*di(3)*zr
            t1zzp = t1zzp + fact2*di(3)*zr
c
C 5TH TERM W/OUT SCREENING
c
            t2xxd = t2xxd + fact1*dk(1)*xr
            t2xxp = t2xxp + fact2*dk(1)*xr
            t2xyd = t2xyd + fact1*dk(1)*yr
            t2xyp = t2xyp + fact2*dk(1)*yr
            t2xzd = t2xzd + fact1*dk(1)*zr
            t2xzp = t2xzp + fact2*dk(1)*zr

            t2yxd = t2yxd + fact1*dk(2)*xr
            t2yxp = t2yxp + fact2*dk(2)*xr
            t2yyd = t2yyd + fact1*dk(2)*yr
            t2yyp = t2yyp + fact2*dk(2)*yr
            t2yzd = t2yzd + fact1*dk(2)*zr
            t2yzp = t2yzp + fact2*dk(2)*zr

            t2zxd = t2zxd + fact1*dk(3)*xr
            t2zxp = t2zxp + fact2*dk(3)*xr
            t2zyd = t2zyd + fact1*dk(3)*yr
            t2zyp = t2zyp + fact2*dk(3)*yr
            t2zzd = t2zzd + fact1*dk(3)*zr
            t2zzp = t2zzp + fact2*dk(3)*zr
c
C 6TH TERM W/OUT SCREENING
c
            !fact1=-rr5*dsc5
            !fact2=-rr5*psc5
            fact1=rr5*dsc5
            fact2=rr5*psc5

            t2xxd = t2xxd + fact1*qk(1)
            t2xxp = t2xxp + fact2*qk(1)
            t2xyd = t2xyd + fact1*qk(4)
            t2xyp = t2xyp + fact2*qk(4)
            t2xzd = t2xzd + fact1*qk(7)
            t2xzp = t2xzp + fact2*qk(7)

            t2yxd = t2yxd + fact1*qk(2)
            t2yxp = t2yxp + fact2*qk(2)
            t2yyd = t2yyd + fact1*qk(5)
            t2yyp = t2yyp + fact2*qk(5)
            t2yzd = t2yzd + fact1*qk(8)
            t2yzp = t2yzp + fact2*qk(8)

            t2zxd = t2zxd + fact1*qk(3)
            t2zxp = t2zxp + fact2*qk(3)
            t2zyd = t2zyd + fact1*qk(6)
            t2zyp = t2zyp + fact2*qk(6)
            t2zzd = t2zzd + fact1*qk(9)
            t2zzp = t2zzp + fact2*qk(9)

            t1xxd = t1xxd - fact1*qi(1)
            t1xxp = t1xxp - fact2*qi(1)
            t1xyd = t1xyd - fact1*qi(4)
            t1xyp = t1xyp - fact2*qi(4)
            t1xzd = t1xzd - fact1*qi(7)
            t1xzp = t1xzp - fact2*qi(7)

            t1yxd = t1yxd - fact1*qi(2)
            t1yxp = t1yxp - fact2*qi(2)
            t1yyd = t1yyd - fact1*qi(5)
            t1yyp = t1yyp - fact2*qi(5)
            t1yzd = t1yzd - fact1*qi(8)
            t1yzp = t1yzp - fact2*qi(8)

            t1zxd = t1zxd - fact1*qi(3)
            t1zxp = t1zxp - fact2*qi(3)
            t1zyd = t1zyd - fact1*qi(6)
            t1zyp = t1zyp - fact2*qi(6)
            t1zzd = t1zzd - fact1*qi(9)
            t1zzp = t1zzp - fact2*qi(9)
c
C 7TH TERM W/OUT SCREENING
c
            !fact1=-rr7*dsc7
            !fact2=-rr7*psc7

            fact1=rr7*dsc7
            fact2=rr7*psc7

            t1xxd = t1xxd + fact1*qir(1)*xr
            t1xxp = t1xxp + fact2*qir(1)*xr
            t1xyd = t1xyd + fact1*qir(1)*yr
            t1xyp = t1xyp + fact2*qir(1)*yr
            t1xzd = t1xzd + fact1*qir(1)*zr
            t1xzp = t1xzp + fact2*qir(1)*zr

            t1yxd = t1yxd + fact1*qir(2)*xr
            t1yxp = t1yxp + fact2*qir(2)*xr
            t1yyd = t1yyd + fact1*qir(2)*yr
            t1yyp = t1yyp + fact2*qir(2)*yr
            t1yzd = t1yzd + fact1*qir(2)*zr
            t1yzp = t1yzp + fact2*qir(2)*zr

            t1zxd = t1zxd + fact1*qir(3)*xr
            t1zxp = t1zxp + fact2*qir(3)*xr
            t1zyd = t1zyd + fact1*qir(3)*yr
            t1zyp = t1zyp + fact2*qir(3)*yr
            t1zzd = t1zzd + fact1*qir(3)*zr
            t1zzp = t1zzp + fact2*qir(3)*zr
c
C 8TH TERM W/OUT SCREENING
c
            t2xxd = t2xxd - fact1*qkr(1)*xr
            t2xxp = t2xxp - fact2*qkr(1)*xr
            t2xyd = t2xyd - fact1*qkr(1)*yr
            t2xyp = t2xyp - fact2*qkr(1)*yr
            t2xzd = t2xzd - fact1*qkr(1)*zr
            t2xzp = t2xzp - fact2*qkr(1)*zr

            t2yxd = t2yxd - fact1*qkr(2)*xr
            t2yxp = t2yxp - fact2*qkr(2)*xr
            t2yyd = t2yyd - fact1*qkr(2)*yr
            t2yyp = t2yyp - fact2*qkr(2)*yr
            t2yzd = t2yzd - fact1*qkr(2)*zr
            t2yzp = t2yzp - fact2*qkr(2)*zr

            t2zxd = t2zxd - fact1*qkr(3)*xr
            t2zxp = t2zxp - fact2*qkr(3)*xr
            t2zyd = t2zyd - fact1*qkr(3)*yr
            t2zyp = t2zyp - fact2*qkr(3)*yr
            t2zzd = t2zzd - fact1*qkr(3)*zr
            t2zzp = t2zzp - fact2*qkr(3)*zr
c
c     account for partially excluded interactions in grad field tensor
c
c
C  TENSOR TERMS CORRESPONDING TO TEMP3 FOR 'fridmp(j)' terms
c
            fact1=-0.5d0*rr3*dscale(kk)
            fact2=-0.5d0*rr3*pscale(kk)
            t2xxd = t2xxd + fact1*ddsc3(1)*(ck*xr + dk(1))
            t2xxp = t2xxp + fact2*ddsc3(1)*(ck*xr + dk(1)) 
            t2xyd = t2xyd + fact1*ddsc3(1)*(ck*yr + dk(2)) 
            t2xyp = t2xyp + fact2*ddsc3(1)*(ck*yr + dk(2)) 
            t2xzd = t2xzd + fact1*ddsc3(1)*(ck*zr + dk(3)) 
            t2xzp = t2xzp + fact2*ddsc3(1)*(ck*zr + dk(3)) 
            t2yxd = t2yxd + fact1*ddsc3(2)*(ck*xr + dk(1)) 
            t2yxp = t2yxp + fact2*ddsc3(2)*(ck*xr + dk(1)) 
            t2yyd = t2yyd + fact1*ddsc3(2)*(ck*yr + dk(2)) 
            t2yyp = t2yyp + fact2*ddsc3(2)*(ck*yr + dk(2)) 
            t2yzd = t2yzd + fact1*ddsc3(2)*(ck*zr + dk(3)) 
            t2yzp = t2yzp + fact2*ddsc3(2)*(ck*zr + dk(3)) 
            t2zxd = t2zxd + fact1*ddsc3(3)*(ck*xr + dk(1)) 
            t2zxp = t2zxp + fact2*ddsc3(3)*(ck*xr + dk(1)) 
            t2zyd = t2zyd + fact1*ddsc3(3)*(ck*yr + dk(2)) 
            t2zyp = t2zyp + fact2*ddsc3(3)*(ck*yr + dk(2)) 
            t2zzd = t2zzd + fact1*ddsc3(3)*(ck*zr + dk(3)) 
            t2zzp = t2zzp + fact2*ddsc3(3)*(ck*zr + dk(3)) 

            fact1=-0.5d0*rr3*dscale(kk)
            fact2=-0.5d0*rr3*pscale(kk)
            t1xxd = t1xxd + fact1*ddsc3(1)*(-ci*xr + di(1)) 
            t1xxp = t1xxp + fact2*ddsc3(1)*(-ci*xr + di(1)) 
            t1xyd = t1xyd + fact1*ddsc3(1)*(-ci*yr + di(2)) 
            t1xyp = t1xyp + fact2*ddsc3(1)*(-ci*yr + di(2)) 
            t1xzd = t1xzd + fact1*ddsc3(1)*(-ci*zr + di(3)) 
            t1xzp = t1xzp + fact2*ddsc3(1)*(-ci*zr + di(3)) 
            t1yxd = t1yxd + fact1*ddsc3(2)*(-ci*xr + di(1)) 
            t1yxp = t1yxp + fact2*ddsc3(2)*(-ci*xr + di(1)) 
            t1yyd = t1yyd + fact1*ddsc3(2)*(-ci*yr + di(2)) 
            t1yyp = t1yyp + fact2*ddsc3(2)*(-ci*yr + di(2)) 
            t1yzd = t1yzd + fact1*ddsc3(2)*(-ci*zr + di(3)) 
            t1yzp = t1yzp + fact2*ddsc3(2)*(-ci*zr + di(3)) 
            t1zxd = t1zxd + fact1*ddsc3(3)*(-ci*xr + di(1)) 
            t1zxp = t1zxp + fact2*ddsc3(3)*(-ci*xr + di(1)) 
            t1zyd = t1zyd + fact1*ddsc3(3)*(-ci*yr + di(2)) 
            t1zyp = t1zyp + fact2*ddsc3(3)*(-ci*yr + di(2)) 
            t1zzd = t1zzd + fact1*ddsc3(3)*(-ci*zr + di(3)) 
            t1zzp = t1zzp + fact2*ddsc3(3)*(-ci*zr + di(3)) 
c         
C  TENSOR TERMS CORRESPONDING TO TEMP5 FOR 'fridmp(j)' terms
c
            fact1=0.5d0*rr5*dscale(kk)*sc(3)
            fact2=0.5d0*rr5*pscale(kk)*sc(3)
            t1xxd = t1xxd + ddsc5(1)*fact1*xr 
            t1xxp = t1xxp + ddsc5(1)*fact2*xr 
            t1xyd = t1xyd + ddsc5(1)*fact1*yr 
            t1xyp = t1xyp + ddsc5(1)*fact2*yr 
            t1xzd = t1xzd + ddsc5(1)*fact1*zr 
            t1xzp = t1xzp + ddsc5(1)*fact2*zr 

            t1yxd = t1yxd + ddsc5(2)*fact1*xr 
            t1yxp = t1yxp + ddsc5(2)*fact2*xr 
            t1yyd = t1yyd + ddsc5(2)*fact1*yr 
            t1yyp = t1yyp + ddsc5(2)*fact2*yr 
            t1yzd = t1yzd + ddsc5(2)*fact1*zr 
            t1yzp = t1yzp + ddsc5(2)*fact2*zr 

            t1zxd = t1zxd + ddsc5(3)*fact1*xr 
            t1zxp = t1zxp + ddsc5(3)*fact2*xr 
            t1zyd = t1zyd + ddsc5(3)*fact1*yr 
            t1zyp = t1zyp + ddsc5(3)*fact2*yr 
            t1zzd = t1zzd + ddsc5(3)*fact1*zr 
            t1zzp = t1zzp + ddsc5(3)*fact2*zr 
       
            fact1=0.5d0*rr5*dscale(kk)*sc(4)
            fact2=0.5d0*rr5*pscale(kk)*sc(4)
            t2xxd = t2xxd + ddsc5(1)*fact1*xr 
            t2xxp = t2xxp + ddsc5(1)*fact2*xr 
            t2xyd = t2xyd + ddsc5(1)*fact1*yr 
            t2xyp = t2xyp + ddsc5(1)*fact2*yr 
            t2xzd = t2xzd + ddsc5(1)*fact1*zr 
            t2xzp = t2xzp + ddsc5(1)*fact2*zr 

            t2yxd = t2yxd + ddsc5(2)*fact1*xr 
            t2yxp = t2yxp + ddsc5(2)*fact2*xr 
            t2yyd = t2yyd + ddsc5(2)*fact1*yr 
            t2yyp = t2yyp + ddsc5(2)*fact2*yr 
            t2yzd = t2yzd + ddsc5(2)*fact1*zr 
            t2yzp = t2yzp + ddsc5(2)*fact2*zr 

            t2zxd = t2zxd + ddsc5(3)*fact1*xr 
            t2zxp = t2zxp + ddsc5(3)*fact2*xr 
            t2zyd = t2zyd + ddsc5(3)*fact1*yr 
            t2zyp = t2zyp + ddsc5(3)*fact2*yr 
            t2zzd = t2zzd + ddsc5(3)*fact1*zr 
            t2zzp = t2zzp + ddsc5(3)*fact2*zr 

            fact1=-rr5*dscale(kk)
            fact2=-rr5*pscale(kk)
            t1xxd = t1xxd + fact1*qir(1)*ddsc5(1)
            t1xxp = t1xxp + fact2*qir(1)*ddsc5(1)
            t1xyd = t1xyd + fact1*qir(2)*ddsc5(1)
            t1xyp = t1xyp + fact2*qir(2)*ddsc5(1)
            t1xzd = t1xzd + fact1*qir(3)*ddsc5(1)
            t1xzp = t1xzp + fact2*qir(3)*ddsc5(1)

            t1yxd = t1yxd + fact1*qir(1)*ddsc5(2)
            t1yxp = t1yxp + fact2*qir(1)*ddsc5(2)
            t1yyd = t1yyd + fact1*qir(2)*ddsc5(2)
            t1yyp = t1yyp + fact2*qir(2)*ddsc5(2)
            t1yzd = t1yzd + fact1*qir(3)*ddsc5(2)
            t1yzp = t1yzp + fact2*qir(3)*ddsc5(2)

            t1zxd = t1zxd + fact1*qir(1)*ddsc5(3)
            t1zxp = t1zxp + fact2*qir(1)*ddsc5(3)
            t1zyd = t1zyd + fact1*qir(2)*ddsc5(3)
            t1zyp = t1zyp + fact2*qir(2)*ddsc5(3)
            t1zzd = t1zzd + fact1*qir(3)*ddsc5(3)
            t1zzp = t1zzp + fact2*qir(3)*ddsc5(3)

            t2xxd = t2xxd - fact1*qkr(1)*ddsc5(1)
            t2xxp = t2xxp - fact2*qkr(1)*ddsc5(1)
            t2xyd = t2xyd - fact1*qkr(2)*ddsc5(1)
            t2xyp = t2xyp - fact2*qkr(2)*ddsc5(1)
            t2xzd = t2xzd - fact1*qkr(3)*ddsc5(1)
            t2xzp = t2xzp - fact2*qkr(3)*ddsc5(1)

            t2yxd = t2yxd - fact1*qkr(1)*ddsc5(2)
            t2yxp = t2yxp - fact2*qkr(1)*ddsc5(2)
            t2yyd = t2yyd - fact1*qkr(2)*ddsc5(2)
            t2yyp = t2yyp - fact2*qkr(2)*ddsc5(2)
            t2yzd = t2yzd - fact1*qkr(3)*ddsc5(2)
            t2yzp = t2yzp - fact2*qkr(3)*ddsc5(2)

            t2zxd = t2zxd - fact1*qkr(1)*ddsc5(3)
            t2zxp = t2zxp - fact2*qkr(1)*ddsc5(3)
            t2zyd = t2zyd - fact1*qkr(2)*ddsc5(3)
            t2zyp = t2zyp - fact2*qkr(2)*ddsc5(3)
            t2zzd = t2zzd - fact1*qkr(3)*ddsc5(3)
            t2zzp = t2zzp - fact2*qkr(3)*ddsc5(3)

            fact1=-0.5d0*rr7*dscale(kk)*sc(6)
            fact2=-0.5d0*rr7*pscale(kk)*sc(6)
            t2xxd = t2xxd + fact1*xr*ddsc7(1)
            t2xxp = t2xxp + fact2*xr*ddsc7(1)
            t2xyd = t2xyd + fact1*yr*ddsc7(1)
            t2xyp = t2xyp + fact2*yr*ddsc7(1)
            t2xzd = t2xzd + fact1*zr*ddsc7(1)
            t2xzp = t2xzp + fact2*zr*ddsc7(1)

            t2yxd = t2yxd + fact1*xr*ddsc7(2)
            t2yxp = t2yxp + fact2*xr*ddsc7(2)
            t2yyd = t2yyd + fact1*yr*ddsc7(2)
            t2yyp = t2yyp + fact2*yr*ddsc7(2)
            t2yzd = t2yzd + fact1*zr*ddsc7(2)
            t2yzp = t2yzp + fact2*zr*ddsc7(2)

            t2zxd = t2zxd + fact1*xr*ddsc7(3)
            t2zxp = t2zxp + fact2*xr*ddsc7(3)
            t2zyd = t2zyd + fact1*yr*ddsc7(3)
            t2zyp = t2zyp + fact2*yr*ddsc7(3)
            t2zzd = t2zzd + fact1*zr*ddsc7(3)
            t2zzp = t2zzp + fact2*zr*ddsc7(3)

            fact1=0.5d0*rr7*dscale(kk)*sc(5)
            fact2=0.5d0*rr7*pscale(kk)*sc(5)
            t1xxd = t1xxd + fact1*xr*ddsc7(1)
            t1xxp = t1xxp + fact2*xr*ddsc7(1)
            t1xyd = t1xyd + fact1*yr*ddsc7(1)
            t1xyp = t1xyp + fact2*yr*ddsc7(1)
            t1xzd = t1xzd + fact1*zr*ddsc7(1)
            t1xzp = t1xzp + fact2*zr*ddsc7(1)

            t1yxd = t1yxd + fact1*xr*ddsc7(2)
            t1yxp = t1yxp + fact2*xr*ddsc7(2)
            t1yyd = t1yyd + fact1*yr*ddsc7(2)
            t1yyp = t1yyp + fact2*yr*ddsc7(2)
            t1yzd = t1yzd + fact1*zr*ddsc7(2)
            t1yzp = t1yzp + fact2*zr*ddsc7(2)

            t1zxd = t1zxd + fact1*xr*ddsc7(3)
            t1zxp = t1zxp + fact2*xr*ddsc7(3)
            t1zyd = t1zyd + fact1*yr*ddsc7(3)
            t1zyp = t1zyp + fact2*yr*ddsc7(3)
            t1zzd = t1zzd + fact1*zr*ddsc7(3)
            t1zzp = t1zzp + fact2*zr*ddsc7(3)
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

             !nlocal = nlocal + 1
             !ilocal(1,nlocal) = i
             !ilocal(2,nlocal) = k

C  NEED LOOK-UP TABLE!!! 

             !dlocal1d(1,nlocal)=t1xxd
             !dlocal1d(2,nlocal)=t1xyd
             !dlocal1d(3,nlocal)=t1xzd
             !dlocal1d(4,nlocal)=t1yxd
             !dlocal1d(5,nlocal)=t1yyd
             !dlocal1d(6,nlocal)=t1yzd
             !dlocal1d(7,nlocal)=t1zxd
             !dlocal1d(8,nlocal)=t1zyd
             !dlocal1d(9,nlocal)=t1zzd
       
             !dlocal2d(1,nlocal)=t2xxd
             !dlocal2d(2,nlocal)=t2xyd
             !dlocal2d(3,nlocal)=t2xzd
             !dlocal2d(4,nlocal)=t2yxd
             !dlocal2d(5,nlocal)=t2yyd
             !dlocal2d(6,nlocal)=t2yzd
             !dlocal2d(7,nlocal)=t2zxd
             !dlocal2d(8,nlocal)=t2zyd
             !dlocal2d(9,nlocal)=t2zzd

       
             !dlocal1p(1,nlocal)=t1xxp
             !dlocal1p(2,nlocal)=t1xyp
             !dlocal1p(3,nlocal)=t1xzp
             !dlocal1p(4,nlocal)=t1yxp
             !dlocal1p(5,nlocal)=t1yyp
             !dlocal1p(6,nlocal)=t1yzp
             !dlocal1p(7,nlocal)=t1zxp
             !dlocal1p(8,nlocal)=t1zyp
             !dlocal1p(9,nlocal)=t1zzp
       
             !dlocal2p(1,nlocal)=t2xxp
             !dlocal2p(2,nlocal)=t2xyp
             !dlocal2p(3,nlocal)=t2xzp
             !dlocal2p(4,nlocal)=t2yxp
             !dlocal2p(5,nlocal)=t2yyp
             !dlocal2p(6,nlocal)=t2yzp
             !dlocal2p(7,nlocal)=t2zxp
             !dlocal2p(8,nlocal)=t2zyp
             !dlocal2p(9,nlocal)=t2zzp



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


C   NOW W/OUT SCREENING

      tau1d(1)=tau1d(1)+ (0.5d0*rr5*dsc5*dixr(1)*xr
     &        +0.5d0*2.0d0*rr5*dsc5*(yr*qi(3)-zr*qi(2))
     &         - rr7*dsc7*rxqir(1)*xr)!*(-1.0d0)
      tau1d(2)=tau1d(2)+(-0.5d0*rr3*dsc3*(-di(3))
     &             +0.5d0*rr5*dsc5*dixr(1)*yr
     &          + 0.5d0*2.0d0*rr5*dsc5*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -rr7*dsc7*rxqir(1)*yr)!*(-1.0d0)
      tau1d(3)=tau1d(3)+ (-0.5d0*rr3*dsc3*di(2)
     &           + 0.5d0*rr5*dsc5*dixr(1)*zr
     &          + 0.5d0*2.0d0*rr5*dsc5*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -rr7*dsc7*rxqir(1)*zr)!*(-1.0d0)


      tau1d(4)= tau1d(4)+ (-0.5d0*rr3*dsc3*di(3)
     &             + 0.5d0*rr5*dsc5*dixr(2)*xr
     &          +0.5d0*2.0d0*rr5*dsc5*( -qir(3) + zr*qi(1)-xr*qi(3))
     &          -rr7*dsc7*rxqir(2)*xr)!*(-1.0d0)
      tau1d(5)= tau1d(5)+ (0.5d0*rr5*dsc5*dixr(2)*yr
     &         +0.5d0*2.0d0*rr5*dsc5*( zr*qi(4)-xr*qi(6))
     &         -rr7*dsc7*rxqir(2)*yr)!*(-1.0d0)
      tau1d(6)=tau1d(6)+ (-0.5d0*rr3*dsc3*(-di(1))
     &               +0.5d0*rr5*dsc5*dixr(2)*zr
     &         +0.5d0*2.0d0*rr5*dsc5*( qir(1) + zr*qi(7)-xr*qi(9))
     &         -rr7*dsc7*rxqir(2)*zr)!*(-1.0d0)

      tau1d(7)=tau1d(7)+ (-0.5d0*rr3*dsc3*(-di(2))
     &         + 0.5d0*rr5*dsc5*dixr(3)*xr
     &         +0.5d0*2.0d0*rr5*dsc5*(qir(2) + xr*qi(2)-yr*qi(1))
     &         -rr7*dsc7*rxqir(3)*xr)!*(-1.0d0)
      tau1d(8)=tau1d(8)+ (-0.5d0*rr3*dsc3*di(1)
     &         + 0.5d0*rr5*dsc5*dixr(3)*yr
     &        +0.5d0*2.0d0*rr5*dsc5*(-qir(1) + xr*qi(5)-yr*qi(4))
     &        -rr7*dsc7*rxqir(3)*yr)!*(-1.0d0)
      tau1d(9)=tau1d(9)+ (0.5d0*rr5*dsc5*dixr(3)*zr
     &        +0.5d0*2.0d0*rr5*dsc5*(xr*qi(8)-yr*qi(7))
     &        -rr7*dsc7*rxqir(3)*zr)!*(-1.0d0)


      tau2d(1)=tau2d(1) + (0.5d0*rr5*dsc5*dkxr(1)*xr
     &        -0.5d0*2.0d0*rr5*dsc5*(yr*qk(3)-zr*qk(2))
     &        -(-rr7*dsc7)*rxqkr(1)*xr)!*(-1.0d0)
      tau2d(2)=tau2d(2)+ (-0.5d0*rr3*dsc3*(-dk(3))
     &           +0.5d0*rr5*dsc5*dkxr(1)*yr
     &       -0.5d0*2.0d0*rr5*dsc5*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-rr7*dsc7)*rxqkr(1)*yr)!*(-1.0d0)
      tau2d(3)=tau2d(3)+ (-0.5d0*rr3*dsc3*dk(2)
     &          + 0.5d0*rr5*dsc5*dkxr(1)*zr
     &        -0.5d0*2.0d0*rr5*dsc5*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-rr7*dsc7)*rxqkr(1)*zr)!*(-1.0d0)
c
c
      tau2d(4)= tau2d(4)+ (-0.5d0*rr3*dsc3*dk(3)
     &                 + 0.5d0*rr5*dsc5*dkxr(2)*xr
     &          -0.5d0*2.0d0*rr5*dsc5*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-rr7*dsc7)*rxqkr(2)*xr)!*(-1.0d0)
      tau2d(5)= tau2d(5)+ (0.5d0*rr5*dsc5*dkxr(2)*yr
     &         -0.5d0*2.0d0*rr5*dsc5*( zr*qk(4)-xr*qk(6))
     &         -(-rr7*dsc7)*rxqkr(2)*yr)!*(-1.0d0)
      tau2d(6)=tau2d(6)+ (-0.5d0*rr3*dsc3*(-dk(1))
     &                  +0.5d0*rr5*dsc5*dkxr(2)*zr
     &         -0.5d0*2.0d0*rr5*dsc5*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-rr7*dsc7)*rxqkr(2)*zr)!*(-1.0d0)

      tau2d(7)=tau2d(7)+ (-0.5d0*rr3*dsc3*(-dk(2))
     &                     +0.5d0*rr5*dsc5*dkxr(3)*xr
     &         -0.5d0*2.0d0*rr5*dsc5*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-rr7*dsc7)*rxqkr(3)*xr)!*(-1.0d0)
      tau2d(8)=tau2d(8)+ (-0.5d0*rr3*dsc3*dk(1)
     &                     + 0.5d0*rr5*dsc5*dkxr(3)*yr
     &        -0.5d0*2.0d0*rr5*dsc5*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-rr7*dsc7)*rxqkr(3)*yr)!*(-1.0d0)
      tau2d(9)=tau2d(9)+(0.5d0*rr5*dsc5*dkxr(3)*zr
     &        -0.5d0*2.0d0*rr5*dsc5*(xr*qk(8)-yr*qk(7))
     &        -(-rr7*dsc7)*rxqkr(3)*zr)!*(-1.0d0)

      tau1p(1)=tau1p(1)+ (0.5d0*rr5*psc5*dixr(1)*xr
     &        +0.5d0*2.0d0*rr5*psc5*(yr*qi(3)-zr*qi(2))
     &         - rr7*psc7*rxqir(1)*xr)!*(-1.0d0)
      tau1p(2)=tau1p(2)+(-0.5d0*rr3*psc3*(-di(3))
     &             +0.5d0*rr5*psc5*dixr(1)*yr
     &          + 0.5d0*2.0d0*rr5*psc5*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -rr7*psc7*rxqir(1)*yr)!*(-1.0d0)
      tau1p(3)=tau1p(3)+ (-0.5d0*rr3*psc3*di(2)
     &           + 0.5d0*rr5*psc5*dixr(1)*zr
     &          + 0.5d0*2.0d0*rr5*psc5*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -rr7*psc7*rxqir(1)*zr)!*(-1.0d0)


      tau1p(4)= tau1p(4)+ (-0.5d0*rr3*psc3*di(3)
     &             + 0.5d0*rr5*psc5*dixr(2)*xr
     &          +0.5d0*2.0d0*rr5*psc5*( -qir(3) + zr*qi(1) - xr*qi(3))
     &          -rr7*psc7*rxqir(2)*xr)!*(-1.0d0)
      tau1p(5)= tau1p(5)+ (0.5d0*rr5*psc5*dixr(2)*yr
     &         +0.5d0*2.0d0*rr5*psc5*( zr*qi(4)-xr*qi(6))
     &         -rr7*psc7*rxqir(2)*yr)!*(-1.0d0)
      tau1p(6)=tau1p(6)+ (-0.5d0*rr3*psc3*(-di(1))
     &               +0.5d0*rr5*psc5*dixr(2)*zr
     &         +0.5d0*2.0d0*rr5*psc5*( qir(1) + zr*qi(7)-xr*qi(9) )
     &         -rr7*psc7*rxqir(2)*zr)!*(-1.0d0)

      tau1p(7)=tau1p(7)+ (-0.5d0*rr3*psc3*(-di(2))
     &         + 0.5d0*rr5*psc5*dixr(3)*xr
     &         +0.5d0*2.0d0*rr5*psc5*(qir(2) + xr*qi(2)-yr*qi(1) )
     &         -rr7*psc7*rxqir(3)*xr)!*(-1.0d0)
      tau1p(8)=tau1p(8)+ (-0.5d0*rr3*psc3*di(1)
     &         + 0.5d0*rr5*psc5*dixr(3)*yr
     &        +0.5d0*2.0d0*rr5*psc5*(-qir(1) + xr*qi(5)-yr*qi(4) )
     &        -rr7*psc7*rxqir(3)*yr)!*(-1.0d0)
      tau1p(9)=tau1p(9)+ (0.5d0*rr5*psc5*dixr(3)*zr
     &        +0.5d0*2.0d0*rr5*psc5*(xr*qi(8)-yr*qi(7))
     &        -rr7*psc7*rxqir(3)*zr)!*(-1.0d0)

      tau2p(1)=tau2p(1) + (0.5d0*rr5*psc5*dkxr(1)*xr
     &        -0.5d0*2.0d0*rr5*psc5*(yr*qk(3)-zr*qk(2))
     &        -(-rr7*psc7)*rxqkr(1)*xr)!*(-1.0d0)
      tau2p(2)=tau2p(2)+ (-0.5d0*rr3*psc3*(-dk(3))
     &           +0.5d0*rr5*psc5*dkxr(1)*yr
     &       -0.5d0*2.0d0*rr5*psc5*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-rr7*psc7)*rxqkr(1)*yr)!*(-1.0d0)
      tau2p(3)=tau2p(3)+ (-0.5d0*rr3*psc3*dk(2)
     &          + 0.5d0*rr5*psc5*dkxr(1)*zr
     &        -0.5d0*2.0d0*rr5*psc5*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-rr7*psc7)*rxqkr(1)*zr)!*(-1.0d0)


      tau2p(4)= tau2p(4)+ (-0.5d0*rr3*psc3*dk(3)
     &                 + 0.5d0*rr5*psc5*dkxr(2)*xr
     &          -0.5d0*2.0d0*rr5*psc5*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-rr7*psc7)*rxqkr(2)*xr)!*(-1.0d0)
      tau2p(5)= tau2p(5)+ (0.5d0*rr5*psc5*dkxr(2)*yr
     &         -0.5d0*2.0d0*rr5*psc5*( zr*qk(4)-xr*qk(6))
     &         -(-rr7*psc7)*rxqkr(2)*yr)!*(-1.0d0)
      tau2p(6)=tau2p(6)+ (-0.5d0*rr3*psc3*(-dk(1))
     &                  +0.5d0*rr5*psc5*dkxr(2)*zr
     &         -0.5d0*2.0d0*rr5*psc5*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-rr7*psc7)*rxqkr(2)*zr)!*(-1.0d0)

      tau2p(7)=tau2p(7)+ (-0.5d0*rr3*psc3*(-dk(2))
     &                     +0.5d0*rr5*psc5*dkxr(3)*xr
     &         -0.5d0*2.0d0*rr5*psc5*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-rr7*psc7)*rxqkr(3)*xr)!*(-1.0d0)
      tau2p(8)=tau2p(8)+ (-0.5d0*rr3*psc3*dk(1)
     &                     + 0.5d0*rr5*psc5*dkxr(3)*yr
     &        -0.5d0*2.0d0*rr5*psc5*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-rr7*psc7)*rxqkr(3)*yr)!*(-1.0d0)
      tau2p(9)=tau2p(9)+(0.5d0*rr5*psc5*dkxr(3)*zr
     &        -0.5d0*2.0d0*rr5*psc5*(xr*qk(8)-yr*qk(7))
     &        -(-rr7*psc7)*rxqkr(3)*zr)!*(-1.0d0)

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


c
c     get the permanent torque with screening
c
c
c     get the permanent torque without screening
c
c
c     handle the case where scaling is used
c
c
c     increment gradient due to force and torque on first site
c
           ! call torque_3b_Perm_new(demi,i,ttm2,ttm2i,frcxi,frcyi,frczi)
       !    call torque3_dEtensor_72015_taunew2_pt1(i,ttm2,ttm2i,frcxi,
     & !     frcyi,frczi,demi,tau1d,tau1p,nlocaltau1,ilocaltau1,
     & !     taulocal1d,taulocal1p,k)

        !   call torque3_dEtensor_72015_taunew2_pt1_vir(i,ttm2,ttm2i,
     &  !frcxi,frcyi,frczi,demi,tau1d,tau1p,nlocaltau1,ilocaltau1,
     &  !taulocal1d,taulocal1p,nlocal,frcztaulocal1d,frcztaulocal1p,
     &  !frcxtaulocal1d,frcxtaulocal1p,frcytaulocal1d,frcytaulocal1p,k)

C LEFT OFF HERE!!!!! NEED TO MODIFY '...taunew3_pt1 & 2..'

!$OMP CRITICAL
             nlocal = nlocal + 1
             elstgrad_tmp(kkk,iter3) = nlocal
             
             dEd1tmp(1,nlocal)=t1xxd
             dEd1tmp(2,nlocal)=t1xyd
             dEd1tmp(3,nlocal)=t1xzd
             dEd1tmp(4,nlocal)=t1yxd
             dEd1tmp(5,nlocal)=t1yyd
             dEd1tmp(6,nlocal)=t1yzd
             dEd1tmp(7,nlocal)=t1zxd
             dEd1tmp(8,nlocal)=t1zyd
             dEd1tmp(9,nlocal)=t1zzd

             dEd2tmp(1,nlocal)=t2xxd
             dEd2tmp(2,nlocal)=t2xyd
             dEd2tmp(3,nlocal)=t2xzd
             dEd2tmp(4,nlocal)=t2yxd
             dEd2tmp(5,nlocal)=t2yyd
             dEd2tmp(6,nlocal)=t2yzd
             dEd2tmp(7,nlocal)=t2zxd
             dEd2tmp(8,nlocal)=t2zyd
             dEd2tmp(9,nlocal)=t2zzd

             dEp1tmp(1,nlocal)=t1xxp
             dEp1tmp(2,nlocal)=t1xyp
             dEp1tmp(3,nlocal)=t1xzp
             dEp1tmp(4,nlocal)=t1yxp
             dEp1tmp(5,nlocal)=t1yyp
             dEp1tmp(6,nlocal)=t1yzp
             dEp1tmp(7,nlocal)=t1zxp
             dEp1tmp(8,nlocal)=t1zyp
             dEp1tmp(9,nlocal)=t1zzp

             dEp2tmp(1,nlocal)=t2xxp
             dEp2tmp(2,nlocal)=t2xyp
             dEp2tmp(3,nlocal)=t2xzp
             dEp2tmp(4,nlocal)=t2yxp
             dEp2tmp(5,nlocal)=t2yyp
             dEp2tmp(6,nlocal)=t2yzp
             dEp2tmp(7,nlocal)=t2zxp
             dEp2tmp(8,nlocal)=t2zyp
             dEp2tmp(9,nlocal)=t2zzp


           call torque3_dEtensor_72015_taunew3_pt1_vir(i,ttm2i,
     &  tau1d,tau1p,nlocaltau1,
     &  nlocal,
     &  k,kkk,iter3)
c
c     increment gradient due to force and torque on second site
c
          ! call torque_3b_Perm_new (demk,k,ttm3,ttm3i,frcxk,frcyk,frczk)
        !   call torque3_dEtensor_72015_taunew2_pt2(i,ttm3,ttm3i,frcxk,
     &  !     frcyk,frczk,demk,tau2d,tau2p,nlocaltau2,ilocaltau2,
     &  !     taulocal2d,taulocal2p,k)

        !   call torque3_dEtensor_72015_taunew2_pt2_vir(i,ttm3,ttm3i,
     &  !frcxk,frcyk,frczk,demk,tau2d,tau2p,nlocaltau2,ilocaltau2,
     &  !taulocal2d,taulocal2p,nlocal,frcztaulocal2d,frcztaulocal2p,
     &  !frcxtaulocal2d,frcxtaulocal2p,frcytaulocal2d,frcytaulocal2p,k)

           call torque3_dEtensor_72015_taunew3_pt2_vir(i,ttm3i,
     &  tau2d,tau2p,nlocaltau2,
     &  nlocal,
     &  k,kkk,iter3)
!$OMP END CRITICAL

c
c     increment the internal virial tensor components
c

           end if
   10       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
           pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
           pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
           pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
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

ccc!$OMP CRITICAL
ccc      tid = 0
ccc!$    tid = omp_get_thread_num ()
ccc      toffset(tid) = toffset0
ccc      toffsettau1(tid)=toffset0tau1
ccc      toffsettau2(tid)=toffset0tau2
cccc      print*,"tid=",tid,"nlocal=",nlocal
ccc      toffset0 = toffset0 + nlocal
ccc      toffset0tau1 = toffset0tau1 + nlocaltau1
ccc      toffset0tau2 = toffset0tau2 + nlocaltau2
ccc      ntpair_dEtmp = toffset0
ccc      ntpair_tau1tmp = toffset0tau1
ccc      ntpair_tau2tmp = toffset0tau2
ccc      !ntpair_dEtmp = nlocal
ccc      !ntpair_tau1tmp = nlocaltau1
ccc      !ntpair_tau2tmp = nlocaltau2
cccc      print*,"taskid",taskid,"thrdnum=",tid,"maxlocal=",maxlocal
cccc      print*,"taskid",taskid,"thrdnum=",tid,"nlocal=",nlocal
cccc      print*,"taskid",taskid,"thrdnum=",tid,"maxlocaltau=",maxlocaltau
cccc      print*,"taskid",taskid,"thrdnum=",tid,"nlocaltau1=",nlocaltau1
cccc      print*,"taskid=",taskid,"ntpair_dEtmp=",ntpair_dEtmp
ccc!$OMP END CRITICAL
ccc         k = toffset(tid)
ccc         do i = 1, nlocal
ccc            m = k + i
ccc           ! m = i
ccc            dEindextmp(1,m) = ilocal(1,i)
ccc            dEindextmp(2,m) = ilocal(2,i)
ccc            do j = 1, 9
ccc               dEd1tmp(j,m) = dlocal1d(j,i)
ccc               dEd2tmp(j,m) = dlocal2d(j,i)
ccc               dEp1tmp(j,m) = dlocal1p(j,i)
ccc               dEp2tmp(j,m) = dlocal2p(j,i)
ccc
ccc          frcztau1d(j,m)= frcztaulocal1d(j,i)
ccc          frcytau1d(j,m)= frcytaulocal1d(j,i)
ccc          frcxtau1d(j,m)= frcxtaulocal1d(j,i)
ccc
ccc          frcztau1p(j,m)= frcztaulocal1p(j,i)
ccc          frcytau1p(j,m)= frcytaulocal1p(j,i)
ccc          frcxtau1p(j,m)= frcxtaulocal1p(j,i)
ccc
ccc          frcztau2d(j,m)= frcztaulocal2d(j,i)
ccc          frcytau2d(j,m)= frcytaulocal2d(j,i)
ccc          frcxtau2d(j,m)= frcxtaulocal2d(j,i)
ccc
ccc          frcztau2p(j,m)= frcztaulocal2p(j,i)
ccc          frcytau2p(j,m)= frcytaulocal2p(j,i)
ccc          frcxtau2p(j,m)= frcxtaulocal2p(j,i)
ccc
ccc            end do
ccc         end do

ccc         k = toffsettau1(tid)
ccc         do i = 1, nlocaltau1
ccc            m = k + i
ccc           ! m = i
ccc            tau1indextmp(1,m) = ilocaltau1(1,i)
ccc            tau1indextmp(2,m) = ilocaltau1(2,i)
ccc            do j = 1, 9
ccc               taud1tmp(j,m) = taulocal1d(j,i)
ccc               taup1tmp(j,m) = taulocal1p(j,i)
ccc            end do
ccc         end do
ccc
ccc         k = toffsettau2(tid)
ccc         do i = 1, nlocaltau2
ccc            m = k + i
ccc           ! m = i
ccc            tau2indextmp(1,m) = ilocaltau2(1,i)
ccc            tau2indextmp(2,m) = ilocaltau2(2,i)
ccc            do j = 1, 9
ccc               taud2tmp(j,m) = taulocal2d(j,i)
ccc               taup2tmp(j,m) = taulocal2p(j,i)
ccc            end do
ccc         end do
ccc
ccc         deallocate (ilocal)
ccc         deallocate (ilocaltau1)
ccc         deallocate (ilocaltau2)
ccc         deallocate (dlocal1d)
ccc         deallocate (dlocal2d)
ccc         deallocate (dlocal1p)
ccc         deallocate (dlocal2p)
ccc         deallocate (taulocal1d)
ccc         deallocate (taulocal2d)
ccc         deallocate (taulocal1p)
ccc         deallocate (taulocal2p)

ccc         deallocate (frcztaulocal1d)
ccc         deallocate (frcytaulocal1d)
ccc         deallocate (frcxtaulocal1d)
ccc         deallocate (frcztaulocal1p)
ccc         deallocate (frcytaulocal1p)
ccc         deallocate (frcxtaulocal1p)
ccc         deallocate (frcztaulocal2d)
ccc         deallocate (frcytaulocal2d)
ccc         deallocate (frcxtaulocal2d)
ccc         deallocate (frcztaulocal2p)
ccc         deallocate (frcytaulocal2p)
ccc         deallocate (frcxtaulocal2p)

!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c

         ntpair_dEtmp = nlocal
         ntpair_tau1tmp = nlocaltau1
         ntpair_tau2tmp = nlocaltau2
       
      do i = 1, n
         do j = 1, 3
c            dem(j,i) = dem(j,i) + demi(j,i) + demk(j,i)
c            demrealt(j,i) = demrealt(j,i) + demi(j,i) + demk(j,i)
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

c        deallocate (mscale)
        deallocate (pscale)
        deallocate (dscale)
        deallocate (fieldt)
        deallocate (fieldtp)
c        deallocate (demi)
c        deallocate (demk)
c        deallocate (toffset)
c        deallocate (toffsettau1)
c        deallocate (toffsettau2)
      return
      end
