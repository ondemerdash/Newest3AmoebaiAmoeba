c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0a  --  conjugate gradient dipole solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0a" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c
c
      subroutine induce0a_3b_scf_vacomp3bclust(npole3b,pnum,uind,uinp,
     & cnt)
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use mpole
      use polar, only: polarity, thole, pdamp
      use polpot
      use potent
      use units
      use uprior
      use uprior3b
c      use pme, only: nfft1, nfft2, nfft3
      use pme
      use chunks
      use fft
      use totfield
c      use tarray
      use neigh2clust
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 udsum,upsum
      real*8 a,ap,b,bp
      real*8 sum,sump
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      logical done
      character*6 mode
      integer npole3b,pnum(*),l1,l3,cnt
      real*8 uind(3,npole3b),uinp(3,npole3b)
c      real*8, allocatable :: minv(:)
c      integer, allocatable :: mindex(:)
c      real*8 M_tot(3*npole3b,3*npole3b)

c      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
c      real*8 qfac3b(nfft1,nfft2,nfft3)
c      real*8 thetai1_3b(4,bsorder,n)
c      real*8 thetai2_3b(4,bsorder,n)
c      real*8 thetai3_3b(4,bsorder,n)
c      integer pmetable3b(n,nchunk)
c      integer igrid3b(3,n)
c      integer*8 planf3b
c      integer*8 planb3b
c      integer iprime3b(maxprime,3)
c      real*8 ffttable3b(maxtable,3)

c         if (poltyp .eq. 'MUTUAL') then
c            if (.not.allocated(tindex))  
c     &         allocate (tindex(2,npole3b*maxelst))
c            if (.not.allocated(tdipdip))
c     &         allocate (tdipdip(6,npole3b*maxelst))
c         end if

c
c
c     zero out the induced dipoles at each site
c
c      do i = 1, npole
      do l1 = 1, npole3b
         do j = 1, 3
c            uind(j,i) = 0.0d0
c            uinp(j,i) = 0.0d0
            uind(j,l1) = 0.0d0
            uinp(j,l1) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
c      allocate (field(3,npole))
c      allocate (fieldp(3,npole))
c      allocate (udir(3,npole))
c      allocate (udirp(3,npole))

      allocate (field(3,npole3b))
      allocate (fieldp(3,npole3b))
      allocate (udir(3,npole3b))
      allocate (udirp(3,npole3b))
c      allocate (minv(6*maxulst*npole3b))
c      allocate (mindex(npole3b))

         allocate (poli(npole3b))
         allocate (rsd(3,npole3b))
         allocate (rsdp(3,npole3b))
         allocate (zrsd(3,npole3b))
         allocate (zrsdp(3,npole3b))
         allocate (conj(3,npole3b))
         allocate (conjp(3,npole3b))
         allocate (vec(3,npole3b))
         allocate (vecp(3,npole3b))

      do l1=1,npole3b
         poli(l1)=0.0d0
         do j=1,3
            field(j,l1)=0.0d0
            fieldp(j,l1)=0.0d0
            uind(j,l1)=0.0d0
            uinp(j,l1)=0.0d0
            udir(j,l1)=0.0d0
            udirp(j,l1)=0.0d0
            rsd(j,l1)=0.0d0
            rsdp(j,l1)=0.0d0
            zrsd(j,l1)=0.0d0
            zrsdp(j,l1)=0.0d0
            conj(j,l1)=0.0d0
            conjp(j,l1)=0.0d0
            vec(j,l1)=0.0d0
            vecp(j,l1)=0.0d0
         end do
      end do

c         do l1 = 1, 3*npole3b
c            do l3 = 1, 3*npole3b
c               M_tot(l1,l3) = 0.0d0
c            end do
c         end do
      
c
c     get the electrostatic field due to permanent multipoles
c
        if(use_ewaldclust) then
         call dfield0c_3b_vacomp3b(npole3b,pnum,field,fieldp,cnt)
        else
         call dfield0b_3b_vacomp3b(npole3b,pnum,field,fieldp,cnt)
        end if
         !call field_noewald_umutual_rl_3b (field,fieldp,M_tot,
     &   !npole3b,pnum)
c
c     set induced dipoles to polarizability times direct field
c
c        print*,"Doing Vacuum embedding"
c      do i = 1, npole
         do l1 =1,npole3b
            i=pnum(l1)
            do j = 1, 3
c            udir(j,i) = polarity(i) * field(j,i)
c            udirp(j,i) = polarity(i) * fieldp(j,i)
c            uind(j,i) = udir(j,i)
c            uinp(j,i) = udirp(j,i)
               udir(j,l1) = polarity(i) * field(j,l1)
               udirp(j,l1) = polarity(i) * fieldp(j,l1)
               uind(j,l1) = udir(j,l1)
               uinp(j,l1) = udirp(j,l1)
            end do
         end do
c      end if 
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c         print*,"poltyp=",poltyp
c
c     estimated induced dipoles from polynomial predictor
c
      !print*,"use_pred=",use_pred

         if (use_pred) then
           if(nualt3b(cnt).eq.maxualt) then
            call ulspred3b(npole3b,cnt)
            !do i = 1, npole
            do l1 =1, npole3b
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt3b(cnt)-1
                   udsum = udsum + bpred3b(k,cnt)*udalt3b(k,j,l1,cnt)
                   upsum = upsum + bpredp3b(k,cnt)*upalt3b(k,j,l1,cnt)
                  end do
                  uind(j,l1) = udsum
                  uinp(j,l1) = upsum
               end do
            end do
           end if
         end if

c
c     perform dynamic allocation of some local arrays
c
c         allocate (poli(npole))
c         allocate (rsd(3,npole))
c         allocate (rsdp(3,npole))
c         allocate (zrsd(3,npole))
c         allocate (zrsdp(3,npole))
c         allocate (conj(3,npole))
c         allocate (conjp(3,npole))
c         allocate (vec(3,npole))
c         allocate (vecp(3,npole))

c         allocate (poli(npole3b))
c         allocate (rsd(3,npole3b))
c         allocate (rsdp(3,npole3b))
c         allocate (zrsd(3,npole3b))
c         allocate (zrsdp(3,npole3b))
c         allocate (conj(3,npole3b))
c         allocate (conjp(3,npole3b))
c         allocate (vec(3,npole3b))
c         allocate (vecp(3,npole3b))

c
c     get the electrostatic field due to induced dipoles
c
          if(use_ewaldclust) then
              if(useboxclust) then
               !call ufield0c_3b_vacompbox3bclust(npole3b,pnum,
     &         !field,fieldp,uind,uinp,cnt)
              else
                call ufield0c_3b_vacomp (npole3b,pnum,field,fieldp,
     &         uind,uinp)
              end if
          else
           call ufield0a_3b_vacomp (npole3b,pnum,field,fieldp,
     &         uind,uinp)
          end if
c
c     set initial conjugate gradient residual and conjugate vector
c

c         do i = 1, npole
         do l1= 1,npole3b
            i=pnum(l1)
            poli(l1) = max(polmin,polarity(i))
            do j = 1, 3
c               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
c     &                       + field(j,i)
c               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
c     &                       + fieldp(j,i)
               rsd(j,l1) = (udir(j,l1)-uind(j,l1))/poli(l1)
     &                       + field(j,l1)
               rsdp(j,l1) = (udirp(j,l1)-uinp(j,l1))/poli(l1)
     &                       + fieldp(j,l1)

            end do
         end do
         mode = 'BUILD'
c         if (use_mlist) then
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c         else
           if(useboxclust) then
            !  call uscale0b_3b_vacomp3btestbox3bclust(npole3b,pnum,mode,
     &      !    rsd,rsdp,zrsd,zrsdp,cnt)
            ! mode = 'APPLY'
            !  call uscale0b_3b_vacomp3btestbox3bclust(npole3b,pnum,mode,
     &      !    rsd,rsdp,zrsd,zrsdp,cnt)
           else
               call uscale0b_3b_vacomp3btest(npole3b,pnum,mode,rsd,rsdp,
     &          zrsd,zrsdp,cnt)
             mode = 'APPLY'
               call uscale0b_3b_vacomp3btest(npole3b,pnum,mode,rsd,rsdp,
     &          zrsd,zrsdp,cnt)
           end if 
c         end if
c         do i = 1, npole
         do l1 =1,npole3b
            i=pnum(l1)
            do j = 1, 3
c               conj(j,i) = zrsd(j,i)
c               conjp(j,i) = zrsdp(j,i)
               conj(j,l1) = zrsd(j,l1)
               conjp(j,l1) = zrsdp(j,l1)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
c         print*,"POLAR-EPS=",poleps

         do while (.not. done)
            iter = iter + 1
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  vec(j,i) = uind(j,i)
c                  vecp(j,i) = uinp(j,i)
c                  uind(j,i) = conj(j,i)
c                  uinp(j,i) = conjp(j,i)

                  vec(j,l1) = uind(j,l1)
                  vecp(j,l1) = uinp(j,l1)
                  uind(j,l1) = conj(j,l1)
                  uinp(j,l1) = conjp(j,l1)
               end do
            end do
          if(use_ewaldclust) then
              if(useboxclust) then
               ! call ufield0c_3b_vacompbox3bclust(npole3b,pnum,
     &         !     field,fieldp,uind,uinp,cnt)
              else
                call ufield0c_3b_vacomp (npole3b,pnum,field,fieldp,
     &               uind,uinp)
              end if
          else
           call ufield0a_3b_vacomp (npole3b,pnum,field,fieldp,
     &          uind,uinp)
          end if
c            do i = 1, npole
            do l1 =1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  uind(j,i) = vec(j,i)
c                  uinp(j,i) = vecp(j,i)
c                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
c                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  uind(j,l1) = vec(j,l1)
                  uinp(j,l1) = vecp(j,l1)
                  vec(j,l1) = conj(j,l1)/poli(l1) - field(j,l1)
                  vecp(j,l1) = conjp(j,l1)/poli(l1) - fieldp(j,l1)

               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  a = a + conj(j,i)*vec(j,i)
c                  ap = ap + conjp(j,i)*vecp(j,i)
c                  sum = sum + rsd(j,i)*zrsd(j,i)
c                  sump = sump + rsdp(j,i)*zrsdp(j,i)

                  a = a + conj(j,l1)*vec(j,l1)
                  ap = ap + conjp(j,l1)*vecp(j,l1)
                  sum = sum + rsd(j,l1)*zrsd(j,l1)
                  sump = sump + rsdp(j,l1)*zrsdp(j,l1)

               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  uind(j,i) = uind(j,i) + a*conj(j,i)
c                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
c                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
c                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                  uind(j,l1) = uind(j,l1) + a*conj(j,l1)
                  uinp(j,l1) = uinp(j,l1) + ap*conjp(j,l1)
                  rsd(j,l1) = rsd(j,l1) - a*vec(j,l1)
                  rsdp(j,l1) = rsdp(j,l1) - ap*vecp(j,l1)
               end do
            end do
c            if (use_mlist) then
c               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            else
           if(useboxclust) then
            ! call uscale0b_3b_vacomp3btestbox3bclust(npole3b,pnum,mode,
     &      !    rsd,rsdp,zrsd,zrsdp,cnt)
           else
             call uscale0b_3b_vacomp3btest(npole3b,pnum,mode,rsd,rsdp,
     &                    zrsd,zrsdp,cnt)
           end if
c            end if
            b = 0.0d0
            bp = 0.0d0
c            do i = 1, npole
            do l1 =1, npole3b
               i=pnum(l1)
               do j = 1, 3
c                  b = b + rsd(j,i)*zrsd(j,i)
c                  bp = bp + rsdp(j,i)*zrsdp(j,i)
                  b = b + rsd(j,l1)*zrsd(j,l1)
                  bp = bp + rsdp(j,l1)*zrsdp(j,l1)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1) 
               do j = 1, 3
c                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
c                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
c                  epsd = epsd + rsd(j,i)*rsd(j,i)
c                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  conj(j,l1) = zrsd(j,l1) + b*conj(j,l1)
                  conjp(j,l1) = zrsdp(j,l1) + bp*conjp(j,l1)
                  epsd = epsd + rsd(j,l1)*rsd(j,l1)
                  epsp = epsp + rsdp(j,l1)*rsdp(j,l1)

               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp)
c            eps = debye * sqrt(eps/dble(npolar))
            eps = debye * sqrt(eps/dble(npole3b))

            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         !deallocate (poli)
         !deallocate (rsd)
         !deallocate (rsdp)
         !deallocate (zrsd)
         !deallocate (zrsdp)
         !deallocate (conj)
         !deallocate (conjp)
         !deallocate (vec)
         !deallocate (vecp)
         !deallocate (minv)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (conj)
         deallocate (conjp)
         deallocate (vec)
         deallocate (vecp)
c         deallocate (minv)
c         deallocate (mindex)
      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
c
c     update the lists of previous induced dipole values
c
      if (use_pred) then
         nualt3b(cnt) = min(nualt3b(cnt)+1,maxualt)
         !do i = 1, npole
         do l1 = 1, npole3b
            i=pnum(l1)
            do j = 1, 3
               do k = nualt3b(cnt), 2, -1
                  udalt3b(k,j,l1,cnt) = udalt3b(k-1,j,l1,cnt)
                  upalt3b(k,j,l1,cnt) = upalt3b(k-1,j,l1,cnt)
               end do
               udalt3b(1,j,l1,cnt) = uind(j,l1)
               upalt3b(1,j,l1,cnt) = uinp(j,l1)
            end do
         end do
      end if

c         if (poltyp .eq. 'MUTUAL') then
c             deallocate (tindex)
c             deallocate (tdipdip)
c         end if

      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dfield0a_3b  --  direct induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dfield0a_3b" computes the direct electrostatic field due to
c     permanent multipole moments via a double loop
c
c
      subroutine dfield0b_3b_vacomp3b(npole3b,pnum,field,fieldp,
     &  cnt)
      use sizes
      use atoms
      use bound
      use cell
      use couple
      use group
      use mpole
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use shunt
      use tarray
      use openmp
      use neigh
      use neigh2clust
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      character*6 mode
      integer npole3b,pnum(*),l1,l3    
      real*8 field(3,npole3b)
      real*8 fieldp(3,npole3b)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: dlocal(:,:)
      integer nlocal,maxlocal
      integer tid,toffset0
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      integer cnt,kkk
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole3b
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)

      nlocal = 0
      toffset0 = 0
      maxlocal = int(dble(npole3b)*dble(maxelst)/dble(nthread))

c      print*,"off2 in dfield0a_3b",off2 
c
c     perform dynamic allocation of some local arrays
c
      !allocate (dscale(n))
      !allocate (pscale(n))
      allocate (dscale(npole3b))
      allocate (pscale(npole3b))
      allocate (toffset(0:nthread-1))
      allocate (uscale(npole3b))
      allocate (fieldt(3,npole3b))
      allocate (fieldtp(3,npole3b))

         do j = 1,npole3b
            dscale(j) = 1.0d0
            pscale(j) = 1.0d0
            uscale(j) = 1.0d0
         end do

c
c     find the electrostatic field due to permanent multipoles
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,thole,
!$OMP& rpole,p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,
!$OMP& d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,nelst3b,
!$OMP& elst3b,cut2,poltyp,ntpair,tindex,tdipdip,
!$OMP& toffset,toffset0,field,fieldp,fieldt,fieldtp,maxlocal,
!$OMP& npole3b,pnum,off2,cnt)
!$OMP& firstprivate(pscale,dscale,uscale,nlocal)
c
c     perform dynamic allocation of some local arrays
c
      if (poltyp .eq. 'MUTUAL') then
         allocate (ilocal(2,maxlocal))
         allocate (dlocal(6,maxlocal))
      end if
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole3b
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     compute the real space portion of the Ewald summation
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)

      do l1 = 1, npole3b
         i=pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
c         do j = i+1, npole
c            dscale(ipole(j)) = 1.0d0
c            pscale(ipole(j)) = 1.0d0
c         end do
c         do j = 1, n12(ii)
c            pscale(i12(j,ii)) = p2scale
c         end do
c         do j = 1, n13(ii)
c            pscale(i13(j,ii)) = p3scale
c         end do
c         do j = 1, n14(ii)
c            pscale(i14(j,ii)) = p4scale
c            do k = 1, np11(ii)
c               if (i14(j,ii) .eq. ip11(k,ii))
c     &            pscale(i14(j,ii)) = p4scale * p41scale
c            end do
c         end do
c         do j = 1, n15(ii)
c            pscale(i15(j,ii)) = p5scale
c         end do
c         do j = 1, np11(ii)
c            dscale(ip11(j,ii)) = d1scale
c         end do
c         do j = 1, np12(ii)
c            dscale(ip12(j,ii)) = d2scale
c         end do
c         do j = 1, np13(ii)
c            dscale(ip13(j,ii)) = d3scale
c         end do
c         do j = 1, np14(ii)
c            dscale(ip14(j,ii)) = d4scale
c         end do
         !do j = 1,npole3b
         !   dscale(j) = 1.0d0
         !   pscale(j) = 1.0d0
         !end do
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

         !do l3=l1+1,npole3b
         do kkk=1,nelst3b(l1,cnt)
            l3=elst3b(kkk,l1,cnt)
            k=pnum(l3)
            kk = ipole(k)
            proceed = .true.
          !  if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               !if (use_bounds)  call image (xr,yr,zr)
               call image(xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz

               if (poltyp .eq. 'MUTUAL') then
                !TERMS CONTAINING BCN COMMENTED OUT BELOW.  NO EWALDKIRSCHTORTE FOR US!
                 ! bcn(1) = bn(1) - (1.0d0-scale3*uscale(kk))*rr3
                 ! bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(kk))*rr5
                  nlocal = nlocal + 1
                  ilocal(1,nlocal) = l1
                  ilocal(2,nlocal) = l3
                  !dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                  !dlocal(2,nlocal) = bcn(2)*xr*yr
                  !dlocal(3,nlocal) = bcn(2)*xr*zr
                  !dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                  !dlocal(5,nlocal) = bcn(2)*yr*zr
                  !dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr

                  dlocal(1,nlocal) = (-rr3 + rr5*xr*xr)*uscale(l3)
                  dlocal(2,nlocal) = (rr5*xr*yr)*uscale(l3)
                  dlocal(3,nlocal) = (rr5*xr*zr)*uscale(l3)
                  dlocal(4,nlocal) = (-rr3 + rr5*yr*yr)*uscale(l3)
                  dlocal(5,nlocal) = (rr5*yr*zr)*uscale(l3)
                  dlocal(6,nlocal) = (-rr3 + rr5*zr*zr)*uscale(l3)
               end if

                  do j = 1, 3
c                     field(j,l1) = field(j,l1) + fid(j)*dscale(kk)
c                     field(j,l3) = field(j,l3) + fkd(j)*dscale(kk)
c                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(kk)
c                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(kk)

c                     field(j,l1) = field(j,l1) + fid(j)*dscale(l3)
c                     field(j,l3) = field(j,l3) + fkd(j)*dscale(l3)
c                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(l3)
c                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(l3)

                     fieldt(j,l1) = fieldt(j,l1) + fid(j)*dscale(l3)
                     fieldt(j,l3) = fieldt(j,l3) + fkd(j)*dscale(l3)
                     fieldtp(j,l1) = fieldtp(j,l1) + fid(j)*pscale(l3)
                     fieldtp(j,l3) = fieldtp(j,l3) + fkd(j)*pscale(l3)
                  end do
               end if
            end if
         end do
         do j=1,npole3b
            pscale(j)=1.0d0
            dscale(j)=1.0d0
            uscale(j)=1.0d0
         end do
      end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
      do i = 1, npole3b
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
c
c     store terms needed later to compute mutual polarization
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = toffset0
      toffset0 = toffset0 + nlocal
      ntpair = toffset0
!$OMP END CRITICAL
      if (poltyp .eq. 'MUTUAL') then
         k = toffset(tid)
         do i = 1, nlocal
            m = k + i
            tindex(1,m) = ilocal(1,i)
            tindex(2,m) = ilocal(2,i)
            do j = 1, 6
               tdipdip(j,m) = dlocal(j,i)
            end do
         end do
         deallocate (ilocal)
         deallocate (dlocal)
      end if
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (toffset)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end

c
c
c
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0a  --  dipole preconditioner via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0a" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a double loop
c
c
      subroutine uscale0b_3b_vacomp3b(npole3b,pnum,mode,rsd,rsdp,zrsd,
     &  zrsdp,cnt)
      use sizes
      use atoms
      use limits
      use mpole
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use neigh
      use neigh2clust
      use usolve
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp
      real*8 pgamma,off2
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
c      real*8 rsd(3,*)
c      real*8 rsdp(3,*)
c      real*8 zrsd(3,*)
c      real*8 zrsdp(3,*)
      character*6 mode
      integer npole3b,pnum(*),l1,l3
c      real*8 minv(*)
      real*8 rsd(3,npole3b)
      real*8 rsdp(3,npole3b)
      real*8 zrsd(3,npole3b)
      real*8 zrsdp(3,npole3b)
      real*8, allocatable :: zrsdt(:,:)
      real*8, allocatable :: zrsdtp(:,:)
c      integer mindex(*),cnt,kkk
      integer cnt,kkk
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then

         allocate (zrsdt(3,npole3b))
         allocate (zrsdtp(3,npole3b))

c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
c         do i = 1, npole
         do l1 = 1, npole3b
            i=pnum(l1)
            poli = udiag * max(polmin,polarity(i))
            do j = 1, 3
c               zrsd(j,i) = poli * rsd(j,i)
c               zrsdp(j,i) = poli * rsdp(j,i)
               zrsd(j,l1) = 0.0d0
               zrsdp(j,l1) = 0.0d0
               zrsdt(j,l1) = poli * rsd(j,l1)
               zrsdtp(j,l1) = poli * rsdp(j,l1)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
         off2 = usolvcut * usolvcut
c         print*,"usolvcut in uscale0a",usolvcut
         j = 0
c         do i = 1, npole-1
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst3b,ulst3b,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp,
!$OMP& npole3b,pnum,cnt)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do l1 = 1, npole3b
            !i=pnum(l1)
            !ii = ipole(i)
c            do k = i+1, npole
            m = mindex(l1)
            !do l3 = l1+1, npole3b
            do kkk=1,nulst3b(l1,cnt)
               l3=ulst3b(kkk,l1,cnt)
            !   k=pnum(l3) 
            !   kk = ipole(k)
            !   xr = x(kk) - x(ii)
            !   yr = y(kk) - y(ii)
            !   zr = z(kk) - z(ii)
            !   call image (xr,yr,zr)
            !   r2 = xr*xr + yr* yr + zr*zr
            !   if (r2 .le. off2) then
                  m1 = minv(m+1)
                  m2 = minv(m+2)
                  m3 = minv(m+3)
                  m4 = minv(m+4)
                  m5 = minv(m+5)
                  m6 = minv(m+6)
                  m = m + 6
                 zrsdt(1,l1)= zrsdt(1,l1) + m1*rsd(1,l3) + m2*rsd(2,l3)
     &                           + m3*rsd(3,l3)
                 zrsdt(2,l1)= zrsdt(2,l1) + m2*rsd(1,l3) + m4*rsd(2,l3)
     &                           + m5*rsd(3,l3)
                 zrsdt(3,l1)= zrsdt(3,l1) + m3*rsd(1,l3) + m5*rsd(2,l3)
     &                           + m6*rsd(3,l3)
                 zrsdt(1,l3)= zrsdt(1,l3) + m1*rsd(1,l1) + m2*rsd(2,l1)
     &                           + m3*rsd(3,l1)
                 zrsdt(2,l3)= zrsdt(2,l3) + m2*rsd(1,l1) + m4*rsd(2,l1)
     &                           + m5*rsd(3,l1)
                 zrsdt(3,l3)= zrsdt(3,l3) + m3*rsd(1,l1) + m5*rsd(2,l1)
     &                           + m6*rsd(3,l1)
                zrsdtp(1,l1)= zrsdtp(1,l1) + m1*rsdp(1,l3)+m2*rsdp(2,l3)
     &                            + m3*rsdp(3,l3)
                zrsdtp(2,l1)= zrsdtp(2,l1) + m2*rsdp(1,l3)+m4*rsdp(2,l3)
     &                            + m5*rsdp(3,l3)
                zrsdtp(3,l1)= zrsdtp(3,l1) + m3*rsdp(1,l3)+m5*rsdp(2,l3)
     &                            + m6*rsdp(3,l3)
               zrsdtp(1,l3)= zrsdtp(1,l3) + m1*rsdp(1,l1) +m2*rsdp(2,l1)
     &                            + m3*rsdp(3,l1)
              zrsdtp(2,l3)= zrsdtp(2,l3) + m2*rsdp(1,l1) +m4*rsdp(2,l1)
     &                            + m5*rsdp(3,l1)
              zrsdtp(3,l3)= zrsdtp(3,l3) + m3*rsdp(1,l1) +m5*rsdp(2,l1)
     &                            + m6*rsdp(3,l1)

            !   end if
            end do
         end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
         do i = 1, npole3b
            do j = 1, 3
               zrsd(j,i) = zrsdt(j,i) + zrsd(j,i)
               zrsdp(j,i) = zrsdtp(j,i) + zrsdp(j,i)
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (zrsdt)
         deallocate (zrsdtp)

c
c     construct off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then

         off2 = usolvcut * usolvcut

         m = 0
         do l1 = 1, npole3b
            mindex(l1) =m
            m = m + 6*nulst3b(l1,cnt)
         end do

c
c     perform dynamic allocation of some local arrays
c
         allocate (dscale(npole3b))
c
c     set array needed to scale connected atom interactions
c
         !do i = 1, n
         do i = 1, npole3b
            dscale(i) = 1.0d0
         end do


c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,
!$OMP& thole,polarity,u1scale,u2scale,u3scale,u4scale,np11,ip11,
!$OMP& np12,ip12,np13,ip13,np14,ip14,nulst3b,ulst3b,mindex,minv,
!$OMP& npole3b,off2,pnum,cnt)
!$OMP& firstprivate (dscale)

c
c     determine the off-diagonal elements of the preconditioner
c
!$OMP DO schedule(guided)
         do l1 =1,npole3b
            i=pnum(l1)
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp(i)
            pti = thole(i)
            poli = polarity(i)
            do j = 1,npole3b
               dscale(j) = 1.0d0
            end do
            do j = 1, np11(ii)
               do kk=1,npole3b
                  if (pnum(kk).eq.ip11(j,ii)) then
                     dscale(kk)=u1scale
                     goto 36
                  end if
               end do
   36             continue
            end do
            do j = 1, np12(ii)
               do kk=1,npole3b
                  if (pnum(kk).eq.ip12(j,ii)) then
                     dscale(kk)=u2scale
                     goto 37
                  end if
               end do
   37             continue
            end do
            do j = 1, np13(ii)
               do kk=1,npole3b
                  if (pnum(kk).eq.ip13(j,ii)) then
                     dscale(kk)=u3scale
                     goto 38
                  end if
               end do
   38             continue
            end do
            do j = 1, np14(ii)
               do kk=1,npole3b
                  if (pnum(kk).eq.ip14(j,ii)) then
                     dscale(kk)=u4scale
                     goto 39
                  end if
               end do
   39             continue
            end do
            m=mindex(l1)
            !do l3 =l1+1,npole3b
            do kkk=1,nulst3b(l1,cnt)
               l3=ulst3b(kkk,l1,cnt)
               k=pnum(l3)
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               !if (r2 .le. off2) then
                  r = sqrt(r2)
c                  scale3 = dscale(kk)
c                  scale5 = dscale(kk)
                  scale3 = dscale(l3)
                  scale5 = dscale(l3)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                     end if
                  end if
                  polik = poli * polarity(k)
                  rr3 = scale3 * polik / (r*r2)
                  rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
                  minv(m+1) = rr5*xr*xr - rr3
                  minv(m+2) = rr5*xr*yr
                  minv(m+3) = rr5*xr*zr
                  minv(m+4) = rr5*yr*yr - rr3
                  minv(m+5) = rr5*yr*zr
                  minv(m+6) = rr5*zr*zr - rr3
                  m = m + 6
               !end if
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, npole3b
               dscale(j) = 1.0d0
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL

c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0a  --  dipole preconditioner via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0a" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a double loop
c
c
      subroutine uscale0b_3b_vacomp3btest(npole3b,pnum,mode,rsd,rsdp,
     &  zrsd,zrsdp,cnt)
      use sizes
      use atoms
      use limits
      use mpole
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use neigh
      use neigh2clust
      use usolve
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp
      real*8 pgamma,off2
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
c      real*8 rsd(3,*)
c      real*8 rsdp(3,*)
c      real*8 zrsd(3,*)
c      real*8 zrsdp(3,*)
      character*6 mode
      integer npole3b,pnum(*),l1,l3
      real*8 rsd(3,npole3b)
      real*8 rsdp(3,npole3b)
      real*8 zrsd(3,npole3b)
      real*8 zrsdp(3,npole3b)
      real*8, allocatable :: zrsdt(:,:)
      real*8, allocatable :: zrsdtp(:,:)
      integer cnt,kkk
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then

         allocate (zrsdt(3,npole3b))
         allocate (zrsdtp(3,npole3b))

c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
c         do i = 1, npole
         do l1 = 1, npole3b
            i=pnum(l1)
            poli = udiag * max(polmin,polarity(i))
            do j = 1, 3
c               zrsd(j,i) = poli * rsd(j,i)
c               zrsdp(j,i) = poli * rsdp(j,i)
               zrsd(j,l1) = 0.0d0
               zrsdp(j,l1) = 0.0d0
               zrsdt(j,l1) = poli * rsd(j,l1)
               zrsdtp(j,l1) = poli * rsdp(j,l1)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
         off2 = usolvcut * usolvcut
c         print*,"usolvcut in uscale0a",usolvcut
         j = 0
c         do i = 1, npole-1
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst3b,ulst3b,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp,
!$OMP& npole3b,pnum,cnt)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do l1 = 1, npole3b
            !i=pnum(l1)
            !ii = ipole(i)
c            do k = i+1, npole
            m = mindex(l1)
            !do l3 = l1+1, npole3b
            do kkk=1,nulst3b(l1,cnt)
               l3=ulst3b(kkk,l1,cnt)
            !   k=pnum(l3) 
            !   kk = ipole(k)
            !   xr = x(kk) - x(ii)
            !   yr = y(kk) - y(ii)
            !   zr = z(kk) - z(ii)
            !   call image (xr,yr,zr)
            !   r2 = xr*xr + yr* yr + zr*zr
            !   if (r2 .le. off2) then
                  m1 = minv(m+1)
                  m2 = minv(m+2)
                  m3 = minv(m+3)
                  m4 = minv(m+4)
                  m5 = minv(m+5)
                  m6 = minv(m+6)
                  m = m + 6
                 zrsdt(1,l1)= zrsdt(1,l1) + m1*rsd(1,l3) + m2*rsd(2,l3)
     &                           + m3*rsd(3,l3)
                 zrsdt(2,l1)= zrsdt(2,l1) + m2*rsd(1,l3) + m4*rsd(2,l3)
     &                           + m5*rsd(3,l3)
                 zrsdt(3,l1)= zrsdt(3,l1) + m3*rsd(1,l3) + m5*rsd(2,l3)
     &                           + m6*rsd(3,l3)
                 zrsdt(1,l3)= zrsdt(1,l3) + m1*rsd(1,l1) + m2*rsd(2,l1)
     &                           + m3*rsd(3,l1)
                 zrsdt(2,l3)= zrsdt(2,l3) + m2*rsd(1,l1) + m4*rsd(2,l1)
     &                           + m5*rsd(3,l1)
                 zrsdt(3,l3)= zrsdt(3,l3) + m3*rsd(1,l1) + m5*rsd(2,l1)
     &                           + m6*rsd(3,l1)
                zrsdtp(1,l1)= zrsdtp(1,l1) + m1*rsdp(1,l3)+m2*rsdp(2,l3)
     &                            + m3*rsdp(3,l3)
                zrsdtp(2,l1)= zrsdtp(2,l1) + m2*rsdp(1,l3)+m4*rsdp(2,l3)
     &                            + m5*rsdp(3,l3)
                zrsdtp(3,l1)= zrsdtp(3,l1) + m3*rsdp(1,l3)+m5*rsdp(2,l3)
     &                            + m6*rsdp(3,l3)
               zrsdtp(1,l3)= zrsdtp(1,l3) + m1*rsdp(1,l1) +m2*rsdp(2,l1)
     &                            + m3*rsdp(3,l1)
              zrsdtp(2,l3)= zrsdtp(2,l3) + m2*rsdp(1,l1) +m4*rsdp(2,l1)
     &                            + m5*rsdp(3,l1)
              zrsdtp(3,l3)= zrsdtp(3,l3) + m3*rsdp(1,l1) +m5*rsdp(2,l1)
     &                            + m6*rsdp(3,l1)

            !   end if
            end do
         end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
         do i = 1, npole3b
            do j = 1, 3
               zrsd(j,i) = zrsdt(j,i) + zrsd(j,i)
               zrsdp(j,i) = zrsdtp(j,i) + zrsdp(j,i)
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (zrsdt)
         deallocate (zrsdtp)

c
c     construct off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then

         off2 = usolvcut * usolvcut

         m = 0
         do l1 = 1, npole3b
            mindex(l1) =m
            m = m + 6*nulst3b(l1,cnt)
         end do

c
c     perform dynamic allocation of some local arrays
c
         allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            dscale(i) = 1.0d0
         end do


c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,
!$OMP& thole,polarity,u1scale,u2scale,u3scale,u4scale,np11,ip11,
!$OMP& np12,ip12,np13,ip13,np14,ip14,nulst3b,ulst3b,mindex,minv,
!$OMP& npole3b,off2,pnum,cnt)
!$OMP& firstprivate (dscale)

c
c     determine the off-diagonal elements of the preconditioner
c
!$OMP DO schedule(guided)
         do l1 =1,npole3b
            i=pnum(l1)
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp(i)
            pti = thole(i)
            poli = polarity(i)
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            m=mindex(l1)
            !do l3 =l1+1,npole3b
            do kkk=1,nulst3b(l1,cnt)
               l3=ulst3b(kkk,l1,cnt)
               k=pnum(l3)
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               !if (r2 .le. off2) then
                  r = sqrt(r2)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                     end if
                  end if
                  polik = poli * polarity(k)
                  rr3 = scale3 * polik / (r*r2)
                  rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
                  minv(m+1) = rr5*xr*xr - rr3
                  minv(m+2) = rr5*xr*yr
                  minv(m+3) = rr5*xr*zr
                  minv(m+4) = rr5*yr*yr - rr3
                  minv(m+5) = rr5*yr*zr
                  minv(m+6) = rr5*zr*zr - rr3
                  m = m + 6
               !end if
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
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
!$OMP END PARALLEL

c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
