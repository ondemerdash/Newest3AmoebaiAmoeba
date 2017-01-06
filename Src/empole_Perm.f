c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole0d  --  Ewald multipole energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole0d" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and a neighbor list
c
c
      subroutine empole0d_3b_Perm_selfeng_bcast
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use polar
      implicit none
      integer i,ii
      real*8 e,ei,f
      real*8 term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      !em = 0.0d0
      !ep = 0.0d0
      !if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec

      em=em+emreal_tmp
c
c     check the sign of multipole components at chiral sites
c
      !call chkpole
c
c     rotate the multipole components into the global frame
c
      !call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      !call induce
c
c     compute the reciprocal space part of the Ewald summation
c
      !call emrecip
c
c     compute the real space portion of the Ewald summation
c
      !call ereal0d
c
c     compute the self-energy portion of the Ewald summation
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
         !uix = uind(1,i)
         !uiy = uind(2,i)
         !uiz = uind(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         !ei = fterm * term * uii / 3.0d0
         em = em + e
         !ep = ep + ei
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ereal0d  --  real space mpole energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal0d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipoles and dipole polarizability
c     using a neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine LoadBalInnerloop_ereal0d_3b_Perm_sentlist2(
     & emrealt)
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use couple
      use energi
      use ewald
      use math
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use paremneigh
      use mpidat
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 e,ei,f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 emtt,eptt
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      character*6 mode
      integer iter3
      real*8 emrealt
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      !allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
       !  pscale(i) = 1.0d0
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
      emtt = 0.0d0
      !eptt = 0.0d0
c
c     set OpenMP directives for the major loop structure
c
c!$OMP PARALLEL default(shared) firstprivate(f)
c!$OMP& private(i,j,k,ii,kk,kkk,e,ei,efix,eifix,bfac,damp,expdamp,
c!$OMP& pdi,pti,pgamma,scale3,scale5,scale7,alsq2,alsq2n,
c!$OMP& exp2a,ralpha,xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,rr9,
c!$OMP& ci,dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz,uix,uiy,uiz,
c!$OMP& ck,dkx,dky,dkz,qkxx,qkxy,qkxz,qkyy,qkyz,qkzz,ukx,uky,ukz,
c!$OMP& bn,sc,gl,sci,gli)
c!$OMP& firstprivate(mscale,pscale)
c!$OMP DO reduction(+:emtt,eptt)
c!$OMP& schedule(guided)

!$OMP PARALLEL default(shared) firstprivate(f)
!$OMP& private(i,j,k,ii,kk,kkk,e,ei,efix,eifix,bfac,damp,expdamp,
!$OMP& pdi,pti,pgamma,scale3,scale5,scale7,alsq2,alsq2n,
!$OMP& exp2a,ralpha,xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,rr9,
!$OMP& ci,dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz,uix,uiy,uiz,
!$OMP& ck,dkx,dky,dkz,qkxx,qkxy,qkxz,qkyy,qkyz,qkzz,ukx,uky,ukz,
!$OMP& bn,sc,gl,sci,gli,iter3)
!$OMP& firstprivate(mscale)
!$OMP DO reduction(+:emtt)
!$OMP& schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      !do i = 1, npole
      do i=start_emreal2(taskid),last_emreal2(taskid)
         iter3=i-start_emreal2(taskid)+1

         ii = ipole(i)
       !  pdi = pdamp(i)
       !  pti = thole(i)
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
         !uix = uind(1,i)
         !uiy = uind(2,i)
         !uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         !   pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         !   pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         !   pscale(i14(j,ii)) = p4scale
         !   do k = 1, np11(ii)
         !       if (i14(j,ii) .eq. ip11(k,ii))
     &   !         pscale(i14(j,ii)) = p4scale * p41scale
         !   end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         !   pscale(i15(j,ii)) = p5scale
         end do
         !do kkk = 1, nelst(i)
         !   k = elst(kkk,i)
         do kkk = 1, nelst_recv(iter3)
            k = elst_recv(kkk,iter3)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
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
            !   ukx = uind(1,k)
            !   uky = uind(2,k)
            !   ukz = uind(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
            !   sci(2) = uix*dkx + dix*ukx + uiy*dky
     &      !               + diy*uky + uiz*dkz + diz*ukz
            !   sci(3) = uix*xr + uiy*yr + uiz*zr
            !   sci(4) = ukx*xr + uky*yr + ukz*zr
            !   sci(7) = qix*ukx + qiy*uky + qiz*ukz
            !   sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
             !  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
             !  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &       !              - sc(3)*sci(4)
             !  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                + gl(3)*bn(3) + gl(4)*bn(4)
           !    ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needeed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
           !    scale3 = pscale(kk)
           !    scale5 = pscale(kk)
           !    scale7 = pscale(kk)
           !    damp = pdi * pdamp(k)
           !    if (damp .ne. 0.0d0) then
           !       pgamma = min(pti,thole(k))
           !       damp = -pgamma * (r/damp)**3
           !       if (damp .gt. -50.0d0) then
           !          expdamp = exp(damp)
           !          scale3 = scale3 * (1.0d0-expdamp)
           !          scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
           !          scale7 = scale7 * (1.0d0-(1.0d0-damp+0.6d0*damp**2)
     &     !                                          *expdamp)
           !       end if
           !    end if
               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
           !    eifix = gli(1)*rr3*(1.0d0-scale3)
     &     !               + gli(2)*rr5*(1.0d0-scale5)
     &     !               + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
               e = e - efix*(1.0d0-mscale(kk))
            !   ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
               e = f * e
            !   ei = 0.5d0 * f * ei
               emtt = emtt + e
            !   eptt = eptt + ei
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         !   pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         !   pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         !   pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         !   pscale(i15(j,ii)) = 1.0d0
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
      !em = em + emtt
      emrealt=emtt
      !ep = ep + eptt
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      !deallocate (pscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emrecip  --  PME recip space multipole energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions and
c     dipole polarizability
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip_Perm
      use sizes
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use potent
      implicit none
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 a(3,3)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      !allocate (fuind(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (fphi(20,npole))
c
c     copy the multipole moments into local storage areas
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
c     compute B-spline coefficients and spatial decomposition
c
c      if (.not. use_polar) then
         call bspline_fill
         call table_fill
c      end if
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp)
      call fftfront
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
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the charge grid
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
      call fphi_mpole (fphi)
c
c     sum over multipoles and increment total multipole energy
c
      e = 0.0d0
      do i = 1, npole
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
         end do
      end do
      e = 0.5d0 * electric * e
      em = em + e
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      !if (use_polar) then
      !   do i = 1, 3
      !      a(1,i) = dble(nfft1) * recip(i,1)
      !      a(2,i) = dble(nfft2) * recip(i,2)
      !      a(3,i) = dble(nfft3) * recip(i,3)
      !   end do
      !   do i = 1, npole
      !      do k = 1, 3
      !         fuind(k,i) = a(k,1)*uind(1,i) + a(k,2)*uind(2,i)
     &!                         + a(k,3)*uind(3,i)
      !      end do
      !   end do
c
c     sum over induced dipoles and increment total induced energy
c
       !  e = 0.0d0
       !  do i = 1, npole
       !     do k = 1, 3
       !        e = e + fuind(k,i)*fphi(k+1,i)
       !     end do
       !  end do
       !  e = 0.5d0 * electric * e
       !  ep = ep + e
      !end if
c
c     perform deallocation of some local arrays
c
      !deallocate (fuind)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
      return
      end
