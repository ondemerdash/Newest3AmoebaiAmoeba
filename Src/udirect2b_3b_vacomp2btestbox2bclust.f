c
c
c
c
      subroutine udirect2b_3b_vacomp2btestbox2bclust(npole3b,pnum,field,
     &  fieldp,cnt)
      use sizes
      use atoms
      use bound
      use boxes
      use ewald
      use math
      use cell
      use units
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
      real*8 rr1,rr2,rr3,rr5,rr7
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
      real*8 erfc,bfac,exp2a
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 bn(0:3),bcn(3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      external erfc

c
c
c     zero out the value of the field at each site
c
      !do i = 1, npole3b
      !   do j = 1, 3
      !      field(j,i) = 0.0d0
      !      fieldp(j,i) = 0.0d0
      !   end do
      !end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switchbox2bclust (mode,cnt)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
      nlocal = 0
      toffset0 = 0
      maxlocal = int(dble(npole3b)*dble(maxelst)/dble(nthread))

c      print*,"off2 in dfield0a_3b",off2 
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (toffset(0:nthread-1))
      allocate (uscale(n))
      allocate (fieldt(3,npole3b))
      allocate (fieldtp(3,npole3b))

         do j = 1,n
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
!$OMP& n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,nelst2b,
!$OMP& elst2b,cut2,poltyp,ntpair,tindex,tdipdip,
!$OMP& toffset,toffset0,field,fieldp,fieldt,fieldtp,maxlocal,
!$OMP& npole3b,pnum,off2,cnt,aewald,aesq2,aesq2n)
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
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
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

         !do l3=l1+1,npole3b
         do kkk=1,nelst2b(l1,cnt)
            l3=elst2b(kkk,l1,cnt)
            k=pnum(l3)
            kk = ipole(k)
            !proceed = .true.
          !  if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            !if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               !if (use_bounds)  call image (xr,yr,zr)
               call image2b(xr,yr,zr,cnt)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rr1 = 1.0d0 / r
                  rr2 = rr1 * rr1
                  rr3 = rr2 * rr1
                  rr5 = rr2 * rr3
                  rr7 = rr2 * rr5
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
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
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
c                  rr3 = scale3 / (r*r2)
c                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
c                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
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
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(kk))*rr7
               fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(kk))*rr7
               fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz


c                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
c     &                        - rr3*dkx + 2.0d0*rr5*qkx
c                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
c     &                        - rr3*dky + 2.0d0*rr5*qky
c                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
c     &                        - rr3*dkz + 2.0d0*rr5*qkz
c                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
c     &                        - rr3*dix - 2.0d0*rr5*qix
c                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
c     &                        - rr3*diy - 2.0d0*rr5*qiy
c                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
c     &                        - rr3*diz - 2.0d0*rr5*qiz

               if (poltyp .eq. 'MUTUAL') then
                !TERMS CONTAINING BCN COMMENTED OUT BELOW.  NO EWALDKIRSCHTORTE FOR US!
                  bcn(1) = bn(1) - (1.0d0-scale3*uscale(kk))*rr3
                  bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(kk))*rr5
                  nlocal = nlocal + 1
                  ilocal(1,nlocal) = l1
                  ilocal(2,nlocal) = l3
                  dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                  dlocal(2,nlocal) = bcn(2)*xr*yr
                  dlocal(3,nlocal) = bcn(2)*xr*zr
                  dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                  dlocal(5,nlocal) = bcn(2)*yr*zr
                  dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr

                  !dlocal(1,nlocal) = (-rr3 + rr5*xr*xr)*uscale(l3)
                  !dlocal(2,nlocal) = (rr5*xr*yr)*uscale(l3)
                  !dlocal(3,nlocal) = (rr5*xr*zr)*uscale(l3)
                  !dlocal(4,nlocal) = (-rr3 + rr5*yr*yr)*uscale(l3)
                  !dlocal(5,nlocal) = (rr5*yr*zr)*uscale(l3)
                  !dlocal(6,nlocal) = (-rr3 + rr5*zr*zr)*uscale(l3)
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

c                     fieldt(j,l1) = fieldt(j,l1) + fid(j)*dscale(l3)
c                     fieldt(j,l3) = fieldt(j,l3) + fkd(j)*dscale(l3)
c                     fieldtp(j,l1) = fieldtp(j,l1) + fid(j)*pscale(l3)
c                     fieldtp(j,l3) = fieldtp(j,l3) + fkd(j)*pscale(l3)

                     fieldt(j,l1) = fieldt(j,l1) + fimd(j)
                     fieldt(j,l3) = fieldt(j,l3) + fkmd(j)
                     fieldtp(j,l1) = fieldtp(j,l1) + fimp(j)
                     fieldtp(j,l3) = fieldtp(j,l3) + fkmp(j)

                  end do
               end if
            !end if
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
            uscale(ip11(j,ii)) = 1.0d0
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = 1.0d0
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = 1.0d0
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = 1.0d0
            dscale(ip14(j,ii)) = 1.0d0
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

c     #################################################################
c     ##                                                             ##
c
      subroutine udirect2a_3b_vacomptestbox2bclust(npole3b,pnum,field,
     &  fieldp,cnt)
      use sizes
      use atoms
      use bound
      use boxes
      use ewald
      use math
      use cell
      use units
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
      real*8 rr1,rr2,rr3,rr5,rr7
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
      real*8 erfc,bfac,exp2a
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 bn(0:3),bcn(3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      external erfc

c
c
c     zero out the value of the field at each site
c
      !do i = 1, npole3b
      !   do j = 1, 3
      !      field(j,i) = 0.0d0
      !      fieldp(j,i) = 0.0d0
      !   end do
      !end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
      nlocal = 0
      toffset0 = 0
      maxlocal = int(dble(npole3b)*dble(maxelst)/dble(nthread))

c      print*,"off2 in dfield0a_3b",off2 
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (toffset(0:nthread-1))
      allocate (uscale(n))
      allocate (fieldt(3,npole3b))
      allocate (fieldtp(3,npole3b))

         do j = 1,n
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
!$OMP& n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,nelst2b,
!$OMP& elst2b,cut2,poltyp,ntpair,tindex,tdipdip,
!$OMP& toffset,toffset0,field,fieldp,fieldt,fieldtp,maxlocal,
!$OMP& npole3b,pnum,off2,cnt,aewald,aesq2,aesq2n)
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

      do l1 = 1, npole3b-1
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
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
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

         do l3=l1+1,npole3b
         !do kkk=1,nelst2b(l1,cnt)
         !   l3=elst2b(kkk,l1,cnt)
            k=pnum(l3)
            kk = ipole(k)
            proceed = .true.
          !  if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            !if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               !if (use_bounds)  call image (xr,yr,zr)
               call image(xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rr1 = 1.0d0 / r
                  rr2 = rr1 * rr1
                  rr3 = rr2 * rr1
                  rr5 = rr2 * rr3
                  rr7 = rr2 * rr5
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
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
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
c                  rr3 = scale3 / (r*r2)
c                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
c                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
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
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(kk))*rr7
               fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(kk))*rr7
               fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz


c                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
c     &                        - rr3*dkx + 2.0d0*rr5*qkx
c                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
c     &                        - rr3*dky + 2.0d0*rr5*qky
c                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
c     &                        - rr3*dkz + 2.0d0*rr5*qkz
c                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
c     &                        - rr3*dix - 2.0d0*rr5*qix
c                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
c     &                        - rr3*diy - 2.0d0*rr5*qiy
c                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
c     &                        - rr3*diz - 2.0d0*rr5*qiz

               if (poltyp .eq. 'MUTUAL') then
                !TERMS CONTAINING BCN COMMENTED OUT BELOW.  NO EWALDKIRSCHTORTE FOR US!
                  bcn(1) = bn(1) - (1.0d0-scale3*uscale(kk))*rr3
                  bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(kk))*rr5
                  nlocal = nlocal + 1
                  ilocal(1,nlocal) = l1
                  ilocal(2,nlocal) = l3
                  dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                  dlocal(2,nlocal) = bcn(2)*xr*yr
                  dlocal(3,nlocal) = bcn(2)*xr*zr
                  dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                  dlocal(5,nlocal) = bcn(2)*yr*zr
                  dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr

                  !dlocal(1,nlocal) = (-rr3 + rr5*xr*xr)*uscale(l3)
                  !dlocal(2,nlocal) = (rr5*xr*yr)*uscale(l3)
                  !dlocal(3,nlocal) = (rr5*xr*zr)*uscale(l3)
                  !dlocal(4,nlocal) = (-rr3 + rr5*yr*yr)*uscale(l3)
                  !dlocal(5,nlocal) = (rr5*yr*zr)*uscale(l3)
                  !dlocal(6,nlocal) = (-rr3 + rr5*zr*zr)*uscale(l3)
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

c                     fieldt(j,l1) = fieldt(j,l1) + fid(j)*dscale(l3)
c                     fieldt(j,l3) = fieldt(j,l3) + fkd(j)*dscale(l3)
c                     fieldtp(j,l1) = fieldtp(j,l1) + fid(j)*pscale(l3)
c                     fieldtp(j,l3) = fieldtp(j,l3) + fkd(j)*pscale(l3)

                     fieldt(j,l1) = fieldt(j,l1) + fimd(j)
                     fieldt(j,l3) = fieldt(j,l3) + fkmd(j)
                     fieldtp(j,l1) = fieldtp(j,l1) + fimp(j)
                     fieldtp(j,l3) = fieldtp(j,l3) + fkmp(j)

                  end do
               end if
            !end if
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
            uscale(ip11(j,ii)) = 1.0d0
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = 1.0d0
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = 1.0d0
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = 1.0d0
            dscale(ip14(j,ii)) = 1.0d0
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
