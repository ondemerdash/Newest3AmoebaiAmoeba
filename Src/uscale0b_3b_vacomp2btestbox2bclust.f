c
c
c
      subroutine uscale0b_3b_vacomp2btestbox2bclust(npole3b,pnum,mode,
     &  rsd,rsdp,zrsd,zrsdp,cnt)
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
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst2b,ulst2b,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp,
!$OMP& npole3b,pnum,cnt)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do l1 = 1, npole3b
            !i=pnum(l1)
            !ii = ipole(i)
c            do k = i+1, npole
            m = mindex(l1)
            !do l3 = l1+1, npole3b
            do kkk=1,nulst2b(l1,cnt)
               l3=ulst2b(kkk,l1,cnt)
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
            m = m + 6*nulst2b(l1,cnt)
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
!$OMP& np12,ip12,np13,ip13,np14,ip14,nulst2b,ulst2b,mindex,minv,
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
            do kkk=1,nulst2b(l1,cnt)
               l3=ulst2b(kkk,l1,cnt)
               k=pnum(l3)
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image2b (xr,yr,zr,cnt)
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
