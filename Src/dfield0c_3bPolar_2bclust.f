      subroutine dfield0c_3b_vacomp2b(npole3b,pnum,field,fieldp,cnt)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar, only: polarity, thole, pdamp
      use neigh2clust
      implicit none
      integer i,j,ii
      real*8 term
      real*8 ucell(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      integer npole3b,pnum(*),cnt,l1

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
c     get the reciprocal space part of the electrostatic field
c
      if(useboxclust) then
       !call udirect1_3b_newbox2bclust (npole3b,pnum,field,cnt)
      else
       call udirect1_3b_new (npole3b,pnum,field)
      end if

      do i = 1, npole3b
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do

c
c     get the real space portion of the electrostatic field
c
        if(useboxclust) then
        ! call udirect2b_3b_vacomp2btestbox2bclust(npole3b,pnum,field,
     &  ! fieldp,cnt)
        else
         call udirect2b_3b_vacomp2btest(npole3b,pnum,field,fieldp,
     &  cnt) 
        end if
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
         i=pnum(l1)
         do j = 1, 3
            field(j,l1) = field(j,l1) + term*rpole(j+1,i)
            fieldp(j,l1) = fieldp(j,l1) + term*rpole(j+1,i)
         end do
      end do

      return
      end


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
      subroutine udirect2b_3b_vacomp2b(npole3b,pnum,field,fieldp,
     &  cnt)
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
         do kkk=1,nelst2b(l1,cnt)
            l3=elst2b(kkk,l1,cnt)
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
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(l3))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(l3))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(l3))*rr7
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
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(l3))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(l3))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(l3))*rr7
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
                  bcn(1) = bn(1) - (1.0d0-scale3*uscale(l3))*rr3
                  bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(l3))*rr5
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c
      subroutine udirect1_3b_new(npole3b,pnum,field)
      use sizes
      use bound
      use boxes
      use ewald
      use math
      use mpole
      use pme
      implicit none
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 field(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
      integer npole3b,pnum(*),l1
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole3b))
      allocate (fmp(10,npole3b))
      allocate (cphi(10,npole3b))
      allocate (fphi(20,npole3b))
c
c     copy multipole moments and coordinates to local storage
c
      !do i = 1, npole
      do l1=1,npole3b
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
c     compute B-spline coefficients and spatial decomposition
c
C   COMMENTED OUT BELOW BECAUSE IT IS CALCULATED PREVIOUSLY ON ALL TASK 
      !call bspline_fill
      !print*,"Completed bspline_fill"
      !call table_fill
      !print*,"Completed table_fill"
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp3b (npole3b,pnum,cmp,fmp)
      !print*,"Completed cmp_to_fmp3b"
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole3b_new (npole3b,pnum,fmp)
      !print*,"Completed grid_mpole3b_new"
      call fftfront
      !print*,"Completed fftfront"
c
c     make the scalar summation over reciprocal lattice
c
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      ntot = nfft1 * nfft2 * nfft3
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
c     perform 3-D FFT backward transform and get field
c
      call fftback
      !print*,"Done w fftback"
      call fphi_mpole3b_new (npole3b,pnum,fphi)
      !print*,"Done w fphi_mpole3b_new"
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi3b (npole3b,pnum,fphi,cphi)
      !print*,"Done w fphi_to_cphi3b"

c
c     increment the field at each multipole site
c
      do i = 1, npole3b
         field(1,i) = field(1,i) - cphi(2,i)
         field(2,i) = field(2,i) - cphi(3,i)
         field(3,i) = field(3,i) - cphi(4,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end

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
      subroutine udirect2b_3b_vacomp2btest(npole3b,pnum,field,fieldp,
     &  cnt)
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

                  !dlocal(1,nlocal) = (-rr3 + rr5*xr*xr)*uscale(kk)
                  !dlocal(2,nlocal) = (rr5*xr*yr)*uscale(kk)
                  !dlocal(3,nlocal) = (rr5*xr*zr)*uscale(kk)
                  !dlocal(4,nlocal) = (-rr3 + rr5*yr*yr)*uscale(kk)
                  !dlocal(5,nlocal) = (rr5*yr*zr)*uscale(kk)
                  !dlocal(6,nlocal) = (-rr3 + rr5*zr*zr)*uscale(kk)
               end if

                  do j = 1, 3
c                     field(j,l1) = field(j,l1) + fid(j)*dscale(kk)
c                     field(j,l3) = field(j,l3) + fkd(j)*dscale(kk)
c                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(kk)
c                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(kk)

c                     field(j,l1) = field(j,l1) + fid(j)*dscale(kk)
c                     field(j,l3) = field(j,l3) + fkd(j)*dscale(kk)
c                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(kk)
c                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(kk)

c                     fieldt(j,l1) = fieldt(j,l1) + fid(j)*dscale(kk)
c                     fieldt(j,l3) = fieldt(j,l3) + fkd(j)*dscale(kk)
c                     fieldtp(j,l1) = fieldtp(j,l1) + fid(j)*pscale(kk)
c                     fieldtp(j,l3) = fieldtp(j,l3) + fkd(j)*pscale(kk)

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
