c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
      subroutine totfieldnpole3b_noewald_umutual_rl_3b_nofield(
     &   M_tot,npole3b,pnum)
      use sizes
      use atoms
      use mpole
      use couple
      use group
      use polgrp
      use polpot
      use polar, only: pdamp, thole
      use shunt
      use cell
      use bound
      implicit none
c      real*8 field(3,*)
      real*8 off3b
c      real*8 fieldp(3,*)
      real*8, allocatable :: dscale_dir(:)
      real*8, allocatable :: dscale_mut(:)
      real*8, allocatable :: pscale(:)
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3_dir,rr5_dir,rr7_dir
      real*8 rr3,rr5
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3_dir,scale5_dir
      real*8 scale7_dir
      real*8 scale3_mut,scale5_mut
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3),M_tot(3*npole3b,3*npole3b)
      real*8 fip(3),fkp(3)
      logical proceed
      character*6 mode
      integer i,ii,j,l1,k,m,kk
      integer l2,l3,k1,k2,i1,i2
      integer pnum(*),npole3b
      real*8 Txx,Txy,Txz,Tyx
      real*8 Tyy,Tyz,Tzx,Tzy,Tzz
      
c      allocate (dscale_dir(npole3b))
      allocate (dscale_mut(npole3b))
c      allocate (pscale(npole3b))
      !off2=1.0d12
c      print*,"off2 in Amatrix routine",off2
c
c     compute the direct induced dipole moment at each atom
c
      do l1 = 1, npole3b-1
         i = pnum(l1)
         i2 = 3*(l1-1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
c         ci = rpole(1,i)
c         dix = rpole(2,i)
c         diy = rpole(3,i)
c         diz = rpole(4,i)
c         qixx = rpole(5,i)
c         qixy = rpole(6,i)
c         qixz = rpole(7,i)
c         qiyy = rpole(9,i)
c         qiyz = rpole(10,i)
c         qizz = rpole(13,i)
         do j = 1,npole3b
c            dscale_dir(j) = 1.0d0
c            pscale(j) = 1.0d0
            dscale_mut(j) = 1.0d0
         end do
c         do j = 1, n12(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i12(j,ii)) then
c                  pscale(kk)=p2scale
c                  goto 31
c               end if
c            end do
c   31             continue
c         end do
c         do j = 1, n13(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i13(j,ii)) then
c                  pscale(kk)=p3scale
c                  goto 32
c               end if
c            end do
c   32             continue
c         end do
c         do j = 1, n14(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i14(j,ii)) then
c                  pscale(kk)=p4scale
c                  goto 33
c               end if
c            end do
c   33             continue
c            do k = 1, np11(ii)
c               do kk=1,npole3b
c                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
c     &              (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(kk) = p4scale * p41scale
c                   goto 34
c                 end if
c               end do
c   34             continue
c            end do
c         end do
c         do j = 1, n15(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i15(j,ii)) then
c                  pscale(kk)=p5scale
c                  goto 35
c               end if
c            end do
c   35             continue
c         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
c                  dscale_dir(kk)=d1scale
                  dscale_mut(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
c                  dscale_dir(kk)=d2scale
                  dscale_mut(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue

         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
c                  dscale_dir(kk)=d3scale
                  dscale_mut(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
c                  dscale_dir(kk)=d4scale
                  dscale_mut(kk)=u4scale
                  goto 39
               end if
            end do
   39             continue
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
c               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               !if (r2 .le. off2) then
                  r = sqrt(r2)
c                  ck = rpole(1,k)
c                  dkx = rpole(2,k)
c                  dky = rpole(3,k)
c                  dkz = rpole(4,k)
c                  qkxx = rpole(5,k)
c                  qkxy = rpole(6,k)
c                  qkxz = rpole(7,k)
c                  qkyy = rpole(9,k)
c                  qkyz = rpole(10,k)
c                  qkzz = rpole(13,k)
c                  scale3_dir = 1.0d0
c                  scale5_dir = 1.0d0
c                  scale7_dir = 1.0d0
                  scale3_mut = dscale_mut(l3)
                  scale5_mut = dscale_mut(l3)

                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
c                        scale3_dir = 1.0d0 - expdamp
c                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
c                        scale7_dir = 1.0d0 - expdamp
c     &                              *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut = scale3_mut * (1.0d0-expdamp)
                        scale5_mut = scale5_mut * (1.0d0-expdamp
     &                                        *(1.0d0-damp))

                     end if
                  end if
c                  rr3_dir = scale3_dir / (r*r2)
c                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
c                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

                  rr3 = scale3_mut / (r*r2)
                  rr5 = 3.0d0 * scale5_mut / (r*r2*r2)


                  Txx = -(-rr3 + xr*xr*rr5)
                  Txy = -(xr*yr*rr5)
                  Txz = -(xr*zr*rr5)
                  Tyx = Txy
                  Tyy = -(-rr3 + yr*yr*rr5)
                  Tyz = -(yr*zr*rr5)
                  Tzx = Txz
                  Tzy = Tyz
                  Tzz = -(-rr3 + zr*zr*rr5)

                  M_tot(i2+1,k2+1) = Txx
                  M_tot(i2+1,k2+2) = Txy
                  M_tot(i2+1,k2+3) = Txz
                  M_tot(i2+2,k2+1) = Tyx
                  M_tot(i2+2,k2+2) = Tyy
                  M_tot(i2+2,k2+3) = Tyz
                  M_tot(i2+3,k2+1) = Tzx
                  M_tot(i2+3,k2+2) = Tzy
                  M_tot(i2+3,k2+3) = Tzz

                  M_tot(k2+1,i2+1) = Txx
                  M_tot(k2+1,i2+2) = Txy
                  M_tot(k2+1,i2+3) = Txz
                  M_tot(k2+2,i2+1) = Tyx
                  M_tot(k2+2,i2+2) = Tyy
                  M_tot(k2+2,i2+3) = Tyz
                  M_tot(k2+3,i2+1) = Tzx
                  M_tot(k2+3,i2+2) = Tzy
                  M_tot(k2+3,i2+3) = Tzz

c                  dir = dix*xr + diy*yr + diz*zr
c                  qix = qixx*xr + qixy*yr + qixz*zr
c                  qiy = qixy*xr + qiyy*yr + qiyz*zr
c                  qiz = qixz*xr + qiyz*yr + qizz*zr
c                  qir = qix*xr + qiy*yr + qiz*zr
c                  dkr = dkx*xr + dky*yr + dkz*zr
c                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
c                  qky = qkxy*xr + qkyy*yr + qkyz*zr
c                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c                  qkr = qkx*xr + qky*yr + qkz*zr
c                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
c     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
c                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
c     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
c                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
c     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
c                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
c     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
c                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
c     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
c                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
c     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
c                  do j = 1, 3
c                     field(j,l1) = field(j,l1) + fid(j)*dscale_dir(l3)
c                     field(j,l3) = field(j,l3) + fkd(j)*dscale_dir(l3)
c                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(l3)
c                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(l3)
c                  end do
               !end if
            end if
         end do
      end do

c      if (use_replica) then
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         i2 = 3*(l1-1)
c         ii = ipole(i)
c         pdi = pdamp(i)
c         pti = thole(i)
c         ci = rpole(1,i)
c         dix = rpole(2,i)
c         diy = rpole(3,i)
c         diz = rpole(4,i)
c         qixx = rpole(5,i)
c         qixy = rpole(6,i)
c         qixz = rpole(7,i)
c         qiyy = rpole(9,i)
c         qiyz = rpole(10,i)
c         qizz = rpole(13,i)
c         do j = 1,npole3b
c            dscale_dir(j) = 1.0d0
c            pscale(j) = 1.0d0
c            dscale_mut(j) = 1.0d0
c         end do
c         do j = 1, n12(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i12(j,ii)) then
c                  pscale(kk)=p2scale
c                  goto 41
c               end if
c            end do
c   41             continue
c         end do
c         do j = 1, n13(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i13(j,ii)) then
c                  pscale(kk)=p3scale
c                  goto 42
c               end if
c            end do
c   42             continue
c         end do
c         do j = 1, n14(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i14(j,ii)) then
c                  pscale(kk)=p4scale
c                  goto 43
c               end if
c            end do
c   43             continue
c            do k = 1, np11(ii)
c               do kk=1,npole3b
c                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
c     &              (pnum(kk).eq.ip11(k,ii)) ) then
c                   pscale(kk) = p4scale * p41scale
c                   goto 44
c                 end if
c               end do
c   44             continue
c            end do
c         end do
c         do j = 1, n15(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.i15(j,ii)) then
c                  pscale(kk)=p5scale
c                  goto 45
c               end if
c            end do
c   45             continue
c         end do
c         do j = 1, np11(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip11(j,ii)) then
c                  dscale_dir(kk)=d1scale
c                  dscale_mut(kk)=u1scale
c                  goto 46
c               end if
c            end do
c   46             continue
c         end do
c         do j = 1, np12(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip12(j,ii)) then
c                  dscale_dir(kk)=d2scale
c                  dscale_mut(kk)=u2scale
c                  goto 47
c               end if
c            end do
c   47             continue
c
c         end do
c         do j = 1, np13(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip13(j,ii)) then
c                  dscale_dir(kk)=d3scale
c                  dscale_mut(kk)=u3scale
c                  goto 48
c               end if
c            end do
c   48             continue
c         end do
c         do j = 1, np14(ii)
c            do kk=1,npole3b
c               if (pnum(kk).eq.ip14(j,ii)) then
c                  dscale_dir(kk)=d4scale
c                  dscale_mut(kk)=u4scale
c                  goto 49
c               end if
c            end do
c   49             continue
c         end do
c         do l3 = l1, npole3b
c            k = pnum(l3)
c            k2 = 3*(l3-1)
c            kk = ipole(k)
c            ck = rpole(1,k)
c            dkx = rpole(2,k)
c            dky = rpole(3,k)
c            dkz = rpole(4,k)
c            qkxx = rpole(5,k)
c            qkxy = rpole(6,k)
c            qkxz = rpole(7,k)
c            qkyy = rpole(9,k)
c            qkyz = rpole(10,k)
c            qkzz = rpole(13,k)
c            do m = 1,ncell
c               xr = x(kk) - x(ii)
c               yr = y(kk) - y(ii)
c               zr = z(kk) - z(ii)
c               call imager (xr,yr,zr,m)
c               r2 = xr*xr + yr* yr + zr*zr
c               if (r2 .le. off2) then
c                  r = sqrt(r2)
c                  scale3_dir = 1.0d0
c                  scale5_dir = 1.0d0
c                  scale7_dir = 1.0d0
c                  scale3_mut = 1.0d0 
c                  scale5_mut = 1.0d0
c
c                  damp = pdi * pdamp(k)
c                  if (damp .ne. 0.0d0) then
c                     pgamma = min(pti,thole(k))
c                     damp = -pgamma * (r/damp)**3
c                     if (damp .gt. -50.0d0) then
c                        expdamp = exp(damp)
c                        scale3_dir = 1.0d0 - expdamp
c                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
c                        scale7_dir = 1.0d0 - expdamp
c     &                              *(1.0d0-damp+0.6d0*damp**2)
c                        scale3_mut =  (1.0d0-expdamp)
c                        scale5_mut =  (1.0d0-expdamp
c     &                                        *(1.0d0-damp))
c
c                     end if
c                  end if
c                  rr3_dir = scale3_dir / (r*r2)
c                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
c                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)
c
c                  rr3 = scale3_mut / (r*r2)
c                  rr5 = 3.0d0 * scale5_mut / (r*r2*r2)
c
c                  if (use_polymer .and. r2 .le. polycut2) then
c                   rr3 = rr3*dscale_mut(l3)
c                   rr5 = rr5*dscale_mut(l3)
c                  end if
c                  Txx = -(-rr3 + xr*xr*rr5)
c                  Txy = -(xr*yr*rr5)
c                  Txz = -(xr*zr*rr5)
c                  Tyx = Txy
c                  Tyy = -(-rr3 + yr*yr*rr5)
c                  Tyz = -(yr*zr*rr5)
c                  Tzx = Txz
c                  Tzy = Tyz
c                  Tzz = -(-rr3 + zr*zr*rr5)
c
c                  if (ii .ne. kk) then
c                  M_tot(i2+1,k2+1) = Txx + M_tot(i2+1,k2+1)
c                  M_tot(i2+1,k2+2) = Txy + M_tot(i2+1,k2+2)
c                  M_tot(i2+1,k2+3) = Txz + M_tot(i2+1,k2+3)
c                  M_tot(i2+2,k2+1) = Tyx + M_tot(i2+2,k2+1)
c                  M_tot(i2+2,k2+2) = Tyy + M_tot(i2+2,k2+2)
c                  M_tot(i2+2,k2+3) = Tyz + M_tot(i2+2,k2+3)
c                  M_tot(i2+3,k2+1) = Tzx + M_tot(i2+3,k2+1)
c                  M_tot(i2+3,k2+2) = Tzy + M_tot(i2+3,k2+2)
c                  M_tot(i2+3,k2+3) = Tzz + M_tot(i2+3,k2+3)
c
c                  M_tot(k2+1,i2+1) = Txx + M_tot(k2+1,i2+1)
c                  M_tot(k2+1,i2+2) = Txy + M_tot(k2+1,i2+2)
c                  M_tot(k2+1,i2+3) = Txz + M_tot(k2+1,i2+3)
c                  M_tot(k2+2,i2+1) = Tyx + M_tot(k2+2,i2+1)
c                  M_tot(k2+2,i2+2) = Tyy + M_tot(k2+2,i2+2)
c                  M_tot(k2+2,i2+3) = Tyz + M_tot(k2+2,i2+3)
c                  M_tot(k2+3,i2+1) = Tzx + M_tot(k2+3,i2+1)
c                  M_tot(k2+3,i2+2) = Tzy + M_tot(k2+3,i2+2)
c                  M_tot(k2+3,i2+3) = Tzz + M_tot(k2+3,i2+3)
c                  end if
c                  dir = dix*xr + diy*yr + diz*zr
c                  qix = qixx*xr + qixy*yr + qixz*zr
c                  qiy = qixy*xr + qiyy*yr + qiyz*zr
c                  qiz = qixz*xr + qiyz*yr + qizz*zr
c                  qir = qix*xr + qiy*yr + qiz*zr
c                  dkr = dkx*xr + dky*yr + dkz*zr
c                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
c                  qky = qkxy*xr + qkyy*yr + qkyz*zr
c                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c                  qkr = qkx*xr + qky*yr + qkz*zr
c                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
c     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
c                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
c     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
c                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
c     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
c                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
c     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
c                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
c     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
c                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
c     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
c                  do j = 1, 3
c                     fip(j) = fid(j)
c                     fkp(j) = fkd(j)
c                  end do
c                  if (use_polymer .and. r2 .le. polycut2) then
c                        do j = 1, 3
c                           fid(j) = fid(j) * dscale_dir(l3)
c                           fip(j) = fip(j) * pscale(l3)
c                           fkd(j) = fkd(j) * dscale_dir(l3)
c                           fkp(j) = fkp(j) * pscale(l3)
c                        end do
c                  end if
c                  do j = 1, 3
c                     field(j,l1) = field(j,l1) + fid(j)
c                     fieldp(j,l1) = fieldp(j,l1) + fip(j)
c                     if (ii .ne. kk) then
c                           field(j,l3) = field(j,l3) + fkd(j)
c                           fieldp(j,l3) = fieldp(j,l3) + fkp(j)
c                     end if
c                  end do
c               end if
c            end do
c         end do
c      end do
c
c      end if

c      deallocate (dscale_dir)
      deallocate (dscale_mut)
c      deallocate (pscale)
      return
      end



