c
c
      subroutine empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(
     &   npole3b,pnum,eptemp,deptemp2,virtemp)
      use bound
      use sizes
      use atoms
      use chgpot
      use couple
      use group
      use mpole
      use polar, only: polarity, thole, pdamp
      use polgrp
      use polpot
      use usage
      use cell
      use boxes
      use shunt
      use totfield
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,fgrp,gfd
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale3i
      real*8 scale5,scale5i
      real*8 scale7,scale7i
      real*8 temp3,temp5,temp7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz,qiyx,qiyy,qiyz,qizx,qizy,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz,qkyx,qkyy,qkyz,qkzx,qkzy,qkzz
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 fdir(3)
      real*8 dixuk(3),dkxui(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qirx,qiry,qirz,qkrx,qkry,qkrz
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqkxz,rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 gli(7),glip(7)
      real*8 sc(10),sci(8),scip(8)
      real*8 gf(7),gfi(6),gti(6)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 eptemp,deptemp2(3,npole),virtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3,k1
      logical proceed,usei,usek,doi
      character*6 mode
      real*8, allocatable :: deptemp(:,:)
      eptemp = 0.0d0
      eptemp_userep =0.0d0
c
c   Zero out temporary gradient of polarization energy
c
      allocate (deptemp(3,npole3b))

      do l1 = 1, npole3b
        do j = 1, 3
          deptemp(j,l1) = 0.0d0
        end do
      end do

      do i=1,npole
         do j = 1,3
            deptemp2(j,i)=0.0d0
c            uinp2(j,i)=0.0d0
         end do
      end do
c
c   Zero out temporary virial
c
      do i=1,3
         do j=1,3
           virtemp(i,j)=0.0d0
         end do
      end do

      mode = 'MPOLE'
      call switch (mode)

c      print*,"off2,use_replica = ",off2,use_replica
       
      call induce0a_3b_PolelecOnly_totfieldnpole3b(npole3b,pnum,uind,
     &     uinp)

c      if(npole3b.eq.3) then
c       do l1=1,npole3b
c         i=pnum(l1)
c         print*,"npole3b i TotFielduindx",npole3b,i,uind(1,l1)
c         print*,"npole3b i TotFielduindy",npole3b,i,uind(2,l1)
c         print*,"npole3b i TotFielduindz",npole3b,i,uind(3,l1)
c       end do
c       do l1=1,npole3b
c         i=pnum(l1)
c         print*,"npole3b i smallfielduinpx",npole3b,i,uinp(1,l1)
c         print*,"npole3b i smallfielduinpy",npole3b,i,uinp(2,l1)
c         print*,"npole3b i smallfielduinpz",npole3b,i,uinp(3,l1)
c       end do
c
c      end if

c       print*,"after induce"
      call empole1a1_3b_Polar_totfieldnpole3b(npole3b,pnum,uind,
     &     uinp,eptemp,deptemp,virtemp)


       call empole1a2_3b_Polar_totfieldnpole3bmod_dol1_2(npole3b,pnum,
     &     uinp,deptemp2,virtemp)

      

      do l1 =1,npole3b
          i =pnum(l1) 
          deptemp2(1,i)=deptemp(1,l1)+deptemp2(1,i)
          deptemp2(2,i)=deptemp(2,l1)+deptemp2(2,i)
          deptemp2(3,i)=deptemp(3,l1)+deptemp2(3,i)
      end do
      deallocate(deptemp)

      return
      end


      subroutine induce0a_3b_PolelecOnly_totfieldnpole3b(npole3b,pnum,
     &           uind,uinp)
      use sizes
      use atoms
      use polar, only: polarity, thole, pdamp
      use totfield
      use polpot
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 field(3,npole3b)
      real*8 fieldp(3,npole3b)
      character*6 mode
      integer l1,l2,l3,k1,k2,i1,i2
      real*8 M_tot(3*npole3b,3*npole3b)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      integer npole3b,pnum(*)

c        print*,"npole3b=",npole3b
        do l1=1,npole3b
           do j=1,3
              field(j,l1)=0.0d0
              fieldp(j,l1)=0.0d0
           end do
        end do

         do l1 = 1, 3*npole3b
            do l3 = 1, 3*npole3b
               M_tot(l1,l3) = 0.0d0
            end do
         end do

         do l1 = 1, npole3b
            do j = 1, 3
               uind(j,l1) = 0.0d0
               uinp(j,l1) = 0.0d0
            end do
         end do

       
c       call totfield_noewald_umutual_rl_3b (field,fieldp,M_tot,
c     &   npole3b,pnum)

c       call dfield0a_totfieldnpole3b (field,npole3b,pnum)

c       call dfield0a_totfieldnpole3bmod (field,npole3b,pnum)

       call totfieldnpole3b_noewald_umutual_rl_3b (field,fieldp,
     &   M_tot,npole3b,pnum)
c       print*,"poltyp=",poltyp 
      if (poltyp .eq. 'MUTUAL') then 
         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               i1 = 3*(l1-1)+j
               M_tot(i1,i1) = M_tot(i1,i1)+1.0d0/polarity(i)
            end do
c             print*,"polarity(i)",polarity(i)
         end do

         call invert(3*npole3b,M_tot)

         do l1 = 1, npole3b
            i = pnum(l1)
            do i1 = 1, 3
               i2 = 3*(l1-1) + i1
               do l3 = 1, npole3b
                  k = pnum(l3)
                  k2 = 3*(l3-1)
c                  uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*field(1,l3)+
c     &                                      M_tot(i2,k2+2)*field(2,l3)+
c     &                                      M_tot(i2,k2+3)*field(3,l3)
               uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*fieldnpole(1,k)+
     &                                  M_tot(i2,k2+2)*fieldnpole(2,k)+
     &                                  M_tot(i2,k2+3)*fieldnpole(3,k)
                  uinp(i1,l1)=uinp(i1,l1)+ M_tot(i2,k2+1)*fieldp(1,l3)+
     &                                     M_tot(i2,k2+2)*fieldp(2,l3)+
     &                                     M_tot(i2,k2+3)*fieldp(3,l3)
c               uinp(i1,l1)=uinp(i1,l1)+M_tot(i2,k2+1)*fieldpnpole(1,k)+
c     &                                 M_tot(i2,k2+2)*fieldpnpole(2,k)+
c     &                                 M_tot(i2,k2+3)*fieldpnpole(3,k)
               end do
            end do
         end do
      else
         do l1 = 1, npole3b
            i = pnum(l1)
            do j=1,3
               uind(j,l1)=polarity(i)*fieldnpole(j,i)
               uinp(j,l1)=polarity(i)*fieldp(j,l1)
            end do
         end do        
      end if
       
      return
      end
c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
      subroutine totfieldnpole3b_noewald_umutual_rl_3b(field,fieldp,
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
      real*8 field(3,*)
      real*8 off3b
      real*8 fieldp(3,*)
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
      
      allocate (dscale_dir(npole3b))
      allocate (dscale_mut(npole3b))
      allocate (pscale(npole3b))

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
         do j = 1,npole3b
            dscale_dir(j) = 1.0d0
            pscale(j) = 1.0d0
            dscale_mut(j) = 1.0d0
         end do
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
                  dscale_dir(kk)=d1scale
                  dscale_mut(kk)=u1scale
                  goto 36
               end if
            end do
   36             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale_dir(kk)=d2scale
                  dscale_mut(kk)=u2scale
                  goto 37
               end if
            end do
   37             continue

         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale_dir(kk)=d3scale
                  dscale_mut(kk)=u3scale
                  goto 38
               end if
            end do
   38             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale_dir(kk)=d4scale
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
               !call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               !if (r2 .le. off2) then
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
                  scale3_dir = 1.0d0
                  scale5_dir = 1.0d0
                  scale7_dir = 1.0d0
                  scale3_mut = dscale_mut(l3)
                  scale5_mut = dscale_mut(l3)

                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3_dir = 1.0d0 - expdamp
                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7_dir = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut = scale3_mut * (1.0d0-expdamp)
                        scale5_mut = scale5_mut * (1.0d0-expdamp
     &                                        *(1.0d0-damp))

                     end if
                  end if
                  rr3_dir = scale3_dir / (r*r2)
                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

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
                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
                  do j = 1, 3
                     field(j,l1) = field(j,l1) + fid(j)*dscale_dir(l3)
                     field(j,l3) = field(j,l3) + fkd(j)*dscale_dir(l3)
                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(l3)
                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(l3)
                  end do
               !end if
            end if
         end do
      end do
      use_replica =.false.
      if (use_replica) then

      do l1 = 1, npole3b
         i = pnum(l1)
         i2 = 3*(l1-1)
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
         do j = 1,npole3b
            dscale_dir(j) = 1.0d0
            pscale(j) = 1.0d0
            dscale_mut(j) = 1.0d0
         end do
         do j = 1, n12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i12(j,ii)) then
                  pscale(kk)=p2scale
                  goto 41
               end if
            end do
   41             continue
         end do
         do j = 1, n13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i13(j,ii)) then
                  pscale(kk)=p3scale
                  goto 42
               end if
            end do
   42             continue
         end do
         do j = 1, n14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i14(j,ii)) then
                  pscale(kk)=p4scale
                  goto 43
               end if
            end do
   43             continue
            do k = 1, np11(ii)
               do kk=1,npole3b
                 if( (i14(j,ii) .eq. ip11(k,ii)).and.
     &              (pnum(kk).eq.ip11(k,ii)) ) then
                   pscale(kk) = p4scale * p41scale
                   goto 44
                 end if
               end do
   44             continue
            end do
         end do
         do j = 1, n15(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.i15(j,ii)) then
                  pscale(kk)=p5scale
                  goto 45
               end if
            end do
   45             continue
         end do
         do j = 1, np11(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip11(j,ii)) then
                  dscale_dir(kk)=d1scale
                  dscale_mut(kk)=u1scale
                  goto 46
               end if
            end do
   46             continue
         end do
         do j = 1, np12(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip12(j,ii)) then
                  dscale_dir(kk)=d2scale
                  dscale_mut(kk)=u2scale
                  goto 47
               end if
            end do
   47             continue

         end do
         do j = 1, np13(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip13(j,ii)) then
                  dscale_dir(kk)=d3scale
                  dscale_mut(kk)=u3scale
                  goto 48
               end if
            end do
   48             continue
         end do
         do j = 1, np14(ii)
            do kk=1,npole3b
               if (pnum(kk).eq.ip14(j,ii)) then
                  dscale_dir(kk)=d4scale
                  dscale_mut(kk)=u4scale
                  goto 49
               end if
            end do
   49             continue
         end do
         do l3 = l1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
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
            do m = 1,ncell
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call imager (xr,yr,zr,m)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  scale3_dir = 1.0d0
                  scale5_dir = 1.0d0
                  scale7_dir = 1.0d0
                  scale3_mut = 1.0d0 
                  scale5_mut = 1.0d0

                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3_dir = 1.0d0 - expdamp
                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7_dir = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut =  (1.0d0-expdamp)
                        scale5_mut =  (1.0d0-expdamp
     &                                        *(1.0d0-damp))

                     end if
                  end if
                  rr3_dir = scale3_dir / (r*r2)
                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

                  rr3 = scale3_mut / (r*r2)
                  rr5 = 3.0d0 * scale5_mut / (r*r2*r2)

                  if (use_polymer .and. r2 .le. polycut2) then
                   rr3 = rr3*dscale_mut(l3)
                   rr5 = rr5*dscale_mut(l3)
                  end if
                  Txx = -(-rr3 + xr*xr*rr5)
                  Txy = -(xr*yr*rr5)
                  Txz = -(xr*zr*rr5)
                  Tyx = Txy
                  Tyy = -(-rr3 + yr*yr*rr5)
                  Tyz = -(yr*zr*rr5)
                  Tzx = Txz
                  Tzy = Tyz
                  Tzz = -(-rr3 + zr*zr*rr5)

                  if (ii .ne. kk) then
                  M_tot(i2+1,k2+1) = Txx + M_tot(i2+1,k2+1)
                  M_tot(i2+1,k2+2) = Txy + M_tot(i2+1,k2+2)
                  M_tot(i2+1,k2+3) = Txz + M_tot(i2+1,k2+3)
                  M_tot(i2+2,k2+1) = Tyx + M_tot(i2+2,k2+1)
                  M_tot(i2+2,k2+2) = Tyy + M_tot(i2+2,k2+2)
                  M_tot(i2+2,k2+3) = Tyz + M_tot(i2+2,k2+3)
                  M_tot(i2+3,k2+1) = Tzx + M_tot(i2+3,k2+1)
                  M_tot(i2+3,k2+2) = Tzy + M_tot(i2+3,k2+2)
                  M_tot(i2+3,k2+3) = Tzz + M_tot(i2+3,k2+3)

                  M_tot(k2+1,i2+1) = Txx + M_tot(k2+1,i2+1)
                  M_tot(k2+1,i2+2) = Txy + M_tot(k2+1,i2+2)
                  M_tot(k2+1,i2+3) = Txz + M_tot(k2+1,i2+3)
                  M_tot(k2+2,i2+1) = Tyx + M_tot(k2+2,i2+1)
                  M_tot(k2+2,i2+2) = Tyy + M_tot(k2+2,i2+2)
                  M_tot(k2+2,i2+3) = Tyz + M_tot(k2+2,i2+3)
                  M_tot(k2+3,i2+1) = Tzx + M_tot(k2+3,i2+1)
                  M_tot(k2+3,i2+2) = Tzy + M_tot(k2+3,i2+2)
                  M_tot(k2+3,i2+3) = Tzz + M_tot(k2+3,i2+3)
                  end if
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
                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
                  do j = 1, 3
                     fip(j) = fid(j)
                     fkp(j) = fkd(j)
                  end do
                  if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale_dir(l3)
                           fip(j) = fip(j) * pscale(l3)
                           fkd(j) = fkd(j) * dscale_dir(l3)
                           fkp(j) = fkp(j) * pscale(l3)
                        end do
                  end if
                  do j = 1, 3
                     field(j,l1) = field(j,l1) + fid(j)
                     fieldp(j,l1) = fieldp(j,l1) + fip(j)
                     if (ii .ne. kk) then
                           field(j,l3) = field(j,l3) + fkd(j)
                           fieldp(j,l3) = fieldp(j,l3) + fkp(j)
                     end if
                  end do
               end if
            end do
         end do
      end do

      end if

      deallocate (dscale_dir)
      deallocate (dscale_mut)
      deallocate (pscale)
      return
      end



