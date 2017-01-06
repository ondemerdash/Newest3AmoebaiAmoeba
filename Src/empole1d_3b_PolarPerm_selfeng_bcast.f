c
c
      subroutine empole1d_3b_PolarPerm_selfeng_bcast
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
      use polar
      use polpot
      use virial
      use mplpot
      use deriv3b
      implicit none
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
      real*8 em_recip

c     set the energy unit conversion factor
c
      f = electric / dielec

      em=em+emreal_tmp
      ep3b=ep3b+ep3b_recip

      do i=1,npole
         do j=1,3
            dem(j,i)=dem(j,i)+demreal_tmp(j,i)
            dep3b(j,i)=dep3b(j,i)+dep3b_recip(j,i)
         end do 
      end do 
      do i=1,3
        do j=1,3
           virep3b(j,i)=virep3b(j,i)+virep3b_recip(j,i)
        end do
      end do

c      print*,"After ereal1c_3b_Perm em real",em-em_recip
c
c     compute the Ewald self-energy term over all the atoms
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e

         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         ei = fterm * term * uii / 3.0d0
         ep3b = ep3b+ei
      end do

      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,i)+uinp(1,i))
         uiy = 0.5d0 * (uind(2,i)+uinp(2,i))
         uiz = 0.5d0 * (uind(3,i)+uinp(3,i))
         trqi(1) = term * (diy*uiz-diz*uiy)
         trqi(2) = term * (diz*uix-dix*uiz)
         trqi(3) = term * (dix*uiy-diy*uix)
         call torque_3b_new_npole(dep3b,i,
     &      trq,trqi,frcx,frcy,frcz)
      end do

c
c     compute the cell dipole boundary correction term
c

      if (boundary .eq. 'VACUUM') then
        trqi(1) = 0.0d0
        trqi(2) = 0.0d0
        trqi(3) = 0.0d0
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xu = 0.0d0
        yu = 0.0d0
        zu = 0.0d0
        xup = 0.0d0
        yup = 0.0d0
        zup = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
          yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
          zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
        end do
        term = (2.0d0/3.0d0) * f * (pi/volbox)
        em = em + term*(xd*xd+yd*yd+zd*zd)
        do i = 1, npole
          ii = ipole(i)
          dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
          dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
          dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
        end do
        xdfield = -2.0d0 * term * xd
        ydfield = -2.0d0 * term * yd
        zdfield = -2.0d0 * term * zd
        do i = 1, npole
          trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
          trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
          trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
c          call torque_3b_Perm (i,trq,trqi,frcx,frcy,frcz)
c          call torque_3b_Perm_new (demk,
c     &          k,ttm3,ttm3i,frcxk,frcyk,frczk)
          call torque_3b_Perm_new (dem,
     &          i,trq,trqi,frcx,frcy,frcz)

        end do
c
c     boundary correction to virial due to overall cell dipole
c
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xq = 0.0d0
        yq = 0.0d0
        zq = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i)
          yd = yd + rpole(3,i)
          zd = zd + rpole(4,i)
          xq = xq + rpole(1,i)*x(ii)
          yq = yq + rpole(1,i)*y(ii)
          zq = zq + rpole(1,i)*z(ii)
        end do
        xv = xq * (xd+0.5d0*(xu+xup))
        yv = yq * (yd+0.5d0*(yu+yup))
        zv = zq * (zd+0.5d0*(zu+zup))
        vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
        vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
        vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
        vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
        vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
        vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
        vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
        vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
        vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
        vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
        if (poltyp .eq. 'DIRECT') then
          vterm = term * (xu*xup+yu*yup+zu*zup)
          vir(1,1) = vir(1,1) + vterm
          vir(2,2) = vir(2,2) + vterm
          vir(3,3) = vir(3,3) + vterm
        end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      return
      end
c
