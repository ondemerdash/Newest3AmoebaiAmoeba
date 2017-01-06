c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque3  --  convert torque to force components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque3" takes the torque values on a single site defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame; this
c     version also returns the individual atomic components
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
      subroutine torque3_dEtensor_72015_taunew3_pt1_vir2(
     &  i,trq2,taud,taup,
     &  k,kkk,iter3)
      use sizes
      use atoms
      use deriv
      use mpole
      use dEtensor
      use dEtensor2
      use dEtensor3
      implicit none
      integer i,j
      integer ia,ib,ic,id
      real*8 du,dv,dw,random
      real*8 usiz,vsiz,wsiz
      real*8 rsiz,ssiz
      real*8 t1siz,t2siz
      real*8 uvsiz,uwsiz,vwsiz
      real*8 ursiz,ussiz
      real*8 vssiz,wssiz
      real*8 uvcos,uwcos,vwcos
      real*8 urcos,uscos
      real*8 vscos,wscos
      real*8 ut1cos,ut2cos
      real*8 uvsin,uwsin,vwsin
      real*8 ursin,ussin
      real*8 vssin,wssin
      real*8 ut1sin,ut2sin
      real*8 dphidu,dphidv,dphidw
      real*8 dphidr,dphids
      real*8 trq1(3),trq2(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
c      real*8 demi(3,*)
c      real*8 depi(3,*)
      character*8 axetyp
      real*8 taud(9),taup(9)
c      integer nlocaltau1,ilocaltau1(2,*),ikflag,k,m
c      real*8 taulocal1d(9,*),taulocal1p(9,*)
      integer nlocaltau1,ikflag,k,m
      real*8 dphidu_uindx,dphidu_uindy,dphidu_uindz,dphidv_uindx
      real*8 dphidv_uindy
      real*8 dphidv_uindz,dphidw_uindx,dphidw_uindy,dphidw_uindz
      real*8 dphidu_uinpx,dphidu_uinpy,dphidu_uinpz,dphidv_uinpx
      real*8 dphidv_uinpy
      real*8 dphidv_uinpz,dphidw_uinpx,dphidw_uinpy,dphidw_uinpz
      real*8 duxx_d,duxy_d,duxz_d,duyx_d,duyy_d,duyz_d,duzx_d,duzy_d
      real*8 duzz_d
      real*8 duxx_p,duxy_p,duxz_p,duyx_p,duyy_p,duyz_p,duzx_p,duzy_p
      real*8 duzz_p
      real*8 dvxx_d,dvxy_d,dvxz_d,dvyx_d,dvyy_d,dvyz_d,dvzx_d,dvzy_d
      real*8 dvzz_d
      real*8 dvxx_p,dvxy_p,dvxz_p,dvyx_p,dvyy_p,dvyz_p,dvzx_p,dvzy_p
      real*8 dvzz_p
      logical done_ia,done_ib,done_ic
      integer nlocal,iter3,kkk
c      real*8 frcztaulocal1d(9,*),frcztaulocal1p(9,*)
c      real*8 frcxtaulocal1d(9,*),frcxtaulocal1p(9,*)
c      real*8 frcytaulocal1d(9,*),frcytaulocal1p(9,*)
c
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0d0
         frcx(j) = 0.0d0
         frcy(j) = 0.0d0
      end do
c
c     get the local frame type and the frame-defining atoms
c
      ia = zaxis(i)
      ib = ipole(i)
      ic = xaxis(i)
      id = yaxis(i)
      axetyp = polaxe(i)
      if (axetyp .eq. 'None')  return
c
c     construct the three rotation axes for the local frame
c
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
      else
         v(1) = random ()
         v(2) = random ()
         v(3) = random ()
      end if
      if (axetyp.eq.'Z-Bisect' .or. axetyp.eq.'3-Fold') then
         w(1) = x(id) - x(ib)
         w(2) = y(id) - y(ib)
         w(3) = z(id) - z(ib)
      else
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
      end if
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      if (axetyp .eq. 'Z-Bisect') then
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         ssiz = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
            s(j) = s(j) / ssiz
         end do
      end if
c
c     find the perpendicular and angle for each pair of axes
c
      uv(1) = v(2)*u(3) - v(3)*u(2)
      uv(2) = v(3)*u(1) - v(1)*u(3)
      uv(3) = v(1)*u(2) - v(2)*u(1)
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
      vwsiz = sqrt(vw(1)*vw(1) + vw(2)*vw(2) + vw(3)*vw(3))
      do j = 1, 3
         uv(j) = uv(j) / uvsiz
         uw(j) = uw(j) / uwsiz
         vw(j) = vw(j) / vwsiz
      end do
      if (axetyp .eq. 'Z-Bisect') then
         ur(1) = r(2)*u(3) - r(3)*u(2)
         ur(2) = r(3)*u(1) - r(1)*u(3)
         ur(3) = r(1)*u(2) - r(2)*u(1)
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
         wssiz = sqrt(ws(1)*ws(1) + ws(2)*ws(2) + ws(3)*ws(3))
         do j = 1, 3
            ur(j) = ur(j) / ursiz
            us(j) = us(j) / ussiz
            vs(j) = vs(j) / vssiz
            ws(j) = ws(j) / wssiz
         end do
      end if
c
c     get sine and cosine of angles between the rotation axes
c
      uvcos = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      uvsin = sqrt(1.0d0 - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0d0 - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0d0 - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0d0 - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0d0 - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0d0 - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0d0 - wscos*wscos)
      end if
c
c     compute the projection of v and w onto the ru-plane
c
      if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
           t1(j) = v(j) - s(j)*vscos
           t2(j) = w(j) - s(j)*wscos
         end do
         t1siz = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         t2siz = sqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
         do j = 1, 3
           t1(j) = t1(j) / t1siz
           t2(j) = t2(j) / t2siz
         end do
         ut1cos = u(1)*t1(1) + u(2)*t1(2) + u(3)*t1(3)
         ut1sin = sqrt(1.0d0 - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0d0 - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
c
c     force distribution for the Z-Only local coordinate method
c
c      if (axetyp .eq. 'Z-Only') then
c
c     force distribution for the Z-then-X local coordinate method
c
c      else if (axetyp .eq. 'Z-then-X') then
c
c     force distribution for the Bisector local coordinate method
c
c      else if (axetyp .eq. 'Bisector') then
c
c     force distribution for the Z-Bisect local coordinate method
c
c      else if (axetyp .eq. 'Z-Bisect') then
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
c      else if (axetyp .eq. '3-Fold') then
c      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      !print*,"!!!IN NEW2TAU TENSOR TORQ"
      dphidu = -trq2(1)*u(1) - trq2(2)*u(2) - trq2(3)*u(3)
      dphidv = -trq2(1)*v(1) - trq2(2)*v(2) - trq2(3)*v(3)
      dphidw = -trq2(1)*w(1) - trq2(2)*w(2) - trq2(3)*w(3)

      dphidu_uindx = -taud(1)*u(1) -taud(4)*u(2) -taud(7)*u(3)
      dphidu_uindy = -taud(2)*u(1) -taud(5)*u(2) -taud(8)*u(3)
      dphidu_uindz = -taud(3)*u(1) -taud(6)*u(2) -taud(9)*u(3)

      dphidv_uindx = -taud(1)*v(1) -taud(4)*v(2) -taud(7)*v(3)
      dphidv_uindy = -taud(2)*v(1) -taud(5)*v(2) -taud(8)*v(3)
      dphidv_uindz = -taud(3)*v(1) -taud(6)*v(2) -taud(9)*v(3)

      dphidw_uindx = -taud(1)*w(1) -taud(4)*w(2) -taud(7)*w(3)
      dphidw_uindy = -taud(2)*w(1) -taud(5)*w(2) -taud(8)*w(3)
      dphidw_uindz = -taud(3)*w(1) -taud(6)*w(2) -taud(9)*w(3)

      dphidu_uinpx = -taup(1)*u(1) -taup(4)*u(2) -taup(7)*u(3)
      dphidu_uinpy = -taup(2)*u(1) -taup(5)*u(2) -taup(8)*u(3)
      dphidu_uinpz = -taup(3)*u(1) -taup(6)*u(2) -taup(9)*u(3)

      dphidv_uinpx = -taup(1)*v(1) -taup(4)*v(2) -taup(7)*v(3)
      dphidv_uinpy = -taup(2)*v(1) -taup(5)*v(2) -taup(8)*v(3)
      dphidv_uinpz = -taup(3)*v(1) -taup(6)*v(2) -taup(9)*v(3)

      dphidw_uinpx = -taup(1)*w(1) -taup(4)*w(2) -taup(7)*w(3)
      dphidw_uinpy = -taup(2)*w(1) -taup(5)*w(2) -taup(8)*w(3)
      dphidw_uinpz = -taup(3)*w(1) -taup(6)*w(2) -taup(9)*w(3)



      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq2(1)*r(1) - trq2(2)*r(2) - trq2(3)*r(3)
         dphids = -trq2(1)*s(1) - trq2(2)*s(2) - trq2(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ib) = depi(j,ib) - du
c            frcz(j) = frcz(j) + du

         end do
c
c     force distribution for the Z-then-X local coordinate method
c     WATER AXETYP
      else if (axetyp .eq. 'Z-then-X') then
c         do j = 1, 3
c            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
c            dv = -uv(j)*dphidu/(vsiz*uvsin)
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv
c            frcz(j) = frcz(j) + du
c            frcx(j) = frcx(j) + dv           
c         end do

         duxx_d=uv(1)*dphidv_uindx/(usiz*uvsin)+ uw(1)*dphidw_uindx/usiz
         duxy_d=uv(1)*dphidv_uindy/(usiz*uvsin)+ uw(1)*dphidw_uindy/usiz
         duxz_d=uv(1)*dphidv_uindz/(usiz*uvsin)+ uw(1)*dphidw_uindz/usiz

         duyx_d=uv(2)*dphidv_uindx/(usiz*uvsin)+ uw(2)*dphidw_uindx/usiz
         duyy_d=uv(2)*dphidv_uindy/(usiz*uvsin)+ uw(2)*dphidw_uindy/usiz
         duyz_d=uv(2)*dphidv_uindz/(usiz*uvsin)+ uw(2)*dphidw_uindz/usiz

         duzx_d=uv(3)*dphidv_uindx/(usiz*uvsin)+ uw(3)*dphidw_uindx/usiz
         duzy_d=uv(3)*dphidv_uindy/(usiz*uvsin)+ uw(3)*dphidw_uindy/usiz
         duzz_d=uv(3)*dphidv_uindz/(usiz*uvsin)+ uw(3)*dphidw_uindz/usiz

         duxx_p=uv(1)*dphidv_uinpx/(usiz*uvsin)+ uw(1)*dphidw_uinpx/usiz
         duxy_p=uv(1)*dphidv_uinpy/(usiz*uvsin)+ uw(1)*dphidw_uinpy/usiz
         duxz_p=uv(1)*dphidv_uinpz/(usiz*uvsin)+ uw(1)*dphidw_uinpz/usiz

         duyx_p=uv(2)*dphidv_uinpx/(usiz*uvsin)+ uw(2)*dphidw_uinpx/usiz
         duyy_p=uv(2)*dphidv_uinpy/(usiz*uvsin)+ uw(2)*dphidw_uinpy/usiz
         duyz_p=uv(2)*dphidv_uinpz/(usiz*uvsin)+ uw(2)*dphidw_uinpz/usiz

         duzx_p=uv(3)*dphidv_uinpx/(usiz*uvsin)+ uw(3)*dphidw_uinpx/usiz
         duzy_p=uv(3)*dphidv_uinpy/(usiz*uvsin)+ uw(3)*dphidw_uinpy/usiz
         duzz_p=uv(3)*dphidv_uinpz/(usiz*uvsin)+ uw(3)*dphidw_uinpz/usiz


         dvxx_d = -uv(1)*dphidu_uindx/(vsiz*uvsin)
         dvxy_d = -uv(1)*dphidu_uindy/(vsiz*uvsin)
         dvxz_d = -uv(1)*dphidu_uindz/(vsiz*uvsin)

         dvyx_d = -uv(2)*dphidu_uindx/(vsiz*uvsin) 
         dvyy_d = -uv(2)*dphidu_uindy/(vsiz*uvsin)
         dvyz_d = -uv(2)*dphidu_uindz/(vsiz*uvsin)

         dvzx_d = -uv(3)*dphidu_uindx/(vsiz*uvsin)
         dvzy_d = -uv(3)*dphidu_uindy/(vsiz*uvsin)
         dvzz_d = -uv(3)*dphidu_uindz/(vsiz*uvsin)

         dvxx_p = -uv(1)*dphidu_uinpx/(vsiz*uvsin) 
         dvxy_p = -uv(1)*dphidu_uinpy/(vsiz*uvsin)
         dvxz_p = -uv(1)*dphidu_uinpz/(vsiz*uvsin)

         dvyx_p = -uv(2)*dphidu_uinpx/(vsiz*uvsin)
         dvyy_p = -uv(2)*dphidu_uinpy/(vsiz*uvsin)
         dvyz_p = -uv(2)*dphidu_uinpz/(vsiz*uvsin)

         dvzx_p = -uv(3)*dphidu_uinpx/(vsiz*uvsin)
         dvzy_p = -uv(3)*dphidu_uinpy/(vsiz*uvsin)
         dvzz_p = -uv(3)*dphidu_uinpz/(vsiz*uvsin)

c           done_ia=.false.
c           done_ic=.false.
c           done_ib=.false.
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv
      

c                nlocaltau1=nlocaltau1+1
                elsttau1_index_tmp(1,kkk,iter3)=ia
                !elsttau1_tmp(k,ia)=nlocaltau1
c                elsttau1(1,kkk,i)=nlocaltau1

c                ilocaltau1(1,nlocaltau1)=ia
c                ilocaltau1(2,nlocaltau1)=k
c                taulocal1d(1,nlocaltau1)=duxx_d
c                taulocal1d(2,nlocaltau1)=duxy_d
c                taulocal1d(3,nlocaltau1)=duxz_d
c                taulocal1d(4,nlocaltau1)=duyx_d
c                taulocal1d(5,nlocaltau1)=duyy_d
c                taulocal1d(6,nlocaltau1)=duyz_d
c                taulocal1d(7,nlocaltau1)=duzx_d
c                taulocal1d(8,nlocaltau1)=duzy_d
c                taulocal1d(9,nlocaltau1)=duzz_d
c                taulocal1p(1,nlocaltau1)=duxx_p
c                taulocal1p(2,nlocaltau1)=duxy_p
c                taulocal1p(3,nlocaltau1)=duxz_p
c                taulocal1p(4,nlocaltau1)=duyx_p
c                taulocal1p(5,nlocaltau1)=duyy_p
c                taulocal1p(6,nlocaltau1)=duyz_p
c                taulocal1p(7,nlocaltau1)=duzx_p
c                taulocal1p(8,nlocaltau1)=duzy_p
c                taulocal1p(9,nlocaltau1)=duzz_p

                taud1_3tmp(1,1,kkk,iter3)=duxx_d
                taud1_3tmp(2,1,kkk,iter3)=duxy_d
                taud1_3tmp(3,1,kkk,iter3)=duxz_d
                taud1_3tmp(4,1,kkk,iter3)=duyx_d
                taud1_3tmp(5,1,kkk,iter3)=duyy_d
                taud1_3tmp(6,1,kkk,iter3)=duyz_d
                taud1_3tmp(7,1,kkk,iter3)=duzx_d
                taud1_3tmp(8,1,kkk,iter3)=duzy_d
                taud1_3tmp(9,1,kkk,iter3)=duzz_d
                taup1_3tmp(1,1,kkk,iter3)=duxx_p
                taup1_3tmp(2,1,kkk,iter3)=duxy_p
                taup1_3tmp(3,1,kkk,iter3)=duxz_p
                taup1_3tmp(4,1,kkk,iter3)=duyx_p
                taup1_3tmp(5,1,kkk,iter3)=duyy_p
                taup1_3tmp(6,1,kkk,iter3)=duyz_p
                taup1_3tmp(7,1,kkk,iter3)=duzx_p
                taup1_3tmp(8,1,kkk,iter3)=duzy_p
                taup1_3tmp(9,1,kkk,iter3)=duzz_p

c            frcz(j) = frcz(j) + du
c            frcx(j) = frcx(j) + dv           
c         end do

c                frcztaulocal1d(1,nlocal)=duxx_d
c                frcztaulocal1d(2,nlocal)=duxy_d
c                frcztaulocal1d(3,nlocal)=duxz_d
c                frcztaulocal1d(4,nlocal)=duyx_d
c                frcztaulocal1d(5,nlocal)=duyy_d
c                frcztaulocal1d(6,nlocal)=duyz_d
c                frcztaulocal1d(7,nlocal)=duzx_d
c                frcztaulocal1d(8,nlocal)=duzy_d
c                frcztaulocal1d(9,nlocal)=duzz_d
c                frcztaulocal1p(1,nlocal)=duxx_p
c                frcztaulocal1p(2,nlocal)=duxy_p
c                frcztaulocal1p(3,nlocal)=duxz_p
c                frcztaulocal1p(4,nlocal)=duyx_p
c                frcztaulocal1p(5,nlocal)=duyy_p
c                frcztaulocal1p(6,nlocal)=duyz_p
c                frcztaulocal1p(7,nlocal)=duzx_p
c                frcztaulocal1p(8,nlocal)=duzy_p
c                frcztaulocal1p(9,nlocal)=duzz_p

                frcztau1dtot_3tmp(1,kkk,iter3)=duxx_d
                frcztau1dtot_3tmp(2,kkk,iter3)=duxy_d
                frcztau1dtot_3tmp(3,kkk,iter3)=duxz_d
                frcztau1dtot_3tmp(4,kkk,iter3)=duyx_d
                frcztau1dtot_3tmp(5,kkk,iter3)=duyy_d
                frcztau1dtot_3tmp(6,kkk,iter3)=duyz_d
                frcztau1dtot_3tmp(7,kkk,iter3)=duzx_d
                frcztau1dtot_3tmp(8,kkk,iter3)=duzy_d
                frcztau1dtot_3tmp(9,kkk,iter3)=duzz_d
                frcztau1ptot_3tmp(1,kkk,iter3)=duxx_p
                frcztau1ptot_3tmp(2,kkk,iter3)=duxy_p
                frcztau1ptot_3tmp(3,kkk,iter3)=duxz_p
                frcztau1ptot_3tmp(4,kkk,iter3)=duyx_p
                frcztau1ptot_3tmp(5,kkk,iter3)=duyy_p
                frcztau1ptot_3tmp(6,kkk,iter3)=duyz_p
                frcztau1ptot_3tmp(7,kkk,iter3)=duzx_p
                frcztau1ptot_3tmp(8,kkk,iter3)=duzy_p
                frcztau1ptot_3tmp(9,kkk,iter3)=duzz_p


c                frcxtaulocal1d(1,nlocal)=dvxx_d
c                frcxtaulocal1d(2,nlocal)=dvxy_d
c                frcxtaulocal1d(3,nlocal)=dvxz_d
c                frcxtaulocal1d(4,nlocal)=dvyx_d
c                frcxtaulocal1d(5,nlocal)=dvyy_d
c                frcxtaulocal1d(6,nlocal)=dvyz_d
c                frcxtaulocal1d(7,nlocal)=dvzx_d
c                frcxtaulocal1d(8,nlocal)=dvzy_d
c                frcxtaulocal1d(9,nlocal)=dvzz_d
c                frcxtaulocal1p(1,nlocal)=dvxx_p
c                frcxtaulocal1p(2,nlocal)=dvxy_p
c                frcxtaulocal1p(3,nlocal)=dvxz_p
c                frcxtaulocal1p(4,nlocal)=dvyx_p
c                frcxtaulocal1p(5,nlocal)=dvyy_p
c                frcxtaulocal1p(6,nlocal)=dvyz_p
c                frcxtaulocal1p(7,nlocal)=dvzx_p
c                frcxtaulocal1p(8,nlocal)=dvzy_p
c                frcxtaulocal1p(9,nlocal)=dvzz_p

                frcxtau1dtot_3tmp(1,kkk,iter3)=dvxx_d
                frcxtau1dtot_3tmp(2,kkk,iter3)=dvxy_d
                frcxtau1dtot_3tmp(3,kkk,iter3)=dvxz_d
                frcxtau1dtot_3tmp(4,kkk,iter3)=dvyx_d
                frcxtau1dtot_3tmp(5,kkk,iter3)=dvyy_d
                frcxtau1dtot_3tmp(6,kkk,iter3)=dvyz_d
                frcxtau1dtot_3tmp(7,kkk,iter3)=dvzx_d
                frcxtau1dtot_3tmp(8,kkk,iter3)=dvzy_d
                frcxtau1dtot_3tmp(9,kkk,iter3)=dvzz_d
                frcxtau1ptot_3tmp(1,kkk,iter3)=dvxx_p
                frcxtau1ptot_3tmp(2,kkk,iter3)=dvxy_p
                frcxtau1ptot_3tmp(3,kkk,iter3)=dvxz_p
                frcxtau1ptot_3tmp(4,kkk,iter3)=dvyx_p
                frcxtau1ptot_3tmp(5,kkk,iter3)=dvyy_p
                frcxtau1ptot_3tmp(6,kkk,iter3)=dvyz_p
                frcxtau1ptot_3tmp(7,kkk,iter3)=dvzx_p
                frcxtau1ptot_3tmp(8,kkk,iter3)=dvzy_p
                frcxtau1ptot_3tmp(9,kkk,iter3)=dvzz_p


c                frcytaulocal1d(1,nlocal)=0.0d0
c                frcytaulocal1d(2,nlocal)=0.0d0
c                frcytaulocal1d(3,nlocal)=0.0d0
c                frcytaulocal1d(4,nlocal)=0.0d0
c                frcytaulocal1d(5,nlocal)=0.0d0
c                frcytaulocal1d(6,nlocal)=0.0d0
c                frcytaulocal1d(7,nlocal)=0.0d0
c                frcytaulocal1d(8,nlocal)=0.0d0
c                frcytaulocal1d(9,nlocal)=0.0d0
c                frcytaulocal1p(1,nlocal)=0.0d0
c                frcytaulocal1p(2,nlocal)=0.0d0
c                frcytaulocal1p(3,nlocal)=0.0d0
c                frcytaulocal1p(4,nlocal)=0.0d0
c                frcytaulocal1p(5,nlocal)=0.0d0
c                frcytaulocal1p(6,nlocal)=0.0d0
c                frcytaulocal1p(7,nlocal)=0.0d0
c                frcytaulocal1p(8,nlocal)=0.0d0
c                frcytaulocal1p(9,nlocal)=0.0d0

                frcytau1dtot_3tmp(1,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(2,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(3,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(4,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(5,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(6,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(7,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(8,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(9,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(1,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(2,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(3,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(4,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(5,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(6,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(7,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(8,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(9,kkk,iter3)=0.0d0

c                nlocaltau1=nlocaltau1+1
c                ilocaltau1(1,nlocaltau1)=ic
c                ilocaltau1(2,nlocaltau1)=k
                elsttau1_index_tmp(2,kkk,iter3)=ic
                !elsttau1_tmp(k,ia)=nlocaltau1
c                elsttau1(2,kkk,i)=nlocaltau1

c                taulocal1d(1,nlocaltau1)=dvxx_d
c                taulocal1d(2,nlocaltau1)=dvxy_d
c                taulocal1d(3,nlocaltau1)=dvxz_d
c                taulocal1d(4,nlocaltau1)=dvyx_d
c                taulocal1d(5,nlocaltau1)=dvyy_d
c                taulocal1d(6,nlocaltau1)=dvyz_d
c                taulocal1d(7,nlocaltau1)=dvzx_d
c                taulocal1d(8,nlocaltau1)=dvzy_d
c                taulocal1d(9,nlocaltau1)=dvzz_d
c                taulocal1p(1,nlocaltau1)=dvxx_p
c                taulocal1p(2,nlocaltau1)=dvxy_p
c                taulocal1p(3,nlocaltau1)=dvxz_p
c                taulocal1p(4,nlocaltau1)=dvyx_p
c                taulocal1p(5,nlocaltau1)=dvyy_p
c                taulocal1p(6,nlocaltau1)=dvyz_p
c                taulocal1p(7,nlocaltau1)=dvzx_p
c                taulocal1p(8,nlocaltau1)=dvzy_p
c                taulocal1p(9,nlocaltau1)=dvzz_p

                taud1_3tmp(1,2,kkk,iter3)=dvxx_d
                taud1_3tmp(2,2,kkk,iter3)=dvxy_d
                taud1_3tmp(3,2,kkk,iter3)=dvxz_d
                taud1_3tmp(4,2,kkk,iter3)=dvyx_d
                taud1_3tmp(5,2,kkk,iter3)=dvyy_d
                taud1_3tmp(6,2,kkk,iter3)=dvyz_d
                taud1_3tmp(7,2,kkk,iter3)=dvzx_d
                taud1_3tmp(8,2,kkk,iter3)=dvzy_d
                taud1_3tmp(9,2,kkk,iter3)=dvzz_d
                taup1_3tmp(1,2,kkk,iter3)=dvxx_p
                taup1_3tmp(2,2,kkk,iter3)=dvxy_p
                taup1_3tmp(3,2,kkk,iter3)=dvxz_p
                taup1_3tmp(4,2,kkk,iter3)=dvyx_p
                taup1_3tmp(5,2,kkk,iter3)=dvyy_p
                taup1_3tmp(6,2,kkk,iter3)=dvyz_p
                taup1_3tmp(7,2,kkk,iter3)=dvzx_p
                taup1_3tmp(8,2,kkk,iter3)=dvzy_p
                taup1_3tmp(9,2,kkk,iter3)=dvzz_p


c                nlocaltau1=nlocaltau1+1
c                ilocaltau1(1,nlocaltau1)=ib
c                ilocaltau1(2,nlocaltau1)=k
                elsttau1_index_tmp(3,kkk,iter3)=ib
                !elsttau1_tmp(k,ia)=nlocaltau1
c                elsttau1(3,kkk,i)=nlocaltau1

c                taulocal1d(1,nlocaltau1)= -duxx_d-dvxx_d
c                taulocal1d(2,nlocaltau1)= -duxy_d-dvxy_d
c                taulocal1d(3,nlocaltau1)= -duxz_d-dvxz_d
c                taulocal1d(4,nlocaltau1)= -duyx_d-dvyx_d
c                taulocal1d(5,nlocaltau1)= -duyy_d-dvyy_d
c                taulocal1d(6,nlocaltau1)= -duyz_d-dvyz_d
c                taulocal1d(7,nlocaltau1)= -duzx_d-dvzx_d
c                taulocal1d(8,nlocaltau1)= -duzy_d-dvzy_d
c                taulocal1d(9,nlocaltau1)= -duzz_d-dvzz_d
c                taulocal1p(1,nlocaltau1)= -duxx_p-dvxx_p
c                taulocal1p(2,nlocaltau1)= -duxy_p-dvxy_p
c                taulocal1p(3,nlocaltau1)= -duxz_p-dvxz_p
c                taulocal1p(4,nlocaltau1)= -duyx_p-dvyx_p
c                taulocal1p(5,nlocaltau1)= -duyy_p-dvyy_p
c                taulocal1p(6,nlocaltau1)= -duyz_p-dvyz_p
c                taulocal1p(7,nlocaltau1)= -duzx_p-dvzx_p
c                taulocal1p(8,nlocaltau1)= -duzy_p-dvzy_p
c                taulocal1p(9,nlocaltau1)= -duzz_p-dvzz_p

                taud1_3tmp(1,3,kkk,iter3)= -duxx_d-dvxx_d
                taud1_3tmp(2,3,kkk,iter3)= -duxy_d-dvxy_d
                taud1_3tmp(3,3,kkk,iter3)= -duxz_d-dvxz_d
                taud1_3tmp(4,3,kkk,iter3)= -duyx_d-dvyx_d
                taud1_3tmp(5,3,kkk,iter3)= -duyy_d-dvyy_d
                taud1_3tmp(6,3,kkk,iter3)= -duyz_d-dvyz_d
                taud1_3tmp(7,3,kkk,iter3)= -duzx_d-dvzx_d
                taud1_3tmp(8,3,kkk,iter3)= -duzy_d-dvzy_d
                taud1_3tmp(9,3,kkk,iter3)= -duzz_d-dvzz_d
                taup1_3tmp(1,3,kkk,iter3)= -duxx_p-dvxx_p
                taup1_3tmp(2,3,kkk,iter3)= -duxy_p-dvxy_p
                taup1_3tmp(3,3,kkk,iter3)= -duxz_p-dvxz_p
                taup1_3tmp(4,3,kkk,iter3)= -duyx_p-dvyx_p
                taup1_3tmp(5,3,kkk,iter3)= -duyy_p-dvyy_p
                taup1_3tmp(6,3,kkk,iter3)= -duyz_p-dvyz_p
                taup1_3tmp(7,3,kkk,iter3)= -duzx_p-dvzx_p
                taup1_3tmp(8,3,kkk,iter3)= -duzy_p-dvzy_p
                taup1_3tmp(9,3,kkk,iter3)= -duzz_p-dvzz_p

                elsttau1_index_tmp(4,kkk,iter3)=0

c
c     force distribution for the Bisector local coordinate method
c     WATER AXETYP 
      else if (axetyp .eq. 'Bisector') then
c         do j = 1, 3
c            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
c            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv
c            frcz(j) = frcz(j) + du
c            frcx(j) = frcx(j) + dv
c         end do

        duxx_d=uv(1)*dphidv_uindx/(usiz*uvsin) 
     &  + 0.5d0*uw(1)*dphidw_uindx/usiz 
        duxy_d=uv(1)*dphidv_uindy/(usiz*uvsin) 
     &  + 0.5d0*uw(1)*dphidw_uindy/usiz
        duxz_d=uv(1)*dphidv_uindz/(usiz*uvsin) 
     &  + 0.5d0*uw(1)*dphidw_uindz/usiz

        duyx_d=uv(2)*dphidv_uindx/(usiz*uvsin) 
     &  + 0.5d0*uw(2)*dphidw_uindx/usiz   
        duyy_d=uv(2)*dphidv_uindy/(usiz*uvsin) 
     &  + 0.5d0*uw(2)*dphidw_uindy/usiz
        duyz_d=uv(2)*dphidv_uindz/(usiz*uvsin) 
     &  + 0.5d0*uw(2)*dphidw_uindz/usiz

        duzx_d=uv(3)*dphidv_uindx/(usiz*uvsin) 
     &   + 0.5d0*uw(3)*dphidw_uindx/usiz
        duzy_d=uv(3)*dphidv_uindy/(usiz*uvsin) 
     &   + 0.5d0*uw(3)*dphidw_uindy/usiz
        duzz_d=uv(3)*dphidv_uindz/(usiz*uvsin) 
     &   + 0.5d0*uw(3)*dphidw_uindz/usiz

        duxx_p=uv(1)*dphidv_uinpx/(usiz*uvsin) 
     &    + 0.5d0*uw(1)*dphidw_uinpx/usiz
        duxy_p=uv(1)*dphidv_uinpy/(usiz*uvsin) 
     &    + 0.5d0*uw(1)*dphidw_uinpy/usiz
        duxz_p=uv(1)*dphidv_uinpz/(usiz*uvsin) 
     &    + 0.5d0*uw(1)*dphidw_uinpz/usiz

        duyx_p=uv(2)*dphidv_uinpx/(usiz*uvsin) 
     &   + 0.5d0*uw(2)*dphidw_uinpx/usiz
        duyy_p=uv(2)*dphidv_uinpy/(usiz*uvsin) 
     &   + 0.5d0*uw(2)*dphidw_uinpy/usiz
        duyz_p=uv(2)*dphidv_uinpz/(usiz*uvsin) 
     &   + 0.5d0*uw(2)*dphidw_uinpz/usiz

        duzx_p=uv(3)*dphidv_uinpx/(usiz*uvsin) 
     &  + 0.5d0*uw(3)*dphidw_uinpx/usiz
        duzy_p=uv(3)*dphidv_uinpy/(usiz*uvsin) 
     &  + 0.5d0*uw(3)*dphidw_uinpy/usiz
        duzz_p=uv(3)*dphidv_uinpz/(usiz*uvsin) 
     &  + 0.5d0*uw(3)*dphidw_uinpz/usiz

        dvxx_d=-uv(1)*dphidu_uindx/(vsiz*uvsin)
     &  + 0.5d0*vw(1)*dphidw_uindx/vsiz
        dvxy_d=-uv(1)*dphidu_uindy/(vsiz*uvsin)
     &  + 0.5d0*vw(1)*dphidw_uindy/vsiz
        dvxz_d=-uv(1)*dphidu_uindz/(vsiz*uvsin)
     &  + 0.5d0*vw(1)*dphidw_uindz/vsiz

        dvyx_d=-uv(2)*dphidu_uindx/(vsiz*uvsin)
     &  + 0.5d0*vw(2)*dphidw_uindx/vsiz
        dvyy_d=-uv(2)*dphidu_uindy/(vsiz*uvsin)
     &  + 0.5d0*vw(2)*dphidw_uindy/vsiz
        dvyz_d=-uv(2)*dphidu_uindz/(vsiz*uvsin)
     &  + 0.5d0*vw(2)*dphidw_uindz/vsiz

        dvzx_d=-uv(3)*dphidu_uindx/(vsiz*uvsin)
     &  + 0.5d0*vw(3)*dphidw_uindx/vsiz
        dvzy_d=-uv(3)*dphidu_uindy/(vsiz*uvsin)
     &  + 0.5d0*vw(3)*dphidw_uindy/vsiz
        dvzz_d=-uv(3)*dphidu_uindz/(vsiz*uvsin)
     &  + 0.5d0*vw(3)*dphidw_uindz/vsiz

        dvxx_p=-uv(1)*dphidu_uinpx/(vsiz*uvsin)
     &  + 0.5d0*vw(1)*dphidw_uinpx/vsiz
        dvxy_p=-uv(1)*dphidu_uinpy/(vsiz*uvsin)
     &  + 0.5d0*vw(1)*dphidw_uinpy/vsiz
        dvxz_p=-uv(1)*dphidu_uinpz/(vsiz*uvsin)
     &  + 0.5d0*vw(1)*dphidw_uinpz/vsiz

        dvyx_p=-uv(2)*dphidu_uinpx/(vsiz*uvsin)
     &  + 0.5d0*vw(2)*dphidw_uinpx/vsiz
        dvyy_p=-uv(2)*dphidu_uinpy/(vsiz*uvsin)
     &  + 0.5d0*vw(2)*dphidw_uinpy/vsiz
        dvyz_p=-uv(2)*dphidu_uinpz/(vsiz*uvsin)
     &  + 0.5d0*vw(2)*dphidw_uinpz/vsiz

        dvzx_p=-uv(3)*dphidu_uinpx/(vsiz*uvsin)
     &    + 0.5d0*vw(3)*dphidw_uinpx/vsiz
        dvzy_p=-uv(3)*dphidu_uinpy/(vsiz*uvsin)
     &    + 0.5d0*vw(3)*dphidw_uinpy/vsiz
        dvzz_p=-uv(3)*dphidu_uinpz/(vsiz*uvsin)
     &    + 0.5d0*vw(3)*dphidw_uinpz/vsiz
c           done_ia=.false.
c           done_ic=.false.
c           done_ib=.false.
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv

c                nlocaltau1=nlocaltau1+1
c                ilocaltau1(1,nlocaltau1)=ia
c                ilocaltau1(2,nlocaltau1)=k
                elsttau1_index_tmp(1,kkk,iter3)=ia
                !elsttau1_tmp(k,ia)=nlocaltau1
c                elsttau1(1,kkk,i)=nlocaltau1

c                taulocal1d(1,nlocaltau1)=duxx_d
c                taulocal1d(2,nlocaltau1)=duxy_d
c                taulocal1d(3,nlocaltau1)=duxz_d
c                taulocal1d(4,nlocaltau1)=duyx_d
c                taulocal1d(5,nlocaltau1)=duyy_d
c                taulocal1d(6,nlocaltau1)=duyz_d
c                taulocal1d(7,nlocaltau1)=duzx_d
c                taulocal1d(8,nlocaltau1)=duzy_d
c                taulocal1d(9,nlocaltau1)=duzz_d
c
c                taulocal1p(1,nlocaltau1)=duxx_p
c                taulocal1p(2,nlocaltau1)=duxy_p
c                taulocal1p(3,nlocaltau1)=duxz_p
c                taulocal1p(4,nlocaltau1)=duyx_p
c                taulocal1p(5,nlocaltau1)=duyy_p
c                taulocal1p(6,nlocaltau1)=duyz_p
c                taulocal1p(7,nlocaltau1)=duzx_p
c                taulocal1p(8,nlocaltau1)=duzy_p
c                taulocal1p(9,nlocaltau1)=duzz_p

                taud1_3tmp(1,1,kkk,iter3)=duxx_d
                taud1_3tmp(2,1,kkk,iter3)=duxy_d
                taud1_3tmp(3,1,kkk,iter3)=duxz_d
                taud1_3tmp(4,1,kkk,iter3)=duyx_d
                taud1_3tmp(5,1,kkk,iter3)=duyy_d
                taud1_3tmp(6,1,kkk,iter3)=duyz_d
                taud1_3tmp(7,1,kkk,iter3)=duzx_d
                taud1_3tmp(8,1,kkk,iter3)=duzy_d
                taud1_3tmp(9,1,kkk,iter3)=duzz_d

                taup1_3tmp(1,1,kkk,iter3)=duxx_p
                taup1_3tmp(2,1,kkk,iter3)=duxy_p
                taup1_3tmp(3,1,kkk,iter3)=duxz_p
                taup1_3tmp(4,1,kkk,iter3)=duyx_p
                taup1_3tmp(5,1,kkk,iter3)=duyy_p
                taup1_3tmp(6,1,kkk,iter3)=duyz_p
                taup1_3tmp(7,1,kkk,iter3)=duzx_p
                taup1_3tmp(8,1,kkk,iter3)=duzy_p
                taup1_3tmp(9,1,kkk,iter3)=duzz_p


c                nlocaltau1=nlocaltau1+1
c                ilocaltau1(1,nlocaltau1)=ic
c                ilocaltau1(2,nlocaltau1)=k
                elsttau1_index_tmp(2,kkk,iter3)=ic
                !elsttau1_tmp(k,ia)=nlocaltau1
c                elsttau1(2,kkk,i)=nlocaltau1

c                taulocal1d(1,nlocaltau1)=dvxx_d
c                taulocal1d(2,nlocaltau1)=dvxy_d
c                taulocal1d(3,nlocaltau1)=dvxz_d
c                taulocal1d(4,nlocaltau1)=dvyx_d
c                taulocal1d(5,nlocaltau1)=dvyy_d
c                taulocal1d(6,nlocaltau1)=dvyz_d
c                taulocal1d(7,nlocaltau1)=dvzx_d
c                taulocal1d(8,nlocaltau1)=dvzy_d
c                taulocal1d(9,nlocaltau1)=dvzz_d
c
c                taulocal1p(1,nlocaltau1)=dvxx_p
c                taulocal1p(2,nlocaltau1)=dvxy_p
c                taulocal1p(3,nlocaltau1)=dvxz_p
c                taulocal1p(4,nlocaltau1)=dvyx_p
c                taulocal1p(5,nlocaltau1)=dvyy_p
c                taulocal1p(6,nlocaltau1)=dvyz_p
c                taulocal1p(7,nlocaltau1)=dvzx_p
c                taulocal1p(8,nlocaltau1)=dvzy_p
c                taulocal1p(9,nlocaltau1)=dvzz_p

                taud1_3tmp(1,2,kkk,iter3)=dvxx_d
                taud1_3tmp(2,2,kkk,iter3)=dvxy_d
                taud1_3tmp(3,2,kkk,iter3)=dvxz_d
                taud1_3tmp(4,2,kkk,iter3)=dvyx_d
                taud1_3tmp(5,2,kkk,iter3)=dvyy_d
                taud1_3tmp(6,2,kkk,iter3)=dvyz_d
                taud1_3tmp(7,2,kkk,iter3)=dvzx_d
                taud1_3tmp(8,2,kkk,iter3)=dvzy_d
                taud1_3tmp(9,2,kkk,iter3)=dvzz_d

                taup1_3tmp(1,2,kkk,iter3)=dvxx_p
                taup1_3tmp(2,2,kkk,iter3)=dvxy_p
                taup1_3tmp(3,2,kkk,iter3)=dvxz_p
                taup1_3tmp(4,2,kkk,iter3)=dvyx_p
                taup1_3tmp(5,2,kkk,iter3)=dvyy_p
                taup1_3tmp(6,2,kkk,iter3)=dvyz_p
                taup1_3tmp(7,2,kkk,iter3)=dvzx_p
                taup1_3tmp(8,2,kkk,iter3)=dvzy_p
                taup1_3tmp(9,2,kkk,iter3)=dvzz_p

c                nlocaltau1=nlocaltau1+1
c                ilocaltau1(1,nlocaltau1)=ib
c                ilocaltau1(2,nlocaltau1)=k
                elsttau1_index_tmp(3,kkk,iter3)=ib
                !elsttau1_tmp(k,ia)=nlocaltau1
c                elsttau1(3,kkk,i)=nlocaltau1

c                taulocal1d(1,nlocaltau1)=-duxx_d-dvxx_d
c                taulocal1d(2,nlocaltau1)=-duxy_d-dvxy_d
c                taulocal1d(3,nlocaltau1)=-duxz_d-dvxz_d
c                taulocal1d(4,nlocaltau1)=-duyx_d-dvyx_d
c                taulocal1d(5,nlocaltau1)=-duyy_d-dvyy_d
c                taulocal1d(6,nlocaltau1)=-duyz_d-dvyz_d
c                taulocal1d(7,nlocaltau1)=-duzx_d-dvzx_d
c                taulocal1d(8,nlocaltau1)=-duzy_d-dvzy_d
c                taulocal1d(9,nlocaltau1)=-duzz_d-dvzz_d
c
c                taulocal1p(1,nlocaltau1)=-duxx_p-dvxx_p
c                taulocal1p(2,nlocaltau1)=-duxy_p-dvxy_p
c                taulocal1p(3,nlocaltau1)=-duxz_p-dvxz_p
c                taulocal1p(4,nlocaltau1)=-duyx_p-dvyx_p
c                taulocal1p(5,nlocaltau1)=-duyy_p-dvyy_p
c                taulocal1p(6,nlocaltau1)=-duyz_p-dvyz_p
c                taulocal1p(7,nlocaltau1)=-duzx_p-dvzx_p
c                taulocal1p(8,nlocaltau1)=-duzy_p-dvzy_p
c                taulocal1p(9,nlocaltau1)=-duzz_p-dvzz_p

                taud1_3tmp(1,3,kkk,iter3)=-duxx_d-dvxx_d
                taud1_3tmp(2,3,kkk,iter3)=-duxy_d-dvxy_d
                taud1_3tmp(3,3,kkk,iter3)=-duxz_d-dvxz_d
                taud1_3tmp(4,3,kkk,iter3)=-duyx_d-dvyx_d
                taud1_3tmp(5,3,kkk,iter3)=-duyy_d-dvyy_d
                taud1_3tmp(6,3,kkk,iter3)=-duyz_d-dvyz_d
                taud1_3tmp(7,3,kkk,iter3)=-duzx_d-dvzx_d
                taud1_3tmp(8,3,kkk,iter3)=-duzy_d-dvzy_d
                taud1_3tmp(9,3,kkk,iter3)=-duzz_d-dvzz_d

                taup1_3tmp(1,3,kkk,iter3)=-duxx_p-dvxx_p
                taup1_3tmp(2,3,kkk,iter3)=-duxy_p-dvxy_p
                taup1_3tmp(3,3,kkk,iter3)=-duxz_p-dvxz_p
                taup1_3tmp(4,3,kkk,iter3)=-duyx_p-dvyx_p
                taup1_3tmp(5,3,kkk,iter3)=-duyy_p-dvyy_p
                taup1_3tmp(6,3,kkk,iter3)=-duyz_p-dvyz_p
                taup1_3tmp(7,3,kkk,iter3)=-duzx_p-dvzx_p
                taup1_3tmp(8,3,kkk,iter3)=-duzy_p-dvzy_p
                taup1_3tmp(9,3,kkk,iter3)=-duzz_p-dvzz_p


c                frcztaulocal1d(1,nlocal)=duxx_d
c                frcztaulocal1d(2,nlocal)=duxy_d
c                frcztaulocal1d(3,nlocal)=duxz_d
c                frcztaulocal1d(4,nlocal)=duyx_d
c                frcztaulocal1d(5,nlocal)=duyy_d
c                frcztaulocal1d(6,nlocal)=duyz_d
c                frcztaulocal1d(7,nlocal)=duzx_d
c                frcztaulocal1d(8,nlocal)=duzy_d
c                frcztaulocal1d(9,nlocal)=duzz_d
c                frcztaulocal1p(1,nlocal)=duxx_p
c                frcztaulocal1p(2,nlocal)=duxy_p
c                frcztaulocal1p(3,nlocal)=duxz_p
c                frcztaulocal1p(4,nlocal)=duyx_p
c                frcztaulocal1p(5,nlocal)=duyy_p
c                frcztaulocal1p(6,nlocal)=duyz_p
c                frcztaulocal1p(7,nlocal)=duzx_p
c                frcztaulocal1p(8,nlocal)=duzy_p
c                frcztaulocal1p(9,nlocal)=duzz_p

                frcztau1dtot_3tmp(1,kkk,iter3)=duxx_d
                frcztau1dtot_3tmp(2,kkk,iter3)=duxy_d
                frcztau1dtot_3tmp(3,kkk,iter3)=duxz_d
                frcztau1dtot_3tmp(4,kkk,iter3)=duyx_d
                frcztau1dtot_3tmp(5,kkk,iter3)=duyy_d
                frcztau1dtot_3tmp(6,kkk,iter3)=duyz_d
                frcztau1dtot_3tmp(7,kkk,iter3)=duzx_d
                frcztau1dtot_3tmp(8,kkk,iter3)=duzy_d
                frcztau1dtot_3tmp(9,kkk,iter3)=duzz_d
                frcztau1ptot_3tmp(1,kkk,iter3)=duxx_p
                frcztau1ptot_3tmp(2,kkk,iter3)=duxy_p
                frcztau1ptot_3tmp(3,kkk,iter3)=duxz_p
                frcztau1ptot_3tmp(4,kkk,iter3)=duyx_p
                frcztau1ptot_3tmp(5,kkk,iter3)=duyy_p
                frcztau1ptot_3tmp(6,kkk,iter3)=duyz_p
                frcztau1ptot_3tmp(7,kkk,iter3)=duzx_p
                frcztau1ptot_3tmp(8,kkk,iter3)=duzy_p
                frcztau1ptot_3tmp(9,kkk,iter3)=duzz_p


c                frcxtaulocal1d(1,nlocal)=dvxx_d
c                frcxtaulocal1d(2,nlocal)=dvxy_d
c                frcxtaulocal1d(3,nlocal)=dvxz_d
c                frcxtaulocal1d(4,nlocal)=dvyx_d
c                frcxtaulocal1d(5,nlocal)=dvyy_d
c                frcxtaulocal1d(6,nlocal)=dvyz_d
c                frcxtaulocal1d(7,nlocal)=dvzx_d
c                frcxtaulocal1d(8,nlocal)=dvzy_d
c                frcxtaulocal1d(9,nlocal)=dvzz_d
c                frcxtaulocal1p(1,nlocal)=dvxx_p
c                frcxtaulocal1p(2,nlocal)=dvxy_p
c                frcxtaulocal1p(3,nlocal)=dvxz_p
c                frcxtaulocal1p(4,nlocal)=dvyx_p
c                frcxtaulocal1p(5,nlocal)=dvyy_p
c                frcxtaulocal1p(6,nlocal)=dvyz_p
c                frcxtaulocal1p(7,nlocal)=dvzx_p
c                frcxtaulocal1p(8,nlocal)=dvzy_p
c                frcxtaulocal1p(9,nlocal)=dvzz_p

                frcxtau1dtot_3tmp(1,kkk,iter3)=dvxx_d
                frcxtau1dtot_3tmp(2,kkk,iter3)=dvxy_d
                frcxtau1dtot_3tmp(3,kkk,iter3)=dvxz_d
                frcxtau1dtot_3tmp(4,kkk,iter3)=dvyx_d
                frcxtau1dtot_3tmp(5,kkk,iter3)=dvyy_d
                frcxtau1dtot_3tmp(6,kkk,iter3)=dvyz_d
                frcxtau1dtot_3tmp(7,kkk,iter3)=dvzx_d
                frcxtau1dtot_3tmp(8,kkk,iter3)=dvzy_d
                frcxtau1dtot_3tmp(9,kkk,iter3)=dvzz_d
                frcxtau1ptot_3tmp(1,kkk,iter3)=dvxx_p
                frcxtau1ptot_3tmp(2,kkk,iter3)=dvxy_p
                frcxtau1ptot_3tmp(3,kkk,iter3)=dvxz_p
                frcxtau1ptot_3tmp(4,kkk,iter3)=dvyx_p
                frcxtau1ptot_3tmp(5,kkk,iter3)=dvyy_p
                frcxtau1ptot_3tmp(6,kkk,iter3)=dvyz_p
                frcxtau1ptot_3tmp(7,kkk,iter3)=dvzx_p
                frcxtau1ptot_3tmp(8,kkk,iter3)=dvzy_p
                frcxtau1ptot_3tmp(9,kkk,iter3)=dvzz_p
                
c                frcytaulocal1d(1,nlocal)=0.0d0
c                frcytaulocal1d(2,nlocal)=0.0d0
c                frcytaulocal1d(3,nlocal)=0.0d0
c                frcytaulocal1d(4,nlocal)=0.0d0
c                frcytaulocal1d(5,nlocal)=0.0d0
c                frcytaulocal1d(6,nlocal)=0.0d0
c                frcytaulocal1d(7,nlocal)=0.0d0
c                frcytaulocal1d(8,nlocal)=0.0d0
c                frcytaulocal1d(9,nlocal)=0.0d0
c                frcytaulocal1p(1,nlocal)=0.0d0
c                frcytaulocal1p(2,nlocal)=0.0d0
c                frcytaulocal1p(3,nlocal)=0.0d0
c                frcytaulocal1p(4,nlocal)=0.0d0
c                frcytaulocal1p(5,nlocal)=0.0d0
c                frcytaulocal1p(6,nlocal)=0.0d0
c                frcytaulocal1p(7,nlocal)=0.0d0
c                frcytaulocal1p(8,nlocal)=0.0d0
c                frcytaulocal1p(9,nlocal)=0.0d0

                frcytau1dtot_3tmp(1,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(2,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(3,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(4,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(5,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(6,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(7,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(8,kkk,iter3)=0.0d0
                frcytau1dtot_3tmp(9,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(1,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(2,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(3,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(4,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(5,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(6,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(7,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(8,kkk,iter3)=0.0d0
                frcytau1ptot_3tmp(9,kkk,iter3)=0.0d0

                elsttau1_index_tmp(4,kkk,iter3)=0

c
c     force distribution for the Z-Bisect local coordinate method
c
      else if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
            dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &              / (vsiz*(ut1sin+ut2sin))
            dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &              / (wsiz*(ut1sin+ut2sin))
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,id) = depi(j,id) + dw
c            depi(j,ib) = depi(j,ib) - du - dv - dw
c            frcz(j) = frcz(j) + du
c            frcx(j) = frcx(j) + dv
c            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
      else if (axetyp .eq. '3-Fold') then
         do j = 1, 3
            du = uw(j)*dphidw/(usiz*uwsin)
     &              + uv(j)*dphidv/(usiz*uvsin)
     &              - uw(j)*dphidu/(usiz*uwsin)
     &              - uv(j)*dphidu/(usiz*uvsin)
            dv = vw(j)*dphidw/(vsiz*vwsin)
     &              - uv(j)*dphidu/(vsiz*uvsin)
     &              - vw(j)*dphidv/(vsiz*vwsin)
     &              + uv(j)*dphidv/(vsiz*uvsin)
            dw = -uw(j)*dphidu/(wsiz*uwsin)
     &              - vw(j)*dphidv/(wsiz*vwsin)
     &              + uw(j)*dphidw/(wsiz*uwsin)
     &              + vw(j)*dphidw/(wsiz*vwsin)
            du = du / 3.0d0
            dv = dv / 3.0d0
            dw = dw / 3.0d0
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,id) = depi(j,id) + dw
c            depi(j,ib) = depi(j,ib) - du - dv - dw
c            frcz(j) = frcz(j) + du
c            frcx(j) = frcx(j) + dv
c            frcy(j) = frcy(j) + dw
         end do
      end if
      return
      end
c
c
