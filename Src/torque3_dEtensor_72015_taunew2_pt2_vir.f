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
      subroutine torque3_dEtensor_72015_taunew2_pt2_vir(i,trq1,trq2,
     &  frcx,frcy,frcz,demi,taud,taup,nlocaltau2,ilocaltau2,
     &  taulocal2d,taulocal2p,
     &  nlocal,frcztaulocal2d,frcztaulocal2p,
     &  frcxtaulocal2d,frcxtaulocal2p,frcytaulocal2d,frcytaulocal2p,k)
      use sizes
      use atoms
      use deriv
      use mpole
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
      real*8 demi(3,*)
c      real*8 depi(3,*)
      character*8 axetyp
      real*8 taud(9),taup(9)
      integer nlocaltau2,ilocaltau2(2,*),ikflag,k,m
      real*8 taulocal2d(9,*)
      real*8 taulocal2p(9,*)
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
      integer nlocal
      real*8 frcztaulocal2d(9,*),frcztaulocal2p(9,*)
      real*8 frcxtaulocal2d(9,*),frcxtaulocal2p(9,*)
      real*8 frcytaulocal2d(9,*),frcytaulocal2p(9,*)
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
      ia = zaxis(k)
      ib = ipole(k)
      ic = xaxis(k)
      id = yaxis(k)
      axetyp = polaxe(k)
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
      dphidu = -trq1(1)*u(1) - trq1(2)*u(2) - trq1(3)*u(3)
      dphidv = -trq1(1)*v(1) - trq1(2)*v(2) - trq1(3)*v(3)
      dphidw = -trq1(1)*w(1) - trq1(2)*w(2) - trq1(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq1(1)*r(1) - trq1(2)*r(2) - trq1(3)*r(3)
         dphids = -trq1(1)*s(1) - trq1(2)*s(2) - trq1(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            demi(j,ia) = demi(j,ia) + du
            demi(j,ib) = demi(j,ib) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            demi(j,ia) = demi(j,ia) + du
            demi(j,ic) = demi(j,ic) + dv
            demi(j,ib) = demi(j,ib) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
            demi(j,ia) = demi(j,ia) + du
            demi(j,ic) = demi(j,ic) + dv
            demi(j,ib) = demi(j,ib) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
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
            demi(j,ia) = demi(j,ia) + du
            demi(j,ic) = demi(j,ic) + dv
            demi(j,id) = demi(j,id) + dw
            demi(j,ib) = demi(j,ib) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
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
            demi(j,ia) = demi(j,ia) + du
            demi(j,ic) = demi(j,ic) + dv
            demi(j,id) = demi(j,id) + dw
            demi(j,ib) = demi(j,ib) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
      end if
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

         dvxx_d=-uv(1)*dphidu_uindx/(vsiz*uvsin)
         dvxy_d=-uv(1)*dphidu_uindy/(vsiz*uvsin)
         dvxz_d=-uv(1)*dphidu_uindz/(vsiz*uvsin)

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
      

                nlocaltau2=nlocaltau2+1
                ilocaltau2(1,nlocaltau2)=i
                ilocaltau2(2,nlocaltau2)=ia
                taulocal2d(1,nlocaltau2)=duxx_d
                taulocal2d(2,nlocaltau2)=duxy_d
                taulocal2d(3,nlocaltau2)=duxz_d
                taulocal2d(4,nlocaltau2)=duyx_d
                taulocal2d(5,nlocaltau2)=duyy_d
                taulocal2d(6,nlocaltau2)=duyz_d
                taulocal2d(7,nlocaltau2)=duzx_d
                taulocal2d(8,nlocaltau2)=duzy_d
                taulocal2d(9,nlocaltau2)=duzz_d

                taulocal2p(1,nlocaltau2)=duxx_p
                taulocal2p(2,nlocaltau2)=duxy_p
                taulocal2p(3,nlocaltau2)=duxz_p
                taulocal2p(4,nlocaltau2)=duyx_p
                taulocal2p(5,nlocaltau2)=duyy_p
                taulocal2p(6,nlocaltau2)=duyz_p
                taulocal2p(7,nlocaltau2)=duzx_p
                taulocal2p(8,nlocaltau2)=duzy_p
                taulocal2p(9,nlocaltau2)=duzz_p

                frcztaulocal2d(1,nlocal)=duxx_d
                frcztaulocal2d(2,nlocal)=duxy_d
                frcztaulocal2d(3,nlocal)=duxz_d
                frcztaulocal2d(4,nlocal)=duyx_d
                frcztaulocal2d(5,nlocal)=duyy_d
                frcztaulocal2d(6,nlocal)=duyz_d
                frcztaulocal2d(7,nlocal)=duzx_d
                frcztaulocal2d(8,nlocal)=duzy_d
                frcztaulocal2d(9,nlocal)=duzz_d
                frcztaulocal2p(1,nlocal)=duxx_p
                frcztaulocal2p(2,nlocal)=duxy_p
                frcztaulocal2p(3,nlocal)=duxz_p
                frcztaulocal2p(4,nlocal)=duyx_p
                frcztaulocal2p(5,nlocal)=duyy_p
                frcztaulocal2p(6,nlocal)=duyz_p
                frcztaulocal2p(7,nlocal)=duzx_p
                frcztaulocal2p(8,nlocal)=duzy_p
                frcztaulocal2p(9,nlocal)=duzz_p

                frcxtaulocal2d(1,nlocal)=dvxx_d
                frcxtaulocal2d(2,nlocal)=dvxy_d
                frcxtaulocal2d(3,nlocal)=dvxz_d
                frcxtaulocal2d(4,nlocal)=dvyx_d
                frcxtaulocal2d(5,nlocal)=dvyy_d
                frcxtaulocal2d(6,nlocal)=dvyz_d
                frcxtaulocal2d(7,nlocal)=dvzx_d
                frcxtaulocal2d(8,nlocal)=dvzy_d
                frcxtaulocal2d(9,nlocal)=dvzz_d
                frcxtaulocal2p(1,nlocal)=dvxx_p
                frcxtaulocal2p(2,nlocal)=dvxy_p
                frcxtaulocal2p(3,nlocal)=dvxz_p
                frcxtaulocal2p(4,nlocal)=dvyx_p
                frcxtaulocal2p(5,nlocal)=dvyy_p
                frcxtaulocal2p(6,nlocal)=dvyz_p
                frcxtaulocal2p(7,nlocal)=dvzx_p
                frcxtaulocal2p(8,nlocal)=dvzy_p
                frcxtaulocal2p(9,nlocal)=dvzz_p

                frcytaulocal2d(1,nlocal)=0.0d0
                frcytaulocal2d(2,nlocal)=0.0d0
                frcytaulocal2d(3,nlocal)=0.0d0
                frcytaulocal2d(4,nlocal)=0.0d0
                frcytaulocal2d(5,nlocal)=0.0d0
                frcytaulocal2d(6,nlocal)=0.0d0
                frcytaulocal2d(7,nlocal)=0.0d0
                frcytaulocal2d(8,nlocal)=0.0d0
                frcytaulocal2d(9,nlocal)=0.0d0
                frcytaulocal2p(1,nlocal)=0.0d0
                frcytaulocal2p(2,nlocal)=0.0d0
                frcytaulocal2p(3,nlocal)=0.0d0
                frcytaulocal2p(4,nlocal)=0.0d0
                frcytaulocal2p(5,nlocal)=0.0d0
                frcytaulocal2p(6,nlocal)=0.0d0
                frcytaulocal2p(7,nlocal)=0.0d0
                frcytaulocal2p(8,nlocal)=0.0d0
                frcytaulocal2p(9,nlocal)=0.0d0
 
                nlocaltau2=nlocaltau2+1
                ilocaltau2(1,nlocaltau2)=i
                ilocaltau2(2,nlocaltau2)=ic
                taulocal2d(1,nlocaltau2)=dvxx_d
                taulocal2d(2,nlocaltau2)=dvxy_d
                taulocal2d(3,nlocaltau2)=dvxz_d
                taulocal2d(4,nlocaltau2)=dvyx_d
                taulocal2d(5,nlocaltau2)=dvyy_d
                taulocal2d(6,nlocaltau2)=dvyz_d
                taulocal2d(7,nlocaltau2)=dvzx_d
                taulocal2d(8,nlocaltau2)=dvzy_d
                taulocal2d(9,nlocaltau2)=dvzz_d

                taulocal2p(1,nlocaltau2)=dvxx_p
                taulocal2p(2,nlocaltau2)=dvxy_p
                taulocal2p(3,nlocaltau2)=dvxz_p
                taulocal2p(4,nlocaltau2)=dvyx_p
                taulocal2p(5,nlocaltau2)=dvyy_p
                taulocal2p(6,nlocaltau2)=dvyz_p
                taulocal2p(7,nlocaltau2)=dvzx_p
                taulocal2p(8,nlocaltau2)=dvzy_p
                taulocal2p(9,nlocaltau2)=dvzz_p


                nlocaltau2=nlocaltau2+1
                ilocaltau2(1,nlocaltau2)=i
                ilocaltau2(2,nlocaltau2)=ib
                taulocal2d(1,nlocaltau2)=-duxx_d-dvxx_d
                taulocal2d(2,nlocaltau2)=-duxy_d-dvxy_d
                taulocal2d(3,nlocaltau2)=-duxz_d-dvxz_d
                taulocal2d(4,nlocaltau2)=-duyx_d-dvyx_d
                taulocal2d(5,nlocaltau2)=-duyy_d-dvyy_d
                taulocal2d(6,nlocaltau2)=-duyz_d-dvyz_d
                taulocal2d(7,nlocaltau2)=-duzx_d-dvzx_d
                taulocal2d(8,nlocaltau2)=-duzy_d-dvzy_d
                taulocal2d(9,nlocaltau2)=-duzz_d-dvzz_d

                taulocal2p(1,nlocaltau2)=-duxx_p-dvxx_p
                taulocal2p(2,nlocaltau2)=-duxy_p-dvxy_p
                taulocal2p(3,nlocaltau2)=-duxz_p-dvxz_p
                taulocal2p(4,nlocaltau2)=-duyx_p-dvyx_p
                taulocal2p(5,nlocaltau2)=-duyy_p-dvyy_p
                taulocal2p(6,nlocaltau2)=-duyz_p-dvyz_p
                taulocal2p(7,nlocaltau2)=-duzx_p-dvzx_p
                taulocal2p(8,nlocaltau2)=-duzy_p-dvzy_p
                taulocal2p(9,nlocaltau2)=-duzz_p-dvzz_p
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
     &          + 0.5d0*uw(1)*dphidw_uindx/usiz
        duxy_d=uv(1)*dphidv_uindy/(usiz*uvsin) 
     &          + 0.5d0*uw(1)*dphidw_uindy/usiz
        duxz_d=uv(1)*dphidv_uindz/(usiz*uvsin) 
     &          + 0.5d0*uw(1)*dphidw_uindz/usiz

        duyx_d=uv(2)*dphidv_uindx/(usiz*uvsin) 
     &          + 0.5d0*uw(2)*dphidw_uindx/usiz
        duyy_d=uv(2)*dphidv_uindy/(usiz*uvsin) 
     &          + 0.5d0*uw(2)*dphidw_uindy/usiz
        duyz_d=uv(2)*dphidv_uindz/(usiz*uvsin) 
     &          + 0.5d0*uw(2)*dphidw_uindz/usiz

        duzx_d=uv(3)*dphidv_uindx/(usiz*uvsin) 
     &          + 0.5d0*uw(3)*dphidw_uindx/usiz
        duzy_d=uv(3)*dphidv_uindy/(usiz*uvsin) 
     &          + 0.5d0*uw(3)*dphidw_uindy/usiz
        duzz_d=uv(3)*dphidv_uindz/(usiz*uvsin) 
     &          + 0.5d0*uw(3)*dphidw_uindz/usiz

        duxx_p=uv(1)*dphidv_uinpx/(usiz*uvsin) 
     &          + 0.5d0*uw(1)*dphidw_uinpx/usiz
        duxy_p=uv(1)*dphidv_uinpy/(usiz*uvsin) 
     &          + 0.5d0*uw(1)*dphidw_uinpy/usiz
        duxz_p=uv(1)*dphidv_uinpz/(usiz*uvsin) 
     &          + 0.5d0*uw(1)*dphidw_uinpz/usiz

        duyx_p=uv(2)*dphidv_uinpx/(usiz*uvsin) 
     &          + 0.5d0*uw(2)*dphidw_uinpx/usiz
        duyy_p=uv(2)*dphidv_uinpy/(usiz*uvsin) 
     &          + 0.5d0*uw(2)*dphidw_uinpy/usiz
        duyz_p=uv(2)*dphidv_uinpz/(usiz*uvsin) 
     &          + 0.5d0*uw(2)*dphidw_uinpz/usiz

        duzx_p=uv(3)*dphidv_uinpx/(usiz*uvsin) 
     &          + 0.5d0*uw(3)*dphidw_uinpx/usiz
        duzy_p=uv(3)*dphidv_uinpy/(usiz*uvsin) 
     &          + 0.5d0*uw(3)*dphidw_uinpy/usiz
        duzz_p=uv(3)*dphidv_uinpz/(usiz*uvsin) 
     &          + 0.5d0*uw(3)*dphidw_uinpz/usiz

        dvxx_d=-uv(1)*dphidu_uindx/(vsiz*uvsin)
     &          + 0.5d0*vw(1)*dphidw_uindx/vsiz
        dvxy_d=-uv(1)*dphidu_uindy/(vsiz*uvsin)
     &          + 0.5d0*vw(1)*dphidw_uindy/vsiz
        dvxz_d=-uv(1)*dphidu_uindz/(vsiz*uvsin)
     &          + 0.5d0*vw(1)*dphidw_uindz/vsiz

        dvyx_d=-uv(2)*dphidu_uindx/(vsiz*uvsin)
     &          + 0.5d0*vw(2)*dphidw_uindx/vsiz
        dvyy_d=-uv(2)*dphidu_uindy/(vsiz*uvsin)
     &          + 0.5d0*vw(2)*dphidw_uindy/vsiz
        dvyz_d=-uv(2)*dphidu_uindz/(vsiz*uvsin)
     &          + 0.5d0*vw(2)*dphidw_uindz/vsiz

        dvzx_d=-uv(3)*dphidu_uindx/(vsiz*uvsin)
     &          + 0.5d0*vw(3)*dphidw_uindx/vsiz
        dvzy_d=-uv(3)*dphidu_uindy/(vsiz*uvsin)
     &          + 0.5d0*vw(3)*dphidw_uindy/vsiz
        dvzz_d=-uv(3)*dphidu_uindz/(vsiz*uvsin)
     &          + 0.5d0*vw(3)*dphidw_uindz/vsiz

        dvxx_p=-uv(1)*dphidu_uinpx/(vsiz*uvsin)
     &          + 0.5d0*vw(1)*dphidw_uinpx/vsiz
        dvxy_p=-uv(1)*dphidu_uinpy/(vsiz*uvsin)
     &          + 0.5d0*vw(1)*dphidw_uinpy/vsiz
        dvxz_p=-uv(1)*dphidu_uinpz/(vsiz*uvsin)
     &          + 0.5d0*vw(1)*dphidw_uinpz/vsiz

        dvyx_p=-uv(2)*dphidu_uinpx/(vsiz*uvsin)
     &          + 0.5d0*vw(2)*dphidw_uinpx/vsiz
        dvyy_p=-uv(2)*dphidu_uinpy/(vsiz*uvsin)
     &          + 0.5d0*vw(2)*dphidw_uinpy/vsiz
        dvyz_p=-uv(2)*dphidu_uinpz/(vsiz*uvsin)
     &          + 0.5d0*vw(2)*dphidw_uinpz/vsiz

        dvzx_p=-uv(3)*dphidu_uinpx/(vsiz*uvsin)
     &          + 0.5d0*vw(3)*dphidw_uinpx/vsiz
        dvzy_p=-uv(3)*dphidu_uinpy/(vsiz*uvsin)
     &          + 0.5d0*vw(3)*dphidw_uinpy/vsiz
        dvzz_p=-uv(3)*dphidu_uinpz/(vsiz*uvsin)
     &          + 0.5d0*vw(3)*dphidw_uinpz/vsiz

c           done_ia=.false.
c           done_ic=.false.
c           done_ib=.false.
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv

                nlocaltau2=nlocaltau2+1
                ilocaltau2(1,nlocaltau2)=i
                ilocaltau2(2,nlocaltau2)=ia
                taulocal2d(1,nlocaltau2)=duxx_d
                taulocal2d(2,nlocaltau2)=duxy_d
                taulocal2d(3,nlocaltau2)=duxz_d
                taulocal2d(4,nlocaltau2)=duyx_d
                taulocal2d(5,nlocaltau2)=duyy_d
                taulocal2d(6,nlocaltau2)=duyz_d
                taulocal2d(7,nlocaltau2)=duzx_d
                taulocal2d(8,nlocaltau2)=duzy_d
                taulocal2d(9,nlocaltau2)=duzz_d

                taulocal2p(1,nlocaltau2)=duxx_p
                taulocal2p(2,nlocaltau2)=duxy_p
                taulocal2p(3,nlocaltau2)=duxz_p
                taulocal2p(4,nlocaltau2)=duyx_p
                taulocal2p(5,nlocaltau2)=duyy_p
                taulocal2p(6,nlocaltau2)=duyz_p
                taulocal2p(7,nlocaltau2)=duzx_p
                taulocal2p(8,nlocaltau2)=duzy_p
                taulocal2p(9,nlocaltau2)=duzz_p


                nlocaltau2=nlocaltau2+1
                ilocaltau2(1,nlocaltau2)=i
                ilocaltau2(2,nlocaltau2)=ic
                taulocal2d(1,nlocaltau2)=dvxx_d
                taulocal2d(2,nlocaltau2)=dvxy_d
                taulocal2d(3,nlocaltau2)=dvxz_d
                taulocal2d(4,nlocaltau2)=dvyx_d
                taulocal2d(5,nlocaltau2)=dvyy_d
                taulocal2d(6,nlocaltau2)=dvyz_d
                taulocal2d(7,nlocaltau2)=dvzx_d
                taulocal2d(8,nlocaltau2)=dvzy_d
                taulocal2d(9,nlocaltau2)=dvzz_d

                taulocal2p(1,nlocaltau2)=dvxx_p
                taulocal2p(2,nlocaltau2)=dvxy_p
                taulocal2p(3,nlocaltau2)=dvxz_p
                taulocal2p(4,nlocaltau2)=dvyx_p
                taulocal2p(5,nlocaltau2)=dvyy_p
                taulocal2p(6,nlocaltau2)=dvyz_p
                taulocal2p(7,nlocaltau2)=dvzx_p
                taulocal2p(8,nlocaltau2)=dvzy_p
                taulocal2p(9,nlocaltau2)=dvzz_p


                nlocaltau2=nlocaltau2+1
                ilocaltau2(1,nlocaltau2)=i
                ilocaltau2(2,nlocaltau2)=ib
                taulocal2d(1,nlocaltau2)=-duxx_d-dvxx_d
                taulocal2d(2,nlocaltau2)=-duxy_d-dvxy_d
                taulocal2d(3,nlocaltau2)=-duxz_d-dvxz_d
                taulocal2d(4,nlocaltau2)=-duyx_d-dvyx_d
                taulocal2d(5,nlocaltau2)=-duyy_d-dvyy_d
                taulocal2d(6,nlocaltau2)=-duyz_d-dvyz_d
                taulocal2d(7,nlocaltau2)=-duzx_d-dvzx_d
                taulocal2d(8,nlocaltau2)=-duzy_d-dvzy_d
                taulocal2d(9,nlocaltau2)=-duzz_d-dvzz_d

                taulocal2p(1,nlocaltau2)=-duxx_p-dvxx_p
                taulocal2p(2,nlocaltau2)=-duxy_p-dvxy_p
                taulocal2p(3,nlocaltau2)=-duxz_p-dvxz_p
                taulocal2p(4,nlocaltau2)=-duyx_p-dvyx_p
                taulocal2p(5,nlocaltau2)=-duyy_p-dvyy_p
                taulocal2p(6,nlocaltau2)=-duyz_p-dvyz_p
                taulocal2p(7,nlocaltau2)=-duzx_p-dvzx_p
                taulocal2p(8,nlocaltau2)=-duzy_p-dvzy_p
                taulocal2p(9,nlocaltau2)=-duzz_p-dvzz_p

                frcztaulocal2d(1,nlocal)=duxx_d
                frcztaulocal2d(2,nlocal)=duxy_d
                frcztaulocal2d(3,nlocal)=duxz_d
                frcztaulocal2d(4,nlocal)=duyx_d
                frcztaulocal2d(5,nlocal)=duyy_d
                frcztaulocal2d(6,nlocal)=duyz_d
                frcztaulocal2d(7,nlocal)=duzx_d
                frcztaulocal2d(8,nlocal)=duzy_d
                frcztaulocal2d(9,nlocal)=duzz_d
                frcztaulocal2p(1,nlocal)=duxx_p
                frcztaulocal2p(2,nlocal)=duxy_p
                frcztaulocal2p(3,nlocal)=duxz_p
                frcztaulocal2p(4,nlocal)=duyx_p
                frcztaulocal2p(5,nlocal)=duyy_p
                frcztaulocal2p(6,nlocal)=duyz_p
                frcztaulocal2p(7,nlocal)=duzx_p
                frcztaulocal2p(8,nlocal)=duzy_p
                frcztaulocal2p(9,nlocal)=duzz_p

                frcxtaulocal2d(1,nlocal)=dvxx_d
                frcxtaulocal2d(2,nlocal)=dvxy_d
                frcxtaulocal2d(3,nlocal)=dvxz_d
                frcxtaulocal2d(4,nlocal)=dvyx_d
                frcxtaulocal2d(5,nlocal)=dvyy_d
                frcxtaulocal2d(6,nlocal)=dvyz_d
                frcxtaulocal2d(7,nlocal)=dvzx_d
                frcxtaulocal2d(8,nlocal)=dvzy_d
                frcxtaulocal2d(9,nlocal)=dvzz_d
                frcxtaulocal2p(1,nlocal)=dvxx_p
                frcxtaulocal2p(2,nlocal)=dvxy_p
                frcxtaulocal2p(3,nlocal)=dvxz_p
                frcxtaulocal2p(4,nlocal)=dvyx_p
                frcxtaulocal2p(5,nlocal)=dvyy_p
                frcxtaulocal2p(6,nlocal)=dvyz_p
                frcxtaulocal2p(7,nlocal)=dvzx_p
                frcxtaulocal2p(8,nlocal)=dvzy_p
                frcxtaulocal2p(9,nlocal)=dvzz_p

                frcytaulocal2d(1,nlocal)=0.0d0
                frcytaulocal2d(2,nlocal)=0.0d0
                frcytaulocal2d(3,nlocal)=0.0d0
                frcytaulocal2d(4,nlocal)=0.0d0
                frcytaulocal2d(5,nlocal)=0.0d0
                frcytaulocal2d(6,nlocal)=0.0d0
                frcytaulocal2d(7,nlocal)=0.0d0
                frcytaulocal2d(8,nlocal)=0.0d0
                frcytaulocal2d(9,nlocal)=0.0d0
                frcytaulocal2p(1,nlocal)=0.0d0
                frcytaulocal2p(2,nlocal)=0.0d0
                frcytaulocal2p(3,nlocal)=0.0d0
                frcytaulocal2p(4,nlocal)=0.0d0
                frcytaulocal2p(5,nlocal)=0.0d0
                frcytaulocal2p(6,nlocal)=0.0d0
                frcytaulocal2p(7,nlocal)=0.0d0
                frcytaulocal2p(8,nlocal)=0.0d0
                frcytaulocal2p(9,nlocal)=0.0d0


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
