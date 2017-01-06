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
      subroutine torque3_dEtensorMod_pt2(i,trq1,trq2,frcx,frcy,frcz,
     &          demi,taud,taup,nlocal,ilocal,dlocal1d,
     &          dlocal1p,dlocal2d,dlocal2p,k)
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
      integer nlocal,ilocal(2,*),ikflag,k,m
      real*8 dlocal1d(9,*),dlocal1p(9,*),dlocal2d(9,*),dlocal2p(9,*)
      real*8 dphidux_d,dphiduy_d,dphiduz_d,dphidvx_d,dphidvy_d
      real*8 dphidvz_d,dphidwx_d,dphidwy_d,dphidwz_d
      real*8 dphidux_p,dphiduy_p,dphiduz_p,dphidvx_p,dphidvy_p
      real*8 dphidvz_p,dphidwx_p,dphidwy_p,dphidwz_p
      real*8 duxx_d,duxy_d,duxz_d,duyx_d,duyy_d,duyz_d,duzx_d,duzy_d
      real*8 duzz_d
      real*8 duxx_p,duxy_p,duxz_p,duyx_p,duyy_p,duyz_p,duzx_p,duzy_p
      real*8 duzz_p
      real*8 dvxx_d,dvxy_d,dvxz_d,dvyx_d,dvyy_d,dvyz_d,dvzx_d,dvzy_d
      real*8 dvzz_d
      real*8 dvxx_p,dvxy_p,dvxz_p,dvyx_p,dvyy_p,dvyz_p,dvzx_p,dvzy_p
      real*8 dvzz_p
      logical done_ia,done_ib,done_ic
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
      print*,"!!!IN PT2 NEW TENSOR TORQ"
      dphidu = -trq2(1)*u(1) - trq2(2)*u(2) - trq2(3)*u(3)
      dphidv = -trq2(1)*v(1) - trq2(2)*v(2) - trq2(3)*v(3)
      dphidw = -trq2(1)*w(1) - trq2(2)*w(2) - trq2(3)*w(3)

      dphidux_d= -(u(1)*taud(1) + u(2)*taud(2) + u(3)*taud(3))
      dphiduy_d= -(u(1)*taud(4) + u(2)*taud(5) + u(3)*taud(6))
      dphiduz_d= -(u(1)*taud(7) + u(2)*taud(8) + u(3)*taud(9))

      dphidvx_d= -(v(1)*taud(1) + v(2)*taud(2) + v(3)*taud(3))
      dphidvy_d= -(v(1)*taud(4) + v(2)*taud(5) + v(3)*taud(6))
      dphidvz_d= -(v(1)*taud(7) + v(2)*taud(8) + v(3)*taud(9))

      dphidwx_d= -(w(1)*taud(1) + w(2)*taud(2) + w(3)*taud(3))
      dphidwy_d= -(w(1)*taud(4) + w(2)*taud(5) + w(3)*taud(6))
      dphidwz_d= -(w(1)*taud(7) + w(2)*taud(8) + w(3)*taud(9))

      dphidux_p= -(u(1)*taup(1) + u(2)*taup(2) + u(3)*taup(3))
      dphiduy_p= -(u(1)*taup(4) + u(2)*taup(5) + u(3)*taup(6))
      dphiduz_p= -(u(1)*taup(7) + u(2)*taup(8) + u(3)*taup(9))

      dphidvx_p= -(v(1)*taup(1) + v(2)*taup(2) + v(3)*taup(3))
      dphidvy_p= -(v(1)*taup(4) + v(2)*taup(5) + v(3)*taup(6))
      dphidvz_p= -(v(1)*taup(7) + v(2)*taup(8) + v(3)*taup(9))

      dphidwx_p= -(w(1)*taup(1) + w(2)*taup(2) + w(3)*taup(3))
      dphidwy_p= -(w(1)*taup(4) + w(2)*taup(5) + w(3)*taup(6))
      dphidwz_p= -(w(1)*taup(7) + w(2)*taup(8) + w(3)*taup(9))


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

         duxx_d = uv(1)*dphidvx_d/(usiz*uvsin) + uw(1)*dphidwx_d/usiz
         duxy_d = uv(1)*dphidvy_d/(usiz*uvsin) + uw(1)*dphidwy_d/usiz
         duxz_d = uv(1)*dphidvz_d/(usiz*uvsin) + uw(1)*dphidwz_d/usiz

         duyx_d = uv(2)*dphidvx_d/(usiz*uvsin) + uw(2)*dphidwx_d/usiz
         duyy_d = uv(2)*dphidvy_d/(usiz*uvsin) + uw(2)*dphidwy_d/usiz
         duyz_d = uv(2)*dphidvz_d/(usiz*uvsin) + uw(2)*dphidwz_d/usiz

         duzx_d = uv(3)*dphidvx_d/(usiz*uvsin) + uw(3)*dphidwx_d/usiz
         duzy_d = uv(3)*dphidvy_d/(usiz*uvsin) + uw(3)*dphidwy_d/usiz
         duzz_d = uv(3)*dphidvz_d/(usiz*uvsin) + uw(3)*dphidwz_d/usiz

         duxx_p = uv(1)*dphidvx_p/(usiz*uvsin) + uw(1)*dphidwx_p/usiz
         duxy_p = uv(1)*dphidvy_p/(usiz*uvsin) + uw(1)*dphidwy_p/usiz
         duxz_p = uv(1)*dphidvz_p/(usiz*uvsin) + uw(1)*dphidwz_p/usiz

         duyx_p = uv(2)*dphidvx_p/(usiz*uvsin) + uw(2)*dphidwx_p/usiz
         duyy_p = uv(2)*dphidvy_p/(usiz*uvsin) + uw(2)*dphidwy_p/usiz
         duyz_p = uv(2)*dphidvz_p/(usiz*uvsin) + uw(2)*dphidwz_p/usiz

         duzx_p = uv(3)*dphidvx_p/(usiz*uvsin) + uw(3)*dphidwx_p/usiz
         duzy_p = uv(3)*dphidvy_p/(usiz*uvsin) + uw(3)*dphidwy_p/usiz
         duzz_p = uv(3)*dphidvz_p/(usiz*uvsin) + uw(3)*dphidwz_p/usiz

         dvxx_d = -uv(1)*dphidux_d/(vsiz*uvsin)
         dvxy_d = -uv(1)*dphiduy_d/(vsiz*uvsin)
         dvxz_d = -uv(1)*dphiduz_d/(vsiz*uvsin)

         dvyx_d = -uv(2)*dphidux_d/(vsiz*uvsin) 
         dvyy_d = -uv(2)*dphiduy_d/(vsiz*uvsin)
         dvyz_d = -uv(2)*dphiduz_d/(vsiz*uvsin)

         dvzx_d = -uv(3)*dphidux_d/(vsiz*uvsin)
         dvzy_d = -uv(3)*dphiduy_d/(vsiz*uvsin)
         dvzz_d = -uv(3)*dphiduz_d/(vsiz*uvsin)

         dvxx_p = -uv(1)*dphidux_p/(vsiz*uvsin) 
         dvxy_p = -uv(1)*dphiduy_p/(vsiz*uvsin)
         dvxz_p = -uv(1)*dphiduz_p/(vsiz*uvsin)

         dvyx_p = -uv(2)*dphidux_p/(vsiz*uvsin)
         dvyy_p = -uv(2)*dphiduy_p/(vsiz*uvsin)
         dvyz_p = -uv(2)*dphiduz_p/(vsiz*uvsin)

         dvzx_p = -uv(3)*dphidux_p/(vsiz*uvsin)
         dvzy_p = -uv(3)*dphiduy_p/(vsiz*uvsin)
         dvzz_p = -uv(3)*dphiduz_p/(vsiz*uvsin)

           done_ia=.false.
           done_ic=.false.
           done_ib=.false.
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv
      
           do m=1,nlocal
              if((ilocal(1,m).eq.i).and.(ilocal(2,m).eq.ia))then
                done_ia=.true.
                dlocal2d(1,m)=dlocal2d(1,m)+duxx_d
                dlocal2d(2,m)=dlocal2d(2,m)+duxy_d
                dlocal2d(3,m)=dlocal2d(3,m)+duxz_d
                dlocal2d(4,m)=dlocal2d(4,m)+duyx_d
                dlocal2d(5,m)=dlocal2d(5,m)+duyy_d
                dlocal2d(6,m)=dlocal2d(6,m)+duyz_d
                dlocal2d(7,m)=dlocal2d(7,m)+duzx_d
                dlocal2d(8,m)=dlocal2d(8,m)+duzy_d
                dlocal2d(9,m)=dlocal2d(9,m)+duzz_d
              
                dlocal2p(1,m)=dlocal2p(1,m)+duxx_p
                dlocal2p(2,m)=dlocal2p(2,m)+duxy_p
                dlocal2p(3,m)=dlocal2p(3,m)+duxz_p
                dlocal2p(4,m)=dlocal2p(4,m)+duyx_p
                dlocal2p(5,m)=dlocal2p(5,m)+duyy_p
                dlocal2p(6,m)=dlocal2p(6,m)+duyz_p
                dlocal2p(7,m)=dlocal2p(7,m)+duzx_p
                dlocal2p(8,m)=dlocal2p(8,m)+duzy_p
                dlocal2p(9,m)=dlocal2p(9,m)+duzz_p
             else if((ilocal(1,m).eq.ia).and.(ilocal(2,m).eq.i))then
                done_ia=.true.
                dlocal1d(1,m)=dlocal1d(1,m)+duxx_d
                dlocal1d(2,m)=dlocal1d(2,m)+duxy_d
                dlocal1d(3,m)=dlocal1d(3,m)+duxz_d
                dlocal1d(4,m)=dlocal1d(4,m)+duyx_d
                dlocal1d(5,m)=dlocal1d(5,m)+duyy_d
                dlocal1d(6,m)=dlocal1d(6,m)+duyz_d
                dlocal1d(7,m)=dlocal1d(7,m)+duzx_d
                dlocal1d(8,m)=dlocal1d(8,m)+duzy_d
                dlocal1d(9,m)=dlocal1d(9,m)+duzz_d

                dlocal1p(1,m)=dlocal1p(1,m)+duxx_p
                dlocal1p(2,m)=dlocal1p(2,m)+duxy_p
                dlocal1p(3,m)=dlocal1p(3,m)+duxz_p
                dlocal1p(4,m)=dlocal1p(4,m)+duyx_p
                dlocal1p(5,m)=dlocal1p(5,m)+duyy_p
                dlocal1p(6,m)=dlocal1p(6,m)+duyz_p
                dlocal1p(7,m)=dlocal1p(7,m)+duzx_p
                dlocal1p(8,m)=dlocal1p(8,m)+duzy_p
                dlocal1p(9,m)=dlocal1p(9,m)+duzz_p
              else if((ilocal(1,m).eq.i).and.(ilocal(2,m).eq.ic)) then
                done_ic=.true.
                dlocal2d(1,m)=dlocal2d(1,m)+dvxx_d
                dlocal2d(2,m)=dlocal2d(2,m)+dvxy_d
                dlocal2d(3,m)=dlocal2d(3,m)+dvxz_d
                dlocal2d(4,m)=dlocal2d(4,m)+dvyx_d
                dlocal2d(5,m)=dlocal2d(5,m)+dvyy_d
                dlocal2d(6,m)=dlocal2d(6,m)+dvyz_d
                dlocal2d(7,m)=dlocal2d(7,m)+dvzx_d
                dlocal2d(8,m)=dlocal2d(8,m)+dvzy_d
                dlocal2d(9,m)=dlocal2d(9,m)+dvzz_d

                dlocal2p(1,m)=dlocal2p(1,m)+dvxx_p
                dlocal2p(2,m)=dlocal2p(2,m)+dvxy_p
                dlocal2p(3,m)=dlocal2p(3,m)+dvxz_p
                dlocal2p(4,m)=dlocal2p(4,m)+dvyx_p
                dlocal2p(5,m)=dlocal2p(5,m)+dvyy_p
                dlocal2p(6,m)=dlocal2p(6,m)+dvyz_p
                dlocal2p(7,m)=dlocal2p(7,m)+dvzx_p
                dlocal2p(8,m)=dlocal2p(8,m)+dvzy_p
                dlocal2p(9,m)=dlocal2p(9,m)+dvzz_p
              else if((ilocal(1,m).eq.ic).and.(ilocal(2,m).eq.i)) then
                done_ic=.true.
                dlocal1d(1,m)=dlocal1d(1,m)+dvxx_d
                dlocal1d(2,m)=dlocal1d(2,m)+dvxy_d
                dlocal1d(3,m)=dlocal1d(3,m)+dvxz_d
                dlocal1d(4,m)=dlocal1d(4,m)+dvyx_d
                dlocal1d(5,m)=dlocal1d(5,m)+dvyy_d
                dlocal1d(6,m)=dlocal1d(6,m)+dvyz_d
                dlocal1d(7,m)=dlocal1d(7,m)+dvzx_d
                dlocal1d(8,m)=dlocal1d(8,m)+dvzy_d
                dlocal1d(9,m)=dlocal1d(9,m)+dvzz_d

                dlocal1p(1,m)=dlocal1p(1,m)+dvxx_p
                dlocal1p(2,m)=dlocal1p(2,m)+dvxy_p
                dlocal1p(3,m)=dlocal1p(3,m)+dvxz_p
                dlocal1p(4,m)=dlocal1p(4,m)+dvyx_p
                dlocal1p(5,m)=dlocal1p(5,m)+dvyy_p
                dlocal1p(6,m)=dlocal1p(6,m)+dvyz_p
                dlocal1p(7,m)=dlocal1p(7,m)+dvzx_p
                dlocal1p(8,m)=dlocal1p(8,m)+dvzy_p
                dlocal1p(9,m)=dlocal1p(9,m)+dvzz_p

              else if((ilocal(1,m).eq.i).and.(ilocal(2,m).eq.ib)) then
                done_ib=.true.
                dlocal2d(1,m)=dlocal2d(1,m)-duxx_d-dvxx_d
                dlocal2d(2,m)=dlocal2d(2,m)-duxy_d-dvxy_d
                dlocal2d(3,m)=dlocal2d(3,m)-duxz_d-dvxz_d
                dlocal2d(4,m)=dlocal2d(4,m)-duyx_d-dvyx_d
                dlocal2d(5,m)=dlocal2d(5,m)-duyy_d-dvyy_d
                dlocal2d(6,m)=dlocal2d(6,m)-duyz_d-dvyz_d
                dlocal2d(7,m)=dlocal2d(7,m)-duzx_d-dvzx_d
                dlocal2d(8,m)=dlocal2d(8,m)-duzy_d-dvzy_d
                dlocal2d(9,m)=dlocal2d(9,m)-duzz_d-dvzz_d

                dlocal2p(1,m)=dlocal2p(1,m)-duxx_p-dvxx_p
                dlocal2p(2,m)=dlocal2p(2,m)-duxy_p-dvxy_p
                dlocal2p(3,m)=dlocal2p(3,m)-duxz_p-dvxz_p
                dlocal2p(4,m)=dlocal2p(4,m)-duyx_p-dvyx_p
                dlocal2p(5,m)=dlocal2p(5,m)-duyy_p-dvyy_p
                dlocal2p(6,m)=dlocal2p(6,m)-duyz_p-dvyz_p
                dlocal2p(7,m)=dlocal2p(7,m)-duzx_p-dvzx_p
                dlocal2p(8,m)=dlocal2p(8,m)-duzy_p-dvzy_p
                dlocal2p(9,m)=dlocal2p(9,m)-duzz_p-dvzz_p
              else if((ilocal(1,m).eq.ib).and.(ilocal(2,m).eq.i)) then
                done_ib=.true.
                dlocal1d(1,m)=dlocal1d(1,m)-duxx_d-dvxx_d
                dlocal1d(2,m)=dlocal1d(2,m)-duxy_d-dvxy_d
                dlocal1d(3,m)=dlocal1d(3,m)-duxz_d-dvxz_d
                dlocal1d(4,m)=dlocal1d(4,m)-duyx_d-dvyx_d
                dlocal1d(5,m)=dlocal1d(5,m)-duyy_d-dvyy_d
                dlocal1d(6,m)=dlocal1d(6,m)-duyz_d-dvyz_d
                dlocal1d(7,m)=dlocal1d(7,m)-duzx_d-dvzx_d
                dlocal1d(8,m)=dlocal1d(8,m)-duzy_d-dvzy_d
                dlocal1d(9,m)=dlocal1d(9,m)-duzz_d-dvzz_d

                dlocal1p(1,m)=dlocal1p(1,m)-duxx_p-dvxx_p
                dlocal1p(2,m)=dlocal1p(2,m)-duxy_p-dvxy_p
                dlocal1p(3,m)=dlocal1p(3,m)-duxz_p-dvxz_p
                dlocal1p(4,m)=dlocal1p(4,m)-duyx_p-dvyx_p
                dlocal1p(5,m)=dlocal1p(5,m)-duyy_p-dvyy_p
                dlocal1p(6,m)=dlocal1p(6,m)-duyz_p-dvyz_p
                dlocal1p(7,m)=dlocal1p(7,m)-duzx_p-dvzx_p
                dlocal1p(8,m)=dlocal1p(8,m)-duzy_p-dvzy_p
                dlocal1p(9,m)=dlocal1p(9,m)-duzz_p-dvzz_p

              end if
              if( (done_ia.eq..true.) .and. (done_ib.eq..true.)
     &           .and. (done_ic.eq..true.) ) then
                goto 31
              end if
           end do  
   31  continue

              if(done_ia.eq..false.) then
                if(ia.le.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=ia
                ilocal(2,nlocal)=i
                dlocal1d(1,nlocal)=dlocal1d(1,nlocal)+duxx_d
                dlocal1d(2,nlocal)=dlocal1d(2,nlocal)+duxy_d
                dlocal1d(3,nlocal)=dlocal1d(3,nlocal)+duxz_d
                dlocal1d(4,nlocal)=dlocal1d(4,nlocal)+duyx_d
                dlocal1d(5,nlocal)=dlocal1d(5,nlocal)+duyy_d
                dlocal1d(6,nlocal)=dlocal1d(6,nlocal)+duyz_d
                dlocal1d(7,nlocal)=dlocal1d(7,nlocal)+duzx_d
                dlocal1d(8,nlocal)=dlocal1d(8,nlocal)+duzy_d
                dlocal1d(9,nlocal)=dlocal1d(9,nlocal)+duzz_d

                dlocal1p(1,nlocal)=dlocal1p(1,nlocal)+duxx_p
                dlocal1p(2,nlocal)=dlocal1p(2,nlocal)+duxy_p
                dlocal1p(3,nlocal)=dlocal1p(3,nlocal)+duxz_p
                dlocal1p(4,nlocal)=dlocal1p(4,nlocal)+duyx_p
                dlocal1p(5,nlocal)=dlocal1p(5,nlocal)+duyy_p
                dlocal1p(6,nlocal)=dlocal1p(6,nlocal)+duyz_p
                dlocal1p(7,nlocal)=dlocal1p(7,nlocal)+duzx_p
                dlocal1p(8,nlocal)=dlocal1p(8,nlocal)+duzy_p
                dlocal1p(9,nlocal)=dlocal1p(9,nlocal)+duzz_p
                   if(ia.eq.i) then
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+duxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+duxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+duxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+duyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+duyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+duyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+duzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+duzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+duzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+duxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+duxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+duxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+duyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+duyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+duyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+duzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+duzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+duzz_p
                   end if
                else if(ia.gt.i) then 
                nlocal=nlocal+1
                ilocal(1,nlocal)=i
                ilocal(2,nlocal)=ia
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+duxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+duxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+duxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+duyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+duyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+duyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+duzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+duzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+duzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+duxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+duxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+duxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+duyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+duyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+duyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+duzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+duzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+duzz_p
                end if
              else if(done_ic.eq..false.) then
                if(ic.le.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=ic
                ilocal(2,nlocal)=i
                dlocal1d(1,nlocal)=dlocal1d(1,nlocal)+dvxx_d
                dlocal1d(2,nlocal)=dlocal1d(2,nlocal)+dvxy_d
                dlocal1d(3,nlocal)=dlocal1d(3,nlocal)+dvxz_d
                dlocal1d(4,nlocal)=dlocal1d(4,nlocal)+dvyx_d
                dlocal1d(5,nlocal)=dlocal1d(5,nlocal)+dvyy_d
                dlocal1d(6,nlocal)=dlocal1d(6,nlocal)+dvyz_d
                dlocal1d(7,nlocal)=dlocal1d(7,nlocal)+dvzx_d
                dlocal1d(8,nlocal)=dlocal1d(8,nlocal)+dvzy_d
                dlocal1d(9,nlocal)=dlocal1d(9,nlocal)+dvzz_d

                dlocal1p(1,nlocal)=dlocal1p(1,nlocal)+dvxx_p
                dlocal1p(2,nlocal)=dlocal1p(2,nlocal)+dvxy_p
                dlocal1p(3,nlocal)=dlocal1p(3,nlocal)+dvxz_p
                dlocal1p(4,nlocal)=dlocal1p(4,nlocal)+dvyx_p
                dlocal1p(5,nlocal)=dlocal1p(5,nlocal)+dvyy_p
                dlocal1p(6,nlocal)=dlocal1p(6,nlocal)+dvyz_p
                dlocal1p(7,nlocal)=dlocal1p(7,nlocal)+dvzx_p
                dlocal1p(8,nlocal)=dlocal1p(8,nlocal)+dvzy_p
                dlocal1p(9,nlocal)=dlocal1p(9,nlocal)+dvzz_p
                   if(ic.eq.i) then
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+dvzz_p
                   end if
                else if(ic.gt.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=i
                ilocal(2,nlocal)=ic
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+dvzz_p
                end if
              else if(done_ib.eq..false.) then
                if(ib.le.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=ib
                ilocal(2,nlocal)=i
                done_ib=.true.
                dlocal1d(1,nlocal)=dlocal1d(1,nlocal)-duxx_d-dvxx_d
                dlocal1d(2,nlocal)=dlocal1d(2,nlocal)-duxy_d-dvxy_d
                dlocal1d(3,nlocal)=dlocal1d(3,nlocal)-duxz_d-dvxz_d
                dlocal1d(4,nlocal)=dlocal1d(4,nlocal)-duyx_d-dvyx_d
                dlocal1d(5,nlocal)=dlocal1d(5,nlocal)-duyy_d-dvyy_d
                dlocal1d(6,nlocal)=dlocal1d(6,nlocal)-duyz_d-dvyz_d
                dlocal1d(7,nlocal)=dlocal1d(7,nlocal)-duzx_d-dvzx_d
                dlocal1d(8,nlocal)=dlocal1d(8,nlocal)-duzy_d-dvzy_d
                dlocal1d(9,nlocal)=dlocal1d(9,nlocal)-duzz_d-dvzz_d

                dlocal1p(1,nlocal)=dlocal1p(1,nlocal)-duxx_p-dvxx_p
                dlocal1p(2,nlocal)=dlocal1p(2,nlocal)-duxy_p-dvxy_p
                dlocal1p(3,nlocal)=dlocal1p(3,nlocal)-duxz_p-dvxz_p
                dlocal1p(4,nlocal)=dlocal1p(4,nlocal)-duyx_p-dvyx_p
                dlocal1p(5,nlocal)=dlocal1p(5,nlocal)-duyy_p-dvyy_p
                dlocal1p(6,nlocal)=dlocal1p(6,nlocal)-duyz_p-dvyz_p
                dlocal1p(7,nlocal)=dlocal1p(7,nlocal)-duzx_p-dvzx_p
                dlocal1p(8,nlocal)=dlocal1p(8,nlocal)-duzy_p-dvzy_p
                dlocal1p(9,nlocal)=dlocal1p(9,nlocal)-duzz_p-dvzz_p
                   if(ib.eq.i) then
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)-duxx_d-dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)-duxy_d-dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)-duxz_d-dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)-duyx_d-dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)-duyy_d-dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)-duyz_d-dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)-duzx_d-dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)-duzy_d-dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)-duzz_d-dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)-duxx_p-dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)-duxy_p-dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)-duxz_p-dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)-duyx_p-dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)-duyy_p-dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)-duyz_p-dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)-duzx_p-dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)-duzy_p-dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)-duzz_p-dvzz_p
                   end if
                else if(ib.gt.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=i
                ilocal(2,nlocal)=ib
                done_ib=.true.
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)-duxx_d-dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)-duxy_d-dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)-duxz_d-dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)-duyx_d-dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)-duyy_d-dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)-duyz_d-dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)-duzx_d-dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)-duzy_d-dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)-duzz_d-dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)-duxx_p-dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)-duxy_p-dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)-duxz_p-dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)-duyx_p-dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)-duyy_p-dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)-duyz_p-dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)-duzx_p-dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)-duzy_p-dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)-duzz_p-dvzz_p
                end if
              end if
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

        duxx_d=uv(1)*dphidvx_d/(usiz*uvsin) + 0.5d0*uw(1)*dphidwx_d/usiz
        duxy_d=uv(1)*dphidvy_d/(usiz*uvsin) + 0.5d0*uw(1)*dphidwy_d/usiz
        duxz_d=uv(1)*dphidvz_d/(usiz*uvsin) + 0.5d0*uw(1)*dphidwz_d/usiz

        duyx_d=uv(2)*dphidvx_d/(usiz*uvsin) + 0.5d0*uw(2)*dphidwx_d/usiz
        duyy_d=uv(2)*dphidvy_d/(usiz*uvsin) + 0.5d0*uw(2)*dphidwy_d/usiz
        duyz_d=uv(2)*dphidvz_d/(usiz*uvsin) + 0.5d0*uw(2)*dphidwz_d/usiz

        duzx_d=uv(3)*dphidvx_d/(usiz*uvsin) + 0.5d0*uw(3)*dphidwx_d/usiz
        duzy_d=uv(3)*dphidvy_d/(usiz*uvsin) + 0.5d0*uw(3)*dphidwy_d/usiz
        duzz_d=uv(3)*dphidvz_d/(usiz*uvsin) + 0.5d0*uw(3)*dphidwz_d/usiz

        duxx_p=uv(1)*dphidvx_p/(usiz*uvsin) + 0.5d0*uw(1)*dphidwx_p/usiz
        duxy_p=uv(1)*dphidvy_p/(usiz*uvsin) + 0.5d0*uw(1)*dphidwy_p/usiz
        duxz_p=uv(1)*dphidvz_p/(usiz*uvsin) + 0.5d0*uw(1)*dphidwz_p/usiz

        duyx_p=uv(2)*dphidvx_p/(usiz*uvsin) + 0.5d0*uw(2)*dphidwx_p/usiz
        duyy_p=uv(2)*dphidvy_p/(usiz*uvsin) + 0.5d0*uw(2)*dphidwy_p/usiz
        duyz_p=uv(2)*dphidvz_p/(usiz*uvsin) + 0.5d0*uw(2)*dphidwz_p/usiz

        duzx_p=uv(3)*dphidvx_p/(usiz*uvsin) + 0.5d0*uw(3)*dphidwx_p/usiz
        duzy_p=uv(3)*dphidvy_p/(usiz*uvsin) + 0.5d0*uw(3)*dphidwy_p/usiz
        duzz_p=uv(3)*dphidvz_p/(usiz*uvsin) + 0.5d0*uw(3)*dphidwz_p/usiz

        dvxx_d=-uv(1)*dphidux_d/(vsiz*uvsin)+ 0.5d0*vw(1)*dphidwx_d/vsiz
        dvxy_d=-uv(1)*dphiduy_d/(vsiz*uvsin)+ 0.5d0*vw(1)*dphidwy_d/vsiz
        dvxz_d=-uv(1)*dphiduz_d/(vsiz*uvsin)+ 0.5d0*vw(1)*dphidwz_d/vsiz

        dvyx_d=-uv(2)*dphidux_d/(vsiz*uvsin)+ 0.5d0*vw(2)*dphidwx_d/vsiz
        dvyy_d=-uv(2)*dphiduy_d/(vsiz*uvsin)+ 0.5d0*vw(2)*dphidwy_d/vsiz
        dvyz_d=-uv(2)*dphiduz_d/(vsiz*uvsin)+ 0.5d0*vw(2)*dphidwz_d/vsiz

        dvzx_d=-uv(3)*dphidux_d/(vsiz*uvsin)+ 0.5d0*vw(3)*dphidwx_d/vsiz
        dvzy_d=-uv(3)*dphiduy_d/(vsiz*uvsin)+ 0.5d0*vw(3)*dphidwy_d/vsiz
        dvzz_d=-uv(3)*dphiduz_d/(vsiz*uvsin)+ 0.5d0*vw(3)*dphidwz_d/vsiz

        dvxx_p=-uv(1)*dphidux_p/(vsiz*uvsin)+ 0.5d0*vw(1)*dphidwx_p/vsiz
        dvxy_p=-uv(1)*dphiduy_p/(vsiz*uvsin)+ 0.5d0*vw(1)*dphidwy_p/vsiz
        dvxz_p=-uv(1)*dphiduz_p/(vsiz*uvsin)+ 0.5d0*vw(1)*dphidwz_p/vsiz

        dvyx_p=-uv(2)*dphidux_p/(vsiz*uvsin)+ 0.5d0*vw(2)*dphidwx_p/vsiz
        dvyy_p=-uv(2)*dphiduy_p/(vsiz*uvsin)+ 0.5d0*vw(2)*dphidwy_p/vsiz
        dvyz_p=-uv(2)*dphiduz_p/(vsiz*uvsin)+ 0.5d0*vw(2)*dphidwz_p/vsiz

        dvzx_p=-uv(3)*dphidux_p/(vsiz*uvsin)+ 0.5d0*vw(3)*dphidwx_p/vsiz
        dvzy_p=-uv(3)*dphiduy_p/(vsiz*uvsin)+ 0.5d0*vw(3)*dphidwy_p/vsiz
        dvzz_p=-uv(3)*dphiduz_p/(vsiz*uvsin)+ 0.5d0*vw(3)*dphidwz_p/vsiz

           done_ia=.false.
           done_ic=.false.
           done_ib=.false.
c            depi(j,ia) = depi(j,ia) + du
c            depi(j,ic) = depi(j,ic) + dv
c            depi(j,ib) = depi(j,ib) - du - dv

           do m=1,nlocal
              if( (ilocal(1,m).eq.i) .and. (ilocal(2,m).eq.ia) ) then
                done_ia=.true.
                dlocal2d(1,m)=dlocal2d(1,m)+duxx_d
                dlocal2d(2,m)=dlocal2d(2,m)+duxy_d
                dlocal2d(3,m)=dlocal2d(3,m)+duxz_d
                dlocal2d(4,m)=dlocal2d(4,m)+duyx_d
                dlocal2d(5,m)=dlocal2d(5,m)+duyy_d
                dlocal2d(6,m)=dlocal2d(6,m)+duyz_d
                dlocal2d(7,m)=dlocal2d(7,m)+duzx_d
                dlocal2d(8,m)=dlocal2d(8,m)+duzy_d
                dlocal2d(9,m)=dlocal2d(9,m)+duzz_d

                dlocal2p(1,m)=dlocal2p(1,m)+duxx_p
                dlocal2p(2,m)=dlocal2p(2,m)+duxy_p
                dlocal2p(3,m)=dlocal2p(3,m)+duxz_p
                dlocal2p(4,m)=dlocal2p(4,m)+duyx_p
                dlocal2p(5,m)=dlocal2p(5,m)+duyy_p
                dlocal2p(6,m)=dlocal2p(6,m)+duyz_p
                dlocal2p(7,m)=dlocal2p(7,m)+duzx_p
                dlocal2p(8,m)=dlocal2p(8,m)+duzy_p
                dlocal2p(9,m)=dlocal2p(9,m)+duzz_p
              else if((ilocal(1,m).eq.ia).and.(ilocal(2,m).eq.i)) then
                done_ia=.true.
                dlocal1d(1,m)=dlocal1d(1,m)+duxx_d
                dlocal1d(2,m)=dlocal1d(2,m)+duxy_d
                dlocal1d(3,m)=dlocal1d(3,m)+duxz_d
                dlocal1d(4,m)=dlocal1d(4,m)+duyx_d
                dlocal1d(5,m)=dlocal1d(5,m)+duyy_d
                dlocal1d(6,m)=dlocal1d(6,m)+duyz_d
                dlocal1d(7,m)=dlocal1d(7,m)+duzx_d
                dlocal1d(8,m)=dlocal1d(8,m)+duzy_d
                dlocal1d(9,m)=dlocal1d(9,m)+duzz_d

                dlocal1p(1,m)=dlocal1p(1,m)+duxx_p
                dlocal1p(2,m)=dlocal1p(2,m)+duxy_p
                dlocal1p(3,m)=dlocal1p(3,m)+duxz_p
                dlocal1p(4,m)=dlocal1p(4,m)+duyx_p
                dlocal1p(5,m)=dlocal1p(5,m)+duyy_p
                dlocal1p(6,m)=dlocal1p(6,m)+duyz_p
                dlocal1p(7,m)=dlocal1p(7,m)+duzx_p
                dlocal1p(8,m)=dlocal1p(8,m)+duzy_p
                dlocal1p(9,m)=dlocal1p(9,m)+duzz_p

              else if((ilocal(1,m).eq.i).and.(ilocal(2,m).eq.ic)) then
                done_ic=.true.
                dlocal2d(1,m)=dlocal2d(1,m)+dvxx_d
                dlocal2d(2,m)=dlocal2d(2,m)+dvxy_d
                dlocal2d(3,m)=dlocal2d(3,m)+dvxz_d
                dlocal2d(4,m)=dlocal2d(4,m)+dvyx_d
                dlocal2d(5,m)=dlocal2d(5,m)+dvyy_d
                dlocal2d(6,m)=dlocal2d(6,m)+dvyz_d
                dlocal2d(7,m)=dlocal2d(7,m)+dvzx_d
                dlocal2d(8,m)=dlocal2d(8,m)+dvzy_d
                dlocal2d(9,m)=dlocal2d(9,m)+dvzz_d

                dlocal2p(1,m)=dlocal2p(1,m)+dvxx_p
                dlocal2p(2,m)=dlocal2p(2,m)+dvxy_p
                dlocal2p(3,m)=dlocal2p(3,m)+dvxz_p
                dlocal2p(4,m)=dlocal2p(4,m)+dvyx_p
                dlocal2p(5,m)=dlocal2p(5,m)+dvyy_p
                dlocal2p(6,m)=dlocal2p(6,m)+dvyz_p
                dlocal2p(7,m)=dlocal2p(7,m)+dvzx_p
                dlocal2p(8,m)=dlocal2p(8,m)+dvzy_p
                dlocal2p(9,m)=dlocal2p(9,m)+dvzz_p
              else if((ilocal(1,m).eq.ic).and.(ilocal(2,m).eq.i)) then
                done_ic=.true.
                dlocal1d(1,m)=dlocal1d(1,m)+dvxx_d
                dlocal1d(2,m)=dlocal1d(2,m)+dvxy_d
                dlocal1d(3,m)=dlocal1d(3,m)+dvxz_d
                dlocal1d(4,m)=dlocal1d(4,m)+dvyx_d
                dlocal1d(5,m)=dlocal1d(5,m)+dvyy_d
                dlocal1d(6,m)=dlocal1d(6,m)+dvyz_d
                dlocal1d(7,m)=dlocal1d(7,m)+dvzx_d
                dlocal1d(8,m)=dlocal1d(8,m)+dvzy_d
                dlocal1d(9,m)=dlocal1d(9,m)+dvzz_d

                dlocal1p(1,m)=dlocal1p(1,m)+dvxx_p
                dlocal1p(2,m)=dlocal1p(2,m)+dvxy_p
                dlocal1p(3,m)=dlocal1p(3,m)+dvxz_p
                dlocal1p(4,m)=dlocal1p(4,m)+dvyx_p
                dlocal1p(5,m)=dlocal1p(5,m)+dvyy_p
                dlocal1p(6,m)=dlocal1p(6,m)+dvyz_p
                dlocal1p(7,m)=dlocal1p(7,m)+dvzx_p
                dlocal1p(8,m)=dlocal1p(8,m)+dvzy_p
                dlocal1p(9,m)=dlocal1p(9,m)+dvzz_p

              else if((ilocal(1,m).eq.i).and.(ilocal(2,m).eq.ib)) then
                done_ib=.true.
                dlocal2d(1,m)=dlocal2d(1,m)-duxx_d-dvxx_d
                dlocal2d(2,m)=dlocal2d(2,m)-duxy_d-dvxy_d
                dlocal2d(3,m)=dlocal2d(3,m)-duxz_d-dvxz_d
                dlocal2d(4,m)=dlocal2d(4,m)-duyx_d-dvyx_d
                dlocal2d(5,m)=dlocal2d(5,m)-duyy_d-dvyy_d
                dlocal2d(6,m)=dlocal2d(6,m)-duyz_d-dvyz_d
                dlocal2d(7,m)=dlocal2d(7,m)-duzx_d-dvzx_d
                dlocal2d(8,m)=dlocal2d(8,m)-duzy_d-dvzy_d
                dlocal2d(9,m)=dlocal2d(9,m)-duzz_d-dvzz_d

                dlocal2p(1,m)=dlocal2p(1,m)-duxx_p-dvxx_p
                dlocal2p(2,m)=dlocal2p(2,m)-duxy_p-dvxy_p
                dlocal2p(3,m)=dlocal2p(3,m)-duxz_p-dvxz_p
                dlocal2p(4,m)=dlocal2p(4,m)-duyx_p-dvyx_p
                dlocal2p(5,m)=dlocal2p(5,m)-duyy_p-dvyy_p
                dlocal2p(6,m)=dlocal2p(6,m)-duyz_p-dvyz_p
                dlocal2p(7,m)=dlocal2p(7,m)-duzx_p-dvzx_p
                dlocal2p(8,m)=dlocal2p(8,m)-duzy_p-dvzy_p
                dlocal2p(9,m)=dlocal2p(9,m)-duzz_p-dvzz_p
              else if((ilocal(1,m).eq.ib).and.(ilocal(2,m).eq.i)) then
                done_ib=.true.
                dlocal1d(1,m)=dlocal1d(1,m)-duxx_d-dvxx_d
                dlocal1d(2,m)=dlocal1d(2,m)-duxy_d-dvxy_d
                dlocal1d(3,m)=dlocal1d(3,m)-duxz_d-dvxz_d
                dlocal1d(4,m)=dlocal1d(4,m)-duyx_d-dvyx_d
                dlocal1d(5,m)=dlocal1d(5,m)-duyy_d-dvyy_d
                dlocal1d(6,m)=dlocal1d(6,m)-duyz_d-dvyz_d
                dlocal1d(7,m)=dlocal1d(7,m)-duzx_d-dvzx_d
                dlocal1d(8,m)=dlocal1d(8,m)-duzy_d-dvzy_d
                dlocal1d(9,m)=dlocal1d(9,m)-duzz_d-dvzz_d

                dlocal1p(1,m)=dlocal1p(1,m)-duxx_p-dvxx_p
                dlocal1p(2,m)=dlocal1p(2,m)-duxy_p-dvxy_p
                dlocal1p(3,m)=dlocal1p(3,m)-duxz_p-dvxz_p
                dlocal1p(4,m)=dlocal1p(4,m)-duyx_p-dvyx_p
                dlocal1p(5,m)=dlocal1p(5,m)-duyy_p-dvyy_p
                dlocal1p(6,m)=dlocal1p(6,m)-duyz_p-dvyz_p
                dlocal1p(7,m)=dlocal1p(7,m)-duzx_p-dvzx_p
                dlocal1p(8,m)=dlocal1p(8,m)-duzy_p-dvzy_p
                dlocal1p(9,m)=dlocal1p(9,m)-duzz_p-dvzz_p

              end if
              if( (done_ia.eq..true.) .and. (done_ib.eq..true.)
     &           .and. (done_ic.eq..true.) ) then
                goto 33
              end if
           end do
   33  continue
              if(done_ia.eq..false.) then
                if(ia.le.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=ia
                ilocal(2,nlocal)=i
                dlocal1d(1,nlocal)=dlocal1d(1,nlocal)+duxx_d
                dlocal1d(2,nlocal)=dlocal1d(2,nlocal)+duxy_d
                dlocal1d(3,nlocal)=dlocal1d(3,nlocal)+duxz_d
                dlocal1d(4,nlocal)=dlocal1d(4,nlocal)+duyx_d
                dlocal1d(5,nlocal)=dlocal1d(5,nlocal)+duyy_d
                dlocal1d(6,nlocal)=dlocal1d(6,nlocal)+duyz_d
                dlocal1d(7,nlocal)=dlocal1d(7,nlocal)+duzx_d
                dlocal1d(8,nlocal)=dlocal1d(8,nlocal)+duzy_d
                dlocal1d(9,nlocal)=dlocal1d(9,nlocal)+duzz_d

                dlocal1p(1,nlocal)=dlocal1p(1,nlocal)+duxx_p
                dlocal1p(2,nlocal)=dlocal1p(2,nlocal)+duxy_p
                dlocal1p(3,nlocal)=dlocal1p(3,nlocal)+duxz_p
                dlocal1p(4,nlocal)=dlocal1p(4,nlocal)+duyx_p
                dlocal1p(5,nlocal)=dlocal1p(5,nlocal)+duyy_p
                dlocal1p(6,nlocal)=dlocal1p(6,nlocal)+duyz_p
                dlocal1p(7,nlocal)=dlocal1p(7,nlocal)+duzx_p
                dlocal1p(8,nlocal)=dlocal1p(8,nlocal)+duzy_p
                dlocal1p(9,nlocal)=dlocal1p(9,nlocal)+duzz_p
                   if(ia.eq.i) then
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+duxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+duxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+duxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+duyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+duyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+duyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+duzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+duzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+duzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+duxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+duxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+duxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+duyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+duyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+duyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+duzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+duzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+duzz_p

                   end if
                else if(ia.gt.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=i
                ilocal(2,nlocal)=ia
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+duxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+duxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+duxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+duyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+duyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+duyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+duzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+duzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+duzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+duxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+duxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+duxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+duyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+duyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+duyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+duzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+duzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+duzz_p
                end if
              else if(done_ic.eq..false.) then
                if (ic.le.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=ic
                ilocal(2,nlocal)=i
                dlocal1d(1,nlocal)=dlocal1d(1,nlocal)+dvxx_d
                dlocal1d(2,nlocal)=dlocal1d(2,nlocal)+dvxy_d
                dlocal1d(3,nlocal)=dlocal1d(3,nlocal)+dvxz_d
                dlocal1d(4,nlocal)=dlocal1d(4,nlocal)+dvyx_d
                dlocal1d(5,nlocal)=dlocal1d(5,nlocal)+dvyy_d
                dlocal1d(6,nlocal)=dlocal1d(6,nlocal)+dvyz_d
                dlocal1d(7,nlocal)=dlocal1d(7,nlocal)+dvzx_d
                dlocal1d(8,nlocal)=dlocal1d(8,nlocal)+dvzy_d
                dlocal1d(9,nlocal)=dlocal1d(9,nlocal)+dvzz_d

                dlocal1p(1,nlocal)=dlocal1p(1,nlocal)+dvxx_p
                dlocal1p(2,nlocal)=dlocal1p(2,nlocal)+dvxy_p
                dlocal1p(3,nlocal)=dlocal1p(3,nlocal)+dvxz_p
                dlocal1p(4,nlocal)=dlocal1p(4,nlocal)+dvyx_p
                dlocal1p(5,nlocal)=dlocal1p(5,nlocal)+dvyy_p
                dlocal1p(6,nlocal)=dlocal1p(6,nlocal)+dvyz_p
                dlocal1p(7,nlocal)=dlocal1p(7,nlocal)+dvzx_p
                dlocal1p(8,nlocal)=dlocal1p(8,nlocal)+dvzy_p
                dlocal1p(9,nlocal)=dlocal1p(9,nlocal)+dvzz_p
                     if(ic.eq.i) then
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+dvzz_p
                     end if
                else if(ic.gt.i) then
                 nlocal=nlocal+1
                ilocal(1,nlocal)=i
                ilocal(2,nlocal)=ic
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)+dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)+dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)+dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)+dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)+dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)+dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)+dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)+dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)+dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)+dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)+dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)+dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)+dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)+dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)+dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)+dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)+dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)+dvzz_p
                end if
              else if(done_ib.eq..false.) then
                if(ib.le.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=ib
                ilocal(2,nlocal)=i
                done_ib=.true.
                dlocal1d(1,nlocal)=dlocal1d(1,nlocal)-duxx_d-dvxx_d
                dlocal1d(2,nlocal)=dlocal1d(2,nlocal)-duxy_d-dvxy_d
                dlocal1d(3,nlocal)=dlocal1d(3,nlocal)-duxz_d-dvxz_d
                dlocal1d(4,nlocal)=dlocal1d(4,nlocal)-duyx_d-dvyx_d
                dlocal1d(5,nlocal)=dlocal1d(5,nlocal)-duyy_d-dvyy_d
                dlocal1d(6,nlocal)=dlocal1d(6,nlocal)-duyz_d-dvyz_d
                dlocal1d(7,nlocal)=dlocal1d(7,nlocal)-duzx_d-dvzx_d
                dlocal1d(8,nlocal)=dlocal1d(8,nlocal)-duzy_d-dvzy_d
                dlocal1d(9,nlocal)=dlocal1d(9,nlocal)-duzz_d-dvzz_d

                dlocal1p(1,nlocal)=dlocal1p(1,nlocal)-duxx_p-dvxx_p
                dlocal1p(2,nlocal)=dlocal1p(2,nlocal)-duxy_p-dvxy_p
                dlocal1p(3,nlocal)=dlocal1p(3,nlocal)-duxz_p-dvxz_p
                dlocal1p(4,nlocal)=dlocal1p(4,nlocal)-duyx_p-dvyx_p
                dlocal1p(5,nlocal)=dlocal1p(5,nlocal)-duyy_p-dvyy_p
                dlocal1p(6,nlocal)=dlocal1p(6,nlocal)-duyz_p-dvyz_p
                dlocal1p(7,nlocal)=dlocal1p(7,nlocal)-duzx_p-dvzx_p
                dlocal1p(8,nlocal)=dlocal1p(8,nlocal)-duzy_p-dvzy_p
                dlocal1p(9,nlocal)=dlocal1p(9,nlocal)-duzz_p-dvzz_p
                   if(ib.eq.i) then
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)-duxx_d-dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)-duxy_d-dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)-duxz_d-dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)-duyx_d-dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)-duyy_d-dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)-duyz_d-dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)-duzx_d-dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)-duzy_d-dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)-duzz_d-dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)-duxx_p-dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)-duxy_p-dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)-duxz_p-dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)-duyx_p-dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)-duyy_p-dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)-duyz_p-dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)-duzx_p-dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)-duzy_p-dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)-duzz_p-dvzz_p
                   end if
                else if (ib.gt.i) then
                nlocal=nlocal+1
                ilocal(1,nlocal)=i
                ilocal(2,nlocal)=ib
                done_ib=.true.
                dlocal2d(1,nlocal)=dlocal2d(1,nlocal)-duxx_d-dvxx_d
                dlocal2d(2,nlocal)=dlocal2d(2,nlocal)-duxy_d-dvxy_d
                dlocal2d(3,nlocal)=dlocal2d(3,nlocal)-duxz_d-dvxz_d
                dlocal2d(4,nlocal)=dlocal2d(4,nlocal)-duyx_d-dvyx_d
                dlocal2d(5,nlocal)=dlocal2d(5,nlocal)-duyy_d-dvyy_d
                dlocal2d(6,nlocal)=dlocal2d(6,nlocal)-duyz_d-dvyz_d
                dlocal2d(7,nlocal)=dlocal2d(7,nlocal)-duzx_d-dvzx_d
                dlocal2d(8,nlocal)=dlocal2d(8,nlocal)-duzy_d-dvzy_d
                dlocal2d(9,nlocal)=dlocal2d(9,nlocal)-duzz_d-dvzz_d

                dlocal2p(1,nlocal)=dlocal2p(1,nlocal)-duxx_p-dvxx_p
                dlocal2p(2,nlocal)=dlocal2p(2,nlocal)-duxy_p-dvxy_p
                dlocal2p(3,nlocal)=dlocal2p(3,nlocal)-duxz_p-dvxz_p
                dlocal2p(4,nlocal)=dlocal2p(4,nlocal)-duyx_p-dvyx_p
                dlocal2p(5,nlocal)=dlocal2p(5,nlocal)-duyy_p-dvyy_p
                dlocal2p(6,nlocal)=dlocal2p(6,nlocal)-duyz_p-dvyz_p
                dlocal2p(7,nlocal)=dlocal2p(7,nlocal)-duzx_p-dvzx_p
                dlocal2p(8,nlocal)=dlocal2p(8,nlocal)-duzy_p-dvzy_p
                dlocal2p(9,nlocal)=dlocal2p(9,nlocal)-duzz_p-dvzz_p
                end if
              end if


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
