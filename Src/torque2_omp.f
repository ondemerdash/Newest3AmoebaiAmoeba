c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque2  --  convert all site torques to forces  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque2" takes the torque values on all sites defined by
c     local coordinate frames and finds the total Cartesian force
c     components on original sites and sites specifying local frames
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
      subroutine torque2_omp (trq,derivs)
      use sizes
      use atoms
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
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
      real*8 trq(3,*)
      real*8 derivs(3,*)
      character*8 axetyp
      real*8, allocatable :: derivst(:,:)

      allocate(derivst(3,npole))
c
c
c     get the local frame type and the frame-defining atoms
c
      do i=1,npole
         derivst(1,i) = derivs(1,i)
         derivst(2,i) = derivs(2,i)
         derivst(3,i) = derivs(3,i)
      end do
   
      !print*,"OMP torque!! Use the Schwarz, Luke!"

!$OMP PARALLEL default(private) shared(npole,zaxis,ipole,xaxis,yaxis,
!$OMP& polaxe,x,y,z,trq,derivst,derivs)
!$OMP DO reduction(+:derivst)
!$OMP& schedule(guided)
      do i = 1, npole
         ia = zaxis(i)
         ib = ipole(i)
         ic = xaxis(i)
         id = yaxis(i)
         axetyp = polaxe(i)
         if (axetyp .eq. 'None')  goto 10
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
         dphidu = -trq(1,i)*u(1) - trq(2,i)*u(2) - trq(3,i)*u(3)
         dphidv = -trq(1,i)*v(1) - trq(2,i)*v(2) - trq(3,i)*v(3)
         dphidw = -trq(1,i)*w(1) - trq(2,i)*w(2) - trq(3,i)*w(3)
         if (axetyp .eq. 'Z-Bisect') then
            dphidr = -trq(1,i)*r(1) - trq(2,i)*r(2) - trq(3,i)*r(3)
            dphids = -trq(1,i)*s(1) - trq(2,i)*s(2) - trq(3,i)*s(3)
         end if
c
c     force distribution for the Z-Only local coordinate method
c
         if (axetyp .eq. 'Z-Only') then
            do j = 1, 3
               du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
               derivst(j,ia) = derivst(j,ia) + du
               derivst(j,ib) = derivst(j,ib) - du
            end do
c
c     force distribution for the Z-then-X local coordinate method
c
         else if (axetyp .eq. 'Z-then-X') then
            do j = 1, 3
               du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
               dv = -uv(j)*dphidu/(vsiz*uvsin)
               derivst(j,ia) = derivst(j,ia) + du
               derivst(j,ic) = derivst(j,ic) + dv
               derivst(j,ib) = derivst(j,ib) - du - dv
            end do
c
c     force distribution for the Bisector local coordinate method
c
         else if (axetyp .eq. 'Bisector') then
            do j = 1, 3
               du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
               dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
               derivst(j,ia) = derivst(j,ia) + du
               derivst(j,ic) = derivst(j,ic) + dv
               derivst(j,ib) = derivst(j,ib) - du - dv
            end do
c
c     force distribution for the Z-Bisect local coordinate method
c
         else if (axetyp .eq. 'Z-Bisect') then
            do j = 1, 3
               du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
               dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &                 / (vsiz*(ut1sin+ut2sin))
               dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &                 / (wsiz*(ut1sin+ut2sin))
               derivst(j,ia) = derivst(j,ia) + du
               derivst(j,ic) = derivst(j,ic) + dv
               derivst(j,id) = derivst(j,id) + dw
               derivst(j,ib) = derivst(j,ib) - du - dv - dw
            end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
         else if (axetyp .eq. '3-Fold') then
            do j = 1, 3
               du = uw(j)*dphidw/(usiz*uwsin)
     &                 + uv(j)*dphidv/(usiz*uvsin)
     &                 - uw(j)*dphidu/(usiz*uwsin)
     &                 - uv(j)*dphidu/(usiz*uvsin)
               dv = vw(j)*dphidw/(vsiz*vwsin)
     &                 - uv(j)*dphidu/(vsiz*uvsin)
     &                 - vw(j)*dphidv/(vsiz*vwsin)
     &                 + uv(j)*dphidv/(vsiz*uvsin)
               dw = -uw(j)*dphidu/(wsiz*uwsin)
     &                 - vw(j)*dphidv/(wsiz*vwsin)
     &                 + uw(j)*dphidw/(wsiz*uwsin)
     &                 + vw(j)*dphidw/(wsiz*vwsin)
               du = du / 3.0d0
               dv = dv / 3.0d0
               dw = dw / 3.0d0
               derivst(j,ia) = derivst(j,ia) + du
               derivst(j,ic) = derivst(j,ic) + dv
               derivst(j,id) = derivst(j,id) + dw
               derivst(j,ib) = derivst(j,ib) - du - dv - dw
            end do
         end if
   10    continue
      end do
!$OMP END DO

!$OMP DO
      do i=1,npole
         derivs(1,i)=derivst(1,i)
         derivs(2,i)=derivst(2,i)
         derivs(3,i)=derivst(3,i)
      end do
!$OMP END DO

!$OMP END PARALLEL

      deallocate(derivst)
      return
      end

