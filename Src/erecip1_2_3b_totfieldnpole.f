
      subroutine erecip1_2_3b_totfieldnpole(npole3b,pnum,uind,uinp,
     & eptemp,deptemp,virtemp)
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'boxes.i'
c      include 'chgpot.i'
c      include 'ewald.i'
c      include 'ewreg2.i'
c      include 'math.i'
c      include 'mpole.i'
c      include 'polar2.i'
c      include 'polpot.i'
c      include 'units.i'
c      include 'virial.i'
      use sizes
      use atoms
      use chgpot
      use mpole
      use polar, only: polarity, thole, pdamp
      use boxes
      use ewald
      use ewreg3bpolz
      use math 
      implicit none
      integer i,j,k,l,l1
      integer ii,m1,m2
      integer jmin
      integer kmin
      integer lmin
      real*8 e,ei,etot,f,cut
      real*8 expterm,term,eterm
      real*8 uterm,vterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 t1,t2,t3,t4
      real*8 ukp,t3p,t4p
      real*8 de,det1,det2
      real*8 dei,det1i,det2i
      real*8 wterm(3,3)
      real*8 t5(3,3),t6(3,3)
      real*8 t5u(3,3),t5p(3,3)
      real*8 t6u(3,3),t6p(3,3)
      real*8 qt(3,3),dt(3,3)
      real*8 dtu(3,3),dtp(3,3)
c      real*8 ckr(maxatm),skr(maxatm)
c      real*8 cjk(maxatm),sjk(maxatm)
c      real*8 s1(maxatm),s2(maxatm)
c      real*8 s3(maxatm),s4(maxatm)
c      real*8 s3p(maxatm),s4p(maxatm)
c      real*8 cm(maxatm),dm(3,maxatm)
c      real*8 qm(9,maxatm),um(3,maxatm)
c      real*8 trq(3,maxatm),trqi(3,maxatm)
c      real*8 dkx(maxatm),qkx(maxatm)
c      real*8 dky(maxatm),qky(maxatm)
c      real*8 dkz(maxatm),qkz(maxatm)
c      integer npole3b,pnum(*),pnum2(*)
      integer npole3b,pnum(*)
      !real*8 ckr(npole3b),skr(npole3b)
      !real*8 cjk(npole3b),sjk(npole3b)
      !real*8 s1(npole3b),s2(npole3b)
      !real*8 s3(npole3b),s4(npole3b)
      !real*8 s3p(npole3b),s4p(npole3b)
      !real*8 cm(npole3b),dm(3,npole3b)
      !real*8 qm(9,npole3b),um(3,npole3b)
c      real*8 trq(3,npole3b),trqi(3,npole3b)
c      real*8 trq(3,npole),trqi(3,npole)
      !real*8 trqi(3,npole3b)
      !real*8 dkx(npole3b),qkx(npole3b)
      !real*8 dky(npole3b),qky(npole3b)
      !real*8 dkz(npole3b),qkz(npole3b)
      !real*8 ejc(npole3b,0:maxvec)
      !real*8 ejs(npole3b,0:maxvec)
      !real*8 ekc(npole3b,-maxvec:maxvec)
      !real*8 eks(npole3b,-maxvec:maxvec)
      !real*8 elc(npole3b,-maxvec:maxvec)
      !real*8 els(npole3b,-maxvec:maxvec)
      real*8 uind(3,*),eptemp,deptemp(3,*),virtemp(3,3)
      real*8 frecip,uinp(3,*)
      real*8 t1_3b,t2_3b
      logical doi
      integer k1
      real*8, allocatable :: ckr(:)
      real*8, allocatable :: skr(:)
      real*8, allocatable :: cjk(:)
      real*8, allocatable :: sjk(:)
      real*8, allocatable :: s1(:)
      real*8, allocatable :: s2(:)
      real*8, allocatable :: s3(:)
      real*8, allocatable :: s4(:)
      real*8, allocatable :: s3p(:)
      real*8, allocatable :: s4p(:)
      real*8, allocatable :: cm(:)
      real*8, allocatable :: dm(:,:)
      real*8, allocatable :: qm(:,:)
      real*8, allocatable :: um(:,:)
      real*8, allocatable :: ejc(:,:)
      real*8, allocatable :: ejs(:,:)
      real*8, allocatable :: ekc(:,:)
      real*8, allocatable :: eks(:,:)
      real*8, allocatable :: elc(:,:)
      real*8, allocatable :: els(:,:)

      real*8, allocatable :: trqi(:,:)
      real*8, allocatable :: dkx(:)
      real*8, allocatable :: dky(:)
      real*8, allocatable :: dkz(:)
      real*8, allocatable :: qkx(:)
      real*8, allocatable :: qky(:)
      real*8, allocatable :: qkz(:)

      real*8, allocatable :: ckr3b(:)
      real*8, allocatable :: skr3b(:)
      real*8, allocatable :: cjk3b(:)
      real*8, allocatable :: sjk3b(:)
      real*8, allocatable :: ejc3b(:,:)
      real*8, allocatable :: ejs3b(:,:)
      real*8, allocatable :: ekc3b(:,:)
      real*8, allocatable :: eks3b(:,:)
      real*8, allocatable :: elc3b(:,:)
      real*8, allocatable :: els3b(:,:)

      allocate (ckr(npole))
      allocate (skr(npole))
      allocate (cjk(npole))
      allocate (sjk(npole))
      allocate (cm(npole))
      allocate (dm(3,npole))
      allocate (qm(9,npole))
      allocate (um(3,npole3b))
      allocate (ejc(npole,0:maxvec))
      allocate (ejs(npole,0:maxvec))
      allocate (ekc(npole,-maxvec:maxvec))
      allocate (eks(npole,-maxvec:maxvec))
      allocate (elc(npole,-maxvec:maxvec))
      allocate (els(npole,-maxvec:maxvec))

      allocate (ckr3b(npole3b))
      allocate (skr3b(npole3b))
      allocate (cjk3b(npole3b))
      allocate (sjk3b(npole3b))
      allocate (ejc3b(npole3b,0:maxvec))
      allocate (ejs3b(npole3b,0:maxvec))
      allocate (ekc3b(npole3b,-maxvec:maxvec))
      allocate (eks3b(npole3b,-maxvec:maxvec))
      allocate (elc3b(npole3b,-maxvec:maxvec))
      allocate (els3b(npole3b,-maxvec:maxvec))


      allocate (s1(npole))
      allocate (s2(npole))
      allocate (s3(npole3b))
      allocate (s4(npole3b))
      allocate (s3p(npole3b))
      allocate (s4p(npole3b))
         
      allocate (trqi(3,npole))
      allocate (dkx(npole))
      allocate (dky(npole))
      allocate (dkz(npole))
      allocate (qkx(npole))
      allocate (qky(npole))
      allocate (qkz(npole))


c      real*8 aewald3b,aewald3b3b,ewaldcut3b
c      integer jmax3b3b,kmax3b3b,lmax3b3b
c
c     return if the Ewald coefficient is zero
c
c      aewald3b=aewald3b3b

c      if (aewald3b .lt. 1.0d-6)  return

      frecip = 0.5d0
      f = electric / dielec
      term = -0.25d0 / aewald3b**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1


      !print*,"j,k,l in erecip1_3b_totfieldnpole3b",jmax3b,kmax3b,lmax3b

      do l1 = 1, npole3b
         i = pnum(l1)
         um(1,l1) = uind(1,l1)
         um(2,l1) = uind(2,l1)
         um(3,l1) = uind(3,l1)
      end do

      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
c         um(1,l1) = uind(1,l1)
c         um(2,l1) = uind(2,l1)
c         um(3,l1) = uind(3,l1)
         trqi(1,i) = 0.0d0
         trqi(2,i) = 0.0d0
         trqi(3,i) = 0.0d0
      end do

c      do l1 = 1, npole3b
c         i = pnum(l1)
c         trqi(1,i) = 0.0d0
c         trqi(2,i) = 0.0d0
c         trqi(3,i) = 0.0d0
c         trqi(1,l1) = 0.0d0
c         trqi(2,l1) = 0.0d0
c         trqi(3,l1) = 0.0d0
c      end do

      !do l1 = 1, npole3b
      !   i = pnum(l1)
      do i = 1,npole
         l1 = i
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(l1,0) = 1.0d0
         ejs(l1,0) = 0.0d0
         ekc(l1,0) = 1.0d0
         eks(l1,0) = 0.0d0
         elc(l1,0) = 1.0d0
         els(l1,0) = 0.0d0
         ejc(l1,1) = cos(xfr)
         ejs(l1,1) = sin(xfr)
         ekc(l1,1) = cos(yfr)
         eks(l1,1) = sin(yfr)
         elc(l1,1) = cos(zfr)
         els(l1,1) = sin(zfr)
         ekc(l1,-1) = ekc(l1,1)
         eks(l1,-1) = -eks(l1,1)
         elc(l1,-1) = elc(l1,1)
         els(l1,-1) = -els(l1,1)
         do j = 2, jmax3b
            ejc(l1,j) = ejc(l1,j-1)*ejc(l1,1) - ejs(l1,j-1)*ejs(l1,1)
            ejs(l1,j) = ejs(l1,j-1)*ejc(l1,1) + ejc(l1,j-1)*ejs(l1,1)
         end do
         do j = 2, kmax3b
            ekc(l1,j) = ekc(l1,j-1)*ekc(l1,1) - eks(l1,j-1)*eks(l1,1)
            eks(l1,j) = eks(l1,j-1)*ekc(l1,1) + ekc(l1,j-1)*eks(l1,1)
            ekc(l1,-j) = ekc(l1,j)
            eks(l1,-j) = -eks(l1,j)
         end do
         do j = 2, lmax3b
            elc(l1,j) = elc(l1,j-1)*elc(l1,1) - els(l1,j-1)*els(l1,1)
            els(l1,j) = els(l1,j-1)*elc(l1,1) + elc(l1,j-1)*els(l1,1)
            elc(l1,-j) = elc(l1,j)
            els(l1,-j) = -els(l1,j)
         end do
      end do


      do l1 = 1, npole3b
         i = pnum(l1)
      !do i = 1,npole
      !   l1 = i
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc3b(l1,0) = 1.0d0
         ejs3b(l1,0) = 0.0d0
         ekc3b(l1,0) = 1.0d0
         eks3b(l1,0) = 0.0d0
         elc3b(l1,0) = 1.0d0
         els3b(l1,0) = 0.0d0
         ejc3b(l1,1) = cos(xfr)
         ejs3b(l1,1) = sin(xfr)
         ekc3b(l1,1) = cos(yfr)
         eks3b(l1,1) = sin(yfr)
         elc3b(l1,1) = cos(zfr)
         els3b(l1,1) = sin(zfr)
         ekc3b(l1,-1) = ekc3b(l1,1)
         eks3b(l1,-1) = -eks3b(l1,1)
         elc3b(l1,-1) = elc3b(l1,1)
         els3b(l1,-1) = -els3b(l1,1)
         do j = 2, jmax3b
            ejc3b(l1,j) = ejc3b(l1,j-1)*ejc3b(l1,1) 
     &                    - ejs3b(l1,j-1)*ejs3b(l1,1)
            ejs3b(l1,j) = ejs3b(l1,j-1)*ejc3b(l1,1) 
     &                    + ejc3b(l1,j-1)*ejs3b(l1,1)
         end do
         do j = 2, kmax3b
            ekc3b(l1,j) = ekc3b(l1,j-1)*ekc3b(l1,1) 
     &                      - eks3b(l1,j-1)*eks3b(l1,1)
            eks3b(l1,j) = eks3b(l1,j-1)*ekc3b(l1,1) 
     &                      + ekc3b(l1,j-1)*eks3b(l1,1)
            ekc3b(l1,-j) = ekc3b(l1,j)
            eks3b(l1,-j) = -eks3b(l1,j)
         end do
         do j = 2, lmax3b
            elc3b(l1,j) = elc3b(l1,j-1)*elc3b(l1,1) 
     &                      - els3b(l1,j-1)*els3b(l1,1)
            els3b(l1,j) = els3b(l1,j-1)*elc3b(l1,1) 
     &                     + elc3b(l1,j-1)*els3b(l1,1)
            elc3b(l1,-j) = elc3b(l1,j)
            els3b(l1,-j) = -els3b(l1,j)
         end do
      end do


      do j = jmin, jmax3b
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax3b
            rk = 2.0d0 * pi * dble(k)
            !do l1 = 1, npole3b
            !   i = pnum(l1)
            do i = 1,npole    
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
c               cjk(l1) = ejc(l1,j)*ekc(l1,k) - ejs(l1,j)*eks(l1,k)
c               sjk(l1) = ejs(l1,j)*ekc(l1,k) + ejc(l1,j)*eks(l1,k)
            end do
            do l1 =1, npole3b
               cjk3b(l1) = ejc3b(l1,j)*ekc3b(l1,k) 
     &                      - ejs3b(l1,j)*eks3b(l1,k)
               sjk3b(l1) = ejs3b(l1,j)*ekc3b(l1,k) 
     &                      + ejc3b(l1,j)*eks3b(l1,k)
            end do
            do l = lmin, lmax3b
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
c               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t1_3b = 0.0d0
                  t2_3b = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  t3p = 0.0d0
                  t4p = 0.0d0
                  do m2 = 1, 3
                     do m1 = 1, 3
                        t5(m1,m2) = 0.0d0
                        t5u(m1,m2) = 0.0d0
                        t5p(m1,m2) = 0.0d0
                        t6(m1,m2) = 0.0d0
                        t6u(m1,m2) = 0.0d0
                        t6p(m1,m2) = 0.0d0
                     end do
                  end do
                  !do l1 = 1, npole3b
                  !   i = pnum(l1)
                  do i = 1,npole
                     doi=.false.
                     do k1=1,npole3b
                       if(pnum(k1).eq.i) then
                         doi=.true.
                         l1=k1
                         goto 21
                       end if
                     end do
   21    continue               
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                    
                     if(doi) then
                 ckr3b(l1)=cjk3b(l1)*elc3b(l1,l) - sjk3b(l1)*els3b(l1,l)
                 skr3b(l1)=sjk3b(l1)*elc3b(l1,l) + cjk3b(l1)*els3b(l1,l)
                     end if

                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     dkx(i) = h3*dm(2,i) - h2*dm(3,i)
                     dky(i) = h1*dm(3,i) - h3*dm(1,i)
                     dkz(i) = h2*dm(1,i) - h1*dm(2,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     qkx(i) = h3*q2 - h2*q3
                     qky(i) = h1*q3 - h3*q1
                     qkz(i) = h2*q1 - h1*q2
                     if(doi) then
                      uk = h1*uind(1,l1) + h2*uind(2,l1) + h3*uind(3,l1)
                     ukp = h1*uinp(1,l1) + h2*uinp(2,l1) + h3*uinp(3,l1)
                     else
                      uk = 0.0d0
                      ukp = 0.0d0       
                     end if
                     s1(i) = (ck-qk)*skr(i) + dk*ckr(i)
                     s2(i) = (ck-qk)*ckr(i) - dk*skr(i)

                     if (doi) then
                     ! s3(l1) = uk * ckr3b(l1)
                     ! s4(l1) = -uk * skr3b(l1)
                     ! s3p(l1) = ukp * ckr3b(l1)
                     ! s4p(l1) = -ukp * skr3b(l1)
                      s3(l1) = uk * ckr(i)
                      s4(l1) = -uk * skr(i)
                      s3p(l1) = ukp * ckr(i)
                      s4p(l1) = -ukp * skr(i)
                      t3 = t3 + s3(l1)
                      t4 = t4 + s4(l1)
                      t3p = t3p + s3p(l1)
                      t4p = t4p + s4p(l1)
                     end if

                     t1 = t1 + s1(i)
                     t2 = t2 + s2(i)

c
c     terms needed for subsequent virial tensor calculation
c
                  qt(1,1)=h1*(h1*qm(1,i)+h2*qm(4,i)+h3*qm(7,i))
                  qt(2,1)=h1*(h1*qm(2,i)+h2*qm(5,i)+h3*qm(8,i))
                  qt(3,1) = h1*(h1*qm(3,i)+h2*qm(6,i)+h3*qm(9,i))
                  qt(1,2) = h2*(h1*qm(1,i)+h2*qm(4,i)+h3*qm(7,i))
                  qt(2,2) = h2*(h1*qm(2,i)+h2*qm(5,i)+h3*qm(8,i))
                  qt(3,2) = h2*(h1*qm(3,i)+h2*qm(6,i)+h3*qm(9,i))
                  qt(1,3) = h3*(h1*qm(1,i)+h2*qm(4,i)+h3*qm(7,i))
                  qt(2,3) = h3*(h1*qm(2,i)+h2*qm(5,i)+h3*qm(8,i))
                  qt(3,3) = h3*(h1*qm(3,i)+h2*qm(6,i)+h3*qm(9,i))
                     dt(1,1) = h1 * dm(1,i)
                     dt(2,1) = h1 * dm(2,i)
                     dt(3,1) = h1 * dm(3,i)
                     dt(1,2) = h2 * dm(1,i)
                     dt(2,2) = h2 * dm(2,i)
                     dt(3,2) = h2 * dm(3,i)
                     dt(1,3) = h3 * dm(1,i)
                     dt(2,3) = h3 * dm(2,i)
                     dt(3,3) = h3 * dm(3,i)
                     if(doi) then
                      dtu(1,1) = h1 * uind(1,l1)
                      dtu(2,1) = h1 * uind(2,l1)
                      dtu(3,1) = h1 * uind(3,l1)
                      dtu(1,2) = h2 * uind(1,l1)
                      dtu(2,2) = h2 * uind(2,l1)
                      dtu(3,2) = h2 * uind(3,l1)
                      dtu(1,3) = h3 * uind(1,l1)
                      dtu(2,3) = h3 * uind(2,l1)
                      dtu(3,3) = h3 * uind(3,l1)
                      dtp(1,1) = h1 * uinp(1,l1)
                      dtp(2,1) = h1 * uinp(2,l1)
                      dtp(3,1) = h1 * uinp(3,l1)
                      dtp(1,2) = h2 * uinp(1,l1)
                      dtp(2,2) = h2 * uinp(2,l1)
                      dtp(3,2) = h2 * uinp(3,l1)
                      dtp(1,3) = h3 * uinp(1,l1)
                      dtp(2,3) = h3 * uinp(2,l1)
                      dtp(3,3) = h3 * uinp(3,l1)
                     else
                      dtu(1,1) = 0.0d0
                      dtu(2,1) = 0.0d0
                      dtu(3,1) = 0.0d0
                      dtu(1,2) = 0.0d0
                      dtu(2,2) = 0.0d0
                      dtu(3,2) = 0.0d0
                      dtu(1,3) = 0.0d0
                      dtu(2,3) = 0.0d0
                      dtu(3,3) = 0.0d0
                      dtp(1,1) = 0.0d0
                      dtp(2,1) = 0.0d0
                      dtp(3,1) = 0.0d0
                      dtp(1,2) = 0.0d0
                      dtp(2,2) = 0.0d0
                      dtp(3,2) = 0.0d0
                      dtp(1,3) = 0.0d0
                      dtp(2,3) = 0.0d0
                      dtp(3,3) = 0.0d0
                     end if
                     do m2 = 1, 3
                        do m1 = 1, 3
                           t5(m1,m2) = t5(m1,m2) - dt(m1,m2)*ckr(i)
     &                                    + 2.0d0*qt(m1,m2)*skr(i)
                           t5u(m1,m2) = t5u(m1,m2) - dtu(m1,m2)*ckr(i)
                           t5p(m1,m2) = t5p(m1,m2) - dtp(m1,m2)*ckr(i)
                           t6(m1,m2) = t6(m1,m2) + dt(m1,m2)*skr(i)
     &                                    + 2.0d0*qt(m1,m2)*ckr(i)
                           t6u(m1,m2) = t6u(m1,m2) + dtu(m1,m2)*skr(i)
                           t6p(m1,m2) = t6p(m1,m2) + dtp(m1,m2)*skr(i)
                        end do
                     end do
                  end do
c
c     get the energy contributions for current reciprocal vector
c
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  ei = expterm * (t1*t3+t2*t4)
                 
                  etot = e + ei
                  eptemp = eptemp + ei

                 ! uterm = expterm * (t1*(t3+t3p) + t3*t3p
     &           !                      + t2*(t4+t4p) + t4*t4p)

                  uterm = expterm * (t1*(t3+t3p) 
     &                                 + t2*(t4+t4p) )

                 ! do m2 = 1, 3
                 !    do m1 = 1, 3

c     NEW wterm w/ just the permanent elec part removed (already done)
                 !       wterm(m1,m2) = 2.0d0 * expterm
     &           !          * (
     &           !             0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
     &           !             + t2*(t6u(m1,m2)+t6p(m1,m2))
     &           !             + (t3+t3p)*t5(m1,m2)
     &           !             + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
     &           !             + (t4+t4p)*t6(m1,m2)
     &           !             + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))
                 !
                 !    end do
                 ! end do


                  do m2 = 1, 3
                     do m1 = 1, 3

c     NEW wterm w/ just the permanent elec part removed (already done)
                        wterm(m1,m2) = 2.0d0 * expterm
     &                     * (
     &                        0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
     &                        + (t3+t3p)*t5(m1,m2)
     &                        !+ t3*t5p(m1,m2) + t3p*t5u(m1,m2)
     &                        + (t4+t4p)*t6(m1,m2)
     &                        ))
                 
                     end do
                  end do

                  wterm(2,1) = 0.5d0 * (wterm(2,1)+wterm(1,2))
                  wterm(3,1) = 0.5d0 * (wterm(3,1)+wterm(1,3))
                  wterm(3,2) = 0.5d0 * (wterm(3,2)+wterm(2,3))
                  wterm(1,2) = wterm(2,1)
                  wterm(1,3) = wterm(3,1)
                  wterm(2,3) = wterm(3,2)
                  vterm = 2.0d0 * uterm * (1.0d0-term*hsq) / hsq

                  virtemp(1,1) = virtemp(1,1) + h1*h1*vterm + wterm(1,1)
     &                            - uterm
                  virtemp(2,1) = virtemp(2,1) + h2*h1*vterm + wterm(2,1)
                  virtemp(3,1) = virtemp(3,1) + h3*h1*vterm + wterm(3,1)
                  virtemp(1,2) = virtemp(1,2) + h1*h2*vterm + wterm(1,2)
                  virtemp(2,2) = virtemp(2,2) + h2*h2*vterm + wterm(2,2)
     &                            - uterm
                  virtemp(3,2) = virtemp(3,2) + h3*h2*vterm + wterm(3,2)
                  virtemp(1,3) = virtemp(1,3) + h1*h3*vterm + wterm(1,3)
                  virtemp(2,3) = virtemp(2,3) + h2*h3*vterm + wterm(2,3)
                  virtemp(3,3) = virtemp(3,3) + h3*h3*vterm + wterm(3,3)
     &                             - uterm

c
c     get the force contributions for current reciprocal vector
c
                  expterm = 2.0d0 * expterm
                 ! do l1 = 1, npole3b
                 !    i = pnum(l1)
                  do i = 1, npole
                     ii = ipole(i)
                     doi=.false.
                     do k1=1,npole3b
                       if(pnum(k1).eq.i) then
                         doi=.true.
                         l1=k1
                         goto 22
                       end if
                     end do
   22    continue

                    ! dei = 0.5d0 * expterm * ((s4(l1)+s4p(l1))*t1
     &              !                         -(s3(l1)+s3p(l1))*t2
     &              !                  +s2(l1)*(t3+t3p)-s1(l1)*(t4+t4p))
 
c                     if (poltyp .eq. 'MUTUAL') then
                    !     dei = dei + 0.5d0 * expterm
     &              !              * (s4p(l1)*t3+s4(l1)*t3p
     &              !                -s3p(l1)*t4-s3(l1)*t4p)
c                     end if

                     if(doi) then
                     dei = 0.5d0 * expterm * ((s4(l1)+s4p(l1))*t1
     &                                       -(s3(l1)+s3p(l1))*t2
     &                                +s2(i)*(t3+t3p)-s1(i)*(t4+t4p))
                     else
                     dei = 0.5d0 * expterm * (
     &                                       
     &                                +s2(i)*(t3+t3p)-s1(i)*(t4+t4p))
                     end if

                    ! det1i = 0.5d0 * expterm * (skr(l1)*(t4+t4p)
     &              !                           -ckr(l1)*(t3+t3p))
                    ! det2i =expterm*(ckr(l1)*(t4+t4p)+skr(l1)*(t3+t3p))

                      
                     det1i = 0.5d0 * expterm * (skr(i)*(t4+t4p)
     &                                         -ckr(i)*(t3+t3p))
                     det2i =expterm*(ckr(i)*(t4+t4p)+skr(i)*(t3+t3p))

                     deptemp(1,i) = deptemp(1,i) + h1*dei
                     deptemp(2,i) = deptemp(2,i) + h2*dei
                     deptemp(3,i) = deptemp(3,i) + h3*dei

                     !deptemp(1,l1) = deptemp(1,l1) + h1*dei
                     !deptemp(2,l1) = deptemp(2,l1) + h2*dei
                     !deptemp(3,l1) = deptemp(3,l1) + h3*dei


c                     trqi(1,i) = trqi(1,i) + dkx(l1)*det1i
c     &                               + qkx(l1)*det2i
c                     trqi(2,i) = trqi(2,i) + dky(l1)*det1i
c     &                               + qky(l1)*det2i
c                     trqi(3,i) = trqi(3,i) + dkz(l1)*det1i
c     &                               + qkz(l1)*det2i

                     trqi(1,i) = trqi(1,i) + dkx(i)*det1i
     &                               + qkx(i)*det2i
                     trqi(2,i) = trqi(2,i) + dky(i)*det1i
     &                               + qky(i)*det2i
                     trqi(3,i) = trqi(3,i) + dkz(i)*det1i
     &                               + qkz(i)*det2i

                  end do
c               end if
            end do
            lmin = -lmax3b
         end do
         kmin = -kmax3b
      end do

c      call torque_ewreg3b (trqi,deptemp,pnum,npole3b)
      ! call torque2_3b (npole3b,pnum,trqi,deptemp) 
       call torque2(trqi,deptemp)

      deallocate (ckr)
      deallocate (skr)
      deallocate (cjk)
      deallocate (sjk)
      deallocate (cm)
      deallocate (dm)
      deallocate (qm)
      deallocate (um)
      deallocate (ejc)
      deallocate (ejs)
      deallocate (ekc)
      deallocate (eks)
      deallocate (elc)
      deallocate (els)

      deallocate (ckr3b)
      deallocate (skr3b)
      deallocate (cjk3b)
      deallocate (sjk3b)
      deallocate (ejc3b)
      deallocate (ejs3b)
      deallocate (ekc3b)
      deallocate (eks3b)
      deallocate (elc3b)
      deallocate (els3b)

      deallocate (s1)
      deallocate (s2)
      deallocate (s3)
      deallocate (s4)
      deallocate (s3p)
      deallocate (s4p)
      deallocate (trqi)
      deallocate (dkx)
      deallocate (dky)
      deallocate (dkz)
      deallocate (qkx)
      deallocate (qky)
      deallocate (qkz)

      return
      end


