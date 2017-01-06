c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine gradrot  --  energy and torsional derivs  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "gradrot" calls subroutines to calculate the potential
c     energy and its torsional first derivatives
c
c
      subroutine gradrot (energy,derivs)
      use sizes
      use atoms
      use deriv
      use domega
      use omega
      use potent
      use rotbnd
      implicit none
      integer i,j,k
      integer base,partner
      real*8 energy,norm
      real*8 xatom,yatom,zatom
      real*8 xdist,ydist,zdist
      real*8 xterm,yterm,zterm
      real*8 derivs(*)
      real*8, allocatable :: g(:,:)
      logical first
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(tesum))  allocate (tesum(nomega))
         if (.not. allocated(teb))  allocate (teb(nomega))
         if (.not. allocated(tea))  allocate (tea(nomega))
         if (.not. allocated(teba))  allocate (teba(nomega))
         if (.not. allocated(teub))  allocate (teub(nomega))
         if (.not. allocated(teaa))  allocate (teaa(nomega))
         if (.not. allocated(teopb))  allocate (teopb(nomega))
         if (.not. allocated(teopd))  allocate (teopd(nomega))
         if (.not. allocated(teid))  allocate (teid(nomega))
         if (.not. allocated(teit))  allocate (teit(nomega))
         if (.not. allocated(tet))  allocate (tet(nomega))
         if (.not. allocated(tept))  allocate (tept(nomega))
         if (.not. allocated(tebt))  allocate (tebt(nomega))
         if (.not. allocated(teat))  allocate (teat(nomega))
         if (.not. allocated(tett))  allocate (tett(nomega))
         if (.not. allocated(tev))  allocate (tev(nomega))
         if (.not. allocated(tec))  allocate (tec(nomega))
         if (.not. allocated(tecd))  allocate (tecd(nomega))
         if (.not. allocated(ted))  allocate (ted(nomega))
         if (.not. allocated(tem))  allocate (tem(nomega))
         if (.not. allocated(tep))  allocate (tep(nomega))
         if (.not. allocated(ter))  allocate (ter(nomega))
         if (.not. allocated(tes))  allocate (tes(nomega))
         if (.not. allocated(telf))  allocate (telf(nomega))
         if (.not. allocated(teg))  allocate (teg(nomega))
         if (.not. allocated(tex))  allocate (tex(nomega))
      end if
c
c     zero out individual components of torsional derivatives
c
      do i = 1, nomega
         teb(i) = 0.0d0
         tea(i) = 0.0d0
         teba(i) = 0.0d0
         teub(i) = 0.0d0
         teaa(i) = 0.0d0
         teopb(i) = 0.0d0
         teopd(i) = 0.0d0
         teid(i) = 0.0d0
         teit(i) = 0.0d0
         tet(i) = 0.0d0
         tept(i) = 0.0d0
         tebt(i) = 0.0d0
         teat(i) = 0.0d0
         tett(i) = 0.0d0
         tev(i) = 0.0d0
         tec(i) = 0.0d0
         tecd(i) = 0.0d0
         ted(i) = 0.0d0
         tem(i) = 0.0d0
         tep(i) = 0.0d0
         ter(i) = 0.0d0
         tes(i) = 0.0d0
         telf(i) = 0.0d0
         teg(i) = 0.0d0
         tex(i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (g(3,n))
c
c     calculate the energy and Cartesian first derivatives
c
      call gradient (energy,g)
c
c     perform deallocation of some local arrays
c
      deallocate (g)
c
c     transform Cartesian derivatives to torsional space
c
      do i = 1, nomega
         base = iomega(1,i)
         partner = iomega(2,i)
         call rotlist (base,partner)
         xdist = x(base) - x(partner)
         ydist = y(base) - y(partner)
         zdist = z(base) - z(partner)
         norm = sqrt(xdist**2 + ydist**2 + zdist**2)
         xdist = xdist / norm
         ydist = ydist / norm
         zdist = zdist / norm
         do j = 1, nrot
            k = rot(j)
            xatom = x(k) - x(base)
            yatom = y(k) - y(base)
            zatom = z(k) - z(base)
            xterm = ydist*zatom - zdist*yatom
            yterm = zdist*xatom - xdist*zatom
            zterm = xdist*yatom - ydist*xatom
            teb(i) = teb(i) + deb(1,k)*xterm + deb(2,k)*yterm
     &                              + deb(3,k)*zterm
            tea(i) = tea(i) + dea(1,k)*xterm + dea(2,k)*yterm
     &                              + dea(3,k)*zterm
            teba(i) = teba(i) + deba(1,k)*xterm + deba(2,k)*yterm
     &                              + deba(3,k)*zterm
            teub(i) = teub(i) + deub(1,k)*xterm + deub(2,k)*yterm
     &                              + deub(3,k)*zterm
            teaa(i) = teaa(i) + deaa(1,k)*xterm + deaa(2,k)*yterm
     &                              + deaa(3,k)*zterm
            teopb(i) = teopb(i) + deopb(1,k)*xterm + deopb(2,k)*yterm
     &                              + deopb(3,k)*zterm
            teopd(i) = teopd(i) + deopd(1,k)*xterm + deopd(2,k)*yterm
     &                              + deopd(3,k)*zterm
            teid(i) = teid(i) + deid(1,k)*xterm + deid(2,k)*yterm
     &                              + deid(3,k)*zterm
            teit(i) = teit(i) + deit(1,k)*xterm + deit(2,k)*yterm
     &                              + deit(3,k)*zterm
            tet(i) = tet(i) + det(1,k)*xterm + det(2,k)*yterm
     &                              + det(3,k)*zterm
            tept(i) = tept(i) + dept(1,k)*xterm + dept(2,k)*yterm
     &                              + dept(3,k)*zterm
            tebt(i) = tebt(i) + debt(1,k)*xterm + debt(2,k)*yterm
     &                              + debt(3,k)*zterm
            teat(i) = teat(i) + deat(1,k)*xterm + deat(2,k)*yterm
     &                              + deat(3,k)*zterm
            tett(i) = tett(i) + dett(1,k)*xterm + dett(2,k)*yterm
     &                              + dett(3,k)*zterm
            tev(i) = tev(i) + dev(1,k)*xterm + dev(2,k)*yterm
     &                              + dev(3,k)*zterm
            tec(i) = tec(i) + dec(1,k)*xterm + dec(2,k)*yterm
     &                              + dec(3,k)*zterm
            tecd(i) = tecd(i) + decd(1,k)*xterm + decd(2,k)*yterm
     &                              + decd(3,k)*zterm
            ted(i) = ted(i) + ded(1,k)*xterm + ded(2,k)*yterm
     &                              + ded(3,k)*zterm
            tem(i) = tem(i) + dem(1,k)*xterm + dem(2,k)*yterm
     &                              + dem(3,k)*zterm
            tep(i) = tep(i) + dep(1,k)*xterm + dep(2,k)*yterm
     &                              + dep(3,k)*zterm
            ter(i) = ter(i) + der(1,k)*xterm + der(2,k)*yterm
     &                              + der(3,k)*zterm
            tes(i) = tes(i) + des(1,k)*xterm + des(2,k)*yterm
     &                              + des(3,k)*zterm
            telf(i) = telf(i) + delf(1,k)*xterm + delf(2,k)*yterm
     &                              + delf(3,k)*zterm
            teg(i) = teg(i) + deg(1,k)*xterm + deg(2,k)*yterm
     &                              + deg(3,k)*zterm
            tex(i) = tex(i) + dex(1,k)*xterm + dex(2,k)*yterm
     &                              + dex(3,k)*zterm
         end do
      end do
c
c     sum up to give the total torsional first derivatives
c
      do i = 1, nomega
         tesum(i) = teb(i) + tea(i) + teba(i) + teub(i) + teaa(i)
     &                 + teopb(i) + teopd(i) + teid(i) + teit(i)
     &                 + tet(i) + tept(i) + tebt(i) + teat(i) + tett(i)
     &                 + tev(i) + tec(i) + tecd(i) + ted(i) + tem(i)
     &                 + tep(i) + ter(i) + tes(i) + telf(i) + teg(i)
     &                 + tex(i)
         derivs(i) = tesum(i)
      end do
      return
      end
