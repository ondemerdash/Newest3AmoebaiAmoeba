c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     via the Beeman multistep recursion formula; uses original
c     coefficients or Bernie Brooks' "Better Beeman" values
c
c     literature references:
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
c     Temperature and Pressure", DCRT Report, NIH, April 1988
c
c
      subroutine beeman_pt1 (istep,dt,dt_x,part1)
      use sizes
      use atomid
      use atoms
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,istep
      real*8 dt,dt_x,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c      real*8, allocatable :: derivs(:,:)
c
c
c     set time values and coefficients for Beeman integration
c
c      factor = dble(bmnmix)
c      dt_x = dt / factor
c      part1 = 0.5d0*factor + 1.0d0
c      part2 = part1 - 2.0d0
c
c     make half-step temperature and pressure corrections
c
      call temper (dt)
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c      allocate (derivs(3,n))
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + (part1*a(j,i)-aalt(j,i))*dt_x
            end do
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i) + v(1,i)*dt
            y(i) = y(i) + v(2,i)*dt
            z(i) = z(i) + v(3,i)*dt
         end if
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)

      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      return 
      end 
c
c     get the potential energy and atomic forces
c

c      call gradient (epot,derivs)

c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
c      do i = 1, n
c         if (use(i)) then
c            do j = 1, 3
c               aalt(j,i) = a(j,i)
c               a(j,i) = -convert * derivs(j,i) / mass(i)
c               v(j,i) = v(j,i) + (part2*a(j,i)+aalt(j,i))*dt_x
c            end do
c         end if
c      end do
c
c     perform deallocation of some local arrays
c
c      deallocate (xold)
c      deallocate (yold)
c      deallocate (zold)
c      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      subroutine beeman_pt2 (istep,dt,energy)
      use sizes
      use atomid
      use atoms
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,istep
      real*8 dt,dt_x,factor
      real*8 etot,eksum,epot,energy
      real*8 temp,pres
      real*8 part1,part2
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c      real*8, allocatable :: derivs(:,:)
c
c
c     set time values and coefficients for Beeman integration
c
c      factor = dble(bmnmix)
c      dt_x = dt / factor
c      part1 = 0.5d0*factor + 1.0d0
c      part2 = part1 - 2.0d0
      epot = energy

      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature and pressure corrections
c
      call temper2 (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)
      return
      end
