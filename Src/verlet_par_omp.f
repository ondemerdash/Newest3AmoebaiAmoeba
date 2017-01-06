c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine verlet  --  Verlet molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "verlet" performs a single molecular dynamics time step
c     via the velocity Verlet multistep recursion formula
c
c
      subroutine verlet_pt1 (istep,dt)
      use sizes
      use atomid
      use atoms
      use freeze
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,istep
      real*8 dt,dt_2
      real*8 etot,epot
      real*8 eksum
      real*8 temp,pres
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c      real*8, allocatable :: derivs(:,:)
      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time

c
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     make half-step temperature and pressure corrections
c
c      call system_clock(t1,clock_rate)
      call temper (dt)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing temper",tot_time

c      print*,"dt in verlet_pt1",dt
c      print*,"dt_2 in verlet_pt1",dt_2
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c      allocate (derivs(3,n))
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Verlet recursion
c
c      call system_clock(t1,clock_rate)

!$OMP PARALLEL DO default(shared) private(i,j)
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i) + v(1,i)*dt
            y(i) = y(i) + v(2,i)*dt
            z(i) = z(i) + v(3,i)*dt
         end if
      end do
!$OMP END PARALLEL DO

c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing Verlet 1/2step v, full step x",tot_time

c      print*,"In verlet_pt1 x(3) y(3) z(3)",x(3),y(3),z(3)
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)

c
c     get the potential energy and atomic forces
c
      return
      end

c      call gradient (epot,derivs)

c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
c      do i = 1, n
c         if (use(i)) then
c            do j = 1, 3
c               a(j,i) = -convert * derivs(j,i) / mass(i)
c               v(j,i) = v(j,i) + a(j,i)*dt_2
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
      subroutine verlet_pt2 (istep,dt,energy)
      use sizes
      use atomid
      use atoms
      use freeze
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,istep
      real*8 dt,dt_2
      real*8 etot,epot,energy
      real*8 eksum
      real*8 temp,pres
      real*8 ekin(3,3)
      real*8 stress(3,3)
c      real*8, allocatable :: xold(:)
c      real*8, allocatable :: yold(:)
c      real*8, allocatable :: zold(:)
c      real*8, allocatable :: derivs(:,:)
      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time

      epot = energy
      if (use_rattle)  call rattle2 (dt)
c      print*,"In verlet_pt2",dt
      
c
c     make full-step temperature and pressure corrections
c
c      call system_clock(t1,clock_rate)
      call temper2 (dt,eksum,ekin,temp)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing temper2",tot_time

c      call system_clock(t1,clock_rate)
      call pressure (dt,epot,ekin,temp,pres,stress)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing pressure",tot_time

c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
c      call system_clock(t1,clock_rate)      
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing mdstat",tot_time

c      call system_clock(t1,clock_rate)
      call mdsave (istep,dt,epot,eksum)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing mdsave",tot_time

c      call system_clock(t1,clock_rate)
      call mdrest (istep)
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing mdrest",tot_time

      return
      end
