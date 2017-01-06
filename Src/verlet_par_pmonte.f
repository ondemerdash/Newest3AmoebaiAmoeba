c
      subroutine verlet_pt2_pmonte (istep,dt,energy)
      use sizes
      use atomid
      use atoms
      use freeze
      use moldyn
      use units
      use usage
      use mpidat
      use boxes
      use virial
      implicit none
      include 'mpif.h'
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
      real*8 factor
      integer ierr

      epot = energy
      
      if(taskid.eq.master) then

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
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
c
c     set isotropic pressure to the average of tensor diagonal
c
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0

      end if
c      call system_clock(t1,clock_rate)
      call mpi_bcast(temp,1,mpi_real8,master,
     &       mpi_comm_world,ierr)
      call mpi_bcast(pres,1,mpi_real8,master,
     &       mpi_comm_world,ierr)

CCC LEFT OFF HERE
      ! call pressure (dt,epot,ekin,temp,pres,stress)
      
       call pmonte_par2(epot,temp,istep)
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
      if(taskid.eq.master) then

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
      end if

c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing mdrest",tot_time

      return
      end
