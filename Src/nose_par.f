c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2011 by Teresa Head-Gordon & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nose  --  Nose-Hoover NPT molecular dynamics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nose" performs a single molecular dynamics time step via
c     a Nose-Hoover extended system isothermal-isobaric algorithm
c
c     literature reference:
c
c     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
c     "Explicit Reversible Integrators for Extended Systems Dynamics",
c     Molecular Physics, 87, 1117-1157 (1996)
c
c     original version written by Teresa Head-Gordon, November 2011
c
c
      subroutine nose_pt1 (istep,dt,press)
      use sizes
      use atomid
      use atoms
      use bath
      use boxes
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,istep
      real*8 dt,dt_2
      real*8 epot,etot
      real*8 eksum,temp
      real*8 pres,press
      real*8 poly,factor
      real*8 term,expterm
      real*8 term2,eterm2
      real*8 e2,e4,e6,e8
      real*8 ekin(3,3)
      real*8 stress(3,3)
c      real*8, allocatable :: derivs(:,:)
c      save press
c
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
      if (istep .eq. 1)  press = atmsph
c
c     update thermostat and barostat values, scale atomic velocities
c
      call hoover (dt,press)
c
c     get half-step velocities via Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     update atomic positions via coupling to barostat
c
      term = vbar * dt_2
      term2 = term * term
      expterm = exp(term)
      eterm2 = expterm * expterm
      e2 = 1.0d0 / 6.0d0
      e4 = e2 / 20.0d0
      e6 = e4 / 42.0d0
      e8 = e6 / 72.0d0
      poly = 1.0d0 + term2*(e2+term2*(e4+term2*(e6+term2*e8)))
      poly = expterm * poly * dt
      do i = 1, n
         if (use(i)) then
            x(i) = x(i)*eterm2 + v(1,i)*poly
            y(i) = y(i)*eterm2 + v(2,i)*poly
            z(i) = z(i)*eterm2 + v(3,i)*poly
         end if
      end do
c
c     constraints under NH-NPT require the ROLL algorithm
c
      if (use_rattle)  call fatal
c
c     update the periodic box size and total volume
c
      xbox = xbox * eterm2
      ybox = ybox * eterm2
      zbox = zbox * eterm2
      call lattice

      return
      end

c
c     perform dynamic allocation of some local arrays
c

c      allocate (derivs(3,n))

c
c     get the potential energy and atomic forces
c
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

      subroutine nose_pt2 (istep,dt,press,energy)
      use sizes
      use atomid
      use atoms
      use bath
      use boxes
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,istep
      real*8 dt,dt_2
      real*8 epot,etot
      real*8 eksum,temp,energy
      real*8 pres,press
      real*8 poly,factor
      real*8 term,expterm
      real*8 term2,eterm2
      real*8 e2,e4,e6,e8
      real*8 ekin(3,3)
      real*8 stress(3,3)
c      real*8, allocatable :: derivs(:,:)
c      save press

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
c      deallocate (derivs)
c
c     constraints under NH-NPT require the ROLL algorithm
c
      epot=energy

      if (use_rattle)  call fatal
c
c     update thermostat and barostat values, scale atomic velocities
c
      call hoover (dt,press)
c
c     set isotropic pressure to the average of tensor diagonal
c
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (-vir(j,i))
         end do
      end do
      press = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     accumulate the kinetic energy and its outer product
c
      call kinetic (eksum,ekin)
c
c     calculate the stress tensor for anisotropic systems
c
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     get the instantaneous temperature from the kinetic energy
c
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      etot = epot + eksum
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
c      call mdstat_propertiesDaveDeTo_OND(istep,dt,etot,epot,eksum,
c     & temp,pres)
      call mdsave (istep,dt,epot,eksum)
c      call mdsave_writeDyns (istep,dt,epot,eksum)
      call mdrest (istep)
      return
      end

      subroutine bcast_lattice
      use mpidat
      use boxes
      use cell
      implicit none
      include 'mpif.h'
      integer ierr
      call mpi_bcast(xcell,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(ycell,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(zcell,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(xcell2,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(ycell2,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(zcell2,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(xbox,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(ybox,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(zbox,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(xbox2,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(ybox2,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(zbox2,1,mpi_real8,master,
     & mpi_comm_world,ierr)

      call mpi_bcast(beta_sin,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(beta_cos,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(gamma_sin,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(gamma_cos,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(beta_term,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(gamma_term,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(box34,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(volbox,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(alpha,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(beta,1,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(gamma,1,mpi_real8,master,
     & mpi_comm_world,ierr)


c    BROADCAST VARIABLES THAT ARE BEING USED BY PARTICLE MESH EWALD,
c    SINCE THE PME ROUTINE IS CALLED BY MPI TASKID=1, NOT TASKID=0. IT IS,
c    OF COURSE, UPDATED BY TASKID=0 "LATTICE", WHICH IS CALLED IN ANY CONSTANT PRESSURE RUN
      call mpi_bcast(lvec,3*3,mpi_real8,master,
     & mpi_comm_world,ierr)
      call mpi_bcast(recip,3*3,mpi_real8,master,
     & mpi_comm_world,ierr)
      return
      end
