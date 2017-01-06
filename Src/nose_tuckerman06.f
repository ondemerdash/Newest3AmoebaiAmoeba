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
c     Corrections to velocity update for NPT per Mark Tuckerman and in reference to:
c     M. E. Tuckerman, J. Alejandre, R. Lopez-Rendon, A. L. Jochim 
c     and G. J. Martyna, "A Liouville-operator derived measure-preserving integrator
c     for molecular dynamics simulations in the isothermal-isobaric ensemble",
c     Journal of Physics A:  Mathematical and General, 39, 5629-5651 (2006)
c
      subroutine nose_tuckerman06 (istep,dt)
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
      real*8, allocatable :: derivs(:,:)
c 'alph' corresponds to the variable alpha in the paper (1+1/n)
      real*8 alph
      parameter (alph=dble(1+1/n))
      save press
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
c     get half-step velocities via Eq. 5.13 of Tuckerman et al., J.Phys. A, 2006
c
      velterm = -alph*vbar*dt_2/2.0d0
      velterm2 = velterm * velterm
      expvelterm = exp(velterm)
      evelterm2=expvelterm * expvelterm
      e2 = 1.0d0 / 6.0d0
      e4 = e2 / 20.0d0
      e6 = e4 / 42.0d0
      e8 = e6 / 72.0d0
      velpoly = 1.0d0 + velterm2*(e2+velterm2*(e4+velterm2*(e6
     &      +velterm2*e8)))
      velpoly = expvelterm * velpoly * dt_2
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
c               v(j,i) = v(j,i) + a(j,i)*dt_2
               v(j,i) = v(j,i)*evelterm2 + a(j,i)*velpoly
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
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using Eq. 5.13
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
c               v(j,i) = v(j,i) + a(j,i)*dt_2
               v(j,i) = v(j,i)*evelterm2 + a(j,i)*velpoly
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     constraints under NH-NPT require the ROLL algorithm
c
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
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine hoover  --  Nose-Hoover thermostat/barostat  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "hoover" applies a combined thermostat and barostat via a
c     Nose-Hoover chain algorithm
c
c
      subroutine hoover (dt,press)
      use sizes
      use atoms
      use bath
      use boxes
      use mdstuf
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer nc,ns
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8
      real*8 ekt,eksum
      real*8 df,odnf,gn1kt
      real*8 press,dpress
      real*8 expterm,scale
      real*8 w(3),ekin(3,3)
c
c
c     find kinetic energy and set an initial scale factor
c
      call kinetic (eksum,ekin)
      ekt = gasconst * kelvin
      nc = 5
      ns = 3
      dtc = dt / dble(nc)
      w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
      w(2) = 1.0d0 - 2.0d0*w(1)
      w(3) = w(1)
      df = dble(nfree)
      odnf = 1.0d0 + 3.0d0/df
      gn1kt = (1.0d0+df) * ekt
      dpress = (press-atmsph) / prescon
      scale = 1.0d0
c
c     use multiple time steps to apply thermostat and barostat
c
      do k = 1, nc
         do j = 1, ns
            dts = w(j) * dtc
            dt2 = 0.5d0 * dts
            dt4 = 0.25d0 * dts
            dt8 = 0.125d0 * dts
c
c     update thermostat and barostat velocities and forces
c
            gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
            vnh(4) = vnh(4) + gnh(4)*dt4
            gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
            expterm = exp(-vnh(4)*dt8)
            vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
            gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
            expterm = exp(-vnh(3)*dt8)
            vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
            gnh(1) = (2.0d0*eksum+qbar*vbar*vbar-gn1kt) / qnh(1)
            expterm = exp(-vnh(2)*dt8)
            vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
            gbar = (2.0d0*eksum*odnf+3.0d0*volbox*dpress) / qbar
            expterm = exp(-vnh(1)*dt8)
            vbar = expterm * (vbar*expterm+gbar*dt4)
c
c     find velocity scale factor and update kinetic energy
c
            expterm = exp(-(vnh(1)+vbar*odnf)*dt2)
            scale = scale * expterm
            eksum = eksum * expterm * expterm
c
c     update barostat and thermostat velocities and forces
c
            gbar = (2.0d0*eksum*odnf+3.0d0*volbox*dpress) / qbar
            expterm = exp(-vnh(1)*dt8)
            vbar = expterm * (vbar*expterm+gbar*dt4)
            gnh(1) = (2.0d0*eksum+qbar*vbar*vbar-gn1kt) / qnh(1)
            expterm = exp(-vnh(2)*dt8)
            vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
            gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
            expterm = exp(-vnh(3)*dt8)
            vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
            gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
            expterm = exp(-vnh(4)*dt8)
            vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
            gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
            vnh(4) = vnh(4) + gnh(4)*dt4
         end do
      end do
c
c     use scale factor to update the atomic velocities
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
         end if
      end do
      return
      end
