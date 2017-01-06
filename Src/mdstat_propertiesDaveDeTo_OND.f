c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdstat  --  compute averages over a trajectory  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdstat" is called at each molecular dynamics time step to
c     form statistics on various average values and fluctuations,
c     and to periodically save the state of the trajectory
c
c
      subroutine mdstat_propertiesDaveDeTo_OND (istep,dt,etot,epot,ekin,
     & temp,pres)
      use sizes
      use atoms
      use bath
      use bound
      use boxes
      use inform
      use inter
      use iounit
      use limits
      use math
      use mdstuf
      use molcul
      use moment
      use units
      use usage
      use warp
      implicit none
      integer istep,modstep
      real*8 dt,temp,pres
      real*8 etot,epot,ekin
      real*8 pico,dens
      real*8 fluctuate,fluctuate2
      real*8 intfluct,intfluct2
      real*8 potfluct,potfluct2
      real*8 kinfluct,kinfluct2
      real*8 tfluct,pfluct,dfluct
      real*8 tfluct2,pfluct2,dfluct2
      real*8 vfluct, vfluct2
      real*8 etot_sum,etot2_sum
      real*8 eint_sum,eint2_sum
      real*8 etot_ave,etot2_ave
      real*8 eint_ave,eint2_ave
      real*8 epot_sum,epot2_sum
      real*8 ekin_sum,ekin2_sum
      real*8 epot_ave,epot2_ave
      real*8 ekin_ave,ekin2_ave
      real*8 temp_sum,temp2_sum
      real*8 temp_ave,temp2_ave
      real*8 pres_sum,pres2_sum
      real*8 pres_ave,pres2_ave
      real*8 dens_sum,dens2_sum
      real*8 dens_ave,dens2_ave
      real*8 vol_sum, vol2_sum
      real*8 vol_ave, vol2_ave
      real*8 h_sum, h2_sum, hv_sum
      real*8 h_ave, h2_ave, hv_ave
      real*8 xdpl_sum, xdpl2_sum
      real*8 ydpl_sum, ydpl2_sum
      real*8 zdpl_sum, zdpl2_sum
      real*8 dpl2_ave
      real*8 heat_cap, h_vap, isothermal_k, alpha_coeff, diel
      save etot_sum,etot2_sum
      save eint_sum,eint2_sum
      save epot_sum,epot2_sum
      save ekin_sum,ekin2_sum
      save temp_sum,temp2_sum
      save pres_sum,pres2_sum
      save dens_sum,dens2_sum
      save vol_sum, vol2_sum
      save h_sum, h2_sum, hv_sum
      save xdpl_sum, xdpl2_sum
      save ydpl_sum, ydpl2_sum
      save zdpl_sum, zdpl2_sum
c
c
c     set number of steps for block averages of properties
c
      modstep = mod(istep,iprint)
c
c     zero out summation variables for new averaging period
c
      if (modstep.eq.1 .or. iprint.eq.1) then
         etot_sum = 0.0d0
         etot2_sum = 0.0d0
         epot_sum = 0.0d0
         epot2_sum = 0.0d0
         ekin_sum = 0.0d0
         ekin2_sum = 0.0d0
         eint_sum = 0.0d0
         eint2_sum = 0.0d0
         temp_sum = 0.0d0
         temp2_sum = 0.0d0
         pres_sum = 0.0d0
         pres2_sum = 0.0d0
         dens_sum = 0.0d0
         dens2_sum = 0.0d0
         vol_sum = 0.0d0
         vol2_sum = 0.0d0
         h_sum = 0.0d0
         h2_sum = 0.0d0
         hv_sum = 0.0d0
         xdpl_sum = 0.d0
         xdpl2_sum = 0.d0
         ydpl_sum = 0.d0
         ydpl2_sum = 0.d0
         zdpl_sum = 0.d0
         zdpl2_sum = 0.d0
      end if
c
c     print energy, temperature and pressure for current step
c
      if (verbose) then
         if (modstep .eq. 1) then
            if (use_bounds .and. integrate.ne.'STOCHASTIC') then
               write (iout,10)
   10          format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                    5x,'E Kinetic',7x,'Temp',7x,'Pres',/)
            else
               write (iout,20)
   20          format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                    5x,'E Kinetic',7x,'Temp',/)
            end if
         end if
         if (use_bounds .and. integrate.ne.'STOCHASTIC') then
            write (iout,30)  istep,etot,epot,ekin,temp,pres
   30       format (i10,3f14.4,2f11.2)
         else
            write (iout,40)  istep,etot,epot,ekin,temp
   40       format (i10,3f14.4,f11.2)
         end if
      end if
c
c     print header for the averages over a group of recent steps
c
      if (modstep .eq. 0) then
         pico = dble(istep) * dt
         write (iout,50)  iprint,istep
   50    format (/,' Average Values for the Last',i6,' Out of',
     &              i9,' Dynamics Steps')
         if (digits .ge. 8) then
            write (iout,60)  pico
   60       format (/,' Simulation Time',5x,f19.8,' Picosecond')
         else if (digits .ge. 6) then
            write (iout,70)  pico
   70       format (/,' Simulation Time',5x,f17.6,' Picosecond')
         else
            write (iout,80)  pico
   80       format (/,' Simulation Time',5x,f15.4,' Picosecond')
         end if
      end if
c
c     compute total energy and fluctuation for recent steps
c
      etot_sum = etot_sum + etot
      etot2_sum = etot2_sum + etot**2
      if (modstep .eq. 0) then
         etot_ave = etot_sum / dble(iprint)
         etot2_ave = etot2_sum / dble(iprint)
         fluctuate2 = etot2_ave - etot_ave**2
         if (fluctuate2 .gt. 0.0d0) then
            fluctuate = sqrt(fluctuate2)
         else
            fluctuate = 0.0d0
         end if
         if (digits .ge. 8) then
            write (iout,90)  etot_ave,fluctuate
   90       format (' Total Energy',8x,f19.8,' Kcal/mole',3x,
     &                 '(+/-',f13.8,')')
         else if (digits .ge. 6) then
            write (iout,100)  etot_ave,fluctuate
  100       format (' Total Energy',8x,f17.6,' Kcal/mole',3x,
     &                 '(+/-',f11.6,')')
         else
            write (iout,110)  etot_ave,fluctuate
  110       format (' Total Energy',8x,f15.4,' Kcal/mole',3x,
     &                 '(+/-',f9.4,')')
         end if
      end if
c
c     compute average potential energy and its fluctuation
c
      epot_sum = epot_sum + epot
      epot2_sum = epot2_sum + epot**2
      if (modstep .eq. 0) then
         epot_ave = epot_sum / dble(iprint)
         epot2_ave = epot2_sum / dble(iprint)
         potfluct2 = epot2_ave - epot_ave**2
         if (potfluct2 .gt. 0.0d0) then
            potfluct = sqrt(potfluct2)
         else
            potfluct = 0.0d0
         end if
         if (digits .ge. 8) then
            write (iout,120)  epot_ave,potfluct
  120       format (' Potential Energy',4x,f19.8,' Kcal/mole',3x,
     &                 '(+/-',f13.8,')')
         else if (digits .ge. 6) then
            write (iout,130)  epot_ave,potfluct
130       format (' Potential Energy',4x,f17.6,' Kcal/mole',3x,
     &                 '(+/-',f11.6,')')
         else
            write (iout,140)  epot_ave,potfluct
  140       format (' Potential Energy',4x,f15.4,' Kcal/mole',3x,
     &                 '(+/-',f9.4,')')
         end if
      end if
c
c     compute average kinetic energy and its fluctuation
c
      ekin_sum = ekin_sum + ekin
      ekin2_sum = ekin2_sum + ekin**2
      if (modstep .eq. 0) then
         ekin_ave = ekin_sum / dble(iprint)
         ekin2_ave = ekin2_sum / dble(iprint)
         kinfluct2 = ekin2_ave - ekin_ave**2
         if (kinfluct2 .gt. 0.0d0) then
            kinfluct = sqrt(kinfluct2)
         else
            kinfluct = 0.0d0
         end if
         if (digits .ge. 8) then
            write (iout,150)  ekin_ave,kinfluct
  150       format (' Kinetic Energy',6x,f19.8,' Kcal/mole',3x,
     &                 '(+/-',f13.8,')')
         else if (digits .ge. 6) then
            write (iout,160)  ekin_ave,kinfluct
  160       format (' Kinetic Energy',6x,f17.6,' Kcal/mole',3x,
     &                 '(+/-',f11.6,')')
         else
            write (iout,170)  ekin_ave,kinfluct
  170       format (' Kinetic Energy',6x,f15.4,' Kcal/mole',3x,
     &                 '(+/-',f9.4,')')
         end if
      end if
c
c     compute average intermolecular energy and its fluctuation
c
      if (nmol.ne.1 .and. nmol.ne.n .and. .not.use_ewald) then
         eint_sum = eint_sum + einter
         eint2_sum = eint2_sum + einter**2
         if (modstep .eq. 0) then
            eint_ave = eint_sum / dble(iprint)
            eint2_ave = eint2_sum / dble(iprint)
            intfluct2 = eint2_ave - eint_ave**2
            if (intfluct2 .gt. 0.0d0) then
               intfluct = sqrt(intfluct2)
            else
               intfluct = 0.0d0
            end if
            if (digits .ge. 8) then
               write (iout,180)  eint_ave,intfluct
  180          format (' Intermolecular',6x,f19.8,' Kcal/mole',3x,
     &                    '(+/-',f13.8,')')
            else if (digits .ge. 6) then
               write (iout,190)  eint_ave,intfluct
  190          format (' Intermolecular',6x,f17.6,' Kcal/mole',3x,
     &                    '(+/-',f11.6,')')
            else
               write (iout,200)  eint_ave,intfluct
  200          format (' Intermolecular',6x,f15.4,' Kcal/mole',3x,
     &                    '(+/-',f9.4,')')
            end if
         end if
      end if
c
c     compute the average temperature and its fluctuation
c
      temp_sum = temp_sum + temp
      temp2_sum = temp2_sum + temp**2
      if (modstep .eq. 0) then
         temp_ave = temp_sum / dble(iprint)
         temp2_ave = temp2_sum / dble(iprint)
         tfluct2 = temp2_ave - temp_ave**2
         if (tfluct2 .gt. 0.0d0) then
            tfluct = sqrt(tfluct2)
         else
            tfluct = 0.0d0
         end if
         if (digits .ge. 8) then
            write (iout,210)  temp_ave,tfluct
  210       format (' Temperature',9x,f19.6,' Kelvin',6x,
     &                 '(+/-',f13.6,')')
         else if (digits .ge. 6) then
            write (iout,220)  temp_ave,tfluct
  220       format (' Temperature',9x,f17.4,' Kelvin',6x,
     &                 '(+/-',f11.4,')')
         else
            write (iout,230)  temp_ave,tfluct
  230       format (' Temperature',9x,f15.2,' Kelvin',6x,
     &                 '(+/-',f9.2,')')
         end if
      end if
c
c     compute the average pressure and its fluctuation
c
      if (use_bounds) then
         pres_sum = pres_sum + pres
         pres2_sum = pres2_sum + pres**2
         if (modstep .eq. 0) then
            pres_ave = pres_sum / dble(iprint)
            pres2_ave = pres2_sum / dble(iprint)
            pfluct2 = pres2_ave - pres_ave**2
            if (pfluct2 .gt. 0.0d0) then
               pfluct = sqrt(pfluct2)
            else
               pfluct = 0.0d0
            end if
            if (digits .ge. 8) then
               write (iout,240)  pres_ave,pfluct
  240          format (' Pressure',12x,f19.6,' Atmosphere',2x,
     &                    '(+/-',f13.6,')')
            else if (digits .ge. 6) then
               write (iout,250)  pres_ave,pfluct
  250          format (' Pressure',12x,f17.4,' Atmosphere',2x,
     &                    '(+/-',f11.4,')')
            else
               write (iout,260)  pres_ave,pfluct
  260          format (' Pressure',12x,f15.2,' Atmosphere',2x,
     &                    '(+/-',f9.2,')')
            end if
         end if
c
c     compute the average density and its fluctuation
c
         dens = (1.0d24/volbox) * (totmass/avogadro)
         dens_sum = dens_sum + dens
         dens2_sum = dens2_sum + dens**2
         if (modstep .eq. 0) then
            dens_ave = dens_sum / dble(iprint)
            dens2_ave = dens2_sum / dble(iprint)
            dfluct2 = dens2_ave - dens_ave**2
            if (dfluct2 .gt. 0.0d0) then
               dfluct = sqrt(dfluct2)
            else
               dfluct = 0.0d0
            end if
            if (digits .ge. 8) then
               write (iout,270)  dens_ave,dfluct
  270          format (' Density',13x,f19.8,' Grams/cc',4x,
     &                    '(+/-',f13.8,')')
            else if (digits .ge. 6) then
               write (iout,280)  dens_ave,dfluct
  280          format (' Density',13x,f17.6,' Grams/cc',4x,
     &                    '(+/-',f11.6,')')
            else
               write (iout,290)  dens_ave,dfluct
  290          format (' Density',13x,f15.4,' Grams/cc',4x,
     &                    '(+/-',f9.4,')')
            end if
         end if
      end if
c
c     compute average volume and its fluctuation
c

      vol_sum = vol_sum + volbox
      vol2_sum = vol2_sum + volbox**2
      if (modstep .eq.0) then
         vol_ave = vol_sum / dble(iprint)
         vol2_ave = vol2_sum / dble(iprint)
         vfluct2 = vol2_ave - vol_ave**2
         if (vfluct2 > 0.0d0) then
            vfluct = sqrt(vfluct2)
         else
            vfluct = 0.0d0
         end if
         if (digits .ge. 8) then
             write (iout,330)  vol_ave,vfluct
  330        format (' Volume',14x,f19.8,' Angs**3',5x,
     &                  '(+/-',f13.8,')')
         else if (digits .ge. 6) then
             write (iout,340)  vol_ave,vfluct
  340        format (' Volume',14x,f17.6,' Angs**3',5x,
     &                  '(+/-',f11.6,')')
         else
             write (iout,350)  vol_ave,vfluct
  350        format (' Volume',14x,f15.4,' Angs**3',5x,
     &                  '(+/-',f9.4,')')
         end if
      end if

c
c     Average enthalpy and enthalpy-volume for themodynamics
c

      h_sum = h_sum + epot + atmsph*volbox/prescon
      h2_sum = h2_sum + (epot + atmsph*volbox/prescon)**2
      hv_sum = hv_sum + (epot + atmsph*volbox/prescon)*volbox

      if (modstep .eq.0) then
         h_ave = h_sum / dble(iprint)
         h2_ave = h2_sum / dble(iprint)
         hv_ave = hv_sum / dble(iprint)

C ADDED MINUS SIGN.  SEE EQ. 19 IN TIP4PEW PAPER AND ALSO TIP4PEW CODE

         h_vap = -epot_ave/dble(nmol) + gasconst*kelvin
     &    - atmsph*vol_ave/prescon
         
C CHANGED 'kelvin' TO AVT TEMP. NEED 'gasconst', not 'boltzmann', I THINK.

         isothermal_k = (vol2_ave - vol_ave**2)
     &                 /(gasconst*temp_ave*vol_ave*prescon)

C  CHANGED KELVIN TO AVG TEMP. ALSO, UNITS OF KCAL/MOL FOR ENTHALPY, SO NEED 'gasconst'

         heat_cap = (h2_ave - h_ave**2) 
     &             /(gasconst*temp_ave**2)

C  I THINK IN THE FOLLOWING FORMULA, WE NEED 'gasconst' INSTEAD OF 'boltzmann',
C  SINCE ENTHALPY IS IN KCAL/MOLE.  ALSO, I THINK WE REALLY WANT 'temp_ave' instead of 
C  'kelvin', PER TIP4PEW PAPER. 

         alpha_coeff = (hv_ave - vol_ave*h_ave)
     &                /(gasconst*temp_ave**2*vol_ave)

         write(iout,360) h_vap
  360    format(' Enthalpy of Vaporization', 5x, f19.8, ' <Kcal/mol>') 


         write(iout,370) isothermal_k
  370    format(' Isothermal Compress.', 9x, f19.8, ' <Atm^-1>') 


         write(iout,380) heat_cap
  380    format(' Heat Capacity', 16x, f19.8, ' <Kcal/(mol*K)>') 


         write(iout,390) alpha_coeff
  390    format(' Thermal Expansioon Coeff', 5x, f19.8, ' <Kelvin^-1>') 
      end if


c
c     Dielectric constant
c
      
c      call moments_onlyDipole

c      xdpl_sum = xdpl_sum + xdpl
c      xdpl2_sum = xdpl2_sum + xdpl**2
c      ydpl_sum = ydpl_sum + ydpl
c      ydpl2_sum = ydpl2_sum + ydpl**2
c      zdpl_sum = zdpl_sum + zdpl
c      zdpl2_sum = zdpl2_sum + zdpl**2
      
c      if (modstep .eq.0) then
c
c         dpl2_ave = xdpl2_sum/dble(iprint) - xdpl_sum**2/dble(iprint)**2
c     &    +  ydpl2_sum/dble(iprint) - ydpl_sum**2/dble(iprint)**2
c     &    +  zdpl2_sum/dble(iprint) - zdpl_sum**2/dble(iprint)**2
c
c         diel = 1 + 4*pi*dpl2_ave
c     &    /(3*vol_ave*boltzmann*kelvin)
c         
c         write(iout, 400) diel
c  400    format(' Dielectric Constant', 6x, f19.8, '<e0>')
c
c      end if


c
c     declare deformation value for potential energy smoothing
c
      if (use_smooth) then
         if (modstep .eq. 0) then
            if (digits .ge. 8) then
               write (iout,300)  deform
  300          format (' Deformation',9x,f19.8,' Sqr Angs')
            else if (digits .ge. 6) then
               write (iout,310)  deform
  310          format (' Deformation',9x,f17.6,' Sqr Angs')
            else
               write (iout,320)  deform
  320          format (' Deformation',9x,f15.4,' Sqr Angs')
            end if
         end if
      end if
      return
      end
