c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular or stochastic dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program runmoments 
      use sizes
      use atoms
      use bath
      use bndstr
      use bound
      use inform
      use iounit
      use keys
      use mdstuf
      use potent
      use solute
      use stodyn
      use usage
      use rigid
      implicit none
      integer i,istep,nstep
      integer mode,next
      real*8 dt,dtdump
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set up the structure and molecular mechanics calculation
c
      print*,"In runmoments before initial getxyz mech"
      call initial
      call getxyz
      call mechanic

      if (use_bounds .and. .not.use_rigid)  call bounds

      call chkpole
      call rotpole

             
      print*,"In runmoments"
c      call getuindfromfile
c
c     initialize the temperature, pressure and coupling baths
c

c
c     check for keywords containing any altered parameters
c
c
c     initialize the simulation length as number of time steps
c
c     get the length of the dynamics time step in picoseconds
c
c
c     enforce bounds on thermostat and barostat coupling times
c
c
c     set the time between trajectory snapshot coordinate dumps
c
c
c     get choice of statistical ensemble for periodic system
c
c
c     use constant energy or temperature for nonperiodic system
c
c
c     initialize any holonomic constraints and setup dynamics
c
c
c     print out a header line for the dynamics computation
c
c
c     integrate equations of motion to take a time step
c
       call moments_onlyDipole_fixedcharge
c
c     perform any final tasks before program exit
c
      call final
      end
