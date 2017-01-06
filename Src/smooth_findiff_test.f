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
      program smooth_findiff_test 
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
      use mpidat
      use moldyn
      use virial
      use deriv3b
      use molcul
      use neigh3b
      use neigh
      use mpole
      use energi
      use polar
      use polgrp
      use couple
      use group
      use files
      use rgddyn
      use atomid
      use units
      use deriv
      implicit none
      include 'mpif.h'
      integer i,istep,nstep,j
      integer mode,next,prov,ierr
      real*8 dt,dtdump,energy,dt_2
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
      integer offset,remainder,start,moli1rmndr
      real*8, allocatable :: derivs(:,:)
      character*6 switchmode
      integer stat(MPI_STATUS_SIZE)
      integer maxn13,maxn14,maxn15
      integer maxp11,maxp12,maxp13,maxp14
      real*8 press
      real*8 part1,part2
      real*8 dt_x,factor
      character*120 dynfile
      integer idyn,freeunit
      real*8 maxwell,speed
      real*8 vec(3),twall,tcpu
      integer s1,s2,s3,s4,s5,totsize1
      integer emreal_comm_size,orig_group,taskidem(4),new_group
      integer, allocatable :: buf1(:)
      real*8 efindif(3),analyticgx(3)
      real*8 efindif_ep3b(3),analyticg_dep3bx(3)
      real*8 analyticgy(3)
      real*8 analyticg_dep3by(3)
      real*8 analyticgz(3)
      real*8 analyticg_dep3bz(3)
      
      real*8 delx


c
c
c     set up the structure and molecular mechanics calculation
c

      call mpi_init_thread(mpi_thread_funneled,prov,ierr)
      call mpi_comm_rank(mpi_comm_world,taskid,ierr)
      call mpi_comm_size(mpi_comm_world,numtasks,ierr)

      master =0
      master1 =1
      if(taskid.eq.master) then
         call settime2(twall,tcpu)
      end if

      call initial
      call getxyz

      if(taskid.eq.master) then
         call mechanic_pme3
c      else if (taskid.eq.master1) then
c         if (.not.allocated(nelst))  allocate (nelst(n))
c         if (.not.allocated(elst))  allocate (elst(maxelst,n))      
      else if (taskid.ne.master) then
        if (.not. allocated(imol))  allocate (imol(2,n))
        if (.not. allocated(kmol))  allocate (kmol(n))
        if (.not. allocated(molmass))  allocate (molmass(n))

      maxn13 = 3 * maxval
      maxn14 = 9 * maxval
      maxn15 = 27 * maxval
      maxp11 = 150
      maxp12 = 50
      maxp13 = 50
      maxp14 = 50
      if (.not. allocated(n13))  allocate (n13(n))
      if (.not. allocated(n14))  allocate (n14(n))
      if (.not. allocated(n15))  allocate (n15(n))
      if (.not. allocated(i13))  allocate (i13(maxn13,n))
      if (.not. allocated(i14))  allocate (i14(maxn14,n))
      if (.not. allocated(i15))  allocate (i15(maxn15,n))
      if (.not. allocated(iuse))  allocate (iuse(n))
      if (.not. allocated(use))  allocate (use(0:n))
      if (.not. allocated(kgrp))  allocate (kgrp(n))
      if (.not. allocated(grplist))  allocate (grplist(n))
      if (.not. allocated(igrp))  allocate (igrp(2,0:maxgrp))
      if (.not. allocated(grpmass))  allocate (grpmass(0:maxgrp))
      if (.not. allocated(wgrp))  allocate (wgrp(0:maxgrp,0:maxgrp))

      if (.not. allocated(ipole))  allocate (ipole(n))
c      if (.not. allocated(polsiz))  allocate (polsiz(n))
c      if (.not. allocated(pollist))  allocate (pollist(n))
      if (.not. allocated(zaxis))  allocate (zaxis(n))
      if (.not. allocated(xaxis))  allocate (xaxis(n))
      if (.not. allocated(yaxis))  allocate (yaxis(n))
c      if (.not. allocated(pole))  allocate (pole(maxpole,n))
      if (.not. allocated(rpole))  allocate (rpole(maxpole,n))
      if (.not. allocated(polaxe))  allocate (polaxe(n))
      if (.not. allocated(np11))  allocate (np11(n))
      if (.not. allocated(np12))  allocate (np12(n))
      if (.not. allocated(np13))  allocate (np13(n))
      if (.not. allocated(np14))  allocate (np14(n))

      if (.not. allocated(polarity))  allocate (polarity(n))
      if (.not. allocated(thole))  allocate (thole(n))
      if (.not. allocated(pdamp))  allocate (pdamp(n))
      if (.not. allocated(ip11))  allocate (ip11(maxp11,n))
      if (.not. allocated(ip12))  allocate (ip12(maxp12,n))
      if (.not. allocated(ip13))  allocate (ip13(maxp13,n))
      if (.not. allocated(ip14))  allocate (ip14(maxp14,n))
      end if

      call bcast_mechanic
      call bcast_ewald

      if(taskid.ne.master) then
         call mechanic_parallel2
      end if
      if(taskid.eq.master1) then
c         if (.not.allocated(nelst))  allocate (nelst(n))
c         if (.not.allocated(elst))  allocate (elst(maxelst,n)) 
c         print*,"Done w/ mpole list alloc on master 1"
          call cutoffs_ewald
      end if
       
      offset=int(nmol/numtasks)
      remainder=mod(nmol,numtasks)
c      offset=int(nmol/(numtasks-1))
c      remainder=mod(nmol,(numtasks-1))

      if(taskid.le.remainder-1) then
c      if(taskid.le.remainder.and.taskid.ne.master) then
        moli1rmndr=numtasks*offset+taskid+1
c        moli1rmndr=(numtasks-1)*offset+taskid
      else if (taskid.gt.remainder-1) then
        moli1rmndr=0
      end if
      start=taskid*offset+1
c      start=(taskid-1)*offset+1


c
c     initialize the temperature, pressure and coupling baths
c
      if(taskid.eq.master) then
c        print*,"Before call to gettime2"
c        call gettime2 (twall,tcpu)
c        print*,"Wall Time=", twall, "CPU Time=",tcpu
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BEEMAN'

        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           string = record(next:120)
           if (keyword(1:11) .eq. 'INTEGRATOR ') then
             call getword (record,integrate,next)
             call upcase (integrate)
           end if
        end do


c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if

c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do

      dt = 0.001d0 * dt
      dt_2 = 0.5d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c

      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,160)
  160          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if

c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if

c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit_pt1_master
      dynfile = filename(1:leng)//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
      end if
      end if

      call mpi_bcast(exist,1,mpi_logical,master,
     & mpi_comm_world,ierr)

      call mpi_bcast(integrate,11,mpi_character,master,
     & mpi_comm_world,ierr)

      call mpi_bcast(nstep,1,mpi_integer,master,
     & mpi_comm_world,ierr)

      if(exist) then
         call bcast_lattice
         call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
         call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
         call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
      end if


      if(exist.eqv..false.) then
         if (integrate .eq. 'RIGIDBODY') then
            if(taskid.eq.master) then
               do i = 1, ngrp
                 speed = maxwell (grpmass(i),kelvin)
                 call ranvec (vec)
                 do j = 1, 3
                  vcm(j,i) = speed * vec(j)
                  wcm(j,i) = 0.0d0
                  lm(j,i) = 0.0d0
                 end do
               end do
               if (nuse .eq. n)  call mdrest (0)
            end if
         else if (integrate .eq. 'RESPA') then
            call gradient_setup
            if(taskid.eq.master) then
              allocate (derivs(3,n))
              call prep_pole
              call gradslow_serial(energy,derivs)
            end if 
            call mpi_bcast(rpole,13*n,mpi_real8,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(listbcast,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
            if(listbcast) then
            call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
            call mpi_bcast(nmollst2,nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
            call mpi_bcast(mollst2,max2blst*nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
            end if
            listbcast=.false.
            
            switchmode = 'MPOLE'
            call switch (switchmode)
            call mpi_barrier(mpi_comm_world,ierr)
            call smooth2b_gradient_polar(start,start+offset-1,
     &        moli1rmndr)
            if(taskid.eq.master) then
              do i=1,3
               do j=1,3
c                 vir(i,j)=vir(i,j)+virep3b2(i,j)+virep3b3(i,j)+
c     &           viremreal(i,j)+virev(i,j)
                 vir(i,j)=vir(i,j)+virep3b(i,j)
               end do
              end do
              energy = energy +ep3b
              do i = 1, n
                do j = 1, 3
                 derivs(j,i) = derivs(j,i) + dep3b(j,i)
                end do
              end do

              do i = 1, n
                if (use(i) .and. mass(i).ne.0.0d0) then
                 speed = maxwell (mass(i),kelvin)
                 call ranvec (vec)
                 do j = 1, 3
                   v(j,i) = speed * vec(j)
                   a(j,i) = -convert * derivs(j,i) / mass(i)
                 end do
                else
                 do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                 end do
                end if
              end do
         
              call gradfast (energy,derivs)
              do i = 1, n
                if (use(i) .and. mass(i).ne.0.0d0) then
                  do j = 1, 3
                   aalt(j,i) = -convert * derivs(j,i) / mass(i)
                  end do
                else
                  do j = 1, 3
                   aalt(j,i) = 0.0d0
                  end do
                end if
              end do
              deallocate (derivs)
              if (nuse .eq. n)  call mdrest (0)
            end if
         else
                       
            call gradient_setup
            call allocPermElec1 
            if(taskid.eq.master) then
              print*,"In verlet/nose/beeman setup"
              print*,"Before alloc3b_nblist verlet/nose/beeman setup"
c              call alloc3b_nblist
              call alloc3b_vdwlist
c              call allocPermElec
              allocate (derivs(3,n))
              call prep_pole
              print*,"listbcast=",listbcast
            else if(taskid .eq. master1) then
c              call allocPermElec1
c               call allocPermElec1_mlist
               call mlist
             print*,"In verlet/nose/beeman setup, after allocPermElec1"
            end if

            call mpi_barrier(mpi_comm_world,ierr)
            
            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(rpole,13*n,mpi_real8,master,
     &       mpi_comm_world,ierr)

            if(listbcast) then
            call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(nmollst2,nmol,mpi_integer,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(mollst2,max2blst2*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
c            print*,"After all 2body list bcast"
            end if
            listbcast = .false.

            if(taskid.eq.master) then
             call gradient_covalent_vdw 
             call gradient_emrecip
            else if(taskid.eq.master1) then
             call gradient_emreal_bcast
            end if
            
            call smooth2_gradient_polar(start,start+offset-1,moli1rmndr)

            call mpi_barrier(mpi_comm_world,ierr)
        call mpi_bcast(emreal_tmp,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
        call mpi_bcast(demreal_tmp,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
        call mpi_bcast(viremreal_tmp,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)

            listsend_mpole=.false.
 
              if(taskid.eq.master) then
c                 call mpi_recv(emreal,1,mpi_real8,master1,4,
c     &           mpi_comm_world,stat,ierr)
c                 print*,"emreal after first receive",emreal
c                 call mpi_recv(demreal,3*n,mpi_real8,master1,5,
c     &           mpi_comm_world,stat,ierr)
c                 call mpi_recv(viremreal,3*3,mpi_real8,master1,6,
c     &           mpi_comm_world,stat,ierr)
c                 print*,"viremreal=",viremreal  
                 call empole1d_3b_Perm_selfeng_bcast 

                 do i=1,3
                   do j=1,3
                    vir(i,j)=vir(i,j)+virep3b(i,j)+viremreal_tmp(i,j)
     &                       + viremrecip(i,j)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett 
     &                   + ev + em+ep3b
               print*,"em ev in erlet/nose/beeman setup",em,ev
               print*,"ep3b in verlet/nose/beeman setup",ep3b
               print*,"TotTriplesin verlet/nose/beeman setup",ntriples
                do i = 1, n
                  do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i)
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i)
                  end do
                 if (use(i) .and. mass(i).ne.0.0d0) then
                   speed = maxwell (mass(i),kelvin)
                   call ranvec (vec)
                   do j = 1, 3
                    v(j,i) = speed * vec(j)
                    a(j,i) = -convert * derivs(j,i) / mass(i)
                    aalt(j,i) = a(j,i)
                   end do
                 else
                  do j = 1, 3
                   v(j,i) = 0.0d0
                   a(j,i) = 0.0d0
                   aalt(j,i) = 0.0d0
                  end do
                 end if
                end do
                 deallocate (derivs)
                 if (nuse .eq. n)  call mdrest (0)
              end if 
         end if
      end if

c      if(taskid.eq.master) then
c        call mdinit_pt2_master
c      end if
c      call mdinit
c

c     print out a header line for the dynamics computation
c
      if(taskid.eq.master) then      
      if (integrate .eq. 'VERLET') then
         write (iout,330)
  330    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,340)
  340    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,350)
  350    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,360)
  360    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,370)
  370    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,380)
  380    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
      end if
c
c     integrate equations of motion to take a time step
c
        if(taskid.eq.master) then
         delx=0.000001d0
        end if


         if (integrate .eq. 'VERLET') then
          do istep = 1, nstep
            
            !call gradient_setup
            call allocPermElec1

            if(taskid.eq.master) then 
c              call alloc3b_nblist 
              !x(7)=x(7)+delx
              !y(7)=y(7)+delx
              z(6)=z(6)+delx
              call gradient_setup
              call alloc3b_vdwlist
c              call allocPermElec
              allocate (derivs(3,n))              
              call prep_pole
            end if


            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(rpole,13*n,mpi_real8,master,
     &       mpi_comm_world,ierr)

            if(taskid .eq. master1) then
c              call allocPermElec1
c              call allocPermElec1_mlist 
               call mlist
            end if          
            call mpi_barrier(mpi_comm_world,ierr)

            if(listbcast) then
            call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(nmollst2,nmol,mpi_integer,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(mollst2,max2blst2*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            end if
            listbcast = .false.

            if(taskid.eq.master) then
             call gradient_covalent_vdw
             call gradient_emrecip
            else if(taskid.eq.master1) then
             call gradient_emreal_bcast
            end if

             call smooth2_gradient_polar(start,start+offset-1,
     &       moli1rmndr)

            call mpi_barrier(mpi_comm_world,ierr)
        call mpi_bcast(emreal_tmp,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
        call mpi_bcast(demreal_tmp,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
        call mpi_bcast(viremreal_tmp,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)

            if(taskid.eq.master) then
                 call empole1d_3b_Perm_selfeng_bcast

                 do i=1,3
                   do j=1,3
                    vir(i,j)=vir(i,j)+virep3b(i,j)+viremreal_tmp(i,j)
     &                       + viremrecip(i,j)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b
              efindif(istep)=energy
              efindif_ep3b(istep)=ep3b
              
              do i = 1, n
                do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i) 
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i)
                end do
               if (use(i)) then
                do j = 1, 3
                 a(j,i) = -convert * derivs(j,i) / mass(i)
                 v(j,i) = v(j,i) + a(j,i)*dt_2
                end do
               end if
              end do
       print*," istep Vdw, PermElec Eng, Pol Eng",istep,ev,em,ep3b
       print*," istep NumTriples",ntriples

              analyticgx(istep)=derivs(1,6)
              analyticg_dep3bx(istep)=dep3b(1,6) 
              analyticgy(istep)=derivs(2,6)
              analyticg_dep3by(istep)=dep3b(2,6)
              analyticgz(istep)=derivs(3,6)
              analyticg_dep3bz(istep)=dep3b(3,6)

              deallocate (derivs)
c              call verlet_pt2(istep,dt,energy)
              if(istep.eq.3) then
                print*,"Finite Diff TotEng, TotGradz",
     & (efindif(3)-efindif(1))/2.0d0/delx,analyticgz(2)
                print*,"Finite Diff PolEng, PolGradz",
     & (efindif_ep3b(3)-efindif_ep3b(1))/2.0d0/delx,analyticg_dep3bz(2)
              end if 
            end if   
          end do 
         end if 

         if (integrate .eq. 'STOCHASTIC') then
          do istep=1,nstep
             call sdstep (istep,dt)
          end do
         end if

         if (integrate .eq. 'BUSSI') then
          do istep=1,nstep
             call bussi (istep,dt)
          end do
         end if

         if (integrate .eq. 'NOSE-HOOVER') then
          do istep=1,nstep
            call gradient_setup
            if(taskid.eq.master) then
               call nose_pt1(istep,dt,press)
               allocate (derivs(3,n))
               call prep_pole
              print*,"listbcast=",listbcast
               call gradient_serial(energy,derivs)
            end if
            call mpi_bcast(x,n,mpi_real8,master,
     &       mpi_comm_world,ierr)
c         print*,"After x mpi_bcast"
            call mpi_bcast(y,n,mpi_real8,master,
     &       mpi_comm_world,ierr)
c         print*,"After y mpi_bcast"
            call mpi_bcast(z,n,mpi_real8,master,
     &       mpi_comm_world,ierr)
         print*,"After x,y,z mpi_bcast in Nose-Hoover"
            call mpi_bcast(rpole,13*n,mpi_real8,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)
c         print*,"listbcast=",listbcast
            if(listbcast) then
            call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(nmollst2,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
            call mpi_bcast(mollst2,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            end if
            listbcast = .false.
            call bcast_lattice
            switchmode = 'MPOLE'
            call switch (switchmode)
            call mpi_barrier(mpi_comm_world,ierr)
            call smooth2b_gradient_polar(start,start+offset-1,
     &       moli1rmndr)
            if(taskid.eq.master) then
              do i=1,3
               do j=1,3
c                 vir(i,j)=vir(i,j)+virep3b2(i,j)+virep3b3(i,j)+
c     &           viremreal(i,j)+virev(i,j)
                 vir(i,j)=vir(i,j)+virep3b(i,j)
               end do
              end do
              energy = energy +ep3b
c              do i = 1, n
c                do j = 1, 3
c                 derivs(j,i) = derivs(j,i) + dep3b(j,i)
c                end do
c              end do
              do i = 1, n
                do j = 1, 3
                 derivs(j,i) = derivs(j,i) + dep3b(j,i)
                end do
               if (use(i)) then
                do j = 1, 3
                 a(j,i) = -convert * derivs(j,i) / mass(i)
                 v(j,i) = v(j,i) + a(j,i)*dt_2
                end do
               end if
              end do
              deallocate (derivs)
              call nose_pt2(istep,dt,press,energy)
            end if
          end do
         end if
         if (integrate .eq. 'GHMC') then
            call ghmcstep (istep,dt)
         end if
         if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep,dt)
         end if
         if (integrate .eq. 'RESPA') then
            call respa (istep,dt)
         end if
         if (integrate .eq. 'BEEMAN') then
          do istep=1,nstep
            call gradient_setup
            if(taskid.eq.master) then
              allocate (derivs(3,n))
              factor = dble(bmnmix)
              dt_x = dt / factor
              part1 = 0.5d0*factor + 1.0d0
              part2 = part1 - 2.0d0
              call beeman_pt1 (istep,dt,dt_x,part1)
              call prep_pole
              call gradient_serial(energy,derivs)
            end if

            call mpi_bcast(x,n,mpi_real8,master,
     &       mpi_comm_world,ierr)
c         print*,"After x mpi_bcast"
            call mpi_bcast(y,n,mpi_real8,master,
     &       mpi_comm_world,ierr)
c         print*,"After y mpi_bcast"
            call mpi_bcast(z,n,mpi_real8,master,
     &       mpi_comm_world,ierr)
c         print*,"After x,y,z mpi_bcast"
            call mpi_bcast(rpole,13*n,mpi_real8,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)

            if(listbcast) then
            call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(nmollst2,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
            call mpi_bcast(mollst2,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            end if

            listbcast = .false.
            switchmode = 'MPOLE'
            call switch (switchmode)
            call mpi_barrier(mpi_comm_world,ierr)
            call smooth2b_gradient_polar(start,start+offset-1,
     &      moli1rmndr)
            if(taskid.eq.master) then
              do i=1,3
               do j=1,3
c                 vir(i,j)=vir(i,j)+virep3b2(i,j)+virep3b3(i,j)+
c     &           viremreal(i,j)+virev(i,j)
                 vir(i,j)=vir(i,j)+virep3b(i,j)
               end do
              end do
              energy = energy +ep3b
              do i = 1, n
                do j = 1, 3
                 derivs(j,i) = derivs(j,i) + dep3b(j,i)
                end do
                if (use(i)) then
                 do j = 1, 3
                  aalt(j,i) = a(j,i)
                  a(j,i) = -convert * derivs(j,i) / mass(i)
                  v(j,i) = v(j,i) + (part2*a(j,i)+aalt(j,i))*dt_x
                 end do
                end if
              end do 
              deallocate (derivs)
              call beeman_pt2 (istep,dt,energy)
            end if
          end do
         end if
c
c     perform any final tasks before program exit
c
      if(taskid.eq.master) then
        print*,"Before call to gettime2"
        call gettime2 (twall,tcpu)
        print*,"Wall Time=", twall, "CPU Time=",tcpu
      end if 
      call final
      call mpi_finalize(ierr)
      end
