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
      program dEtensor_test_nolistsend_separatereduce_funcdecomp 
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
      use boxes
      use totfield
      use cell
      use limits
      use pme
      use chunks
      use aprx
      use dEtensor
      use vdw
      use vdwpot
      use mutant
      implicit none
      include 'mpif.h'
      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time
      integer i,istep,nstep,j
      integer mode,next,prov,ierr
      real*8 dt,dtdump,energy,dt_2
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
      integer offset,remainder,start,moli1rmndr
      real*8, allocatable :: derivs(:,:)
      character*7 switchmode
      integer stat(MPI_STATUS_SIZE)
      integer stats2(MPI_STATUS_SIZE,6)
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
      real*8 efindif(3),analyticg(3)
      real*8 efindif_ep3b(3),analyticg_dep3b(3)
      real*8 delx,delta,vir_analytic(3,3),virep3b_analytic(3,3)
      real*8 volbox_analytic,dedv_vir,dedv_fd,dedv_vir_ep3b,dedv_fd_ep3b
      real*8 analyticgx(3)
      real*8 analyticg_dep3bx(3)
      real*8 analyticgy(3)
      real*8 analyticg_dep3by(3)
      real*8 analyticgz(3)
      real*8 analyticg_dep3bz(3)
      integer atomind,start_emreal
      integer emreal_group,k
      integer, allocatable :: taskidother(:)
      integer, allocatable :: taskidemreal(:)
      logical firstload_emreal,firstload_polz
      logical firstload_emreal_bcast,firstload_polz_bcast
      logical firstload_vdw
      integer offset2,remainder2,start2,last2,last
      real*8 elrc,vlrc
      save firstload_emreal
      data firstload_emreal  / .true. /
      save firstload_polz
      data firstload_polz  / .true. /
      save firstload_emreal_bcast
      data firstload_emreal_bcast  / .true. /
      save firstload_polz_bcast
      data firstload_polz_bcast  / .true. /
      save firstload_vdw
      data firstload_vdw / .true. /
      save press

c
c
c     set up the structure and molecular mechanics calculation
c
c      print*,"Before MPI initialization"
      call mpi_init_thread(mpi_thread_funneled,prov,ierr)
      call mpi_comm_rank(mpi_comm_world,taskid,ierr)
      call mpi_comm_size(mpi_comm_world,numtasks,ierr)
c      print*,"After MPI initialization"
      master =0
      master1 =1
      master2=2
      if(taskid.eq.master) then
         call settime2(twall,tcpu)
      end if
C BROADCAST SOME LOGICALS.
c      call mpi_bcast(firstload_emreal,1,mpi_logical,master,
c     & mpi_comm_world,ierr)
c
c      call mpi_bcast(firstload_vdw,1,mpi_logical,master,
c     & mpi_comm_world,ierr)
c
c      call mpi_bcast(firstload_emreal_bcast,1,mpi_logical,master,
c     & mpi_comm_world,ierr)
c
c      call mpi_bcast(firstload_polz_bcast,1,mpi_logical,master,
c     & mpi_comm_world,ierr)
c      call mpi_bcast(firstload_polz,1,mpi_logical,master,
c     & mpi_comm_world,ierr)
      call initial
      call getxyz

      if(taskid.eq.master) then
         call mechanic_pme4_eastman
         !call mechanic_pme4
         !call kewald 
         call kewald_eastman
         !call kewaldreg3b
         !call get_numtasks_emreal
         numtasks_emreal=int((numtasks-2)/2)
         numtasks_vdw=numtasks_emreal
         numtasks_polz=numtasks_vdw 
         numtasks_vdw2=numtasks_vdw
         numtasks_emreal2=numtasks_emreal

c        if(.not.allocated(dEd1)) allocate (dEd1(9,n*maxelst))
c        if(.not.allocated(dEd2)) allocate (dEd2(9,n*maxelst))
c        if(.not.allocated(dEp1)) allocate (dEp1(9,n*maxelst))
c        if(.not.allocated(dEp2)) allocate (dEp2(9,n*maxelst))
c        if(.not.allocated(dEindex)) allocate (dEindex(2,n*maxelst))

      else if (taskid.ne.master) then
        if (.not. allocated(imol))  allocate (imol(2,n))
        if (.not. allocated(kmol))  allocate (kmol(n))
        if (.not. allocated(molmass))  allocate (molmass(n))

c        if(.not.allocated(dEd1)) allocate (dEd1(9,n*maxelst))
c        if(.not.allocated(dEd2)) allocate (dEd2(9,n*maxelst))
c        if(.not.allocated(dEp1)) allocate (dEp1(9,n*maxelst))
c        if(.not.allocated(dEp2)) allocate (dEp2(9,n*maxelst))
c        if(.not.allocated(dEindex)) allocate (dEindex(2,n*maxelst))

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
      if (.not. allocated(pole))  allocate (pole(maxpole,n))
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

      if (.not. allocated(ivdw)) allocate (ivdw(n))
      if (.not. allocated(jvdw)) allocate (jvdw(n))
      if (.not. allocated(ired)) allocate (ired(n))
      if (.not. allocated(kred)) allocate (kred(n))
      if (.not. allocated(radmin)) allocate (radmin(maxclass,maxclass))
      if (.not. allocated(epsilon))
     &       allocate (epsilon(maxclass,maxclass))
      if (.not. allocated(radmin4))
     &       allocate (radmin4(maxclass,maxclass))
      if (.not. allocated(epsilon4))
     &       allocate (epsilon4(maxclass,maxclass))
      if (.not. allocated(mut)) allocate(mut(n))

      end if


      call bcast_mechanic
      call bcast_ewald
      call bcast_ewald3b
      call bcast_vdw
      !if(.not.longrangepoldir) then
      !  call kewald3b_totfield
      !end if

      !call kewald3b
      !call kewald3b_totfield
c      if(approxmode.eq.'1BODYMODE') then
       !print*,"Doing uind allocs"
       if(.not. allocated(uind)) allocate(uind(3,n))
       if(.not. allocated(uinp)) allocate(uinp(3,n))
         if(.not.allocated(dep3b_recip)) allocate (dep3b_recip(3,n))

C  ALLOCS NECESSARY FOR DEBUGGING MODE ONLY!! REMOVE LATER!!
       if(.not.allocated(demreal_tmp)) allocate(demreal_tmp(3,n))
         if (.not. allocated(fieldnpole))  allocate (fieldnpole(3,n))
         if (.not. allocated(fieldpnpole))  allocate (fieldpnpole(3,n))
         if(.not.allocated(dep3b)) allocate (dep3b(3,n))
         if(.not.allocated(dev)) allocate (dev(3,n))
         if(.not.allocated(dep3bmut)) allocate (dep3bmut(3,n))
        if(.not.allocated(dep3b1)) allocate (dep3b1(3,n))

C DEFINE A NEW COMM FOR THE MPI/OMP REAL-SPACE PERM ELECTROSTATICS
C USING ONLY 'numtasks_emreal' TASKS    
      call mpi_bcast(numtasks_emreal,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(numtasks_vdw,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(numtasks_polz,1,mpi_integer,master,
     &   mpi_comm_world,ierr)

      call mpi_bcast(numtasks_emreal2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(numtasks_vdw2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)

      !print*,"numtasks_emreal after read and bcast",numtasks_emreal
      !allocate (taskidemreal(0:numtasks_emreal-1))
      !allocate (taskidother(0:(numtasks-1-numtasks_emreal)))

      !do j=0,numtasks_emreal-1
      !  taskidemreal(j)=j
      !end do

      !do j=numtasks_emreal,numtasks-1
      !   k=j-numtasks_emreal
      !  taskidother(k)=j
      !end do

      !call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)

      !if(taskid.lt.numtasks_emreal) then
      !  call mpi_group_incl(orig_group,numtasks_emreal,taskidemreal,
     &!        emreal_group,ierr)
         if (.not. allocated(start_emreal2))
     &          allocate (start_emreal2(0:numtasks_emreal-1))
         if (.not. allocated(last_emreal2))
     &          allocate (last_emreal2(0:numtasks_emreal-1))
         if (.not. allocated(maxsize_elst))
     &         allocate (maxsize_elst(0:numtasks_emreal-1))
        if (.not. allocated(start_vdw2))
     &          allocate (start_vdw2(0:numtasks_vdw-1))
         if (.not. allocated(last_vdw2))
     &          allocate (last_vdw2(0:numtasks_vdw-1))
         if (.not. allocated(maxsize_vlst))
     &         allocate (maxsize_vlst(0:numtasks_vdw-1))

      !else
      !  call mpi_group_incl(orig_group,numtasks-numtasks_emreal,
     &!       taskidother,emreal_group,ierr)
      !end if

        !call MPI_COMM_CREATE(MPI_COMM_WORLD, emreal_group,
     &  !          emreal_comm,ierr)

      if(taskid.ne.master) then
         call mechanic_parallel2
c         call kewald
      end if
      if(taskid.eq.master1) then
         !call kewald
         call kewald_eastman
         !call cutoffs_ewald
         if (.not. allocated(fieldnpole))  allocate (fieldnpole(3,n))
         if (.not. allocated(fieldpnpole))  allocate (fieldpnpole(3,n))
      end if

      if(.not.longrangepoldir) then
         if (.not. allocated(fieldnpole))  allocate (fieldnpole(3,n))
         if (.not. allocated(fieldpnpole))  allocate (fieldpnpole(3,n))
      end if

      if((taskid.gt.1).and.(taskid.lt.numtasks_emreal+2)) then
        call cutoffs_ewald
      else if(taskid.ge.numtasks_emreal+2) then
        call cutoffs_vdw
      end if

      if(approxmode.ne.'1BODYMODE') then
        call cutoffs2_serial
      end if

      call mpi_barrier(mpi_comm_world,ierr)
        ! ewtotfieldsmallsmooth2_gradient_polar_1body(
     &  !     start,start+offset-1,moli1rmndr)

      if(itermode.ne.1) then
         offset=int(nmol/numtasks)
         remainder=mod(nmol,numtasks)

         remainder2=mod(offset,itermode)

         if(taskid.le.remainder-1) then
           moli1rmndr=numtasks*offset+taskid+1
         else if (taskid.gt.remainder-1) then
           moli1rmndr=0
         end if
         start=taskid*offset+1
         offset2=offset-remainder2
         last=start+offset2-1
         start2=last+1
         last2=start2+remainder2-1
      else 
         offset=int(nmol/numtasks)
         remainder=mod(nmol,numtasks)

         if(taskid.le.remainder-1) then
           moli1rmndr=numtasks*offset+taskid+1
         else if (taskid.gt.remainder-1) then
           moli1rmndr=0
         end if
         start=taskid*offset+1
      end if

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
            call smooth_gradient_polar(start,start+offset-1,moli1rmndr)
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

            if(taskid.eq.master) then
              call gradient_setup
c              allocate (derivs(3,n))
c              call prep_pole
c              print*,"listbcast=",listbcast
              call allocPermElec1
              call alloc3b_vdwlist_nolist
c              call alloc_erecip
c            else if(taskid.eq.master1) then
c              call alloc_erecip
            end if

c            call mpi_barrier(mpi_comm_world,ierr)

            call alloc_erecip
            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
c            call mpi_bcast(rpole,13*n,mpi_real8,master,
c     &       mpi_comm_world,ierr)
            call prep_pole

           if((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2)))then
c              print*,"In verlet/nose/beeman setup"
c              print*,"Before alloc3b_nblist verlet/nose/beeman setup"
               call vlist_par
              if(firstload_vdw.and.(taskid.eq.numtasks_emreal+2)) then
c                firstload_vdw=.false.
                call vdw_load_balance_sizelist_half
              end if
c             print*,"In verlet/nose/beeman setup,after vdw list shit"
           else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2))then
               call mlist_par
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
                if(firstload_emreal.and.(taskid.eq.master2)) then
c                 firstload_emreal = .false.
                  if(longrangepoldir) then
                   call emreal_load_balance_sizelist
                  else
                   call emreal_load_balance_sizelist_mol
                  end if
                end if
ccc             print*,"In verlet/nose/beeman setup,after emreal list shit"
ccc            else if(taskid.eq.master2
ccc     &                .and.(approxmode.ne.'1BODYMODE')) then
ccc                call mollist2bodyOO6_8_par
ccc              if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
ccc                firstload_polz = .false.
ccc                call ewpolar_load_balance
ccc              end if
           end if

           ! if((approxmode.ne.'1BODYMODE').and.
     &     !           (taskid.ge.numtasks_emreal+2)) then
            if((approxmode.ne.'1BODYMODE').and.
     &                (taskid.eq.numtasks_emreal+2)) then
              if(approxmode.eq.'2BODYMODE') then
              call mollist2bodyOO6_8_par_1list_par_newesttink7
              else
              call mollist2bodyOO6_8_par_1list
              end if
            !  if(firstload_polz.and.(taskid.eq.numtasks_emreal+2)) then
              if(firstload_polz) then
                !call ewpolar_load_balance_half
                call polar_load_balance7_half
              end if
            end if

            call mpi_barrier(mpi_comm_world,ierr)

            call mpi_bcast(listbcast,1,mpi_logical,numtasks_emreal+2,
     &       mpi_comm_world,ierr)
              
           ! if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
            if(listbcast.and.(approxmode.ne.'1BODYMODE')) then
             firstload_polz = .false.
             !call ewbcast_polar_load_balance_half
             call bcast_polar_load_balance_half
            call mpi_bcast(nmollst,nmol,mpi_integer,numtasks_emreal+2,
     &       mpi_comm_world,ierr)
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,
     &       numtasks_emreal+2,mpi_comm_world,ierr)
            end if
            listbcast = .false.

c              call mpi_bcast(numtasks_emreal2,1,mpi_integer,
c     &             master1,mpi_comm_world,ierr)

            if(firstload_emreal) then
               firstload_emreal = .false.
              call mpi_bcast(numtasks_emreal2,1,mpi_integer,
     &             master2,mpi_comm_world,ierr)
               call mpi_bcast(start_emreal2,numtasks_emreal,mpi_integer,
     &             master2,mpi_comm_world,ierr)
               call mpi_bcast(last_emreal2,numtasks_emreal,mpi_integer,
     &             master2,mpi_comm_world,ierr)
c               call mpi_bcast(maxsize_elst,numtasks_emreal,mpi_integer,
c     &             master1,mpi_comm_world,ierr)
c              call sendmlist_load_balanced3_sizelist
            end if
c            listsend_mpole = .false.

c              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
c     &            numtasks_emreal,mpi_comm_world,ierr)

            if(firstload_vdw) then
                firstload_vdw=.false.
              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
     &            numtasks_emreal+2,mpi_comm_world,ierr)
              call mpi_bcast(start_vdw2,numtasks_vdw,mpi_integer,
     &             numtasks_emreal+2,mpi_comm_world,ierr)
               call mpi_bcast(last_vdw2,numtasks_vdw,mpi_integer,
     &             numtasks_emreal+2,mpi_comm_world,ierr)
c               call mpi_bcast(maxsize_vlst,numtasks_vdw,mpi_integer,
c     &             master,mpi_comm_world,ierr)
c               call sendvlist_load_balanced3_sizelist_half
            end if
c            listsend_vdw = .false.

      call grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2
c            if(approxmode.eq.'1BODYMODE') then
            call mpi_bcast(uind,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
            call mpi_bcast(uinp,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
            call mpi_bcast(ep3b,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
c            else
c            call mpi_bcast(fieldnpole,3*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(fieldpnpole,3*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(fphi_totfield,20*npole,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vxx_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vyx_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vzx_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vyy_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vzy_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vzz_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c
c            call mpi_bcast(thetai1,4*bsorder*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(thetai2,4*bsorder*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(thetai3,4*bsorder*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(qfac,nfft1*nfft2*nfft3,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(pmetable,n*nchunk,mpi_integer,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(igrid,3*n,mpi_integer,master1,
c     &       mpi_comm_world,ierr)
c            end if

         !call ewtotfield_gradient_polar(start,start+offset-1,moli1rmndr)
c            if(approxmode.eq.'1BODYMODE') then
             !call ewtotfieldsmooth2_gradient_polar_1body_dEtensor_mpioff
          call ewtotfieldsmooth2_gradient_polar_1body_dEtensor_simult
c            else
c            call ewtotfieldsmooth2_gradient_polar_load_balanced_dEtensor
c            end if
            !  if(approxmode.ne.'1BODYMODE') then
            !    if(approxmode.eq.'2BODYMODE') then
            !      call smooth2b_gradient_polar_load_balanced
            !    else
            !      call smooth2_gradient_polar_load_balanced
            !    end if
            !  end if


c           call mpi_barrier(mpi_comm_world,ierr)

c            listsend_mpole=.false.
            call mpi_bcast(em,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(dem,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(viremrecip,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           if(uzepmedirpolz) then
           call mpi_bcast(dep3b_recip,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(virep3b_recip,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(ep3b_recip,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(virep3b_recip_pt1,3*3,mpi_real8,master,
     &       mpi_comm_world,ierr)
           end if
c            call mpi_waitall(6,reqs,stats2,ierr)

              if(taskid.eq.master) then
                 print*,"In master after irecv em=",em
               if(uzepmedirpolz) then
                 call empole1d_3b_PolarPerm_selfeng_bcast2
               else
                 call empole1d_3b_Perm_selfeng_bcast
               end if

                 if (use_vcorr) then
                  call evcorr1 (elrc,vlrc)
                  ev = ev + elrc
                  vir(1,1) = vir(1,1) + vlrc
                  vir(2,2) = vir(2,2) + vlrc
                  vir(3,3) = vir(3,3) + vlrc
                 end if

                 do i=1,3
                   do j=1,3
                    vir(j,i)=vir(j,i)+virep3b(j,i)+viremreal_tmp(j,i)
     &               + viremrecip(j,i) +virev(j,i)+virep3bmut(j,i)
     &               + virep3b1(j,i)
       print*,"j i virep3b(j,i)",j,i,virep3b(j,i)+virep3bmut(j,i)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b+ep3bmut
               print*,"em ev in verlet/nose/beeman setup",em,ev
               print*,"ep3b in verlet/nose/beeman setup",ep3b+ep3bmut
     &         !        +ep3b_recip
               print*,"eb in verlet/nose/beeman setup",eb
               print*,"ea in verlet/nose/beeman setup",ea
               print*,"eba in verlet/nose/beeman setup",eba
               print*,"eub in verlet/nose/beeman setup",eub
               print*,"eopb in verlet/nose/beeman setup",eopb
               print*,"et in verlet/nose/beeman setup",et
               print*,"ept in verlet/nose/beeman setup",ept
               print*,"ett in verlet/nose/beeman setup",ett

               allocate (derivs(3,n))

                do i = 1, n
                  do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i)
     &                 + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                 + dev(j,i) + dem(j,i)+dep3b(j,i)+dep3bmut(j,i)
     &                 + dep3b1(j,i)
     &        !      +dep3b_recip(j,i) 

                  end do
                 print*,"dep3b_x",i,dep3b(1,i)+dep3bmut(1,i)+dep3b1(1,i)
     &         !          +dep3b_recip(1,i)
                 print*,"dep3b_y",i,dep3b(2,i)+dep3bmut(2,i)+dep3b1(2,i)
     &          !         +dep3b_recip(2,i)
                 print*,"dep3b_z",i,dep3b(3,i)+dep3bmut(3,i)+dep3b1(3,i)
     &          !          +dep3b_recip(3,i)
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

      if(taskid.eq.master) then
        call mdinit_pt2_master
      end if
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
          delta = 0.000001d0 
          delx=0.000001d0
        end if

c      call mpi_bcast(delV,1,mpi_real8,master,
c     & mpi_comm_world,ierr)

         if (integrate .eq. 'VERLET') then
          do istep = 1, nstep

            if(taskid.eq.master) then 
             call verlet_pt1 (istep,dt)
c             call ptest_parallel

             ! x(6)=x(6)+delx
c              y(6)=y(6)+delx
c              z(6)=z(6)+delx
             call gradient_setup
             call allocPermElec1
c             allocate (derivs(3,n))
             !call prep_pole
             call alloc3b_vdwlist_nolist
c             call alloc_erecip
c            else if(taskid.eq.master1) then
c             call alloc_erecip
            end if
            
c            call mpi_barrier(mpi_comm_world,ierr)
            call alloc_erecip

            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
           ! call mpi_bcast(rpole,13*n,mpi_real8,master,
     &     !  mpi_comm_world,ierr)
            call prep_pole
            !call bcast_lattice

           if((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2)))then
c              print*,"In verlet/nose/beeman setup"
c              print*,"Before alloc3b_nblist verlet/nose/beeman setup"
               call vlist_par
              if(firstload_vdw.and.(taskid.eq.numtasks_emreal+2)) then
c                firstload_vdw=.false.
                call vdw_load_balance_sizelist_half
              end if
c             print*,"In verlet/nose/beeman setup,after vdw list shit"
           else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2))then
               call mlist_par
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
                if(firstload_emreal.and.(taskid.eq.master2)) then
c                 firstload_emreal = .false.
                  if(longrangepoldir) then
                   call emreal_load_balance_sizelist
                  else
                   call emreal_load_balance_sizelist_mol
                  end if
                end if
ccc             print*,"In verlet/nose/beeman setup,after emreal list shit"
ccc            else if(taskid.eq.master2
ccc     &                .and.(approxmode.ne.'1BODYMODE')) then
ccc                call mollist2bodyOO6_8_par
ccc              if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
ccc                firstload_polz = .false.
ccc                call ewpolar_load_balance
ccc              end if
           end if

           ! if((approxmode.ne.'1BODYMODE').and.
     &     !           (taskid.ge.numtasks_emreal+2)) then
            if((approxmode.ne.'1BODYMODE').and.
     &                (taskid.eq.numtasks_emreal+2)) then
              if(approxmode.eq.'2BODYMODE') then
              call mollist2bodyOO6_8_par_1list_par_newesttink7
              else
              call mollist2bodyOO6_8_par_1list
              end if
            !  if(firstload_polz.and.(taskid.eq.numtasks_emreal+2)) then
              if(firstload_polz) then
                call polar_load_balance7_half
              !  call ewpolar_load_balance_half
              end if
            end if

            call mpi_barrier(mpi_comm_world,ierr)

            call mpi_bcast(listbcast,1,mpi_logical,numtasks_emreal+2,
     &       mpi_comm_world,ierr)

           ! if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
            if(listbcast.and.(approxmode.ne.'1BODYMODE')) then
             firstload_polz = .false.
            ! call ewbcast_polar_load_balance_half
              call bcast_polar_load_balance_half
            call mpi_bcast(nmollst,nmol,mpi_integer,numtasks_emreal+2,
     &       mpi_comm_world,ierr)
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,
     &       numtasks_emreal+2,mpi_comm_world,ierr)
            end if
            listbcast=.false.

C IF 'listsend_mpole' IS TRUE, REAL-SPACE PERM ELECTROSTATICS NEIGHBORLIST HAS BEEN REBUILT OR BUILT DE NOVO.
C THEREFORE, LIST MUST BE BROADCAST. ALSO BROADCAST VARIABLES TO DO WITH LOAD BALANCING
c              call mpi_bcast(numtasks_emreal2,1,mpi_integer,
c     &             master1,mpi_comm_world,ierr)

            if(firstload_emreal) then
              firstload_emreal = .false.
              call mpi_bcast(numtasks_emreal2,1,mpi_integer,
     &             master2,mpi_comm_world,ierr)
              call mpi_bcast(start_emreal2,numtasks_emreal,mpi_integer,
     &             master2,mpi_comm_world,ierr)
              call mpi_bcast(last_emreal2,numtasks_emreal,mpi_integer,
     &             master2,mpi_comm_world,ierr)
c               call mpi_bcast(maxsize_elst,numtasks_emreal,mpi_integer,
c     &             master1,mpi_comm_world,ierr)
c              call sendmlist_load_balanced3_sizelist
            end if
c            listsend_mpole = .false.

c              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
c     &         numtasks_emreal,mpi_comm_world,ierr)

            if(firstload_vdw) then
              firstload_vdw = .false.
              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
     &         numtasks_emreal+2,mpi_comm_world,ierr)
              call mpi_bcast(start_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
               call mpi_bcast(last_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
c               call mpi_bcast(maxsize_vlst,numtasks_vdw,mpi_integer,
c     &             master,mpi_comm_world,ierr)
c               call sendvlist_load_balanced3_sizelist_half
            end if
c            listsend_vdw = .false.
 


      call grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2

c            if(approxmode.eq.'1BODYMODE') then
            call mpi_bcast(uind,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
            call mpi_bcast(uinp,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
            call mpi_bcast(ep3b,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
c            else
c            call mpi_bcast(fieldnpole,3*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(fieldpnpole,3*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(fphi_totfield,20*npole,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vxx_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vyx_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vzx_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vyy_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vzy_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(vzz_perm,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(thetai1,4*bsorder*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(thetai2,4*bsorder*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(thetai3,4*bsorder*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(qfac,nfft1*nfft2*nfft3,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(pmetable,n*nchunk,mpi_integer,master1,
c     &       mpi_comm_world,ierr)
c            call mpi_bcast(igrid,3*n,mpi_integer,master1,
c     &       mpi_comm_world,ierr)
c            end if
         !call ewtotfield_gradient_polar(start,start+offset-1,moli1rmndr)

c            if(approxmode.eq.'1BODYMODE') then
           ! call ewtotfieldsmooth2_gradient_polar_1body_dEtensor_mpioff
          call ewtotfieldsmooth2_gradient_polar_1body_dEtensor_simult

             ! if(approxmode.ne.'1BODYMODE') then
             !   if(approxmode.eq.'2BODYMODE') then
             !     call smooth2b_gradient_polar_load_balanced
             !   else
             !     call smooth2_gradient_polar_load_balanced
             !   end if
             ! end if

c            call mpi_waitall(6,reqs,stats2,ierr)

c            else
c            call ewtotfieldsmooth2_gradient_polar_load_balanced_dEtensor
c            end if
c            call mpi_barrier(mpi_comm_world,ierr)
           call mpi_bcast(em,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(dem,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(viremrecip,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)

           if(uzepmedirpolz) then
           call mpi_bcast(dep3b_recip,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(virep3b_recip,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)              
           call mpi_bcast(ep3b_recip,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(virep3b_recip_pt1,3*3,mpi_real8,master,
     &       mpi_comm_world,ierr)
           end if

            if(taskid.eq.master) then
               if(uzepmedirpolz) then
                 call empole1d_3b_PolarPerm_selfeng_bcast2
               else          
                 call empole1d_3b_Perm_selfeng_bcast
               end if
                 if (use_vcorr) then
                  call evcorr1 (elrc,vlrc)
                  ev = ev + elrc
                  vir(1,1) = vir(1,1) + vlrc
                  vir(2,2) = vir(2,2) + vlrc
                  vir(3,3) = vir(3,3) + vlrc
                 end if

                 do i=1,3
                   do j=1,3
                    vir(j,i)=vir(j,i)+viremreal_tmp(j,i)+virep3b(j,i)
     &                  + viremrecip(j,i) + virev(j,i)+virep3bmut(j,i)
     &               +virep3b1(j,i)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b+ep3bmut
c      call system_clock(t1,clock_rate)              
             allocate (derivs(3,n))
c              efindif_ep3b(istep)=ep3b+ep3bmut
c              efindif(istep)=energy!-ep3b-ep3bmut-em-ev
c!$OMP PARALLEL DO default(shared) private(i,j)              
              do i = 1, n
                do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i) 
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i) +dep3bmut(j,i)
     &                   + dep3b1(j,i)
                end do
               if (use(i)) then
                do j = 1, 3
                 a(j,i) = -convert * derivs(j,i) / mass(i)
                 v(j,i) = v(j,i) + a(j,i)*dt_2
c                 print*,"istep i j v(j,i)",istep,i,j,v(j,i)
                end do
               end if
              end do
c!$OMP END PARALLEL DO

c              analyticgx(istep)=derivs(1,6)!-(dep3bmut(1,6)
     &              !+dep3b_recip(1,6)+dep3b1(1,6)+dem(1,6)+dev(1,6))
c              analyticg_dep3bx(istep)=dep3b(1,6)+dep3bmut(1,6)
c     &              +dep3b1(1,6)

c              analyticgy(istep)=derivs(2,6)!-(dep3bmut(2,6)
     &            !  +dep3b_recip(2,6)+dep3b1(2,6)+dem(2,6)+dev(2,6))
c              analyticg_dep3by(istep)=dep3b(2,6)+dep3bmut(2,6)
c     &              +dep3b1(2,6)

c              analyticgz(istep)=derivs(3,6)!-(dep3bmut(3,6)
     &            !  +dep3b_recip(3,6)+dep3b1(3,6)+dem(3,6)+dev(3,6))
c              analyticg_dep3bz(istep)=dep3b(3,6)+dep3bmut(3,6)
c     &              +dep3b1(3,6)


c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing Verlet fullstep v",tot_time

       print*," istep Vdw, PermElec, Pol",istep,ev,em,ep3b+ep3bmut

              deallocate (derivs)
               call verlet_pt2(istep,dt,energy)
c              if(istep.eq.3) then
c                print*,"Finite Diff PairwiseAddEng, PairwiseAddGradz",
c     & (efindif(3)-efindif(1))/2.0d0/delx,analyticgy(2)
c                print*,"Finite Diff PolEng, PolGradz",
c     & (efindif_ep3b(3)-efindif_ep3b(1))/2.0d0/delx,analyticg_dep3by(2)
c              end if

c              if(istep.eq.2) then
c                do i=1,3
c                   do j=1,3
c                 vir_analytic(j,i)=vir(j,i)-virep3b(j,i)-virep3bmut(j,i)
c                 virep3b_analytic(j,i)=virep3b(j,i)+virep3bmut(j,i)
c                print*,"Step 2 virep3b",j,i,virep3b(j,i)+virep3bmut(j,i)
c                   end do
c                end do
c                volbox_analytic=volbox
c              end if
c              if(istep.eq.3) then
c      dedv_vir = (vir_analytic(1,1)+vir_analytic(2,2)+vir_analytic(3,3))
c     &         / (3.0d0*volbox)
c              dedv_fd = (efindif(3)-efindif(1))
c     &        / (2.0d0*delta*volbox)
c
c      dedv_vir_ep3b = (virep3b_analytic(1,1)+virep3b_analytic(2,2)
c     &          +virep3b_analytic(3,3))/ (3.0d0*volbox)
c      dedv_fd_ep3b =(efindif_ep3b(3)-efindif_ep3b(1))
c     &         / (2.0d0*delta*volbox)
c
c      write (iout,31)  dedv_vir
c   31 format (/,'PairwsdE/dV(Virial-based):',f27.16,' Kcal/mol/A**3')
c      write (iout,32)  dedv_fd
c   32 format ('PairwsdE/dV(Finite Diff):',f27.16,' Kcal/mol/A**3')
c      write (iout,33)  dedv_vir_ep3b
c   33 format (/,'PolzdE/dV (Virial-based):',f27.16,' Kcal/mol/A**3')
c      write (iout,34)  dedv_fd_ep3b
c   34 format ('PolzdE/dV (Finite Diff) :',f27.16,' Kcal/mol/A**3')
c              end if

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
            if(taskid.eq.master) then
             call nose_pt1(istep,dt,press)

             call gradient_setup
             call allocPermElec1
             call alloc3b_vdwlist_nolist
c             call alloc_erecip
c            else if(taskid.eq.master1) then
c             call alloc_erecip
            end if

            call alloc_erecip

            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call bcast_lattice
            call prep_pole

           if((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2)))then
c              print*,"In verlet/nose/beeman setup"
c              print*,"Before alloc3b_nblist verlet/nose/beeman setup"
               call vlist_par
              if(firstload_vdw.and.(taskid.eq.numtasks_emreal+2)) then
c                firstload_vdw=.false.
                call vdw_load_balance_sizelist_half
              end if
c             print*,"In verlet/nose/beeman setup,after vdw list shit"
           else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2))then
               call mlist_par
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
                if(firstload_emreal.and.(taskid.eq.master2)) then
c                 firstload_emreal = .false.
                  if(longrangepoldir) then
                   call emreal_load_balance_sizelist
                  else
                   call emreal_load_balance_sizelist_mol
                  end if
                end if
ccc             print*,"In verlet/nose/beeman setup,after emreal list shit"
ccc            else if(taskid.eq.master2
ccc     &                .and.(approxmode.ne.'1BODYMODE')) then
ccc                call mollist2bodyOO6_8_par
ccc              if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
ccc                firstload_polz = .false.
ccc                call ewpolar_load_balance
ccc              end if
           end if

           ! if((approxmode.ne.'1BODYMODE').and.
     &     !           (taskid.ge.numtasks_emreal+2)) then
            if((approxmode.ne.'1BODYMODE').and.
     &                (taskid.eq.numtasks_emreal+2)) then
              if(approxmode.eq.'2BODYMODE') then
              call mollist2bodyOO6_8_par_1list_par_newesttink7
              else
              call mollist2bodyOO6_8_par_1list
              end if
            !  if(firstload_polz.and.(taskid.eq.numtasks_emreal+2)) then
              if(firstload_polz) then
              !  call ewpolar_load_balance_half
                call polar_load_balance7_half
              end if
            end if

            call mpi_barrier(mpi_comm_world,ierr)

            call mpi_bcast(listbcast,1,mpi_logical,numtasks_emreal+2,
     &       mpi_comm_world,ierr)

           ! if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
            if(listbcast.and.(approxmode.ne.'1BODYMODE')) then
             firstload_polz = .false.
           !  call ewbcast_polar_load_balance_half
             call bcast_polar_load_balance_half
            call mpi_bcast(nmollst,nmol,mpi_integer,numtasks_emreal+2,
     &       mpi_comm_world,ierr)
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,
     &       numtasks_emreal+2,mpi_comm_world,ierr)
            end if
            listbcast=.false.

            if(firstload_emreal) then
              firstload_emreal = .false.
              call mpi_bcast(numtasks_emreal2,1,mpi_integer,
     &             master2,mpi_comm_world,ierr)
              call mpi_bcast(start_emreal2,numtasks_emreal,mpi_integer,
     &             master2,mpi_comm_world,ierr)
              call mpi_bcast(last_emreal2,numtasks_emreal,mpi_integer,
     &             master2,mpi_comm_world,ierr)
            end if

            if(firstload_vdw) then
              firstload_vdw = .false.
              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
     &         numtasks_emreal+2,mpi_comm_world,ierr)
              call mpi_bcast(start_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
               call mpi_bcast(last_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
            end if

      call grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2

            call mpi_bcast(uind,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
            call mpi_bcast(uinp,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
            call mpi_bcast(ep3b,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)

          !  call ewtotfieldsmooth2_gradient_polar_1body_dEtensor_mpioff
          call ewtotfieldsmooth2_gradient_polar_1body_dEtensor_simult

          !    if(approxmode.ne.'1BODYMODE') then
          !      if(approxmode.eq.'2BODYMODE') then
          !        call smooth2b_gradient_polar_load_balanced
          !      else
          !        call smooth2_gradient_polar_load_balanced
          !      end if
          !    end if

           call mpi_bcast(em,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(dem,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(viremrecip,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)

           if(uzepmedirpolz) then
           call mpi_bcast(dep3b_recip,3*n,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(virep3b_recip,3*3,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(ep3b_recip,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)
           call mpi_bcast(virep3b_recip_pt1,3*3,mpi_real8,master,
     &       mpi_comm_world,ierr)
           end if

            if(taskid.eq.master) then
             !  if(uzepmedirpolz) then
                 call empole1d_3b_PolarPerm_selfeng_bcast2
             !  else
             !    call empole1d_3b_Perm_selfeng_bcast
             !  end if
                 if (use_vcorr) then
                  call evcorr1 (elrc,vlrc)
                  ev = ev + elrc
                  vir(1,1) = vir(1,1) + vlrc
                  vir(2,2) = vir(2,2) + vlrc
                  vir(3,3) = vir(3,3) + vlrc
                 end if

                 do i=1,3
                   do j=1,3
                    vir(j,i)=vir(j,i)+viremreal_tmp(j,i)+virep3b(j,i)
     &                  + viremrecip(j,i) + virev(j,i)+virep3bmut(j,i)
     &                  +virep3b1(j,i)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b+ep3bmut
c      call system_clock(t1,clock_rate)              
             allocate (derivs(3,n))
              !efindif_ep3b(istep)=ep3b+ep3bmut
              !efindif(istep)=energy-ep3b-ep3bmut
c!$OMP PARALLEL DO default(shared) private(i,j)              
              do i = 1, n
                do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i)
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i) +dep3bmut(j,i)
     &                   + dep3b1(j,i)
                end do
               if (use(i)) then
                do j = 1, 3
                 a(j,i) = -convert * derivs(j,i) / mass(i)
                 v(j,i) = v(j,i) + a(j,i)*dt_2
c                 print*,"istep i j v(j,i)",istep,i,j,v(j,i)
                end do
               end if
              end do
c!$OMP END PARALLEL DO

       print*," istep Vdw, PermElec, Pol",istep,ev,em,ep3b+ep3bmut

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
            call smooth_gradient_polar(start,start+offset-1,moli1rmndr)
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
