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
      program uinddynrMPIOMPkBcast_rloadbal_nolistsend_funcdecomp_clust3
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
      use mpidat2
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
      use neigh2clust
      use deriv1bmat
      use usolve
      use tarray
      use polpot
      use deriv2bmat
      use uprior
      use mpidat4
      use mpidat3
      use uind2bmat
      use uind1bmatrix
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
      real*8 efindif_ep3b(3),analyticg_dep3b(3),efindif_em(3)
      real*8 delx,delta,vir_analytic(3,3),virep3b_analytic(3,3)
      real*8 volbox_analytic,dedv_vir,dedv_fd,dedv_vir_ep3b,dedv_fd_ep3b
      real*8 analyticgx(3)
      real*8 analyticg_dep3bx(3),analyticg_demx(3)
      real*8 analyticgy(3)
      real*8 analyticg_dep3by(3),analyticg_demy(3)
      real*8 analyticgz(3)
      real*8 analyticg_dep3bz(3),analyticg_demz(3)
      integer atomind,start_emreal
      integer emreal_group,emreal_group2,emreal_group3,k
      integer, allocatable :: taskidother(:)
      integer, allocatable :: taskidemreal(:)
      integer, allocatable :: taskidvdw(:)
      integer vdw_group
      logical firstload_emreal,firstload_polz
      logical firstload_emreal_bcast,firstload_polz_bcast
      logical firstload_vdw
      integer offset2,remainder2,start2,last2,last,numtasks_tmp
      real*8 elrc,vlrc
      logical setup_pred1b,setup_pred2b,setup_pred3b
      logical done_clustlist_setup,done_comm_create
      integer numtasks_polar_old,numtasks_polar3b_old

      INTEGER NMAX,TMAX,TMAX2,LMAX
      PARAMETER (NMAX=5000,TMAX=5001,TMAX2=5001)

      INTEGER NMAX2
      PARAMETER (NMAX2=5000)

      INTEGER ISTOREV,DTCALCV,ISKIPV,ISTORED,DTCALCD,ISKIPD,NUM
      integer L,M,Np,Nc,Nw,II,JJ,II2,istep0
      integer ISTEP2,TT,DTC,ISTOREM,DTCALCM,ISKIPM
      integer wlist(NMAX2)
      real*8 xcm,ycm,zcm,norm,dx,dy,dz,Rg,T
      real*8 XC(NMAX),YC(NMAX),ZC(NMAX)


      INTEGER COUNTM(TMAX),COUNTV(TMAX),COUNTD(TMAX2)
      INTEGER countmass
      INTEGER RESTARTTIMESTEP,RESTARTFLAG,SAVE_LENGTH                ! restart variables
      INTEGER SAVE_DATA_FREQ,SAVE_REST_FREQ,t_stp                    ! restart variables
      INTEGER O_IDX,H_IDX,LYR_IDX,ATM_IDX                            ! restart variables
      real*8 XO(NMAX,TMAX),VXO(NMAX,TMAX),VXH(NMAX,TMAX)
      real*8 YO(NMAX,TMAX),VYO(NMAX,TMAX),VYH(NMAX,TMAX)
      real*8 ZO(NMAX,TMAX),VZO(NMAX,TMAX),VZH(NMAX,TMAX)
      real*8 MXC,MYC,MZC,deltax,deltay,deltaz

ccc      real*8 MXI(NMAX,TMAX),MYI(NMAX,TMAX),MZI(NMAX,TMAX)
      real*8 MXI(TMAX),MYI(TMAX),MZI(TMAX)

ccc      real*8 MXIC(NMAX),MYIC(NMAX),MZIC(NMAX)

      real*8 MXCT,MYCT,MZCT

      real*8 DR2,MSDO(TMAX)
      real*8 VACFO(TMAX),VACFH(TMAX)
      real*8 VX(NMAX),VY(NMAX),VZ(NMAX),eksum,ekin,temp
      real*8 drsq
      real*8 DCF(TMAX2)

      real*8 tempMSDO,tempVACFO,tempVACFH                       ! temporary local variables for omp
      real*8 tempCOUNTM,tempCOUNTV                              ! temporary local variables for omp
      real*8 privateVACFO,privateVACFH                          ! temporary local variables for omp

      real*8 summass,drsqmin

      integer KK,NGRP2

      real*8 startTime, finishTime, netStartTime, netFinishtime         ! to get timings
      real*8 wall,cpu       ! To get timings

c     Wall time
      INTEGER nb_ticks_initial,nb_ticks_final,nb_ticks_max,nb_ticks_sec
      INTEGER nb_ticks
      REAL elapsed_time  ! real time in seconds

      CHARACTER FILE1*99,FILE2*99,FILE3*99,FILE4*99,FILE5*99,FILE6*99
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                       c
c     These files save the values from the tape, for use to restart the simulation      c
c                                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CHARACTER FILE7*99,FILE8*99,FILE9*99  ! These files save the X, Y and Z coordinates for the Oxygens of the waters.
      CHARACTER FILE10*99,FILE11*99,FILE12*99 ! These files save the X, Y and Z velocities for the Oxygens of the waters.
      CHARACTER FILE13*99,FILE14*99,FILE15*99 ! These files save the X, Y and Z velocities for the Hydrogens of the waters.
      CHARACTER FILE16*99,FILE17*99,FILE18*99 ! These files save the X, Y and Z dipole moments for the 13 layers.
      CHARACTER FILE19*99,FILE24*99           ! This file writes out the file which contains the index upto which the simulation has progressed, in order to restart from that point.
      CHARACTER FILE20*99,FILE21*99,FILE22*99 ! These files write out the COUNTM, COUNTV and COUNTD values so that they can be multiplied into the values of MSD, VACF_O, VACF_H and DCF read in from the restart files.
      CHARACTER FILE23*99 ! This file writes out the summass array.
ccc     call cpu_time(startTime)
ccc     call cpu_time(netStartTime)
ccc     CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
ccc     CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c
      logical exist1

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
      save setup_pred1b      
      data setup_pred1b / .false. /
      save setup_pred2b
      data setup_pred2b / .false. /
      save setup_pred3b
      data setup_pred3b / .false. /
      save done_clustlist_setup
      data done_clustlist_setup / .false. /
      save done_comm_create
      data done_comm_create / .false. /
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
      !if(taskid.eq.master) then
      !   call settime2(twall,tcpu)
      !end if
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

c      do2waterclustlist=.false.
      if(taskid.eq.master) then
c         call mechanic_pme4_eastman
         call mechanic_pme4
         call kewald3b
c         call kewald 
c         call kewald_eastman
         !call kewaldreg3b
         !call get_numtasks_emreal
c         numtasks_emreal=int((numtasks-2)/2)
c         numtasks_vdw=numtasks_emreal
c         numtasks_vdw2=numtasks_vdw
c         numtasks_emreal2=numtasks_emreal
         if(numtasks.ge.13) then
           numtasks_tmp=numtasks-2
           numtasks_vdw=int(numtasks_tmp/11)
           numtasks_emreal=10*numtasks_vdw+mod(numtasks_tmp,11)
         else if ((numtasks .lt.13).and.(numtasks.ge.6)) then
           numtasks_vdw=2
           numtasks_emreal=numtasks-4
         else if ((numtasks .le.5).and.(numtasks.ge.4)) then
           numtasks_vdw=1
           numtasks_emreal=numtasks-3
         end if
         numtasks_vdw2=numtasks_vdw
         numtasks_emreal2=numtasks_emreal

C  ALLOCATIONS FOR CLUSTER VERSION
        if (.not. allocated(xmolold2)) allocate(xmolold2(nmol))
        if (.not. allocated(ymolold2)) allocate(ymolold2(nmol))
        if (.not. allocated(zmolold2)) allocate(zmolold2(nmol))

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

      !call kewald3b
      !call kewald3b_totfield
c      if(approxmode.eq.'1BODYMODE') then
       !print*,"Doing uind allocs"
c       if(.not. allocated(uind)) allocate(uind(3,n))
c       if(.not. allocated(uinp)) allocate(uinp(3,n))
c         if(.not.allocated(dep3b_recip)) allocate (dep3b_recip(3,n))

C  ALLOCS NECESSARY FOR DEBUGGING MODE ONLY!! REMOVE LATER!!
c        if(.not.allocated(demreal_tmp)) allocate(demreal_tmp(3,n))
c         if (.not. allocated(fieldnpole))  allocate (fieldnpole(3,n))
c         if (.not. allocated(fieldpnpole))  allocate (fieldpnpole(3,n))
c         if(.not.allocated(dep3b)) allocate (dep3b(3,n))
c         if(.not.allocated(dev)) allocate (dev(3,n))
c         if(.not.allocated(dem)) allocate (dem(3,n))
      !   if(.not.allocated(dep3bmut)) allocate (dep3bmut(3,n))
      !  if(.not.allocated(dep3b1)) allocate (dep3b1(3,n))

C DEFINE A NEW COMM FOR THE MPI/OMP REAL-SPACE PERM ELECTROSTATICS
C USING ONLY 'numtasks_emreal' TASKS    
      call mpi_bcast(numtasks_emreal,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(numtasks_vdw,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(numtasks_emreal2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(numtasks_vdw2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)

c      print*,"numtasks_emreal after read and bcast",numtasks_emreal


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


      if(taskid.ne.master) then
         call mechanic_parallel2
c         call kewald
c        if((taskid.ne.master1).and.use_ewaldclust) then
          !call kewald3b
          
c        end if
      end if
      if(taskid.eq.master1) then
         call kewald
         !call kewald_eastman
         !call cutoffs_ewald
      end if

      if((taskid.ne.master1).and.use_ewaldclust) then
         call kewald
      end if

      if((taskid.gt.1).and.(taskid.lt.numtasks_emreal+2)) then
        call cutoffs_ewald
      else if(taskid.ge.numtasks_emreal+2) then
        call cutoffs_vdw
      end if

        call cutoffs2_serial

C  ALLOCATIONS FOR CLUSTER VERSION
      if(inputkmeansclust) then
         call getclust
c        if(.not. allocated(ep1bmat)) allocate(ep1bmat(clustcount)) 
c        if(.not. allocated(dep1bmat)) allocate(dep1bmat(3,
c     &    3*maxsizeclust,clustcount))
c        if(.not. allocated(virep1bmat)) allocate(virep1bmat(3,
c     &    3,clustcount))
c        if(.not. allocated(em1bmat)) allocate(em1bmat(clustcount))
c        if(.not. allocated(dem1bmat)) allocate(dem1bmat(3,
c     &    3*maxsizeclust,clustcount))
         if (poltyp .eq. 'MUTUAL') then
            if (.not.allocated(tindex))
     &         allocate (tindex(2,3*3*maxsizeclust*maxelst))
            if (.not.allocated(tdipdip))
     &         allocate (tdipdip(6,3*3*maxsizeclust*maxelst))
            if (.not.allocated(minv))
     &         allocate (minv(6*maxulst*3*3*maxsizeclust))
            if (.not.allocated(mindex))
     &         allocate (mindex(3*3*maxsizeclust))
         end if

       if(taskid.eq.master) then
        if(.not. allocated(ep1bmat)) allocate(ep1bmat(clustcount))
        if(.not. allocated(dep1bmat)) allocate(dep1bmat(3,
     &    3*maxsizeclust,clustcount))

        if(.not. allocated(uind1bmat)) allocate(uind1bmat(3,
     &    3*maxsizeclust,clustcount))
        if(.not. allocated(virep1bmat)) allocate(virep1bmat(3,
     &    3,clustcount))
c        if(.not. allocated(em1bmat)) allocate(em1bmat(clustcount))
c        if(.not. allocated(dem1bmat)) allocate(dem1bmat(3,
c     &    3*maxsizeclust,clustcount))
c        if(.not. allocated(ep2bmat)) 
c     &           allocate(ep2bmat(clustcount,clustcount))
c        if(.not. allocated(dep2bmat)) allocate(dep2bmat(3,
c     &    2*3*maxsizeclust,clustcount,clustcount))
c        if(.not. allocated(virep2bmat)) allocate(virep2bmat(3,
c     &    3,clustcount,clustcount))
c        if(.not. allocated(em2bmat)) 
c     &        allocate(em2bmat(clustcount,clustcount))
c        if(.not. allocated(dem2bmat)) allocate(dem2bmat(3,
c     &    2*3*maxsizeclust,clustcount,clustcount))
       end if
      else
         if (.not. allocated(clust))  allocate (clust(4,nmol))
         if (.not. allocated(sizeclust))  allocate (sizeclust(nmol))
         if (.not. allocated(clust_cm))  allocate (clust_cm(3,nmol))
         if (.not. allocated(distmax))  allocate (distmax(nmol))
      end if

        if(.not.allocated(nmollst3)) allocate(nmollst3(clustcount))
        if(.not.allocated(mollst3)) 
     &         allocate(mollst3(7*clustcount,clustcount))
        if(.not.allocated(save2bfor3b)) 
     &         allocate(save2bfor3b(clustcount,clustcount))
        if(.not.allocated(ind2b_counter))
     &      allocate(ind2b_counter(clustcount,clustcount))


      call mpi_barrier(mpi_comm_world,ierr)

        ! ewtotfieldsmallsmooth2_gradient_polar_1body(
     &  !     start,start+offset-1,moli1rmndr)

c      if(itermode.ne.1) then
c         offset=int(nmol/numtasks)
c         remainder=mod(nmol,numtasks)
c
c         remainder2=mod(offset,itermode)
c
c         if(taskid.le.remainder-1) then
c           moli1rmndr=numtasks*offset+taskid+1
c         else if (taskid.gt.remainder-1) then
c           moli1rmndr=0
c         end if
c         start=taskid*offset+1
c         offset2=offset-remainder2
c         last=start+offset2-1
c         start2=last+1
c         last2=start2+remainder2-1
c      else 
c         offset=int(nmol/numtasks)
c         remainder=mod(nmol,numtasks)
c
c         if(taskid.le.remainder-1) then
c           moli1rmndr=numtasks*offset+taskid+1
c         else if (taskid.gt.remainder-1) then
c           moli1rmndr=0
c         end if
c         start=taskid*offset+1
c      end if

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
         call mpi_bcast(use_pred,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(polpred,4,mpi_character,
     &   master,mpi_comm_world,ierr)

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
c          beeman/nose/verlet setup

            if(taskid.eq.master) then
              call gradient_setup
              call alloc3b_vdwlist_nolist_vac
              call allocPermElec1

              call alloc3b_uind
            end if

              call alloc_erecip

            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call prep_pole
            call bspline_fill_omp
            call table_fill_omp
            call mpi_barrier(mpi_comm_world,ierr)

            if(inputkmeansclust.and.(taskid.eq.master)) then
             !call init1bmatPolarPerm
             call init1bmat
             call init1bmat_uind
            end if

           if((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2)))then
               call vlist_par
              if(firstload_vdw.and.(taskid.eq.numtasks_emreal+2)) then
                call vdw_load_balance_sizelist_half
              end if
           else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2))then
               call mlist_par
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
                if(firstload_emreal.and.(taskid.eq.master2)) then
c                 firstload_emreal = .false.
                   call emreal_load_balance_sizelist
                end if
           else if(taskid.eq.master) then

              if(do2waterclustlist) then
                 call mollist2body_NewBodieskmeans
              end if
              if(firstload_polz.and.do2waterclustlist) then
                firstload_polz = .false.
                !call polar_load_balance7
                !call polar_load_balance_clust
                call polar_load_balance_clust_splitlist_offset
                 call polar_load_balance_clust_splitlist_offset3b_simple
              end if
           !   done_clustlist_setup = .true.        
           end if


              numtasks_polar1=clustcount
              offset=int(clustcount/numtasks_polar1)
              remainder=mod(clustcount,numtasks_polar1)

              if(taskid.le.remainder-1) then
               moli1rmndr=numtasks_polar1*offset+taskid+1
              else if (taskid.gt.remainder-1) then
               moli1rmndr=0
              end if
              start=taskid*offset+1

           allocate (taskidemreal(0:numtasks_polar1-1))
           do j=0,numtasks_polar1-1
            taskidemreal(j)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
            call mpi_group_incl(orig_group,numtasks_polar1,taskidemreal,
     &        emreal_group,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group,
     &            emreal_comm,ierr)
            deallocate(taskidemreal)

            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)
 
            !if(firstload_polz.and.(approxmode.ne.'1BODYMODE')) then
            if(listbcast.and.(approxmode.ne.'1BODYMODE')
     &                    .and.do2waterclustlist) then
             firstload_polz = .false.
           !  call ewbcast_polar_load_balance
           !  call bcast_polar_load_balance_clust
             call bcast_polar_load_balance_clust_splitlist
             call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
             call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
             !if(approxmode.eq.'3BODYMODE') then
               call bcast_polar_load_balance_clust_splitlist3b
               call mpi_bcast(nmollst3,clustcount,mpi_integer,master,
     &              mpi_comm_world,ierr)
               call mpi_bcast(mollst3,7*clustcount*clustcount,
     &              mpi_integer,master,mpi_comm_world,ierr)
               call mpi_bcast(save2bfor3b,clustcount*clustcount,
     &              mpi_logical,master,mpi_comm_world,ierr)
              call mpi_bcast(num2bsave,1,mpi_integer,master,
     &              mpi_comm_world,ierr)
              call mpi_bcast(ind2b_counter,clustcount*clustcount,
     &              mpi_integer,master,mpi_comm_world,ierr)
c               if(.not.allocated(counter2b_ind))
c     &             allocate(counter2b_ind(2,num2bsave))
c              call mpi_bcast(counter2b_ind,num2bsave*2,mpi_integer,
c     &           master,mpi_comm_world,ierr)
             ! listbcast=.false.
             !end if
            end if
c            listbcast = .false.
       if(taskid.eq.master) then

        if(.not. allocated(ep2bmatmod))
     &           allocate(ep2bmatmod(num2bsave))
        if(.not. allocated(dep2bmatmod)) allocate(dep2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
        if(.not. allocated(virep2bmatmod)) allocate(virep2bmatmod(3,
     &          3,num2bsave))
        if(.not. allocated(uind2bmatmod)) allocate(uind2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
             call init2bmatmodPolar
             call init2bmatmodPolar_uind
       end if
c        if(.not. allocated(em2bmatmod))
c     &        allocate(em2bmatmod(num2bsave))
c        if(.not. allocated(dem2bmatmod)) allocate(dem2bmatmod(3,
c     &    2*3*maxsizeclust,num2bsave))

            if((taskid.lt.numtasks_polar1).and.mlistclust) then
              call clust1blists(start,start+offset-1, moli1rmndr) 
              print*,"Completed clust1blists"
              if(use_pred) then
                call initpolarpred1b(start,start+offset-1,moli1rmndr)
                setup_pred1b=.true.
              end if
            end if

           if((taskid.ge.numtasks_polar1).and.
     &             (taskid.lt.(numtasks_polar+numtasks_polar1))
     &             .and.mlistclust) then
              call clust2blists_offset
              print*,"Completed clust2blists_offset"
              if(use_pred) then
                call initpolarpred2b_offset
                setup_pred2b=.true.
              end if
           end if


           if( (taskid.ge.(numtasks_polar+numtasks_polar1)).and.
     &         (taskid.lt.
     &         (numtasks_polar+numtasks_polar1+numtasks_polar3b))
     &         .and.mlistclust) then 
              call clust3blists_offset
              print*,"Completed clust3blists_offset"
              if(use_pred) then
                call initpolarpred3b_offset
                setup_pred3b=.true.
              print*,"Completed initpolarpred3b_offset"
              end if
           end if

            listbcast=.false.           
            call mpi_barrier(mpi_comm_world,ierr)

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
           allocate(taskidvdw(0:numtasks_vdw2))
           taskidvdw(0)=0
           do j=numtasks_emreal+2,numtasks_vdw2+numtasks_emreal+1
             k=j-(numtasks_emreal+2)+1
             taskidvdw(k)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
           call mpi_group_incl(orig_group,numtasks_vdw2+1,taskidvdw,
     &        vdw_group,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,vdw_group,
     &            vdw_comm,ierr)
           deallocate(taskidvdw)


           allocate (taskidemreal(0:numtasks_polar))
           taskidemreal(0)=0
           do j=numtasks_polar1,numtasks_polar+numtasks_polar1-1
              k=j-numtasks_polar1+1
              taskidemreal(k)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
         call mpi_group_incl(orig_group,numtasks_polar+1,taskidemreal,
     &        emreal_group2,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group2,
     &            emreal_comm2,ierr)
            deallocate(taskidemreal)


            allocate (taskidemreal(0:numtasks_polar3b))
            taskidemreal(0)=0
            do j=numtasks_polar+numtasks_polar1,
     &            numtasks_polar+numtasks_polar1+numtasks_polar3b-1
               k=j-(numtasks_polar+numtasks_polar1)+1
               taskidemreal(k)=j
            end do
         call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
        call mpi_group_incl(orig_group,numtasks_polar3b+1,taskidemreal,
     &        emreal_group3,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group3,
     &            emreal_comm3,ierr)
            deallocate(taskidemreal)
            done_comm_create=.true.

          call uindclst_gradientPolar1b2b3bsimult_ireducematcommsmall2(
     &      start,start+offset-1, moli1rmndr)

            call grad_covemrecip_vandr_reduce_totfield_commsmall

            !end if


            if(taskid.eq.master1) then
             call mpi_wait(reqs1,stat,ierr)
             call mpi_wait(reqs2,stat,ierr)
             call mpi_wait(reqs3,stat,ierr)
            else if(taskid.eq.master) then
             call mpi_wait(reqs4,stat,ierr)
             call mpi_wait(reqs5,stat,ierr)
             call mpi_wait(reqs6,stat,ierr)
            end if
             call mpi_wait(reqs13,stat,ierr)
             call mpi_wait(reqs14,stat,ierr)
             call mpi_wait(reqs15,stat,ierr)

             if((taskid.eq.master).or.((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) ) then
             call mpi_wait(reqs16,stat,ierr)
             call mpi_wait(reqs17,stat,ierr)
             call mpi_wait(reqs18,stat,ierr)
             end if

       if(taskid.lt.numtasks_polar1) then
             call mpi_wait(reqs19,stat,ierr)
             call mpi_wait(reqs20,stat,ierr)
             call mpi_wait(reqs21,stat,ierr)
             call mpi_wait(reqs22,stat,ierr)
             call mpi_wait(reqs23,stat,ierr)
             call mpi_wait(reqs24,stat,ierr)
       end if

      if( (taskid.eq.master) .or. ( (taskid.ge.numtasks_polar1)
     &   .and.(taskid.lt.(numtasks_polar+numtasks_polar1))) ) then
             call mpi_wait(reqs25,stat,ierr)
             call mpi_wait(reqs26,stat,ierr)
             call mpi_wait(reqs27,stat,ierr)
             call mpi_wait(reqs28,stat,ierr)
             call mpi_wait(reqs29,stat,ierr)
             call mpi_wait(reqs30,stat,ierr)
      end if

      if( (taskid.eq.master) .or. (
     & (taskid.ge.(numtasks_polar+numtasks_polar1))
     &    .and. (taskid.lt.(numtasks_polar+numtasks_polar1
     &         +numtasks_polar3b))) ) then
             call mpi_wait(reqs31,stat,ierr)
             call mpi_wait(reqs32,stat,ierr)
             call mpi_wait(reqs33,stat,ierr)
      end if
      call uindsum_gradient_Polar1b2b3bsimult_ireducematcommsmall2

              if(taskid.eq.master) then
                 print*,"In master after bcast or irecv em=",em
                 call empole1d_3b_Perm_selfeng_bcast

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
     &               + viremrecip(j,i) +virev(j,i)
       print*,"j i virep3b(j,i)",j,i,virep3b(j,i)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b!+ep3bmut
               print*,"em ev in verlet/nose/beeman setup",em,ev
               print*,"ep3b in verlet/nose/beeman setup",ep3b!+ep3bmut
     &         !        +ep3b_recip
c               print*,"eb in verlet/nose/beeman setup",eb
c               print*,"ea in verlet/nose/beeman setup",ea
c               print*,"eba in verlet/nose/beeman setup",eba
c               print*,"eub in verlet/nose/beeman setup",eub
c               print*,"eopb in verlet/nose/beeman setup",eopb
c               print*,"et in verlet/nose/beeman setup",et
c               print*,"ept in verlet/nose/beeman setup",ept
c               print*,"ett in verlet/nose/beeman setup",ett
               call moments_onlyDipole

               allocate (derivs(3,n))

                do i = 1, n
                  do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i)
     &                 + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                 + dev(j,i) + dem(j,i)+dep3b(j,i)!+dep3bmut(j,i)
     &               !  + dep3b1(j,i)
     &        !      +dep3b_recip(j,i) 

                  end do
                 print*,"dep3b_x",i,dep3b(1,i)!+dep3bmut(1,i)+dep3b1(1,i)
     &         !          +dep3b_recip(1,i)
                 print*,"dep3b_y",i,dep3b(2,i)!+dep3bmut(2,i)+dep3b1(2,i)
     &          !         +dep3b_recip(2,i)
                 print*,"dep3b_z",i,dep3b(3,i)!+dep3bmut(3,i)+dep3b1(3,i)
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
                do i = 1, n
                  print*,"uind3b_x=",uind3b(1,i)
                  print*,"uind3b_y=",uind3b(2,i)
                  print*,"uind3b_z=",uind3b(3,i)
                end do      
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

c
c     set modified group variables for amoeba
c
      NGRP2=N/3

      
c      call mpi_bcast(delV,1,mpi_real8,master,
c     & mpi_comm_world,ierr)
      if(taskid.eq.master) then
         call settime2(twall,tcpu)
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  ISTOREM !length of MSD    
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  DTCALCM !time-step for MSD
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  ISKIPM !interval to calculate MSD
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  ISTOREV !length of VACF
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  DTCALCV !time-step for VACF
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  ISKIPV !interval to calculate VACF
      call nextarg (string,exist1)
      if (exist1)  read (string,*)  ISTORED !length IN NUM TIMESTEPS IF YOU PUT IN 5000, AND HAVE DTCALCD OF 2FS, THEN THIS CORRESPONDS TO 10ps for time correlation func.
      call nextarg (string,exist1)
      if (exist1)  read (string,*)  DTCALCD !time-step for DCF 2fs
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  ISKIPD !interval to calculate DCF
      call nextarg (string,exist1)
      if (exist1)  read (string,*)  SAVE_DATA_FREQ ! Freq to write the dcf and other production arrays to file (in #steps; so 1000) 
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  SAVE_REST_FREQ  ! Freq to write out coordinates and velocities for restart files.
c      ISKIPM=DTCALCM*INT(ISKIPM/DTCALCM)! make multiple of DTCALCM
c      ISKIPV=DTCALCV*INT(ISKIPV/DTCALCV)! make multiple of DTCALCV
c      ISKIPD=DTCALCD*INT(ISKIPD/DTCALCD)! make multiple of DTCALCD
c
c     read file names for writing output
c         
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE1 !name of dipole moment file
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE2 !name of MSD file
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE3 !name of O vacf file
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE4 !name of H vacf file
      call nextarg (string,exist1)
      if (exist1)  read (string,*)  FILE5 !name of dipole correlation file
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE6 !name of Rg file
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE7 !store file for X positions of Oxygens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE8 !store file for Y positions of Oxygens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE9 !store file for Z positions of Oxygens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE10 !store file for X velocities of Oxygens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE11 !store file for Y velocities of Oxygens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE12 !store file for Z velocities of Oxygens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE13 !store file for X velocities of Hydrogens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE14 !store file for Y velocities of Hydrogens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE15 !store file for Z velocities of Hydrogens of waters
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE16 !store file for X dipole moments for layers
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE17 !store file for Y dipole moments for layers
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE18 !store file for Z dipole moments for layers
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE19 !restart file which stores the index of the time step to which the simulation has progressed
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE20 !File with the count values for the msd array.
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE21 !File with the count values for the vacf array.
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE22 !File with the count values for the dcf array.
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE23 !File with the summass vector.
c      call nextarg (string,exist1)
c      if (exist1)  read (string,*)  FILE24 !Test file to check if the right thing is being read in.
c     open output files
c     
c     INQUIRE(FILE=FILE1,EXIST=EXIST)
c        IF (EXIST) THEN
c           OPEN (UNIT=201,FILE=FILE1,STATUS='OLD',
c    &           ACTION='WRITE',POSITION='APPEND')
c        ELSE
c           OPEN (UNIT=201,FILE=FILE1,STATUS='UNKNOWN')
c        ENDIF
c     INQUIRE(FILE=FILE6,EXIST=EXIST)
c        IF (EXIST) THEN
c           OPEN (UNIT=206,FILE=FILE6,STATUS='OLD',
c    &           ACTION='WRITE',POSITION='APPEND')
c        ELSE
c           OPEN (UNIT=206,FILE=FILE6,STATUS='UNKNOWN')
c        ENDIF
c
c     initialize time-correlation functions
c
      print*,"Read in Thz files and parameters"
c      DO DTC=1,ISTOREM+1
c         MSDO(DTC)=0d0
c         COUNTM(DTC)=0d0
c      ENDDO
c      DO DTC=1,ISTOREV+1
c         VACFO(DTC)=0d0
c         VACFH(DTC)=0d0
c         COUNTV(DTC)=0d0
c      ENDDO
      DO DTC=1,ISTORED+1
         DCF(DTC)=0d0
         COUNTD(DTC)=0d0
      ENDDO
c
c     count water atoms 
c
      Nw=0
      do I=1,N
         Nw=Nw+1
         wlist(Nw)=I
      enddo
c
c     X/Y/Z = boxed water coord, XC/YC/ZC = unboxed water coords
c
      DO I=1,Nw,3
         II=wlist(I) ! Oxygen index
         DO J=0,2
            JJ=wlist(I)+J ! all water indices
            XC(JJ)=X(JJ)-DNINT((X(JJ)-X(II))/XBOX)*XBOX
            YC(JJ)=Y(JJ)-DNINT((Y(JJ)-Y(II))/YBOX)*YBOX
            ZC(JJ)=Z(JJ)-DNINT((Z(JJ)-Z(II))/ZBOX)*ZBOX
         ENDDO
      ENDDO
c
c     calc dipole moment
c
              call alloc3b_uind

      MXC=0.0d0
      MYC=0.0d0
      MZC=0.0d0
      DO II=1,npole
         K=ipole(II)
         MXC=MXC+X(K)*rpole(1,ii)+rpole(2,ii)+uind3b(1,ii)
         MYC=MYC+Y(K)*rpole(1,ii)+rpole(3,ii)+uind3b(2,ii)
         MZC=MZC+Z(K)*rpole(1,ii)+rpole(4,ii)+uind3b(3,ii)
      ENDDO
c     
c     save initial coords to tape
c
c      TT=ISTOREM
      TT=ISTORED ! MODIFIED FROM ORIGINAL INITIALIZATION!!!!!

c      DO I=1,Nw,3
c         K=wlist(I) ! Oxygen index
c         XO(I,TT)=XC(K)
c         YO(I,TT)=YC(K)
c         ZO(I,TT)=ZC(K)
c      add velocities to tape         
c         II=(I+2)/3
c         J=Nw/3+II ! I (J) - Hydrogen 1 (2)   
c         L=K+1 ! Hydrogen 1
c         M=K+2 ! Hydrogen 2
c         VXO(II,TT)=V(1,K)
c         VYO(II,TT)=V(2,K)
c         VZO(II,TT)=V(3,K)
c         VXH(II,TT)=V(1,L)
c         VYH(II,TT)=V(2,L)
c         VZH(II,TT)=V(3,L)
c         VXH(J,TT)=V(1,M)
c         VYH(J,TT)=V(2,M)
c         VZH(J,TT)=V(3,M)
c      ENDDO
c      add dipole moments to tape
ccc      DO I=1,N
ccc         MXI(I,TT)=MXIC(I)
ccc         MYI(I,TT)=MYIC(I)
ccc         MZI(I,TT)=MZIC(I)
ccc      ENDDO
      MXI(TT)=MXC
      MYI(TT)=MYC
      MZI(TT)=MZC
      summass=0d0
      countmass=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     If simulation is being restarted, the ISTEP value is read         c
c     in from the restarttimestep.txt file to start from that point     c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         RESTARTFLAG=0
c         RESTARTTIMESTEP=0
        print*,"Done initializing Thz stuff"
c        INQUIRE(FILE=FILE19,EXIST=EXIST)
c        IF (EXIST) THEN
c           OPEN (UNIT=219,FILE=FILE19, STATUS='OLD', ACTION='READ')
c              READ(219,*) RESTARTTIMESTEP
c              RESTARTFLAG=1
c           CLOSE(219)
c        ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                 c
c     All the arrays that have been written out are read back into memory         c
c     These include: MSD, H_vacf, O_vacf, M, M2, dcf, dcf2,                       c
c     X0, Y0, Z0, VXO, VYO, VZO, VXH, VYH, VZH, MXI, MYI, MZI, MXI2, MYI2, MZI2   c
c                                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This reads in the MSD, O_VACF, H_vacf and dcf files and 
c       updates the MSD, O_vacf, H_vacf and dcf arrays, which are zero
c       since the simulation has been restarted.
c
c       INQUIRE(FILE=FILE2,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=202,FILE=FILE2, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=203,FILE=FILE3, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=204,FILE=FILE4, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=205,FILE=FILE5, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=220,FILE=FILE20, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=221,FILE=FILE21, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=222,FILE=FILE22, STATUS='OLD', ACTION='READ')
c             DO I=1,ISTOREV
c                READ(202,*) t_stp,MSDO(I)
c                READ(203,*) t_stp,VACFO(I)
c                READ(204,*) t_stp,VACFH(I)
c                READ(205,*) t_stp,DCF(I)
c                READ(220,*) t_stp,COUNTM(I)
c                READ(221,*) t_stp,COUNTV(I)
c                READ(222,*) t_stp,COUNTD(I)
c             ENDDO
c          CLOSE(202)
c       ENDIF
c
c       This file reads in the summass array.
c       INQUIRE(FILE=FILE23,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=223,FILE=FILE23, STATUS='OLD', ACTION='READ')
c             READ(223,*) summass
c          CLOSE(223)
c       ENDIF
c
c       This multiplies the values of the MSD, VACF_O, VACF_H and DCF read in
c       from the files by the appropriate COUNTM, COUNTV and COUNTD values to
c       get the actual values for those quantities
c
c       INQUIRE(FILE=FILE19,EXIST=EXIST)
c       IF (EXIST) THEN
c          DO I = 1, ISTOREV
c             MSDO(I)=MSDO(I)*COUNTM(I)
c             VACFO(I)=(VACFO(I)*COUNTV(I)*
c    &         (VACFO(1)/COUNTV(1)))
c             VACFH(I)=(VACFH(I)*COUNTV(I)*
c    &         (VACFH(1)/COUNTV(1)))
c             DCF(I)=DCF(I)*COUNTD(I)*summass
c          ENDDO
c       ENDIF

c
c       This reads in the Oxygen X, Y, Z coordinates and velocities file,
c       written out from the previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE7,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=207,FILE=FILE7, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=208,FILE=FILE8, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=209,FILE=FILE9, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=210,FILE=FILE10, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=211,FILE=FILE11, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=212,FILE=FILE12, STATUS='OLD', ACTION='READ')
c             IF(RESTARTTIMESTEP.LT.ISTOREV) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREV
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(207,*) (XO(O_IDX,DTC),O_IDX=1,Nw)
c                READ(208,*) (YO(O_IDX,DTC),O_IDX=1,Nw)
c                READ(209,*) (ZO(O_IDX,DTC),O_IDX=1,Nw)
c                READ(210,*) (VXO(O_IDX,DTC),O_IDX=1,Nw/3)
c                READ(211,*) (VYO(O_IDX,DTC),O_IDX=1,Nw/3)
c                READ(212,*) (VZO(O_IDX,DTC),O_IDX=1,Nw/3)
c             ENDDO
c          CLOSE(207)
c          CLOSE(208)
c          CLOSE(209)
c          CLOSE(210)
c          CLOSE(211)
c          CLOSE(212)
c       ENDIF
c
c       This reads in the Hydrogen X, Y, Z velocities file, written out from the 
c       previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE13,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=213,FILE=FILE13, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=214,FILE=FILE14, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=215,FILE=FILE15, STATUS='OLD', ACTION='READ')
c             IF (RESTARTTIMESTEP.LT.ISTOREV) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREV
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(213,*) (VXH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c                READ(214,*) (VYH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c                READ(215,*) (VZH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c             ENDDO
c          CLOSE(213)
c          CLOSE(214)
c          CLOSE(215)
c       ENDIF
c
c       This reads in the X, Y, Z dipole moments file, written out from the 
c       previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE13,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=213,FILE=FILE13, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=214,FILE=FILE14, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=215,FILE=FILE15, STATUS='OLD', ACTION='READ')
c             IF (RESTARTTIMESTEP.LT.ISTOREV) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREV
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(213,*) (VXH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c                READ(214,*) (VYH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c                READ(215,*) (VZH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c             ENDDO
c          CLOSE(213)
c          CLOSE(214)
c          CLOSE(215)
c       ENDIF
c
c       This reads in the X, Y, Z dipole moments file, written out from the 
c       previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE16,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=216,FILE=FILE16, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=217,FILE=FILE17, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=218,FILE=FILE18, STATUS='OLD', ACTION='READ')
c             IF (RESTARTTIMESTEP.LT.ISTOREM) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREM
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(216,*) (MXI(DTC))
c                READ(217,*) (MYI(DTC))
c                READ(218,*) (MZI(DTC))
c             ENDDO
c          CLOSE(216)
c          CLOSE(217)
c          CLOSE(218)
c       ENDIF
ccc       call cpu_time(finishTime)
ccc       write(70,*) 'Time for initialization =',finishTime-startTime
ccc       CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc       nb_ticks = nb_ticks_final - nb_ticks_initial
ccc       IF (nb_ticks_final < nb_ticks_initial) 
ccc    &           nb_ticks = nb_ticks + nb_ticks_max
ccc       elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc       write(80,*) 'Wall time for initialization =',elapsed_time
c
c
c     This is where the simulation actually starts
c     integrate equations of motion to take a time step
c
      end if

      
         if (integrate .eq. 'VERLET') then
          do istep0 = 1, nstep           
            if(taskid.eq.master) then 
             print*,"In master, Just beginning VERLET"
             call verlet_pt1 (istep0,dt)
             !call ptest_parallel

             ! x(6)=x(6)+delx
c              y(6)=y(6)+delx
c              z(6)=z(6)+delx
             call gradient_setup
             call allocPermElec1

             call alloc3b_vdwlist_nolist_vac
              call alloc3b_uind
            end if
            
             call alloc_erecip

            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call prep_pole
            call bspline_fill_omp
            call table_fill_omp
            call mpi_barrier(mpi_comm_world,ierr)
            !call bcast_lattice
           if((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2)))then
               call vlist_par
              if(firstload_vdw.and.(taskid.eq.numtasks_emreal+2)) then
                call vdw_load_balance_sizelist_half
              end if
           else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2))then
               call mlist_par
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
                if(firstload_emreal.and.(taskid.eq.master2)) then
c                 firstload_emreal = .false.
                   call emreal_load_balance_sizelist
                end if
           !else if(taskid.eq.master.and.(.not.done_clustlist_setup))then
           else if(taskid.eq.master.and.((mod(istep0,10).eq.0)
     &                            .or.exist))then
              if(do2waterclustlist) then
                 call mollist2body_NewBodieskmeans
              !end if
              !if(firstload_polz.and.do2waterclustlist) then
              !  firstload_polz = .false.
                !call polar_load_balance7
                !call polar_load_balance_clust
                call polar_load_balance_clust_splitlist_offset
                 call polar_load_balance_clust_splitlist_offset3b_simple
              end if
              !done_clustlist_setup = .true.
           end if

            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)

            if(listbcast.and.(approxmode.ne.'1BODYMODE')
     &                    .and.do2waterclustlist) then
             firstload_polz = .false.
           !  call ewbcast_polar_load_balance
           !  call bcast_polar_load_balance_clust
             call bcast_polar_load_balance_clust_splitlist

             call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
             call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
             !if(approxmode.eq.'3BODYMODE') then
               call bcast_polar_load_balance_clust_splitlist3b
               call mpi_bcast(nmollst3,clustcount,mpi_integer,master,
     &              mpi_comm_world,ierr)
               call mpi_bcast(mollst3,5*clustcount*clustcount,
     &              mpi_integer,master,mpi_comm_world,ierr)
               call mpi_bcast(save2bfor3b,clustcount*clustcount,
     &              mpi_logical,master,mpi_comm_world,ierr)
              call mpi_bcast(num2bsave,1,mpi_integer,master,
     &              mpi_comm_world,ierr)
              call mpi_bcast(ind2b_counter,clustcount*clustcount,
     &              mpi_integer,master,mpi_comm_world,ierr)
c               if(.not.allocated(counter2b_ind))
c     &             allocate(counter2b_ind(2,num2bsave))
c              call mpi_bcast(counter2b_ind,num2bsave*2,mpi_integer,
c     &           master,mpi_comm_world,ierr)

             !end if
              !listbcast=.false.
            end if

            if(.not.done_comm_create) then
              numtasks_polar1=clustcount
              offset=int(clustcount/numtasks_polar1)
              remainder=mod(clustcount,numtasks_polar1)

              if(taskid.le.remainder-1) then
               moli1rmndr=numtasks_polar1*offset+taskid+1
              else if (taskid.gt.remainder-1) then
               moli1rmndr=0
              end if
              start=taskid*offset+1


           allocate (taskidemreal(0:numtasks_polar1-1))
           do j=0,numtasks_polar1-1
            taskidemreal(j)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
            call mpi_group_incl(orig_group,numtasks_polar1,taskidemreal,
     &        emreal_group,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group,
     &            emreal_comm,ierr)

            deallocate(taskidemreal)

           allocate(taskidvdw(0:numtasks_vdw2))
           taskidvdw(0)=0
           do j=numtasks_emreal+2,numtasks_vdw2+numtasks_emreal+1
             k=j-(numtasks_emreal+2)+1
             taskidvdw(k)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
           call mpi_group_incl(orig_group,numtasks_vdw2+1,taskidvdw,
     &        vdw_group,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,vdw_group,
     &            vdw_comm,ierr)
           deallocate(taskidvdw)
             done_comm_create=.true.
            end if

            if(listbcast) then
              if((istep0.ne.1).and.(taskid.ge.numtasks_polar1).and.
     &            (taskid.lt.(numtasks_polar_old+numtasks_polar1)))then
                call mpi_comm_free(emreal_comm2,ierr)
              end if

           allocate (taskidemreal(0:numtasks_polar))
           taskidemreal(0)=0
           do j=numtasks_polar1,numtasks_polar+numtasks_polar1-1
              k=j-numtasks_polar1+1
              taskidemreal(k)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
         call mpi_group_incl(orig_group,numtasks_polar+1,taskidemreal,
     &        emreal_group2,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group2,
     &            emreal_comm2,ierr)
            deallocate(taskidemreal)

              if((istep0.ne.1).and.
     &           (taskid.ge.(numtasks_polar_old+numtasks_polar1)).and.
     &         (taskid.lt.
     &        (numtasks_polar_old+numtasks_polar1+numtasks_polar3b_old)
     &          )) then
                call mpi_comm_free(emreal_comm3,ierr)
              end if
            allocate (taskidemreal(0:numtasks_polar3b))
            taskidemreal(0)=0
            do j=numtasks_polar+numtasks_polar1,
     &            numtasks_polar+numtasks_polar1+numtasks_polar3b-1
               k=j-(numtasks_polar+numtasks_polar1)+1
               taskidemreal(k)=j
            end do
         call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
        call mpi_group_incl(orig_group,numtasks_polar3b+1,taskidemreal,
     &        emreal_group3,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group3,
     &            emreal_comm3,ierr)
            deallocate(taskidemreal)
            end if

          if(inputkmeansclust.and.(taskid.eq.master)) then
            ! call init1b2bmatPolarPerm
             call init1bmat!PolarPerm
             call init1bmat_uind
                if(.not. allocated(ep2bmatmod))
     &           allocate(ep2bmatmod(num2bsave))
               if(.not. allocated(dep2bmatmod)) allocate(dep2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
            if(.not. allocated(virep2bmatmod)) allocate(virep2bmatmod(3,
     &          3,num2bsave))
        if(.not. allocated(uind2bmatmod)) allocate(uind2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
             call init2bmatmodPolar
             call init2bmatmodPolar_uind
          end if

             if((taskid.lt.numtasks_polar1).and.mlistclust) then
              if(use_pred.and.(.not.setup_pred1b)) then
                call initpolarpred1b(start,start+offset-1,moli1rmndr)
                setup_pred1b=.true.
              end if
             end if

             if((taskid.ge.numtasks_polar1).and.
     &             (taskid.lt.(numtasks_polar+numtasks_polar1))
     &             .and.mlistclust) then
              !if(use_pred.and.(.not.setup_pred2b)) then
              if(use_pred.and.listbcast) then
                call initpolarpred2b_offset
                setup_pred2b=.true.
              end if
             end if

            if( (taskid.ge.(numtasks_polar+numtasks_polar1)).and.
     &         (taskid.lt.
     &         (numtasks_polar+numtasks_polar1+numtasks_polar3b))
     &         .and.mlistclust) then
              !if(use_pred.and.(.not.setup_pred3b)) then
              if(use_pred.and.listbcast) then
                call initpolarpred3b_offset
                setup_pred3b=.true.
            !  print*,"Completed initpolarpred3b_offset"
              end if
            end if

           if( (mod(istep0,5).eq.0 ) .or. exist) then
             if((taskid.lt.numtasks_polar1).and.mlistclust) then
              call clust1blists(start,start+offset-1, moli1rmndr)
c              print*,"Completed clust1blists"
             end if

             if((taskid.ge.numtasks_polar1).and.
     &             (taskid.lt.(numtasks_polar+numtasks_polar1))
     &             .and.mlistclust) then
              call clust2blists_offset
c              print*,"Completed clust2blists_offset"
             end if


            if( (taskid.ge.(numtasks_polar+numtasks_polar1)).and.
     &         (taskid.lt.
     &         (numtasks_polar+numtasks_polar1+numtasks_polar3b))
     &         .and.mlistclust) then
              call clust3blists_offset
c              print*,"Completed clust3blists_offset"
            end if
           end if
             listbcast=.false.
             exist=.false.
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

            if(firstload_vdw) then
              firstload_vdw = .false.
              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
     &         numtasks_emreal+2,mpi_comm_world,ierr)
              call mpi_bcast(start_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
               call mpi_bcast(last_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
            end if
            numtasks_polar_old=numtasks_polar
            numtasks_polar3b_old=numtasks_polar3b

          call uindclst_gradientPolar1b2b3bsimult_ireducematcommsmall2(
     &       start,start+offset-1, moli1rmndr)

            !call grad_cov_v_reduce
c            call grad_cov_v_reduce_commsmall
            call grad_covemrecip_vandr_reduce_totfield_commsmall

            if(taskid.eq.master1) then
             call mpi_wait(reqs1,stat,ierr)
             call mpi_wait(reqs2,stat,ierr)
             call mpi_wait(reqs3,stat,ierr)
            else if(taskid.eq.master) then
             call mpi_wait(reqs4,stat,ierr)
             call mpi_wait(reqs5,stat,ierr)
             call mpi_wait(reqs6,stat,ierr)
            end if
             call mpi_wait(reqs13,stat,ierr)
             call mpi_wait(reqs14,stat,ierr)
             call mpi_wait(reqs15,stat,ierr)
             if((taskid.eq.master).or.((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) ) then
             call mpi_wait(reqs16,stat,ierr)
             call mpi_wait(reqs17,stat,ierr)
             call mpi_wait(reqs18,stat,ierr)
             end if

       if(taskid.lt.numtasks_polar1) then
             call mpi_wait(reqs19,stat,ierr)
             call mpi_wait(reqs20,stat,ierr)
             call mpi_wait(reqs21,stat,ierr)
             call mpi_wait(reqs22,stat,ierr)
             call mpi_wait(reqs23,stat,ierr)
             call mpi_wait(reqs24,stat,ierr)
       end if

      if( (taskid.eq.master) .or. ( (taskid.ge.numtasks_polar1)
     &   .and.(taskid.lt.(numtasks_polar+numtasks_polar1))) ) then
             call mpi_wait(reqs25,stat,ierr)
             call mpi_wait(reqs26,stat,ierr)
             call mpi_wait(reqs27,stat,ierr)
             call mpi_wait(reqs28,stat,ierr)
             call mpi_wait(reqs29,stat,ierr)
             call mpi_wait(reqs30,stat,ierr)
      end if

      if( (taskid.eq.master) .or. (
     & (taskid.ge.(numtasks_polar+numtasks_polar1))
     &    .and. (taskid.lt.(numtasks_polar+numtasks_polar1
     &         +numtasks_polar3b))) ) then
             call mpi_wait(reqs31,stat,ierr)
             call mpi_wait(reqs32,stat,ierr)
             call mpi_wait(reqs33,stat,ierr)
      end if
      call uindsum_gradient_Polar1b2b3bsimult_ireducematcommsmall2 


            if(taskid.eq.master) then
                 call empole1d_3b_Perm_selfeng_bcast
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
     &                  + viremrecip(j,i) + virev(j,i)!+virep3bmut(j,i)
     &             !  +virep3b1(j,i)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b!+ep3bmut
c      call system_clock(t1,clock_rate)              
             allocate (derivs(3,n))
              !efindif_ep3b(istep)=ep3b!+ep3bmut
              !efindif(istep)=energy-ep3b!-em!-ep3bmut-em-ev
c              efindif_em(istep)=em
c!$OMP PARALLEL DO default(shared) private(i,j)              
              do i = 1, n
                do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i) 
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i) !+dep3bmut(j,i)
     &            !       + dep3b1(j,i)
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

c              analyticgx(istep)=derivs(1,6)-dep3b(1,6)-dem(1,6)!-(dep3bmut(1,6)
     &              !+dep3b_recip(1,6)+dep3b1(1,6)+dem(1,6)+dev(1,6))
c              analyticg_dep3bx(istep)=dep3b(1,6)!+dep3bmut(1,6)
     &             ! +dep3b1(1,6)
c              analyticg_demx(istep)=dem(1,6)

c              analyticgy(istep)=derivs(2,6)-dep3b(2,6)-dem(2,6)!-(dep3bmut(2,6)
     &            !  +dep3b_recip(2,6)+dep3b1(2,6)+dem(2,6)+dev(2,6))
c              analyticg_dep3by(istep)=dep3b(2,6)!+dep3bmut(2,6)
     &            !  +dep3b1(2,6)
c              analyticg_demy(istep)=dem(2,6)

c              analyticgz(istep)=derivs(3,6)-dep3b(3,6)-dem(3,6)!-(dep3bmut(3,6)
     &            !  +dep3b_recip(3,6)+dep3b1(3,6)+dem(3,6)+dev(3,6))
c              analyticg_dep3bz(istep)=dep3b(3,6)!+dep3bmut(3,6)
     &           !   +dep3b1(3,6)
c              analyticg_demz(istep)=dem(3,6)


c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*,"Timing Verlet fullstep v",tot_time

       print*," istep Vdw, PermElec, Pol",istep0,ev,em,ep3b!+ep3bmut

              deallocate (derivs)
               call verlet_pt2(istep0,dt,energy)

c              if(istep.eq.3) then
c                print*,"Finite Diff VdwandCovEng, VdwandCovGradz",
c     & (efindif(3)-efindif(1))/2.0d0/delx,analyticgy(2)
c                print*,"Finite Diff PermElecEng, PermElecGradz",
c     & (efindif_em(3)-efindif_em(1))/2.0d0/delx,analyticg_demy(2)
c                print*,"Finite Diff PolEng, PolGradz",
c     & (efindif_ep3b(3)-efindif_ep3b(1))/2.0d0/delx,analyticg_dep3by(2)
c              end if


c              if(istep.eq.2) then
c                do i=1,3
c                   do j=1,3
c                 vir_analytic(j,i)=vir(j,i)-virep3b(j,i)!-virep3bmut(j,i)
c                 virep3b_analytic(j,i)=virep3b(j,i)!+virep3bmut(j,i)
c                print*,"Step 2 virep3b",j,i,virep3b(j,i)!+virep3bmut(j,i)
c                   end do
c                end do
c                volbox_analytic=volbox
c              end if
c              if(istep.eq.3) then
c      dedv_vir = (vir_analytic(1,1)+vir_analytic(2,2)+vir_analytic(3,3))
c     &         / (3.0d0*volbox)
c              dedv_fd = (efindif(3)-efindif(1))
c     &        / (2.0d0*delta*volbox)
cc
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
ccc        call cpu_time(finishTime)
ccc        write(71,*) 'Time for integration step',istep0,
ccc    &               '=',finishTime-startTime
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc        nb_ticks = nb_ticks_final - nb_ticks_initial
ccc        IF (nb_ticks_final < nb_ticks_initial) 
ccc    &            nb_ticks = nb_ticks + nb_ticks_max
ccc        elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc        write(81,*) 'Wall time for integration step',istep0,
ccc    &               ' =',elapsed_time
c     
c      Take coordinates out of the box for diffusion calculations
c
         DO I=1,N
            XC(I)=X(I)-DNINT((X(I)-XC(I))/XBOX)*XBOX
            YC(I)=Y(I)-DNINT((Y(I)-YC(I))/YBOX)*YBOX
            ZC(I)=Z(I)-DNINT((Z(I)-ZC(I))/ZBOX)*ZBOX
         ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     ISTEP is recaliberated if the simulation is being restarted       c
c     ISTEP=current_simulation_ISTEP_value+ISTEP value read in          c
c     from where previous simulation died.                              c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ISTEP=ISTEP0+RESTARTTIMESTEP
        ISTEP=ISTEP0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     The new modified value of ISTEP is written out to use for         c
c     subsequent restarts.                                              c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        IF(MOD(ISTEP,SAVE_REST_FREQ).EQ.0) THEN
c           OPEN (UNIT=219,FILE=FILE19, STATUS='REPLACE')
c              WRITE(219,'(I8)') ISTEP
c           CLOSE(219)
c        ENDIF
c
c        calculate dipole moment
c
ccc        call cpu_time(startTime)
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
         IF(MOD(ISTEP,DTCALCD).EQ.0) THEN
            MXC=0.0d0
            MYC=0.0d0
            MZC=0.0d0
            DO II=1,npole
                K=ipole(II)
                MXC=MXC+X(K)*rpole(1,ii)+rpole(2,ii)+uind3b(1,ii)
                MYC=MYC+Y(K)*rpole(1,ii)+rpole(3,ii)+uind3b(2,ii)
                MZC=MZC+Z(K)*rpole(1,ii)+rpole(4,ii)+uind3b(3,ii)
            ENDDO
c           WRITE(201,'(4ES16.8)') DT*DBLE(ISTEP),MXC,MYC,MZC
         ENDIF
ccc        call cpu_time(finishTime)
ccc        write(74,*) 'Time for dipole moment at step',istep0,
ccc    &               '=',finishTime-startTime
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc        nb_ticks = nb_ticks_final - nb_ticks_initial
ccc        IF (nb_ticks_final < nb_ticks_initial) 
ccc    &               nb_ticks = nb_ticks + nb_ticks_max
ccc        elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc        write(84,*) 'Wall time for dipole moment at step',istep0,
ccc    &               '=',elapsed_time
c          
c        calc msd & vacf and dcf
c        
c         DO I=1,Nw
c            K=wlist(I)
c            VX(K)=V(1,K)
c            VY(K)=V(2,K)
c            VZ(K)=V(3,K)
c         ENDDO

         ISTOREM=ISTORED  ! NEW
         ISKIPD=DTCALCD   ! NEW
         ISKIPM=ISKIPD    ! NEW
         
         DTCALCM=DTCALCD  ! NEW
         IF (MOD(ISTEP,ISKIPM).EQ.0.AND.ISTEP.GE.ISTOREM*DTCALCM) THEN
ccc           call cpu_time(startTime)
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
            DO ISTEP2=ISTEP-ISTOREM*DTCALCM,ISTEP,DTCALCM
               TT=MOD(ISTEP2/DTCALCM-1,ISTOREM)+1
               IF (TT.EQ.0) TT=ISTOREM
               DTC=(ISTEP-ISTEP2)/DTCALCM
c               tempMSDO=MSDO(DTC+1)
c               tempVACFO=VACFO(DTC+1)
c               tempVACFH=VACFH(DTC+1)
c               tempCOUNTM=COUNTM(DTC+1)
c               tempCOUNTV=COUNTV(DTC+1)
c!$OMP PARALLEL default(private) shared(wlist,
c!$OMP& TT,NW,XO,YO,ZO,XC,YC,ZC,VXO,VYO,VZO,VXH,VYH,VZH,VX,VY,VZ,ISTEP,
c!$OMP& ISTEP2)
c!$OMP& shared(tempMSDO,tempVACFO,tempVACFH,tempCOUNTM,tempCOUNTV)
c!$OMP DO reduction(+:tempMSDO,tempVACFO,tempVACFH,tempCOUNTM,tempCOUNTV)
c               DO I=1,Nw,3
c                  K=wlist(I) ! Oxygen
c                  IF (ISTEP2.LT.ISTEP) THEN
c                     DX=XO(I,TT)-XC(K)
c                     DY=YO(I,TT)-YC(K)
c                     DZ=ZO(I,TT)-ZC(K)
c                    write(*,*) I,TT,XO(I,TT)
c                  ELSE
c                     DX=0d0
c                     DY=0d0
c                     DZ=0d0
c                  ENDIF
c                  DR2=DX*DX+DY*DY+DZ*DZ
c                  tempMSDO=tempMSDO+DR2
c                  tempCOUNTM=tempCOUNTM+1
c
c             Calculate vacf
c
c                  II=(I+2)/3
c                  J=Nw/3+II ! I (J) - Hydrogen 1 (2)           
c                  L=K+1 ! Hydrogen 1
c                  M=K+2 ! Hydrogen 2
c                  IF (ISTEP2.LT.ISTEP) THEN
c                     privateVACFO=
c     +               VXO(II,TT)*VX(K)+VYO(II,TT)*VY(K)+VZO(II,TT)*VZ(K)
c                     privateVACFH=
c     +               VXH(II,TT)*VX(L)+VYH(II,TT)*VY(L)+VZH(II,TT)*VZ(L)
c     +               +VXH(J,TT)*VX(M)+VYH(J,TT)*VY(M)+VZH(J,TT)*VZ(M)
c                  ELSE
c                     privateVACFO=
c     +               VX(K)*VX(K)+VY(K)*VY(K)+VZ(K)*VZ(K)
c                     privateVACFH=
c     +               VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L)
c     +               +VX(M)*VX(M)+VY(M)*VY(M)+VZ(M)*VZ(M)
c                  ENDIF
c                  tempVACFO=tempVACFO+privateVACFO
c                  tempVACFH=tempVACFH+privateVACFH
c                  tempCOUNTV=tempCOUNTV+1
c               ENDDO
c!$OMP END DO
c!$OMP END PARALLEL
c               MSDO(DTC+1)=tempMSDO
c               VACFO(DTC+1)=tempVACFO
c               VACFH(DTC+1)=tempVACFH
c               COUNTM(DTC+1)=tempCOUNTM
c               COUNTV(DTC+1)=tempCOUNTV
c
c             Calculate DCF
c
ccc               DO II=1,N
ccc                  DCF(DTC+1)=DCF(DTC+1)+
ccc     &               MXI(II,TT)*MXC+MYI(II,TT)*MYC+MZI(II,TT)*MZC
ccc               ENDDO
               IF (ISTEP2.LT.ISTEP) THEN
                  DCF(DTC+1)=DCF(DTC+1)+MXI(TT)*MXC+MYI(TT)*MYC+
     &                       MZI(TT)*MZC
               ELSE
                  DCF(DTC+1)=DCF(DTC+1)+MXC*MXC+MYC*MYC+
     &                       MZC*MZC
               ENDIF
               COUNTD(DTC+1)=COUNTD(DTC+1)+1
               summass = 0d0
               do i=1,N
                  summass=summass+mass(i)
               enddo
               countmass=countmass+1
            ENDDO
ccc           call cpu_time(finishTime)
ccc           write(75,*) 'Time for calc msd, vacf, dcf at step',
ccc    &                  istep0,'=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(85,*) 'Wall time for calc msd, vacf, dcf',istep0,
ccc    &                  '=',elapsed_time
            IF (MOD(ISTEP,SAVE_DATA_FREQ).EQ.0) THEN
ccc           call cpu_time(startTime)
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c           OPEN(UNIT=202,FILE=FILE2,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=203,FILE=FILE3,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=204,FILE=FILE4,STATUS='REPLACE',ACTION='WRITE')
           OPEN(UNIT=205,FILE=FILE5,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=220,FILE=FILE20,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=221,FILE=FILE21,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=222,FILE=FILE22,STATUS='REPLACE',ACTION='WRITE')
c           DO DTC=0,ISTOREM
c              WRITE(202,'(2ES16.8)') DT*DBLE(DTC*DTCALCM),
c    &           MSDO(DTC+1)/DBLE(COUNTM(DTC+1))
c              WRITE(203,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             ((VACFO(DTC+1)/(DBLE(COUNTV(DTC+1))))
c    &             /(VACFO(1)/DBLE(COUNTV(1))))
c              WRITE(204,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             ((VACFH(DTC+1)/(DBLE(COUNTV(DTC+1))))
c    &             /(VACFH(1)/DBLE(COUNTV(1))))
               print*,"dt=",dt
               print*,"dtc=",dtc
               print*,"dtcald=",dtcalcd
               WRITE(205,'(2ES16.8)') DT*DBLE(DTC*DTCALCD),
     &             DCF(DTC+1)/DBLE(COUNTD(DTC+1))/summass
c              WRITE(220,'(2ES16.8)') DT*DBLE(DTC*DTCALCM),
c    &             DBLE(COUNTM(DTC+1))
c              WRITE(221,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             DBLE(COUNTV(DTC+1))
c              WRITE(222,'(2ES16.8)') DT*DBLE(DTC*DTCALCD),
c    &             DBLE(COUNTD(DTC+1))
c           ENDDO
c           CLOSE(202)
c           CLOSE(203)
c           CLOSE(204)
           CLOSE(205)
c           CLOSE(220)
c           CLOSE(221)
c           CLOSE(222)
c           OPEN(UNIT=223,FILE=FILE23,STATUS='REPLACE',ACTION='WRITE')
c               WRITE(223,'(1ES16.8)') summass
c           CLOSE(223)
ccc           call cpu_time(finishTime)
ccc           write(76,*) 'Time for write msd, vacf, dcf at step',
ccc    &                  istep0,'=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(86,*) 'Wall time for write msd, vacf, dcf',istep0,
ccc    &                  '=',elapsed_time
            ENDIF
         ENDIF
c        
c        add current positions to tape
c        
         IF(MOD(ISTEP,DTCALCM).EQ.0) THEN
            TT=MOD(ISTEP/DTCALCM-1,ISTOREM)+1
c            DO I=1,Nw,3
c               K=wlist(I) ! Oxygen index
c               XO(I,TT)=XC(K)
c               YO(I,TT)=YC(K)
c               ZO(I,TT)=ZC(K)
c              write(*,*) I,TT
c              write(*,*) XO(I,TT),YO(I,TT),ZO(I,TT)
c               II=(I+2)/3
c         Add velocities to tape               
c               J=Nw/3+II ! I (J) - Hydrogen 1 (2) 
c               L=K+1 ! Hydrogen 1
c               M=K+2 ! Hydrogen 2
c               VXO(II,TT)=VX(K)
c               VYO(II,TT)=VY(K)
c               VZO(II,TT)=VZ(K)
c               VXH(II,TT)=VX(L)
c               VYH(II,TT)=VY(L)
c               VZH(II,TT)=VZ(L)
c               VXH(J,TT)=VX(M)
c               VYH(J,TT)=VY(M)
c               VZH(J,TT)=VZ(M)
c            ENDDO
c         Add dipole moments to tape            
ccc            DO I=1,N
ccc               MXI(I,TT)=MXIC(I)
ccc               MYI(I,TT)=MYIC(I)
ccc               MZI(I,TT)=MZIC(I)
ccc            ENDDO
            MXI(TT)=MXC
            MYI(TT)=MYC
            MZI(TT)=MZC
         ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                               c
c     Write out tape: All the tape values are written out, so that they can be read in          c
c     to use for restarting.                                                                    c
c                                                                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c        Store files for the Oxygen positions
c     
c         IF(MOD(ISTEP,SAVE_REST_FREQ).EQ.0) THEN
ccc        call cpu_time(startTime)
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c           OPEN (UNIT=207,FILE=FILE7,STATUS='REPLACE') ! store file for X positions of Oxygens of waters
c           OPEN (UNIT=208,FILE=FILE8,STATUS='REPLACE') ! store file for Y positions of Oxygens of waters
c           OPEN (UNIT=209,FILE=FILE9,STATUS='REPLACE') ! store file for Z positions of Oxygens of waters
c           OPEN (UNIT=210,FILE=FILE10,STATUS='REPLACE') !store file for X velocities of Oxygens of waters
c           OPEN (UNIT=211,FILE=FILE11,STATUS='REPLACE') !store file for Y velocities of Oxygens of waters
c           OPEN (UNIT=212,FILE=FILE12,STATUS='REPLACE') !store file for Z velocities of Oxygens of waters
c           OPEN (UNIT=213,FILE=FILE13,STATUS='REPLACE') !store file for X velocities of Hydrogens of waters
c           OPEN (UNIT=214,FILE=FILE14,STATUS='REPLACE') !store file for Y velocities of Hydrogens of waters
c           OPEN (UNIT=215,FILE=FILE15,STATUS='REPLACE') !store file for Z velocities of Hydrogens of waters
c           OPEN (UNIT=216,FILE=FILE16,STATUS='REPLACE') !store file for X dipole moment
c           OPEN (UNIT=217,FILE=FILE17,STATUS='REPLACE') !store file for y dipole moment
c           OPEN (UNIT=218,FILE=FILE18,STATUS='REPLACE') !store file for Z dipole moment
c              IF (ISTEP.LT.ISTOREV*DTCALCV) THEN
c                 SAVE_LENGTH = ISTEP/DTCALCV
c              ELSE 
c                 SAVE_LENGTH = ISTOREV
c              END IF
c              DO DTC=1,SAVE_LENGTH
c                DO O_IDX=1,Nw
c                  WRITE(207,'(ES16.8,X)',advance='no') XO(O_IDX,DTC)
c                  WRITE(208,'(ES16.8,X)',advance='no') YO(O_IDX,DTC)
c                  WRITE(209,'(ES16.8,X)',advance='no') ZO(O_IDX,DTC)
c                ENDDO
c                DO O_IDX=1,Nw/3
c                  WRITE(210,'(ES16.8,X)',advance='no') VXO(O_IDX,DTC)
c                  WRITE(211,'(ES16.8,X)',advance='no') VYO(O_IDX,DTC)
c                  WRITE(212,'(ES16.8,X)',advance='no') VZO(O_IDX,DTC)
c                ENDDO
c                DO H_IDX=1,2*Nw/3
c                  WRITE(213,'(ES16.8,X)',advance='no') VXH(H_IDX,DTC)
c                  WRITE(214,'(ES16.8,X)',advance='no') VYH(H_IDX,DTC)
c                  WRITE(215,'(ES16.8,X)',advance='no') VZH(H_IDX,DTC)
c                ENDDO
c                WRITE(216,'(ES16.8)') MXI(DTC)
c                WRITE(217,'(ES16.8)') MYI(DTC)
c                WRITE(218,'(ES16.8)') MZI(DTC)
c                WRITE(207,'(/)',advance='no')
c                WRITE(208,'(/)',advance='no')
c                WRITE(209,'(/)',advance='no')
c                WRITE(210,'(/)',advance='no')
c                WRITE(211,'(/)',advance='no')
c                WRITE(212,'(/)',advance='no')
c                WRITE(213,'(/)',advance='no')
c                WRITE(214,'(/)',advance='no')
c                WRITE(215,'(/)',advance='no')
c              ENDDO
c           CLOSE(207)
c           CLOSE(208)
c           CLOSE(209)
c           CLOSE(210)
c           CLOSE(211)
c           CLOSE(212)
c           CLOSE(213)
c           CLOSE(214)
c           CLOSE(215)
c           CLOSE(216)
c           CLOSE(217)
c           CLOSE(218)
ccc           call cpu_time(finishTime)
ccc           write(77,*) 'Time for savefiles at step',istep0,
ccc    &                  '=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(87,*) 'Wall time for savefiles at step =',istep0,
ccc    &                  '=',elapsed_time
c         ENDIF
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
          do istep0=1,nstep
            if(taskid.eq.master) then
             call nose_pt1 (istep0,dt,press)
             !call ptest_parallel

             ! x(6)=x(6)+delx
c              y(6)=y(6)+delx
c              z(6)=z(6)+delx
             call gradient_setup
             call allocPermElec1

             call alloc3b_vdwlist_nolist_vac
              call alloc3b_uind
            end if

             call alloc_erecip

            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call prep_pole
            call bspline_fill_omp
            call table_fill_omp
            call mpi_barrier(mpi_comm_world,ierr)
            call bcast_lattice

           if((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2)))then
               call vlist_par
              if(firstload_vdw.and.(taskid.eq.numtasks_emreal+2)) then
                call vdw_load_balance_sizelist_half
              end if
           else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2))then
               call mlist_par
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
                if(firstload_emreal.and.(taskid.eq.master2)) then
c                 firstload_emreal = .false.
                   call emreal_load_balance_sizelist
                end if
           !else if(taskid.eq.master.and.(.not.done_clustlist_setup))then
           else if(taskid.eq.master.and.((mod(istep0,10).eq.0)
     &                            .or.exist))then
              if(do2waterclustlist) then
                 call mollist2body_NewBodieskmeans
              !end if
              !if(firstload_polz.and.do2waterclustlist) then
              !  firstload_polz = .false.
                !call polar_load_balance7
                !call polar_load_balance_clust
                call polar_load_balance_clust_splitlist_offset
                 call polar_load_balance_clust_splitlist_offset3b_simple
              end if
              !done_clustlist_setup = .true.
           end if

            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)

            if(listbcast.and.(approxmode.ne.'1BODYMODE')
     &                    .and.do2waterclustlist) then
             firstload_polz = .false.
           !  call ewbcast_polar_load_balance
           !  call bcast_polar_load_balance_clust
             call bcast_polar_load_balance_clust_splitlist

             call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
             call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
             !if(approxmode.eq.'3BODYMODE') then
               call bcast_polar_load_balance_clust_splitlist3b
               call mpi_bcast(nmollst3,clustcount,mpi_integer,master,
     &              mpi_comm_world,ierr)
               call mpi_bcast(mollst3,5*clustcount*clustcount,
     &              mpi_integer,master,mpi_comm_world,ierr)
               call mpi_bcast(save2bfor3b,clustcount*clustcount,
     &              mpi_logical,master,mpi_comm_world,ierr)
              call mpi_bcast(num2bsave,1,mpi_integer,master,
     &              mpi_comm_world,ierr)
              call mpi_bcast(ind2b_counter,clustcount*clustcount,
     &              mpi_integer,master,mpi_comm_world,ierr)
c               if(.not.allocated(counter2b_ind))
c     &             allocate(counter2b_ind(2,num2bsave))
c              call mpi_bcast(counter2b_ind,num2bsave*2,mpi_integer,
c     &           master,mpi_comm_world,ierr)

             !end if
              !listbcast=.false.
            end if

            if(.not.done_comm_create) then
              numtasks_polar1=clustcount
              offset=int(clustcount/numtasks_polar1)
              remainder=mod(clustcount,numtasks_polar1)

              if(taskid.le.remainder-1) then
               moli1rmndr=numtasks_polar1*offset+taskid+1
              else if (taskid.gt.remainder-1) then
               moli1rmndr=0
              end if
              start=taskid*offset+1


           allocate (taskidemreal(0:numtasks_polar1-1))
           do j=0,numtasks_polar1-1
            taskidemreal(j)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
            call mpi_group_incl(orig_group,numtasks_polar1,taskidemreal,
     &        emreal_group,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group,
     &            emreal_comm,ierr)

            deallocate(taskidemreal)

           allocate(taskidvdw(0:numtasks_vdw2))
           taskidvdw(0)=0
           do j=numtasks_emreal+2,numtasks_vdw2+numtasks_emreal+1
             k=j-(numtasks_emreal+2)+1
             taskidvdw(k)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
           call mpi_group_incl(orig_group,numtasks_vdw2+1,taskidvdw,
     &        vdw_group,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,vdw_group,
     &            vdw_comm,ierr)
           deallocate(taskidvdw)
             done_comm_create=.true.

           allocate (taskidemreal(0:numtasks_polar))
           taskidemreal(0)=0
           do j=numtasks_polar1,numtasks_polar+numtasks_polar1-1
              k=j-numtasks_polar1+1
              taskidemreal(k)=j
           end do
           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
         call mpi_group_incl(orig_group,numtasks_polar+1,taskidemreal,
     &        emreal_group2,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group2,
     &            emreal_comm2,ierr)
            deallocate(taskidemreal)

            allocate (taskidemreal(0:numtasks_polar3b))
            taskidemreal(0)=0
            do j=numtasks_polar+numtasks_polar1,
     &            numtasks_polar+numtasks_polar1+numtasks_polar3b-1
               k=j-(numtasks_polar+numtasks_polar1)+1
               taskidemreal(k)=j
            end do
         call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
        call mpi_group_incl(orig_group,numtasks_polar3b+1,taskidemreal,
     &        emreal_group3,ierr)
           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group3,
     &            emreal_comm3,ierr)
            deallocate(taskidemreal)
            end if 

c   Commented out below due to error msg of too many communicators
c            if(listbcast) then
c              if((istep0.ne.1).and.(taskid.ge.numtasks_polar1).and.
c     &            (taskid.lt.(numtasks_polar_old+numtasks_polar1)))then
c                call mpi_comm_free(emreal_comm2,ierr)
c              end if
c
c           allocate (taskidemreal(0:numtasks_polar))
c           taskidemreal(0)=0
c           do j=numtasks_polar1,numtasks_polar+numtasks_polar1-1
c              k=j-numtasks_polar1+1
c              taskidemreal(k)=j
c           end do
c           call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
c         call mpi_group_incl(orig_group,numtasks_polar+1,taskidemreal,
c     &        emreal_group2,ierr)
c           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group2,
c     &            emreal_comm2,ierr)
c            deallocate(taskidemreal)
c
c              if((istep0.ne.1).and.
c     &           (taskid.ge.(numtasks_polar_old+numtasks_polar1)).and.
c     &         (taskid.lt.
c     &        (numtasks_polar_old+numtasks_polar1+numtasks_polar3b_old)
c     &          )) then
c                call mpi_comm_free(emreal_comm3,ierr)
c              end if
c            allocate (taskidemreal(0:numtasks_polar3b))
c            taskidemreal(0)=0
c            do j=numtasks_polar+numtasks_polar1,
c     &            numtasks_polar+numtasks_polar1+numtasks_polar3b-1
c               k=j-(numtasks_polar+numtasks_polar1)+1
c               taskidemreal(k)=j
c            end do
c         call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)
c        call mpi_group_incl(orig_group,numtasks_polar3b+1,taskidemreal,
c     &        emreal_group3,ierr)
c           call MPI_COMM_CREATE(MPI_COMM_WORLD,emreal_group3,
c     &            emreal_comm3,ierr)
c            deallocate(taskidemreal)
c            end if

          if(inputkmeansclust.and.(taskid.eq.master)) then
            ! call init1b2bmatPolarPerm
             call init1bmat!PolarPerm
             call init1bmat_uind
                if(.not. allocated(ep2bmatmod))
     &           allocate(ep2bmatmod(num2bsave))
               if(.not. allocated(dep2bmatmod)) allocate(dep2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
            if(.not. allocated(virep2bmatmod)) allocate(virep2bmatmod(3,
     &          3,num2bsave))
        if(.not. allocated(uind2bmatmod)) allocate(uind2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
             call init2bmatmodPolar
             call init2bmatmodPolar_uind

          end if

             if((taskid.lt.numtasks_polar1).and.mlistclust) then
              if(use_pred.and.(.not.setup_pred1b)) then
                call initpolarpred1b(start,start+offset-1,moli1rmndr)
                setup_pred1b=.true.
              end if
             end if 

             if((taskid.ge.numtasks_polar1).and.
     &             (taskid.lt.(numtasks_polar+numtasks_polar1))
     &             .and.mlistclust) then
              !if(use_pred.and.(.not.setup_pred2b)) then
              if(use_pred.and.listbcast) then
                call initpolarpred2b_offset
                setup_pred2b=.true.
              end if
             end if

            if( (taskid.ge.(numtasks_polar+numtasks_polar1)).and.
     &         (taskid.lt.
     &         (numtasks_polar+numtasks_polar1+numtasks_polar3b))
     &         .and.mlistclust) then
              !if(use_pred.and.(.not.setup_pred3b)) then
              if(use_pred.and.listbcast) then
                call initpolarpred3b_offset
                setup_pred3b=.true.
            !  print*,"Completed initpolarpred3b_offset"
              end if
            end if

           if( (mod(istep0,5).eq.0 ) .or. exist) then
             if((taskid.lt.numtasks_polar1).and.mlistclust) then
              call clust1blists(start,start+offset-1, moli1rmndr)
c              print*,"Completed clust1blists"
             end if

             if((taskid.ge.numtasks_polar1).and.
     &             (taskid.lt.(numtasks_polar+numtasks_polar1))
     &             .and.mlistclust) then
              call clust2blists_offset
c              print*,"Completed clust2blists_offset"
             end if


            if( (taskid.ge.(numtasks_polar+numtasks_polar1)).and.
     &         (taskid.lt.
     &         (numtasks_polar+numtasks_polar1+numtasks_polar3b))
     &         .and.mlistclust) then
              call clust3blists_offset
c              print*,"Completed clust3blists_offset"
            end if
           end if
             listbcast=.false.
             exist=.false.
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

            if(firstload_vdw) then
              firstload_vdw = .false.
              call mpi_bcast(numtasks_vdw2,1,mpi_integer,
     &         numtasks_emreal+2,mpi_comm_world,ierr)
              call mpi_bcast(start_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
               call mpi_bcast(last_vdw2,numtasks_vdw,mpi_integer,
     &          numtasks_emreal+2,mpi_comm_world,ierr)
            end if
            numtasks_polar_old=numtasks_polar
            numtasks_polar3b_old=numtasks_polar3b

          call uindclst_gradientPolar1b2b3bsimult_ireducematcommsmall2(
     &       start,start+offset-1, moli1rmndr)

            !call grad_cov_v_reduce
c            call grad_cov_v_reduce_commsmall
            call grad_covemrecip_vandr_reduce_totfield_commsmall

            if(taskid.eq.master1) then
             call mpi_wait(reqs1,stat,ierr)
             call mpi_wait(reqs2,stat,ierr)
             call mpi_wait(reqs3,stat,ierr)
            else if(taskid.eq.master) then
             call mpi_wait(reqs4,stat,ierr)
             call mpi_wait(reqs5,stat,ierr)
             call mpi_wait(reqs6,stat,ierr)
            end if
             call mpi_wait(reqs13,stat,ierr)
             call mpi_wait(reqs14,stat,ierr)
             call mpi_wait(reqs15,stat,ierr)
             if((taskid.eq.master).or.((taskid.ge.numtasks_emreal+2)
     &       .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) ) then
             call mpi_wait(reqs16,stat,ierr)
             call mpi_wait(reqs17,stat,ierr)
             call mpi_wait(reqs18,stat,ierr)
             end if

       if(taskid.lt.numtasks_polar1) then
             call mpi_wait(reqs19,stat,ierr)
             call mpi_wait(reqs20,stat,ierr)
             call mpi_wait(reqs21,stat,ierr)
             call mpi_wait(reqs22,stat,ierr)
             call mpi_wait(reqs23,stat,ierr)
             call mpi_wait(reqs24,stat,ierr)
       end if

      if( (taskid.eq.master) .or. ( (taskid.ge.numtasks_polar1)
     &   .and.(taskid.lt.(numtasks_polar+numtasks_polar1))) ) then
             call mpi_wait(reqs25,stat,ierr)
             call mpi_wait(reqs26,stat,ierr)
             call mpi_wait(reqs27,stat,ierr)
             call mpi_wait(reqs28,stat,ierr)
             call mpi_wait(reqs29,stat,ierr)
             call mpi_wait(reqs30,stat,ierr)
      end if

      if( (taskid.eq.master) .or. (
     & (taskid.ge.(numtasks_polar+numtasks_polar1))
     &    .and. (taskid.lt.(numtasks_polar+numtasks_polar1
     &         +numtasks_polar3b))) ) then
             call mpi_wait(reqs31,stat,ierr)
             call mpi_wait(reqs32,stat,ierr)
             call mpi_wait(reqs33,stat,ierr)
      end if

      call uindsum_gradient_Polar1b2b3bsimult_ireducematcommsmall2


            if(taskid.eq.master) then
                 call empole1d_3b_Perm_selfeng_bcast
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
     &                  + viremrecip(j,i) + virev(j,i)!+virep3bmut(j,i)
     &             !  +virep3b1(j,i)
                   end do
                 end do

                energy = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b!+ep3bmut
c      call system_clock(t1,clock_rate)              
             allocate (derivs(3,n))
              !efindif_ep3b(istep)=ep3b!+ep3bmut
              !efindif(istep)=energy-ep3b!-em!-ep3bmut-em-ev
c              efindif_em(istep)=em
c!$OMP PARALLEL DO default(shared) private(i,j)              
              do i = 1, n
                do j = 1, 3
               derivs(j,i) =deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i)
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i) !+dep3bmut(j,i)
     &            !       + dep3b1(j,i)
                end do
               if (use(i)) then
                do j = 1, 3
                 a(j,i) = -convert * derivs(j,i) / mass(i)
                 v(j,i) = v(j,i) + a(j,i)*dt_2
c                 print*,"istep i j v(j,i)",istep,i,j,v(j,i)
                end do
               end if
              end do

       print*,"Nose_istep Vdw, PermElec, Pol",istep0,ev,em,ep3b!+ep3bmut

              deallocate (derivs)
               call nose_pt2(istep0,dt,press,energy)

         DO I=1,N
            XC(I)=X(I)-DNINT((X(I)-XC(I))/XBOX)*XBOX
            YC(I)=Y(I)-DNINT((Y(I)-YC(I))/YBOX)*YBOX
            ZC(I)=Z(I)-DNINT((Z(I)-ZC(I))/ZBOX)*ZBOX
         ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     ISTEP is recaliberated if the simulation is being restarted       c
c     ISTEP=current_simulation_ISTEP_value+ISTEP value read in          c
c     from where previous simulation died.                              c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ISTEP=ISTEP0+RESTARTTIMESTEP
        ISTEP=ISTEP0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     The new modified value of ISTEP is written out to use for         c
c     subsequent restarts.                                              c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        IF(MOD(ISTEP,SAVE_REST_FREQ).EQ.0) THEN
c           OPEN (UNIT=219,FILE=FILE19, STATUS='REPLACE')
c              WRITE(219,'(I8)') ISTEP
c           CLOSE(219)
c        ENDIF
c
c        calculate dipole moment
c
ccc        call cpu_time(startTime)
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
         IF(MOD(ISTEP,DTCALCD).EQ.0) THEN
            MXC=0.0d0
            MYC=0.0d0
            MZC=0.0d0
            DO II=1,npole
                K=ipole(II)
                MXC=MXC+X(K)*rpole(1,ii)+rpole(2,ii)+uind3b(1,ii)
                MYC=MYC+Y(K)*rpole(1,ii)+rpole(3,ii)+uind3b(2,ii)
                MZC=MZC+Z(K)*rpole(1,ii)+rpole(4,ii)+uind3b(3,ii)
            ENDDO
c           WRITE(201,'(4ES16.8)') DT*DBLE(ISTEP),MXC,MYC,MZC
         ENDIF
ccc        call cpu_time(finishTime)
ccc        write(74,*) 'Time for dipole moment at step',istep0,
ccc    &               '=',finishTime-startTime
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc        nb_ticks = nb_ticks_final - nb_ticks_initial
ccc        IF (nb_ticks_final < nb_ticks_initial) 
ccc    &               nb_ticks = nb_ticks + nb_ticks_max
ccc        elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc        write(84,*) 'Wall time for dipole moment at step',istep0,
ccc    &               '=',elapsed_time
c          
c        calc msd & vacf and dcf
c        
c         DO I=1,Nw
c            K=wlist(I)
c            VX(K)=V(1,K)
c            VY(K)=V(2,K)
c            VZ(K)=V(3,K)
c         ENDDO

         ISTOREM=ISTORED  ! NEW
         ISKIPD=DTCALCD   ! NEW
         ISKIPM=ISKIPD    ! NEW

         DTCALCM=DTCALCD  ! NEW
         IF (MOD(ISTEP,ISKIPM).EQ.0.AND.ISTEP.GE.ISTOREM*DTCALCM) THEN
ccc           call cpu_time(startTime)
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
            DO ISTEP2=ISTEP-ISTOREM*DTCALCM,ISTEP,DTCALCM
               TT=MOD(ISTEP2/DTCALCM-1,ISTOREM)+1
               IF (TT.EQ.0) TT=ISTOREM
               DTC=(ISTEP-ISTEP2)/DTCALCM
c               print*,ISTOREM,ISTORED,ISKIPD,DTCALCD,ISKIPM,DTCALCM
c               print*,ISTEP,ISTEP2,TT,DTC
c               tempMSDO=MSDO(DTC+1)
c               tempVACFO=VACFO(DTC+1)
c               tempVACFH=VACFH(DTC+1)
c               tempCOUNTM=COUNTM(DTC+1)
c               tempCOUNTV=COUNTV(DTC+1)
c!$OMP PARALLEL default(private) shared(wlist,
c!$OMP& TT,NW,XO,YO,ZO,XC,YC,ZC,VXO,VYO,VZO,VXH,VYH,VZH,VX,VY,VZ,ISTEP,
c!$OMP& ISTEP2)
c!$OMP& shared(tempMSDO,tempVACFO,tempVACFH,tempCOUNTM,tempCOUNTV)
c!$OMP DO reduction(+:tempMSDO,tempVACFO,tempVACFH,tempCOUNTM,tempCOUNTV)
c               DO I=1,Nw,3
c                  K=wlist(I) ! Oxygen
c                  IF (ISTEP2.LT.ISTEP) THEN
c                     DX=XO(I,TT)-XC(K)
c                     DY=YO(I,TT)-YC(K)
c                     DZ=ZO(I,TT)-ZC(K)
cc                    write(*,*) I,TT,XO(I,TT)
c                  ELSE
c                     DX=0d0
c                     DY=0d0
c                     DZ=0d0
c                  ENDIF
c                  DR2=DX*DX+DY*DY+DZ*DZ
c                  tempMSDO=tempMSDO+DR2
c                  tempCOUNTM=tempCOUNTM+1
c
c             Calculate vacf
c
c                  II=(I+2)/3
c                  J=Nw/3+II ! I (J) - Hydrogen 1 (2)           
c                  L=K+1 ! Hydrogen 1
c                  M=K+2 ! Hydrogen 2
c                  IF (ISTEP2.LT.ISTEP) THEN
c                     privateVACFO=
c     +               VXO(II,TT)*VX(K)+VYO(II,TT)*VY(K)+VZO(II,TT)*VZ(K)
c                     privateVACFH=
c     +               VXH(II,TT)*VX(L)+VYH(II,TT)*VY(L)+VZH(II,TT)*VZ(L)
c     +               +VXH(J,TT)*VX(M)+VYH(J,TT)*VY(M)+VZH(J,TT)*VZ(M)
c                  ELSE
c                     privateVACFO=
c     +               VX(K)*VX(K)+VY(K)*VY(K)+VZ(K)*VZ(K)
c                     privateVACFH=
c     +               VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L)
c     +               +VX(M)*VX(M)+VY(M)*VY(M)+VZ(M)*VZ(M)
c                  ENDIF
c                  tempVACFO=tempVACFO+privateVACFO
c                  tempVACFH=tempVACFH+privateVACFH
c                  tempCOUNTV=tempCOUNTV+1
c               ENDDO
c!$OMP END DO
c!$OMP END PARALLEL
c               MSDO(DTC+1)=tempMSDO
c               VACFO(DTC+1)=tempVACFO
c               VACFH(DTC+1)=tempVACFH
c               COUNTM(DTC+1)=tempCOUNTM
c               COUNTV(DTC+1)=tempCOUNTV
c
c             Calculate DCF
c
ccc               DO II=1,N
ccc                  DCF(DTC+1)=DCF(DTC+1)+
ccc     &               MXI(II,TT)*MXC+MYI(II,TT)*MYC+MZI(II,TT)*MZC
ccc               ENDDO
               IF (ISTEP2.LT.ISTEP) THEN
                  DCF(DTC+1)=DCF(DTC+1)+MXI(TT)*MXC+MYI(TT)*MYC+
     &                       MZI(TT)*MZC
               ELSE
                  DCF(DTC+1)=DCF(DTC+1)+MXC*MXC+MYC*MYC+
     &                       MZC*MZC
               ENDIF
               COUNTD(DTC+1)=COUNTD(DTC+1)+1
               summass = 0d0
               do i=1,N
                  summass=summass+mass(i)
               enddo
               countmass=countmass+1
            ENDDO
ccc           call cpu_time(finishTime)
ccc           write(75,*) 'Time for calc msd, vacf, dcf at step',
ccc    &                  istep0,'=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(85,*) 'Wall time for calc msd, vacf, dcf',istep0,
ccc    &                  '=',elapsed_time
            IF (MOD(ISTEP,SAVE_DATA_FREQ).EQ.0) THEN
ccc           call cpu_time(startTime)
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c           OPEN(UNIT=202,FILE=FILE2,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=203,FILE=FILE3,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=204,FILE=FILE4,STATUS='REPLACE',ACTION='WRITE')
           OPEN(UNIT=205,FILE=FILE5,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=220,FILE=FILE20,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=221,FILE=FILE21,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=222,FILE=FILE22,STATUS='REPLACE',ACTION='WRITE')
           DO DTC=0,ISTOREM
c              WRITE(202,'(2ES16.8)') DT*DBLE(DTC*DTCALCM),
c    &           MSDO(DTC+1)/DBLE(COUNTM(DTC+1))
c              WRITE(203,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             ((VACFO(DTC+1)/(DBLE(COUNTV(DTC+1))))
c    &             /(VACFO(1)/DBLE(COUNTV(1))))
c              WRITE(204,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             ((VACFH(DTC+1)/(DBLE(COUNTV(DTC+1))))
c    &             /(VACFH(1)/DBLE(COUNTV(1))))
c               print*,"dt=",dt
c               print*,"dtc=",dtc
c               print*,"dtcalcd=",dtcalcd
c               print*,"istep0=",istep0
c               print*,"istep=",istep
c               print*,"istep2=",istep2
c               print*,"tt=",TT
               WRITE(205,'(2ES16.8)') DT*DBLE(DTC*DTCALCD),
     &             DCF(DTC+1)/DBLE(COUNTD(DTC+1))/summass
c              WRITE(220,'(2ES16.8)') DT*DBLE(DTC*DTCALCM),
c    &             DBLE(COUNTM(DTC+1))
c              WRITE(221,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             DBLE(COUNTV(DTC+1))
c              WRITE(222,'(2ES16.8)') DT*DBLE(DTC*DTCALCD),
c    &             DBLE(COUNTD(DTC+1))
           ENDDO
c           CLOSE(202)
c           CLOSE(203)
c           CLOSE(204)
           CLOSE(205)
c           CLOSE(220)
c           CLOSE(221)
c           CLOSE(222)
c           OPEN(UNIT=223,FILE=FILE23,STATUS='REPLACE',ACTION='WRITE')
c               WRITE(223,'(1ES16.8)') summass
c           CLOSE(223)
ccc           call cpu_time(finishTime)
ccc           write(76,*) 'Time for write msd, vacf, dcf at step',
ccc    &                  istep0,'=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(86,*) 'Wall time for write msd, vacf, dcf',istep0,
ccc    &                  '=',elapsed_time
            ENDIF
         ENDIF
c        
c        add current positions to tape
c        
         IF(MOD(ISTEP,DTCALCM).EQ.0) THEN
            TT=MOD(ISTEP/DTCALCM-1,ISTOREM)+1
c            DO I=1,Nw,3
c               K=wlist(I) ! Oxygen index
c               XO(I,TT)=XC(K)
c               YO(I,TT)=YC(K)
c               ZO(I,TT)=ZC(K)
c              write(*,*) I,TT
c              write(*,*) XO(I,TT),YO(I,TT),ZO(I,TT)
c               II=(I+2)/3
c         Add velocities to tape               
c               J=Nw/3+II ! I (J) - Hydrogen 1 (2) 
c               L=K+1 ! Hydrogen 1
c               M=K+2 ! Hydrogen 2
c               VXO(II,TT)=VX(K)
c               VYO(II,TT)=VY(K)
c               VZO(II,TT)=VZ(K)
c               VXH(II,TT)=VX(L)
c               VYH(II,TT)=VY(L)
c               VZH(II,TT)=VZ(L)
c               VXH(J,TT)=VX(M)
c               VYH(J,TT)=VY(M)
c               VZH(J,TT)=VZ(M)
c            ENDDO
c         Add dipole moments to tape            
ccc            DO I=1,N
ccc               MXI(I,TT)=MXIC(I)
ccc               MYI(I,TT)=MYIC(I)
ccc               MZI(I,TT)=MZIC(I)
ccc            ENDDO
            MXI(TT)=MXC
            MYI(TT)=MYC
            MZI(TT)=MZC
         ENDIF

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
