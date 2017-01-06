c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pmonte  --  Monte Carlo barostat trial moves  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pmonte" implements a Monte Carlo barostat via random trial
c     changes in the periodic box volume and shape
c
c     literature references:
c
c     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
c     2nd Edition", Academic Press, San Diego, CA, 2002; Section 5.4
c
c     original version written by Alan Grossfield, January 2004;
c     anisotropic modification implemented by Lee-Ping Wang, Stanford
c     University, March 2013
c
c
      subroutine pmonte_par (epot,temp,istep)
      use sizes
      use atomid
      use atoms
      use bath
      use boxes
      use group
      use math
      use mdstuf
      use molcul
      use units
      use usage
      use mpidat 
      use deriv3b
      use energi
      use neigh3b
      use neigh
      use mpole
      implicit none
      include 'mpif.h'
      integer i,j,k
      integer start,stop
      real*8 epot,temp,term
      real*8 energy,random
      real*8 third,weigh
      real*8 step,scale
      real*8 enew,rnd6
      real*8 xcm,ycm,zcm
      real*8 vold,vfrac
      real*8 cosine,diff
      real*8 xmove,ymove,zmove
      real*8 xboxold,yboxold,zboxold
      real*8 alphaold,betaold,gammaold
      real*8 kt,de,dv,lnv,expterm
      real*8 temp3(3,3)
      real*8 hbox(3,3)
      real*8 ascale(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      integer ierr,istep
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     save the old lattice parameters, box size and coordinates
c
      print*,"Here we are in parallelized pmonte for 3-amoeba"
      if (random() .lt. 1.0d0/dble(voltrial)) then
         xboxold = xbox
         yboxold = ybox
         zboxold = zbox
         alphaold = alpha
         betaold  = beta
         gammaold = gamma
         vold = volbox
         do i = 1, n
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
         end do
         third = 1.0d0 / 3.0d0
c
c     for the isotropic case, change the lattice lengths uniformly
c
         if (.not.anisotrop .or. random().lt.0.5d0) then
            step = volmove * (2.0d0*random()-1.0d0)
            volbox = volbox + step
            scale = (volbox/vold)**third
            xbox = xbox * scale
            ybox = ybox * scale
            zbox = zbox * scale
            call lattice
            if (integrate .eq. 'RIGIDBODY') then
               diff = scale - 1.0d0
               do i = 1, ngrp
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = igrp(1,i)
                  stop = igrp(2,i)
                  do j = start, stop
                     k = kgrp(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xmove = diff * xcm/grpmass(i)
                  ymove = diff * ycm/grpmass(i)
                  zmove = diff * zcm/grpmass(i)
                  do j = start, stop
                     k = kgrp(j)
                     x(k) = x(k) + xmove
                     y(k) = y(k) + ymove
                     z(k) = z(k) + zmove
                  end do
               end do
            else if (volscale .eq. 'MOLECULAR') then
               diff = scale - 1.0d0
               do i = 1, nmol
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = imol(1,i)
                  stop = imol(2,i)
                  do j = start, stop
                     k = kmol(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xmove = diff * xcm/molmass(i)
                  ymove = diff * ycm/molmass(i)
                  zmove = diff * zcm/molmass(i)
                  do j = start, stop
                     k = kmol(j)
                     if (use(k)) then
                        x(k) = x(k) + xmove
                        y(k) = y(k) + ymove
                        z(k) = z(k) + zmove
                     end if
                  end do
               end do
            else
               do i = 1, n
                  if (use(i)) then
                     x(i) = x(i) * scale
                     y(i) = y(i) * scale
                     z(i) = z(i) * scale
                  end if
               end do
            end if
c
c     for anisotropic case alter lattice angles, then scale lengths
c
         else
            rnd6 = 6.0d0*random()
            step  = volmove * (2.0d0*random()-1.0d0)
            scale = (1.0d0+step/vold)**third
            ascale(1,1) = 1.0d0
            ascale(2,2) = 1.0d0
            ascale(3,3) = 1.0d0
            if (monoclinic .or. triclinic) then
               if (rnd6 .lt. 1.0d0) then
                  ascale(1,1) = scale
               else if (rnd6 .lt. 2.0d0) then
                  ascale(2,2) = scale
               else if (rnd6 .lt. 3.0d0) then
                  ascale(3,3) = scale
               else if (rnd6 .lt. 4.0d0) then
                  ascale(1,2) = scale - 1.0d0
                  ascale(2,1) = scale - 1.0d0
               else if (rnd6 .lt. 5.0d0) then
                  ascale(1,3) = scale - 1.0d0
                  ascale(3,1) = scale - 1.0d0
               else
                  ascale(2,3) = scale - 1.0d0
                  ascale(3,2) = scale - 1.0d0
               end if
            else
               if (rnd6 .lt. 2.0d0) then
                  ascale(1,1) = scale
               else if (rnd6 .lt. 4.0d0) then
                  ascale(2,2) = scale
               else
                  ascale(3,3) = scale
               end if
            end if
c
c     modify the current periodic box dimension values
c
            temp3(1,1) = xbox
            temp3(2,1) = 0.0d0
            temp3(3,1) = 0.0d0
            temp3(1,2) = ybox * gamma_cos
            temp3(2,2) = ybox * gamma_sin
            temp3(3,2) = 0.0d0
            temp3(1,3) = zbox * beta_cos
            temp3(2,3) = zbox * beta_term
            temp3(3,3) = zbox * gamma_term
            do i = 1, 3
               do j = 1, 3
                  hbox(j,i) = 0.0d0
                  do k = 1, 3
                     hbox(j,i) = hbox(j,i) + ascale(j,k)*temp3(k,i)
                  end do
               end do
            end do
            xbox = sqrt(hbox(1,1)**2 + hbox(2,1)**2 + hbox(3,1)**2)
            ybox = sqrt(hbox(1,2)**2 + hbox(2,2)**2 + hbox(3,2)**2)
            zbox = sqrt(hbox(1,3)**2 + hbox(2,3)**2 + hbox(3,3)**2)
            if (monoclinic) then
               cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
     &                     + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
               beta = radian * acos(cosine)
            else if (triclinic) then
               cosine = (hbox(1,2)*hbox(1,3) + hbox(2,2)*hbox(2,3)
     &                     + hbox(3,2)*hbox(3,3)) / (ybox*zbox)
               alpha = radian * acos(cosine)
               cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
     &                     + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
               beta = radian * acos(cosine)
               cosine = (hbox(1,1)*hbox(1,2) + hbox(2,1)*hbox(2,2)
     &                     + hbox(3,1)*hbox(3,2)) / (xbox*ybox)
               gamma = radian * acos(cosine)
            end if
c
c     find the new box dimensions and other lattice values
c
            call lattice
            vfrac = vold / volbox
            scale = vfrac**third
            xbox = xbox * scale
            ybox = ybox * scale
            zbox = zbox * scale
            call lattice
c
c     scale the coordinates by groups, molecules or atoms
c
            if (integrate .eq. 'RIGIDBODY') then
               ascale(1,1) = ascale(1,1) - 1.0d0
               ascale(2,2) = ascale(2,2) - 1.0d0
               ascale(3,3) = ascale(3,3) - 1.0d0
               do i = 1, ngrp
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = igrp(1,i)
                  stop = igrp(2,i)
                  do j = start, stop
                     k = kgrp(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xcm = xcm / grpmass(i)
                  ycm = ycm / grpmass(i)
                  zcm = zcm / grpmass(i)
                  xmove = ascale(1,1)*xcm + ascale(1,2)*ycm
     &                       + ascale(1,3)*zcm
                  ymove = ascale(2,1)*xcm + ascale(2,2)*ycm
     &                       + ascale(2,3)*zcm
                  zmove = ascale(3,1)*xcm + ascale(3,2)*ycm
     &                       + ascale(3,3)*zcm
                  do j = start, stop
                     k = kgrp(j)
                     x(k) = x(k) + xmove
                     y(k) = y(k) + ymove
                     z(k) = z(k) + zmove
                  end do
               end do
            else if (volscale .eq. 'MOLECULAR') then
               ascale(1,1) = ascale(1,1) - 1.0d0
               ascale(2,2) = ascale(2,2) - 1.0d0
               ascale(3,3) = ascale(3,3) - 1.0d0
               do i = 1, nmol
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = imol(1,i)
                  stop = imol(2,i)
                  do j = start, stop
                     k = kmol(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xcm = xcm / molmass(i)
                  ycm = ycm / molmass(i)
                  zcm = zcm / molmass(i)
                  xmove = ascale(1,1)*xcm + ascale(1,2)*ycm
     &                       + ascale(1,3)*zcm
                  ymove = ascale(2,1)*xcm + ascale(2,2)*ycm
     &                       + ascale(2,3)*zcm
                  zmove = ascale(3,1)*xcm + ascale(3,2)*ycm
     &                       + ascale(3,3)*zcm
                  do j = start, stop
                     k = kmol(j)
                     if (use(k)) then
                        x(k) = x(k) + xmove
                        y(k) = y(k) + ymove
                        z(k) = z(k) + zmove
                     end if
                  end do
               end do
            else
               do i = 1, n
                  if (use(i)) then
                     x(i) = ascale(1,1)*x(i) + ascale(1,2)*y(i)
     &                         + ascale(1,3)*z(i)
                     y(i) = ascale(2,1)*x(i) + ascale(2,2)*y(i)
     &                         + ascale(2,3)*z(i)
                     z(i) = ascale(3,1)*x(i) + ascale(3,2)*y(i)
     &                         + ascale(3,3)*z(i)
                  end if
               end do
            end if
         end if
c
c     find the energy change and PV term for the trial move
c

c         enew = energy ()
            call mpi_barrier(mpi_comm_world,ierr)

            if(taskid.eq.master) then
C PERFORM 1ST STEP OF VERLET UPDATE
C  CALLS TO 'bounds' AND TO 'replica' IN 'gradient_setup'
             call gradient_setup
C  In 'alloc3b_vdwlist':  ALLOCATE POLARIZATION GRADIENT VECTOR ON
C  MASTER TASK.  INITIALIZE POLARIZATION ENERGY, GRADIENT, AND VIRIAL TO ZERO.
C  LASTLY, BUILD/UPDATE VDW NEIGHBOR LIST AND
C  BUILD/UPDATE POLARIZATION NEIGHBOR LIST
            ! call alloc3b_vdwlist
             call alloc3b_vdwlist_1list
C  PERFORM LOAD BALANCING FOR POLARIZATION BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
             !if(firstload_polz) then
             !  firstload_polz = .false.
             !  call polar_load_balance
             !end if
C  ALLOCATE VECTOR FOR THE SUM OF THE CONTRIBUTIONS TO GRADIENT FROM EACH TYPE OF ENERGETIC INTERACTION
C  'prep_pole' CONTAINS CALLS TO 'chkpole' AND 'rotpole'. PLEASE REFER TO THE COMMENTS IN THOSE FILES FOR DETAILS.
             call prep_pole
C  ALLOCATE REAL-SPACE PERM ELECTROSTATICS GRADIENT VECTOR; INITIALIZE THIS VECTOR AND THE ASSOC. ENERGY 
C  AND VIRIAL TO ZERO
            ! call allocPermElec1
              emreal_tmp=0.0d0
            end if

C ALLOCATE PME GRADIENT VECTOR ON ALL TASKS DUE TO LATER BROADCAST
            !call alloc_erecip
               em = 0.0d0

C BROADCAST COORDINATES AT EACH STEP
            call mpi_bcast(x,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(y,n,mpi_real8,master,
     &      mpi_comm_world,ierr)
            call mpi_bcast(z,n,mpi_real8,master,
     &      mpi_comm_world,ierr)

            if(taskid .eq. master1) then
C  BUILD/UPDATE NEIGHBOR LIST FOR REAL SPACE PERMANENT ELECTROSTATICS
              call mlist
C  PERFORM LOAD BALANCING FOR REAL-SPACE PERM ELEC BASED ON THE NUMBER OF NEIGHBORS OF EACH OUTER LOOP MOLECULE INDEX
              !if(firstload_emreal) then
              !  firstload_emreal = .false.
              !  call emreal_load_balance_sizelist
              !else
                call emreal_sizelist_only
              !end if
            end if

            call mpi_barrier(mpi_comm_world,ierr)
C BROADCAST PERMANENT MULTIPOLES AFTER ROTATION TO GLOBAL COORDINATE FRAME, 
C AS WELL AS LOGICALS THAT ARE TRUE IF THE LISTS HAVE BEEN REBUILT.
            call mpi_bcast(listbcast,1,mpi_logical,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(rpole,13*n,mpi_real8,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(listsend_mpole,1,mpi_logical,master1,
     &       mpi_comm_world,ierr)

            if(listbcast) then
C IF 'listbcast' IS TRUE, POLZ NEIGHBOR LIST HAS BEEN REBUILT OR BUILT DE NOVO.  THREFORE, LIST MUST BE BROADCAST, 
C AS WELL AS VARIABLES ASSOCIATED W/ BETTER LOAD-BALANCED POLZ.
                call bcast_polar_load_balance
            call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
            call mpi_bcast(mollst,max2blst*nmol,mpi_integer,master,
     &       mpi_comm_world,ierr)
          !  call mpi_bcast(nmollst2,nmol,mpi_integer,master,
     &    !  mpi_comm_world,ierr)
          !  call mpi_bcast(mollst2,max2blst2*nmol,mpi_integer,master,
     &    !   mpi_comm_world,ierr)
            end if
            listbcast = .false.

C IF 'listsend_mpole' IS TRUE, REAL-SPACE PERM ELECTROSTATICS NEIGHBORLIST HAS BEEN REBUILT OR BUILT DE NOVO.
C THEREFORE, LIST MUST BE BROADCAST. ALSO BROADCAST VARIABLES TO DO WITH LOAD BALANCING
              call mpi_bcast(numtasks_emreal2,1,mpi_integer,
     &             master1,emreal_comm,ierr)

            if(listsend_mpole.and.taskid.lt.numtasks_emreal) then
              call mpi_bcast(start_emreal2,numtasks_emreal,mpi_integer,
     &             master1,emreal_comm,ierr)
              call mpi_bcast(last_emreal2,numtasks_emreal,mpi_integer,
     &             master1,emreal_comm,ierr)
               call mpi_bcast(maxsize_elst,numtasks_emreal,mpi_integer,
     &             master1,emreal_comm,ierr)
              call sendmlist_load_balanced3_sizelist
            end if
            listsend_mpole = .false.
            if(taskid.eq.master) then
C MASTER TASK PERFORMS ENERGY/GRADIENT/VIRAL CALC OF COVALENT AND VAN DER WAALS.
             !call gradient_covalent_vdw
              call eng_covalent_vdw
            else if (taskid.eq.master1) then
C MASTER1 TASK PERFORMS PARTICLE MESH EWALD AND SENDS IT BACK TO MASTER AND EVERYWHERE ELSE
C VIA BROADCAST
             !call gradient_emrecip
              call emrecip_Perm
            end if

C IF TASKID IS LESS THAN 'numtaks_emreal', THEN CALCULATE A PORTION OF THE ENERGY/GRADIENT/VIRIAL 
C OF REAL-SPACE PERM ELECTROSTATICS AND ACCUMULATE VIA MPI_REDUCE
            if(taskid.le.numtasks_emreal-1) then
             ! call gradient_emreal_reduce_load_balanced
              call eng_emreal_reduce_load_balanced
            end if

C ALL TASKS CALCULATE A PORTION OF ENERGY/GRADIENT/VIRIAL OF POLARIZATION. 
             call smallsmooth2_eng_polar_load_balanced

            call mpi_barrier(mpi_comm_world,ierr)

C BROADCAST PME ENERGY,GRADIENT,VIRIAL
           call mpi_bcast(em,1,mpi_real8,master1,
     &       mpi_comm_world,ierr)

            if(taskid.eq.master) then
C SUM REAL AND RECIPROCAL SPACE PORTIONS OF REAL-SPACE AND COMPUTE SELF-ENERGY TERM.
                 call empole0d_3b_Perm_selfeng_bcast

                enew = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b
       print*," PMONTE Vdw, PermElecEng, PolEng",istep,ev,em,ep3b
            end if

            call mpi_bcast(enew,1,mpi_real8,
     &             master,mpi_comm_world,ierr)

         de = enew - epot
         dv = atmsph * (volbox-vold) / prescon
c
c     set the entropy of mixing term based on system type
c
         if (isothermal) then
            kt = gasconst * kelvin
         else
            kt = gasconst * temp
         end if
         if (integrate .eq. 'RIGIDBODY') then
            lnv = dble(ngrp) * kt * log(volbox/vold)
         else if (volscale .eq. 'MOLECULAR') then
            lnv = dble(nmol) * kt * log(volbox/vold)
         else
            lnv = dble(nuse) * kt * log(volbox/vold)
         end if
c
c     acceptance ratio from energy change, PV and mixing term
c
         term = -(de+dv-lnv) / kt
         expterm = exp(term)
c
c     reject the step, restore old box size and coordinates
c
         if (random() .gt. expterm) then
            xbox = xboxold
            ybox = yboxold
            zbox = zboxold
            call lattice
            do i = 1, n
               x(i) = xold(i)
               y(i) = yold(i)
               z(i) = zold(i)
            end do
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      return
      end
c
