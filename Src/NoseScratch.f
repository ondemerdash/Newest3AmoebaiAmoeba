
         if (integrate .eq. 'NOSE-HOOVER') then
          do istep=1,nstep
            if(taskid.eq.master) then
             call nose_pt1 (istep,dt,press)
             !call ptest_parallel

             ! x(6)=x(6)+delx
c              y(6)=y(6)+delx
c              z(6)=z(6)+delx
             call gradient_setup
             call allocPermElec1

             call alloc3b_vdwlist_nolist_vac
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
           else if(taskid.eq.master.and.((mod(istep,10).eq.0)
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
              if((istep.ne.1).and.(taskid.ge.numtasks_polar1).and.
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

              if((istep.ne.1).and.
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
                if(.not. allocated(ep2bmatmod))
     &           allocate(ep2bmatmod(num2bsave))
               if(.not. allocated(dep2bmatmod)) allocate(dep2bmatmod(3,
     &          2*3*maxsizeclust,num2bsave))
            if(.not. allocated(virep2bmatmod)) allocate(virep2bmatmod(3,
     &          3,num2bsave))
             call init2bmatmodPolar
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

           if( (mod(istep,5).eq.0 ) .or. exist) then
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

          call clust_gradient_Polar1b2b3bsimult_ireducematcommsmall2(
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
      call sum_gradient_Polar1b2b3bsimult_ireducematcommsmall2

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

       print*,"Nose_istep Vdw, PermElec, Pol",istep,ev,em,ep3b!+ep3bmut

              deallocate (derivs)
               call nose_pt2(istep,dt,press,energy)
            end if

          end do
         end if
