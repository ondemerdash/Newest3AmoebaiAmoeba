      subroutine gradient_emreal_reduce_load_balanced_totfield_dEtensor
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      use totfield
      use math
      use ewald
      use dEtensor
      use sizes
      use polar
      use aprx
      use deriv3b
      use chgpot   
      use paremneigh
      implicit none
      include 'mpif.h'
      real*8 emrealt,viremrealt(3,3)
      real*8, allocatable :: demrealt(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6
      real*8, allocatable :: fieldnpolet(:,:)
      real*8, allocatable :: fieldpnpolet(:,:)
      real*8 ucell(3),term,f
      integer k,maxsize_master,maxoffset_master,offset_cum
      integer stat(MPI_STATUS_SIZE),sz,maxlocal_allthread
      integer maxlocal_allthread_tau
c      real*8, allocatable :: dEd1tmp(:,:)
c      real*8, allocatable :: dEd2tmp(:,:)
c      real*8, allocatable :: dEp1tmp(:,:)
c      real*8, allocatable :: dEp2tmp(:,:)
c      integer, allocatable :: dEindextmp(:,:)
c      integer ntpair_dEtmp
      allocate (demrealt(3,npole))
      allocate (fieldnpolet(3,npole))
      allocate (fieldpnpolet(3,npole))

      f = electric / dielec
      
      emrealt=0.0d0
      do i=1,npole
         do j=1,3
            demrealt(j,i)=0.0d0
            fieldnpolet(j,i)=0.0d0
            fieldpnpolet(j,i) = 0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           viremrealt(j,i)=0.0d0
        end do
      end do
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"
             !t1=mpi_Wtime()

      if(taskid.lt.numtasks_emreal2) then

      if(allocated(dEindextmp)) deallocate(dEindextmp)
      if(allocated(dEd1tmp)) deallocate(dEd1tmp)
      if(allocated(dEd2tmp)) deallocate(dEd2tmp)
      if(allocated(dEp1tmp)) deallocate(dEp1tmp)
      if(allocated(dEp2tmp)) deallocate(dEp2tmp)

      if(allocated(taud1tmp)) deallocate(taud1tmp)
      if(allocated(taud2tmp)) deallocate(taud2tmp)
      if(allocated(taup1tmp)) deallocate(taup1tmp)
      if(allocated(taup2tmp)) deallocate(taup2tmp)

      if(allocated(tau1indextmp)) deallocate(tau1indextmp)
      if(allocated(tau2indextmp)) deallocate(tau2indextmp)
      
      sz=last_emreal2(taskid)-start_emreal2(taskid)+1
      
c      maxlocal_allthread=int(dble(sz)*dble(maxsize_elst(taskid)))
         maxlocal_allthread=0
         
           do i=start_emreal2(taskid),last_emreal2(taskid)
              j=i-start_emreal2(taskid)+1
c             print*,"i",i,"j",j,"nelst_recv(i)",nelst_recv(j)
              maxlocal_allthread=maxlocal_allthread+nelst_recv(j)     
           end do
      maxlocal_allthread_tau=3*maxlocal_allthread

      !print*,"taskid=",taskid,"sz=",sz
      !maxlocal_allthread=int(dble(sz)*dble(maxelst))
      

      allocate(dEindextmp(2,maxlocal_allthread))
      allocate(dEd1tmp(9,maxlocal_allthread))
      allocate(dEd2tmp(9,maxlocal_allthread))
      allocate(dEp1tmp(9,maxlocal_allthread))
      allocate(dEp2tmp(9,maxlocal_allthread))

      allocate(taud1tmp(9,maxlocal_allthread_tau))
      allocate(taud2tmp(9,maxlocal_allthread_tau))
      allocate(taup1tmp(9,maxlocal_allthread_tau))
      allocate(taup2tmp(9,maxlocal_allthread_tau))

      allocate(tau1indextmp(2,maxlocal_allthread_tau))
      allocate(tau2indextmp(2,maxlocal_allthread_tau))

      ! call LoadBal_ereal1d_3b_Perm_sentlist2_totfield_dEtensorOmp2( 
     & ! emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet) 
      !call LoadBalInner_thg4omp_mpi_72015torque2 (
     & !emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
      if(npole.le.21000) then
c        call thg5_small(
c     & emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
        call thg6_save_small(
     & emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

      else
c        call thg5(
c     & emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
        call thg6_save(
     & emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

      end if

        if(taskid.ne.master.and.approxmode.ne.'1BODYMODE') then
      print*,"After LoadBal..blah..,before mpi_send"

          call mpi_send(ntpair_dEtmp,1,mpi_integer,master,
     &                 6*taskid,emreal_comm,ierr)
      print*,"After a single mpi_send"

          call mpi_send(dEindextmp,2*ntpair_dEtmp,mpi_integer,master,
     &                 6*taskid+1,emreal_comm,ierr)
      print*,"After dEindex send"
          call mpi_send(dEd1tmp,9*ntpair_dEtmp,mpi_real8,master,
     &                 6*taskid+2,emreal_comm,ierr)
      print*,"After dEd1tmp send"
          call mpi_send(dEd2tmp,9*ntpair_dEtmp,mpi_real8,master,
     &                 6*taskid+3,emreal_comm,ierr)
      print*,"After dEd2tmp send"

          call mpi_send(dEp1tmp,9*ntpair_dEtmp,mpi_real8,master,
     &                 6*taskid+4,emreal_comm,ierr)
          call mpi_send(dEp2tmp,9*ntpair_dEtmp,mpi_real8,master,
     &                 6*taskid+5,emreal_comm,ierr)
       print*,"After all sends"
        end if
       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      end if
        !    call mpi_barrier(emreal_comm,ierr)

         !    t1=mpi_Wtime()
       call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)

       call mpi_reduce(fieldnpolet,fieldnpole,npole*3,mpi_real8,
     &      mpi_sum,master1,mpi_comm_world,ierr)
       call mpi_reduce(fieldpnpolet,fieldpnpole,npole*3,mpi_real8,
     &      mpi_sum,master1,mpi_comm_world,ierr)
c       print*,"PermElec Real After all reduce and send!"


      if(taskid.eq.master1) then
c
c     get the self-energy portion of the electrostatic field
c
       term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c       do i = 1, npole
c         do j = 1, 3
c            fieldnpole(j,i)=fieldnpolet(j,i)
c            fieldpnpole(j,i)=fieldpnpolet(j,i)
c            demreal_tmp(j,i)=demrealt(j,i)
           ! fieldnpole(j,i) = fieldnpole(j,i) +fieldnpole_rcp(j,i)
           ! fieldpnpole(j,i) = fieldpnpole(j,i) +fieldnpole_rcp(j,i)
           ! fieldnpole(j,i) = fieldnpole(j,i) + term*rpole(j+1,i)
           ! fieldpnpole(j,i) = fieldpnpole(j,i) + term*rpole(j+1,i)
c         end do
c       end do


      else if(taskid.eq.master.and.approxmode.ne.'1BODYMODE') then
        if(allocated(dEindextmp_rcv)) deallocate(dEindextmp_rcv)
        if(allocated(dEd1tmp_rcv)) deallocate(dEd1tmp_rcv)
        if(allocated(dEd2tmp_rcv)) deallocate(dEd2tmp_rcv)
        if(allocated(dEp1tmp_rcv)) deallocate(dEd1tmp_rcv)
        if(allocated(dEd2tmp_rcv)) deallocate(dEd2tmp_rcv)

        do i=0,numtasks_emreal2-1
           if(i.eq.0) then
             maxsize_master=maxsize_elst(i)
             maxoffset_master=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           else             
              if(maxsize_elst(i).gt.maxsize_master) then
                maxsize_master=maxsize_elst(i)
              end if
            
              if( (last_emreal2(taskid)-start_emreal2(taskid)+1).gt.
     &           maxoffset_master) then
                 maxoffset_master=last_emreal2(taskid)
     &                      -start_emreal2(taskid)+1
              end if
           end if
        end do

        allocate(dEindextmp_rcv(2,maxelst*maxoffset_master))
        allocate(dEd1tmp_rcv(9,maxelst*maxoffset_master))
        allocate(dEd2tmp_rcv(9,maxelst*maxoffset_master))
        allocate(dEd1tmp_rcv(9,maxelst*maxoffset_master))
        allocate(dEd2tmp_rcv(9,maxelst*maxoffset_master))

        offset_cum=0
        do i=0,numtasks_emreal2-1
            if(i.eq.master) then
               do k=offset_cum+1,offset_cum+ntpair_dEtmp
                  dEindex(1,k)=dEindextmp(1,k)
                  dEindex(2,k)=dEindextmp(2,k)
                  do j=1,9
                    dEd1(j,k)=dEd1tmp(j,k)
                    dEd2(j,k)=dEd2tmp(j,k)
                    dEp1(j,k)=dEp1tmp(j,k)
                    dEp2(j,k)=dEp2tmp(j,k)
                  end do
               end do
               offset_cum=offset_cum+ntpair_dEtmp
            else
              call mpi_recv(ntpair_dEtmp_rcv,1,mpi_integer,i,6*i,
     &          emreal_comm,stat,ierr)
              call mpi_recv(dEindextmp_rcv,2*ntpair_dEtmp_rcv,
     &          mpi_integer,i,6*i+1,emreal_comm,stat,ierr)
              call mpi_recv(dEd1tmp_rcv,9*ntpair_dEtmp_rcv,
     &          mpi_real8,  i,6*i+2,emreal_comm,stat,ierr)  
              call mpi_recv(dEd2tmp_rcv,9*ntpair_dEtmp_rcv,
     &          mpi_real8,  i,6*i+3,emreal_comm,stat,ierr)
              call mpi_recv(dEp1tmp_rcv,9*ntpair_dEtmp_rcv,
     &          mpi_real8,  i,6*i+4,emreal_comm,stat,ierr)
              call mpi_recv(dEp2tmp_rcv,9*ntpair_dEtmp_rcv,
     &          mpi_real8,  i,6*i+5,emreal_comm,stat,ierr)
              
              do k=offset_cum+1,offset_cum+ntpair_dEtmp_rcv
                 dEindex(1,k)=dEindextmp_rcv(1,k)
                 dEindex(2,k)=dEindextmp_rcv(2,k)
                  do j=1,9
                    dEd1(j,k)=dEd1tmp_rcv(j,k)
                    dEd2(j,k)=dEd2tmp_rcv(j,k)
                    dEp1(j,k)=dEp1tmp_rcv(j,k)
                    dEp2(j,k)=dEp2tmp_rcv(j,k)
                  end do
              end do
               offset_cum=offset_cum+ntpair_dEtmp_rcv
            end if           
        end do
        ntpair_dE=offset_cum
        deallocate(dEindextmp_rcv)
        deallocate(dEd1tmp_rcv)
        deallocate(dEd2tmp_rcv)
        deallocate(dEp1tmp_rcv)

      end if

      if(approxmode.ne.'1BODYMODE') then
       deallocate(dEindextmp)
       deallocate(dEd1tmp)
       deallocate(dEd2tmp)
       deallocate(dEp1tmp)
       deallocate(dEp2tmp)
      call mpi_bcast(ntpair_dE,1,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(dEindex,2*ntpair_dE,mpi_integer,master,
     &              mpi_comm_world,ierr)
      call mpi_bcast(dEd1,9*ntpair_dE,mpi_real8,master,
     &              mpi_comm_world,ierr)
      call mpi_bcast(dEd2,9*ntpair_dE,mpi_real8,master,
     &              mpi_comm_world,ierr)
      call mpi_bcast(dEp1,9*ntpair_dE,mpi_real8,master,
     &              mpi_comm_world,ierr)
      call mpi_bcast(dEp2,9*ntpair_dE,mpi_real8,master,
     &              mpi_comm_world,ierr)
      end if  
      
      deallocate(demrealt)
      deallocate(fieldnpolet)
      deallocate(fieldpnpolet)

      call mpi_barrier(mpi_comm_world,ierr)

      if(taskid.eq.master1.and.approxmode.eq.'1BODYMODE') then
        ep3b=0.0d0
        do i=1,npole
          uind(1,i)=polarity(i) * fieldnpole(1,i)
          uind(2,i)=polarity(i) * fieldnpole(2,i)
          uind(3,i)=polarity(i) * fieldnpole(3,i)

          uinp(1,i)=polarity(i) * fieldpnpole(1,i)
          uinp(2,i)=polarity(i) * fieldpnpole(2,i)
          uinp(3,i)=polarity(i) * fieldpnpole(3,i)
c          print*,"uind(1,i)",uind(1,i)
c          print*,"uind(2,i)",uind(2,i)
c          print*,"uind(3,i)",uinp(3,i)
          ep3b=ep3b - 0.5d0*f*( uind(1,i)*fieldpnpole(1,i) +
     &    uind(2,i)*fieldpnpole(2,i) + uind(3,i)*fieldpnpole(3,i))
        end do
        !print*,"ep3b in grad_r_reduce_totfield_dEtensor",ep3b 
      end if
             !t2=mpi_Wtime()
       !print*,"TensPermElec Real taskid",taskid,"time=",t2-t1

      return
      end

