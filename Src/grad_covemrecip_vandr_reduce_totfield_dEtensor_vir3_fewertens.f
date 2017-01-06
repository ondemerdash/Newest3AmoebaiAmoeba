      subroutine grad_covemrecip_vandr_reduce_totfield_dEtensor_vir7 
      use mpole
      use mpidat
      use mpidat2
      use deriv
      use energi
      use virial
      use totfield
      use math
      use ewald
      use dEtensor
      use sizes
      use sizes2
      use polar
      use aprx
      use deriv3b
      use chgpot   
      use paremneigh
      use neigh
      use neigh2
      use atoms
      use vdw
      use pme
      use chunks
      use dEtensor2
      use mpidat3      
      use dEtensor3
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
      real*8 evtmp,virevtmp(3,3)
      real*8, allocatable :: devtmp(:,:)
      real*8, allocatable :: dep3bt1(:,:)
      real*8 virep3bt1(3,3)
      real*8, allocatable :: uindt(:,:)
      real*8, allocatable :: uinpt(:,:)
      integer, allocatable :: ntpair_dE_array(:)
      integer, allocatable :: ntpair_tau1_array(:)
      integer, allocatable :: ntpair_tau2_array(:)
      integer m,m1,m2,ntpair_dEoffset 
      integer ntpair_tau1offset,ntpair_tau2offset
      integer iter,k2,k3,largest_nelst,largest_nelst2
      logical first
      save first
      data first  / .true. /

      

      allocate (demrealt(3,npole))
      allocate (fieldnpolet(3,npole))
      allocate (fieldpnpolet(3,npole))

      allocate (devtmp(3,nvdw))

      f = electric / dielec
      evtmp=0.0d0
      emrealt=0.0d0
      do i=1,n
         do j=1,3
            demrealt(j,i)=0.0d0
            fieldnpolet(j,i)=0.0d0
            fieldpnpolet(j,i) = 0.0d0
            devtmp(j,i)=0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           viremrealt(j,i)=0.0d0
           virevtmp(j,i)=0.0d0
        end do
      end do
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"
             !t1=mpi_Wtime()
      if(taskid.eq.master) then
         call gradient_covalent2
c         print*,"Finished gradient_covalent2 eub=",eub
      else if(taskid.eq.master1) then
        em = 0.0d0

        if (first) then
         first = .false.
         if (.not. allocated(dem))  allocate (dem(3,n))
        end if

        do i = 1, n
           do j=1,3
              dem(j,i) =0.0d0
           end do
        end do
        do i = 1, 3
           do j=1,3
            viremrecip(j,i) =0.0d0
           end do
        end do

        call emrecip1_3b_Perm2
        call mpi_isend(em,1,mpi_real8,master,1,mpi_comm_world,
     &      reqs1,ierr)
        call mpi_isend(dem,3*n,mpi_real8,master,2,mpi_comm_world,
     &      reqs2,ierr)
        call mpi_isend(viremrecip,3*3,mpi_real8,master,3,
     &      mpi_comm_world,reqs3,ierr)

       !  call mpi_send(em,1,mpi_real8,master,1,mpi_comm_world,
     & !      ierr)
       !call mpi_send(dem,3*n,mpi_real8,master,2,mpi_comm_world,
     & !     ierr)
       !call mpi_send(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
     & !     ierr)
 
      ! if(taskid.lt.numtasks_emreal2) then
      else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then

      ! if(allocated(dEindextmp)) deallocate(dEindextmp)
      if(allocated(dEd1_3tmp)) deallocate(dEd1_3tmp)
      if(allocated(dEd2_3tmp)) deallocate(dEd2_3tmp)
      if(allocated(dEp1_3tmp)) deallocate(dEp1_3tmp)
      if(allocated(dEp2_3tmp)) deallocate(dEp2_3tmp)

      if(allocated(taud1_3tmp)) deallocate(taud1_3tmp)
      if(allocated(taud2_3tmp)) deallocate(taud2_3tmp)
      if(allocated(taup1_3tmp)) deallocate(taup1_3tmp)
      if(allocated(taup2_3tmp)) deallocate(taup2_3tmp)

      !if(allocated(tau1indextmp)) deallocate(tau1indextmp)
      !if(allocated(tau2indextmp)) deallocate(tau2indextmp)

      if(allocated(frcztau1dtot_3tmp)) deallocate (frcztau1dtot_3tmp)
      if(allocated(frcytau1dtot_3tmp)) deallocate (frcytau1dtot_3tmp)
      if(allocated(frcxtau1dtot_3tmp)) deallocate (frcxtau1dtot_3tmp)
      if(allocated(frcztau1ptot_3tmp)) deallocate (frcztau1ptot_3tmp)
      if(allocated(frcytau1ptot_3tmp)) deallocate (frcytau1ptot_3tmp)
      if(allocated(frcxtau1ptot_3tmp)) deallocate (frcxtau1ptot_3tmp)
      if(allocated(frcztau2dtot_3tmp)) deallocate (frcztau2dtot_3tmp)
      if(allocated(frcytau2dtot_3tmp)) deallocate (frcytau2dtot_3tmp)
      if(allocated(frcxtau2dtot_3tmp)) deallocate (frcxtau2dtot_3tmp)
      if(allocated(frcztau2ptot_3tmp)) deallocate (frcztau2ptot_3tmp)
      if(allocated(frcytau2ptot_3tmp)) deallocate (frcytau2ptot_3tmp)
      if(allocated(frcxtau2ptot_3tmp)) deallocate (frcxtau2ptot_3tmp)

      if(allocated(elsttau1_index_tmp)) deallocate(elsttau1_index_tmp)
      if(allocated(elsttau2_index_tmp)) deallocate(elsttau2_index_tmp)

      
      sz=last_emreal2(taskid-2)-start_emreal2(taskid-2)+1
      
       !maxlocal_allthread=int(dble(sz)*dble(maxsize_elst(taskid)))
         maxlocal_allthread=0
         largest_nelst=0
           do i=start_emreal2(taskid-2),last_emreal2(taskid-2)
              !j=i-start_emreal2(taskid)+1
       !       print*,"i",i,"j",j,"nelst_recv(i)",nelst_recv(j)
              !maxlocal_allthread=maxlocal_allthread+nelst_recv(j)     
              maxlocal_allthread=maxlocal_allthread+nelst2(i)
              if(nelst2(i).gt.largest_nelst) then
                largest_nelst=nelst2(i)
              end if
           end do
      maxlocal_allthread_tau=3*maxlocal_allthread

      !print*,"taskid=",taskid,"sz=",sz
      !maxlocal_allthread=int(dble(sz)*dble(maxelst))
      

      allocate (elsttau1_index_tmp(4,largest_nelst,sz))
      allocate (elsttau2_index_tmp(4,largest_nelst,sz))

      allocate(dEd1_3tmp(9,largest_nelst,sz))
      allocate(dEd2_3tmp(9,largest_nelst,sz))
      allocate(dEp1_3tmp(9,largest_nelst,sz))
      allocate(dEp2_3tmp(9,largest_nelst,sz))

      allocate(taud1_3tmp(9,4,largest_nelst,sz))
      allocate(taud2_3tmp(9,4,largest_nelst,sz))
      allocate(taup1_3tmp(9,4,largest_nelst,sz))
      allocate(taup2_3tmp(9,4,largest_nelst,sz))


      allocate (frcztau1dtot_3tmp(9,largest_nelst,sz))
      allocate (frcytau1dtot_3tmp(9,largest_nelst,sz))
      allocate (frcxtau1dtot_3tmp(9,largest_nelst,sz))
      allocate (frcztau1ptot_3tmp(9,largest_nelst,sz))
      allocate (frcytau1ptot_3tmp(9,largest_nelst,sz))
      allocate (frcxtau1ptot_3tmp(9,largest_nelst,sz))
      allocate (frcztau2dtot_3tmp(9,largest_nelst,sz))
      allocate (frcytau2dtot_3tmp(9,largest_nelst,sz))
      allocate (frcxtau2dtot_3tmp(9,largest_nelst,sz))
      allocate (frcztau2ptot_3tmp(9,largest_nelst,sz))
      allocate (frcytau2ptot_3tmp(9,largest_nelst,sz))
      allocate (frcxtau2ptot_3tmp(9,largest_nelst,sz))

           

           call thg5_small_vir_nolistsend_Nn_noewtens2(
     &      fieldnpolet,fieldpnpolet)

        call mpi_isend(dEd1_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,2,mpi_comm_world,reqs22,ierr)
        call mpi_isend(dEd2_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,3,mpi_comm_world,reqs23,ierr)
        call mpi_isend(dEp1_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,4,mpi_comm_world,reqs24,ierr)
        call mpi_isend(dEp2_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,5,mpi_comm_world,reqs25,ierr)

        call mpi_isend(taud1_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,8,mpi_comm_world,reqs28,ierr)
        call mpi_isend(taud2_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,9,mpi_comm_world,reqs29,ierr)
       call mpi_isend(taup1_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,10,mpi_comm_world,reqs30,ierr)
       call mpi_isend(taup2_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,11,mpi_comm_world,reqs31,ierr)


        call mpi_isend(frcztau1dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,12,mpi_comm_world,reqs32,ierr)
        call mpi_isend(frcytau1dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,13,mpi_comm_world,reqs33,ierr)
        call mpi_isend(frcxtau1dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,14,mpi_comm_world,reqs34,ierr)
        call mpi_isend(frcztau1ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,15,mpi_comm_world,reqs35,ierr)
        call mpi_isend(frcytau1ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,16,mpi_comm_world,reqs36,ierr)
        call mpi_isend(frcxtau1ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,17,mpi_comm_world,reqs37,ierr)

        call mpi_isend(frcztau2dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,18,mpi_comm_world,reqs38,ierr)
        call mpi_isend(frcytau2dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,19,mpi_comm_world,reqs39,ierr)
        call mpi_isend(frcxtau2dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,20,mpi_comm_world,reqs40,ierr)
        call mpi_isend(frcztau2ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,21,mpi_comm_world,reqs41,ierr)
        call mpi_isend(frcytau2ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,22,mpi_comm_world,reqs42,ierr)
        call mpi_isend(frcxtau2ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,23,mpi_comm_world,reqs43,ierr)

        call mpi_isend(largest_nelst,1,mpi_integer,
     &       numtasks-1,24,mpi_comm_world,reqs44,ierr)
        call mpi_isend(elsttau1_index_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,26,mpi_comm_world,reqs46,ierr)
        call mpi_isend(elsttau2_index_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,28,mpi_comm_world,reqs48,ierr)

ccc      LEFT OFF HERE;ALSDJF;KJ;KJ

       call LoadBalInnerloop_ereal1d_3b_Perm_nolistsend(
     & emrealt,viremrealt,demrealt)

       !print*,"PermElec Real!"


        ! else
        !    call thg5_vir(
     &  !  emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

        ! end if

       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      else if((taskid.ge.numtasks_emreal+2)
     &         .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) then

          !  call loadbal_ehal1c_half(
     &    !    evtmp,devtmp,virevtmp)
            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)
       !print*,"VDW!"

      end if
      !       call mpi_barrier(mpi_comm_world,ierr)

         !    t1=mpi_Wtime()
       !emrealt=emrealt+evtmp
       !do i=1,n
       !   do j=1,3
       !      demrealt(j,i)=demrealt(j,i)+devtmp(j,i)
       !   end do
       !end do
       !do i=1,3
       !  do j=1,3
       !     viremrealt(j,i)=viremrealt(j,i)+virevtmp(j,i)
       !  end do
       !end do

       !call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       call mpi_ireduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs13,ierr)
       call mpi_ireduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs14,ierr)
       call mpi_ireduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs15,ierr)

c       print*,"PermElec Real After all reduce and send!"

       !call mpi_reduce(evtmp,ev,1,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(devtmp,dev,nvdw*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(virevtmp,virev,3*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       call mpi_ireduce(evtmp,ev,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs16,ierr)
       call mpi_ireduce(devtmp,dev,nvdw*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs17,ierr)
       call mpi_ireduce(virevtmp,virev,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs18,ierr)

       call mpi_reduce(fieldnpolet,fieldnpole,npole*3,mpi_real8,
     &      mpi_sum,master1,mpi_comm_world,ierr)
       call mpi_reduce(fieldpnpolet,fieldpnpole,npole*3,mpi_real8,
     &      mpi_sum,master1,mpi_comm_world,ierr)

      deallocate(demrealt)
      deallocate(fieldnpolet)
      deallocate(fieldpnpolet)
      deallocate(devtmp)

      if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then
       ! call mpi_wait(reqs21,stat,ierr)
        call mpi_wait(reqs22,stat,ierr)
        call mpi_wait(reqs23,stat,ierr)
        call mpi_wait(reqs24,stat,ierr)
        call mpi_wait(reqs25,stat,ierr)
       ! call mpi_wait(reqs26,stat,ierr)
       ! call mpi_wait(reqs27,stat,ierr)
        call mpi_wait(reqs28,stat,ierr)
        call mpi_wait(reqs29,stat,ierr)
        call mpi_wait(reqs30,stat,ierr)
        call mpi_wait(reqs31,stat,ierr)
        call mpi_wait(reqs32,stat,ierr)
        call mpi_wait(reqs33,stat,ierr)
        call mpi_wait(reqs34,stat,ierr)
        call mpi_wait(reqs35,stat,ierr)
        call mpi_wait(reqs36,stat,ierr)
        call mpi_wait(reqs37,stat,ierr)
        call mpi_wait(reqs38,stat,ierr)
        call mpi_wait(reqs39,stat,ierr)
        call mpi_wait(reqs40,stat,ierr)
        call mpi_wait(reqs41,stat,ierr)
        call mpi_wait(reqs42,stat,ierr)
        call mpi_wait(reqs43,stat,ierr)
        call mpi_wait(reqs44,stat,ierr)
       ! call mpi_wait(reqs45,stat,ierr)
        call mpi_wait(reqs46,stat,ierr)
       ! call mpi_wait(reqs47,stat,ierr)
        call mpi_wait(reqs48,stat,ierr)
      !  call mpi_wait(reqs49,stat,ierr)
      end if

         

      !if(taskid.eq.master1) then

      ! get the self-energy portion of the electrostatic field
 
      !  term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      !  do i = 1, npole
      !    do j = 1, 3
      !       fieldnpole(j,i)=fieldnpolet(j,i)
      !       fieldpnpole(j,i)=fieldpnpolet(j,i)
      !       demreal_tmp(j,i)=demrealt(j,i)
           ! fieldnpole(j,i) = fieldnpole(j,i) +fieldnpole_rcp(j,i)
           ! fieldpnpole(j,i) = fieldpnpole(j,i) +fieldnpole_rcp(j,i)
           ! fieldnpole(j,i) = fieldnpole(j,i) + term*rpole(j+1,i)
           ! fieldpnpole(j,i) = fieldpnpole(j,i) + term*rpole(j+1,i)
      !    end do
      !  end do

      ! end if
      if(taskid.eq.master) then
      ! call mpi_recv(em,1,mpi_real8,master1,1,mpi_comm_world,stat,ierr)
      ! call mpi_recv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
      !&     stat,ierr)
      ! call mpi_recv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
      !&     stat,ierr)
      !print*,"Before irecv"
      call mpi_irecv(em,1,mpi_real8,master1,1,mpi_comm_world,
     &    reqs4,ierr)
      call mpi_irecv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
     &    reqs5,ierr)
      call mpi_irecv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
     &    reqs6,ierr)
       print*,"After irecv"

      end if


         largest_nelst2=0

         do i=1,npole
           if(nelst2(i).gt.largest_nelst2) then
             largest_nelst2=nelst2(i)
           end if
         end do

      if(taskid.eq.numtasks-1) then



         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)

         if(allocated(dEd1_3)) deallocate(dEd1_3)
         if(allocated(dEd2_3)) deallocate(dEd2_3)
         if(allocated(dEp1_3)) deallocate(dEp1_3)
         if(allocated(dEp2_3)) deallocate(dEp2_3)

         if(allocated(taud1_3)) deallocate(taud1_3)
         if(allocated(taud2_3)) deallocate(taud2_3)
         if(allocated(taup1_3)) deallocate(taup1_3)
         if(allocated(taup2_3)) deallocate(taup2_3)

         if(allocated(frcztau1dtot_3)) deallocate(frcztau1dtot_3)
         if(allocated(frcytau1dtot_3)) deallocate(frcytau1dtot_3)
         if(allocated(frcxtau1dtot_3)) deallocate(frcxtau1dtot_3)
         if(allocated(frcztau1ptot_3)) deallocate(frcztau1ptot_3)
         if(allocated(frcytau1ptot_3)) deallocate(frcytau1ptot_3)
         if(allocated(frcxtau1ptot_3)) deallocate(frcxtau1ptot_3)
         if(allocated(frcztau2dtot_3)) deallocate(frcztau2dtot_3)
         if(allocated(frcytau2dtot_3)) deallocate(frcytau2dtot_3)
         if(allocated(frcxtau2dtot_3)) deallocate(frcxtau2dtot_3)
         if(allocated(frcztau2ptot_3)) deallocate(frcztau2ptot_3)
         if(allocated(frcytau2ptot_3)) deallocate(frcytau2ptot_3)
         if(allocated(frcxtau2ptot_3)) deallocate(frcxtau2ptot_3)

         allocate(elsttau1_index(4,largest_nelst2,npole))
         allocate(elsttau2_index(4,largest_nelst2,npole))

         allocate(dEd1_3(9,largest_nelst2,npole))
         allocate(dEd2_3(9,largest_nelst2,npole))
         allocate(dEp1_3(9,largest_nelst2,npole))
         allocate(dEp2_3(9,largest_nelst2,npole))

         allocate(taud1_3(9,4,largest_nelst2,npole))
         allocate(taup1_3(9,4,largest_nelst2,npole))
         allocate(taud2_3(9,4,largest_nelst2,npole))
         allocate(taup2_3(9,4,largest_nelst2,npole))

         allocate(frcztau1dtot_3(9,largest_nelst2,npole))
         allocate(frcytau1dtot_3(9,largest_nelst2,npole))
         allocate(frcxtau1dtot_3(9,largest_nelst2,npole))
         allocate(frcztau1ptot_3(9,largest_nelst2,npole))
         allocate(frcytau1ptot_3(9,largest_nelst2,npole))
         allocate(frcxtau1ptot_3(9,largest_nelst2,npole))
         allocate(frcztau2dtot_3(9,largest_nelst2,npole))
         allocate(frcytau2dtot_3(9,largest_nelst2,npole))
         allocate(frcxtau2dtot_3(9,largest_nelst2,npole))
         allocate(frcztau2ptot_3(9,largest_nelst2,npole))
         allocate(frcytau2ptot_3(9,largest_nelst2,npole))
         allocate(frcxtau2ptot_3(9,largest_nelst2,npole))


         do k=2,numtasks_emreal2+1

           call mpi_recv(largest_nelst,1,mpi_integer,k,24,
     &     mpi_comm_world,stat,ierr)

           sz=last_emreal2(k-2)-start_emreal2(k-2)+1

           allocate(dEd1_3tmp_rcv(9,largest_nelst,sz))
           allocate(dEd2_3tmp_rcv(9,largest_nelst,sz))
           allocate(dEp1_3tmp_rcv(9,largest_nelst,sz))
           allocate(dEp2_3tmp_rcv(9,largest_nelst,sz))

           allocate(frcztau1d_3rcv(9,largest_nelst,sz))
           allocate(frcytau1d_3rcv(9,largest_nelst,sz))
           allocate(frcxtau1d_3rcv(9,largest_nelst,sz))
           allocate(frcztau1p_3rcv(9,largest_nelst,sz))
           allocate(frcytau1p_3rcv(9,largest_nelst,sz))
           allocate(frcxtau1p_3rcv(9,largest_nelst,sz))

           allocate(frcztau2d_3rcv(9,largest_nelst,sz))
           allocate(frcytau2d_3rcv(9,largest_nelst,sz))
           allocate(frcxtau2d_3rcv(9,largest_nelst,sz))
           allocate(frcztau2p_3rcv(9,largest_nelst,sz))
           allocate(frcytau2p_3rcv(9,largest_nelst,sz))
           allocate(frcxtau2p_3rcv(9,largest_nelst,sz))

           allocate(taud1_3tmp_rcv(9,4,largest_nelst,sz))
           allocate(taud2_3tmp_rcv(9,4,largest_nelst,sz))
           allocate(taup1_3tmp_rcv(9,4,largest_nelst,sz))
           allocate(taup2_3tmp_rcv(9,4,largest_nelst,sz))

           allocate(elsttau1_index_rcv(4,largest_nelst,sz))
           allocate(elsttau2_index_rcv(4,largest_nelst,sz))
           
c         NEED TO CORRECT THE TAGS
           call mpi_recv(dEd1_3tmp_rcv,9*largest_nelst*sz,mpi_real8,k,2,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEd2_3tmp_rcv,9*largest_nelst*sz,mpi_real8,k,3,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEp1_3tmp_rcv,9*largest_nelst*sz,mpi_real8,k,4,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEp2_3tmp_rcv,9*largest_nelst*sz,mpi_real8,k,5,
     &     mpi_comm_world,stat,ierr)

        call mpi_recv(taud1_3tmp_rcv,9*4*largest_nelst*sz,mpi_real8,k,8,
     &     mpi_comm_world,stat,ierr)
        call mpi_recv(taud2_3tmp_rcv,9*4*largest_nelst*sz,mpi_real8,k,9,
     &     mpi_comm_world,stat,ierr)
       call mpi_recv(taup1_3tmp_rcv,9*4*largest_nelst*sz,mpi_real8,k,10,
     &     mpi_comm_world,stat,ierr)
       call mpi_recv(taup2_3tmp_rcv,9*4*largest_nelst*sz,mpi_real8,k,11,
     &     mpi_comm_world,stat,ierr)

           call mpi_recv(frcztau1d_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         12,mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau1d_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         13,mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau1d_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         14,mpi_comm_world,stat,ierr)
           call mpi_recv(frcztau1p_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         15,mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau1p_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         16,mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau1p_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         17,mpi_comm_world,stat,ierr)

           call mpi_recv(frcztau2d_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         18,mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau2d_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         19,mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau2d_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         20,mpi_comm_world,stat,ierr)
           call mpi_recv(frcztau2p_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         21,mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau2p_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         22,mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau2p_3rcv,9*largest_nelst*sz,mpi_real8,k,
     &         23,mpi_comm_world,stat,ierr)

           call mpi_recv(elsttau1_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,26,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau2_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,28,mpi_comm_world,stat,ierr)

c           do i=1,m
c              k2=ntpair_dEoffset+i
           do i=start_emreal2(k-2),last_emreal2(k-2)
              iter=i-start_emreal2(k-2)+1
              do k2=1,nelst2(i)
                do j=1,9
                  dEd1_3(j,k2,i)=dEd1_3tmp_rcv(j,k2,iter)
                  dEd2_3(j,k2,i)=dEd2_3tmp_rcv(j,k2,iter)
                  dEp1_3(j,k2,i)=dEp1_3tmp_rcv(j,k2,iter)
                  dEp2_3(j,k2,i)=dEp2_3tmp_rcv(j,k2,iter)

                  frcztau1dtot_3(j,k2,i)=frcztau1d_3rcv(j,k2,iter)
                  frcytau1dtot_3(j,k2,i)=frcytau1d_3rcv(j,k2,iter)
                  frcxtau1dtot_3(j,k2,i)=frcxtau1d_3rcv(j,k2,iter)
                  frcztau1ptot_3(j,k2,i)=frcztau1p_3rcv(j,k2,iter)
                  frcytau1ptot_3(j,k2,i)=frcytau1p_3rcv(j,k2,iter)
                  frcxtau1ptot_3(j,k2,i)=frcxtau1p_3rcv(j,k2,iter)

                  frcztau2dtot_3(j,k2,i)=frcztau2d_3rcv(j,k2,iter)
                  frcytau2dtot_3(j,k2,i)=frcytau2d_3rcv(j,k2,iter)
                  frcxtau2dtot_3(j,k2,i)=frcxtau2d_3rcv(j,k2,iter)
                  frcztau2ptot_3(j,k2,i)=frcztau2p_3rcv(j,k2,iter)
                  frcytau2ptot_3(j,k2,i)=frcytau2p_3rcv(j,k2,iter)
                  frcxtau2ptot_3(j,k2,i)=frcxtau2p_3rcv(j,k2,iter)
                end do
                
                do k3=1,4
                  elsttau1_index(k3,k2,i)=elsttau1_index_rcv(k3,k2,iter)
                  elsttau2_index(k3,k2,i)=elsttau2_index_rcv(k3,k2,iter)
                   do j=1,9
                      taud1_3(j,k3,k2,i)=taud1_3tmp_rcv(j,k3,k2,iter)
                      taud2_3(j,k3,k2,i)=taud2_3tmp_rcv(j,k3,k2,iter)
                      taup1_3(j,k3,k2,i)=taup1_3tmp_rcv(j,k3,k2,iter)
                      taup2_3(j,k3,k2,i)=taup2_3tmp_rcv(j,k3,k2,iter)
                   end do
                end do
              end do
           end do

           deallocate(dEd1_3tmp_rcv)
           deallocate(dEd2_3tmp_rcv)
           deallocate(dEp1_3tmp_rcv)
           deallocate(dEp2_3tmp_rcv)

           deallocate(frcztau1d_3rcv)
           deallocate(frcytau1d_3rcv)
           deallocate(frcxtau1d_3rcv)
           deallocate(frcztau1p_3rcv)
           deallocate(frcytau1p_3rcv)
           deallocate(frcxtau1p_3rcv)
           deallocate(frcztau2d_3rcv)
           deallocate(frcytau2d_3rcv)
           deallocate(frcxtau2d_3rcv)
           deallocate(frcztau2p_3rcv)
           deallocate(frcytau2p_3rcv)
           deallocate(frcxtau2p_3rcv)

           deallocate(elsttau1_index_rcv)
           deallocate(elsttau2_index_rcv)

           deallocate(taud1_3tmp_rcv)
           deallocate(taud2_3tmp_rcv)
           deallocate(taup1_3tmp_rcv)
           deallocate(taup2_3tmp_rcv)
         end do

      end if

      ! if(allocated(dEindextmp)) deallocate(dEindextmp)
      if(allocated(dEd1_3tmp)) deallocate(dEd1_3tmp)
      if(allocated(dEd2_3tmp)) deallocate(dEd2_3tmp)
      if(allocated(dEp1_3tmp)) deallocate(dEp1_3tmp)
      if(allocated(dEp2_3tmp)) deallocate(dEp2_3tmp)

      if(allocated(taud1_3tmp)) deallocate(taud1_3tmp)
      if(allocated(taud2_3tmp)) deallocate(taud2_3tmp)
      if(allocated(taup1_3tmp)) deallocate(taup1_3tmp)
      if(allocated(taup2_3tmp)) deallocate(taup2_3tmp)

      !if(allocated(tau1indextmp)) deallocate(tau1indextmp)
      !if(allocated(tau2indextmp)) deallocate(tau2indextmp)

      if(allocated(frcztau1dtot_3tmp)) deallocate (frcztau1dtot_3tmp)
      if(allocated(frcytau1dtot_3tmp)) deallocate (frcytau1dtot_3tmp)
      if(allocated(frcxtau1dtot_3tmp)) deallocate (frcxtau1dtot_3tmp)
      if(allocated(frcztau1ptot_3tmp)) deallocate (frcztau1ptot_3tmp)
      if(allocated(frcytau1ptot_3tmp)) deallocate (frcytau1ptot_3tmp)
      if(allocated(frcxtau1ptot_3tmp)) deallocate (frcxtau1ptot_3tmp)
      if(allocated(frcztau2dtot_3tmp)) deallocate (frcztau2dtot_3tmp)
      if(allocated(frcytau2dtot_3tmp)) deallocate (frcytau2dtot_3tmp)
      if(allocated(frcxtau2dtot_3tmp)) deallocate (frcxtau2dtot_3tmp)
      if(allocated(frcztau2ptot_3tmp)) deallocate (frcztau2ptot_3tmp)
      if(allocated(frcytau2ptot_3tmp)) deallocate (frcytau2ptot_3tmp)
      if(allocated(frcxtau2ptot_3tmp)) deallocate (frcxtau2ptot_3tmp)

      if(allocated(elsttau1_index_tmp)) deallocate(elsttau1_index_tmp)
      if(allocated(elsttau2_index_tmp)) deallocate(elsttau2_index_tmp)

         largest_nelst2=0

         do i=1,npole
           if(nelst2(i).gt.largest_nelst2) then
             largest_nelst2=nelst2(i)
           end if 
         end do


      if(taskid.ne.numtasks-1) then   
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)

         if(allocated(dEd1_3)) deallocate(dEd1_3)
         if(allocated(dEd2_3)) deallocate(dEd2_3)
         if(allocated(dEp1_3)) deallocate(dEp1_3)
         if(allocated(dEp2_3)) deallocate(dEp2_3)

         if(allocated(taud1_3)) deallocate(taud1_3)
         if(allocated(taud2_3)) deallocate(taud2_3)
         if(allocated(taup1_3)) deallocate(taup1_3)
         if(allocated(taup2_3)) deallocate(taup2_3)

         if(allocated(frcztau1dtot_3)) deallocate(frcztau1dtot_3)
         if(allocated(frcytau1dtot_3)) deallocate(frcytau1dtot_3)
         if(allocated(frcxtau1dtot_3)) deallocate(frcxtau1dtot_3)
         if(allocated(frcztau1ptot_3)) deallocate(frcztau1ptot_3)
         if(allocated(frcytau1ptot_3)) deallocate(frcytau1ptot_3)
         if(allocated(frcxtau1ptot_3)) deallocate(frcxtau1ptot_3)
         if(allocated(frcztau2dtot_3)) deallocate(frcztau2dtot_3)
         if(allocated(frcytau2dtot_3)) deallocate(frcytau2dtot_3)
         if(allocated(frcxtau2dtot_3)) deallocate(frcxtau2dtot_3)
         if(allocated(frcztau2ptot_3)) deallocate(frcztau2ptot_3)
         if(allocated(frcytau2ptot_3)) deallocate(frcytau2ptot_3)
         if(allocated(frcxtau2ptot_3)) deallocate(frcxtau2ptot_3)

         allocate(elsttau1_index(4,largest_nelst2,npole))
         allocate(elsttau2_index(4,largest_nelst2,npole))

         allocate(dEd1_3(9,largest_nelst2,npole))
         allocate(dEd2_3(9,largest_nelst2,npole))
         allocate(dEp1_3(9,largest_nelst2,npole))
         allocate(dEp2_3(9,largest_nelst2,npole))

         allocate(taud1_3(9,4,largest_nelst2,npole))
         allocate(taup1_3(9,4,largest_nelst2,npole))
         allocate(taud2_3(9,4,largest_nelst2,npole))
         allocate(taup2_3(9,4,largest_nelst2,npole))

         allocate(frcztau1dtot_3(9,largest_nelst2,npole))
         allocate(frcytau1dtot_3(9,largest_nelst2,npole))
         allocate(frcxtau1dtot_3(9,largest_nelst2,npole))
         allocate(frcztau1ptot_3(9,largest_nelst2,npole))
         allocate(frcytau1ptot_3(9,largest_nelst2,npole))
         allocate(frcxtau1ptot_3(9,largest_nelst2,npole))
         allocate(frcztau2dtot_3(9,largest_nelst2,npole))
         allocate(frcytau2dtot_3(9,largest_nelst2,npole))
         allocate(frcxtau2dtot_3(9,largest_nelst2,npole))
         allocate(frcztau2ptot_3(9,largest_nelst2,npole))
         allocate(frcytau2ptot_3(9,largest_nelst2,npole))
         allocate(frcxtau2ptot_3(9,largest_nelst2,npole))
      end if


      call mpi_bcast(dEd1_3,9*largest_nelst2*npole,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dEd2_3,9*largest_nelst2*npole,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)      
      call mpi_bcast(dEp1_3,9*largest_nelst2*npole,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dEp2_3,9*largest_nelst2*npole,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(taud1_3,9*4*largest_nelst2*npole,mpi_real8,
     &    numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(taud2_3,9*4*largest_nelst2*npole,mpi_real8,
     &    numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(taup1_3,9*4*largest_nelst2*npole,mpi_real8,
     &    numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(taup2_3,9*4*largest_nelst2*npole,mpi_real8,
     &    numtasks-1,mpi_comm_world,ierr)

      call mpi_bcast(frcztau1dtot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcytau1dtot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcxtau1dtot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcztau1ptot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcytau1ptot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcxtau1ptot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcztau2dtot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcytau2dtot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcxtau2dtot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcztau2ptot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcytau2ptot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(frcxtau2ptot_3,9*largest_nelst2*npole,mpi_real8,
     &     numtasks-1,mpi_comm_world,ierr)

      call mpi_bcast(elsttau1_index,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau2_index,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)

      call mpi_bcast(fieldnpole,3*npole,mpi_real8,master1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(fieldpnpole,3*npole,mpi_real8,master1,
     &   mpi_comm_world,ierr)
      return
      end

