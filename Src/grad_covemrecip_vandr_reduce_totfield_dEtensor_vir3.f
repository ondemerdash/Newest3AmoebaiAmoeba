      subroutine grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3 
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
      if(allocated(dEd1tmp)) deallocate(dEd1tmp)
      if(allocated(dEd2tmp)) deallocate(dEd2tmp)
      if(allocated(dEp1tmp)) deallocate(dEp1tmp)
      if(allocated(dEp2tmp)) deallocate(dEp2tmp)

      if(allocated(taud1tmp)) deallocate(taud1tmp)
      if(allocated(taud2tmp)) deallocate(taud2tmp)
      if(allocated(taup1tmp)) deallocate(taup1tmp)
      if(allocated(taup2tmp)) deallocate(taup2tmp)

      !if(allocated(tau1indextmp)) deallocate(tau1indextmp)
      !if(allocated(tau2indextmp)) deallocate(tau2indextmp)

      if(allocated(frcztau1d)) deallocate (frcztau1d)
      if(allocated(frcytau1d)) deallocate (frcytau1d)
      if(allocated(frcxtau1d)) deallocate (frcxtau1d)
      if(allocated(frcztau1p)) deallocate (frcztau1p)
      if(allocated(frcytau1p)) deallocate (frcytau1p)
      if(allocated(frcxtau1p)) deallocate (frcxtau1p)
      if(allocated(frcztau2d)) deallocate (frcztau2d)
      if(allocated(frcytau2d)) deallocate (frcytau2d)
      if(allocated(frcxtau2d)) deallocate (frcxtau2d)
      if(allocated(frcztau2p)) deallocate (frcztau2p)
      if(allocated(frcytau2p)) deallocate (frcytau2p)
      if(allocated(frcxtau2p)) deallocate (frcxtau2p)

      if(allocated(elstgrad_tmp)) deallocate(elstgrad_tmp)
      if(allocated(elsttau1_index_tmp)) deallocate(elsttau1_index_tmp)
      if(allocated(elsttau1_tmp)) deallocate(elsttau1_tmp)
      if(allocated(elsttau2_index_tmp)) deallocate(elsttau2_index_tmp)
      if(allocated(elsttau2_tmp)) deallocate(elsttau2_tmp)

      
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
      

      allocate (elstgrad_tmp(largest_nelst,sz))
      allocate (elsttau1_index_tmp(4,largest_nelst,sz))
      allocate (elsttau1_tmp(4,largest_nelst,sz))
      allocate (elsttau2_index_tmp(4,largest_nelst,sz))
      allocate (elsttau2_tmp(4,largest_nelst,sz))

      ! allocate(dEindextmp(2,maxlocal_allthread))
      allocate(dEd1tmp(9,maxlocal_allthread))
      allocate(dEd2tmp(9,maxlocal_allthread))
      allocate(dEp1tmp(9,maxlocal_allthread))
      allocate(dEp2tmp(9,maxlocal_allthread))

      allocate(taud1tmp(9,maxlocal_allthread_tau))
      allocate(taud2tmp(9,maxlocal_allthread_tau))
      allocate(taup1tmp(9,maxlocal_allthread_tau))
      allocate(taup2tmp(9,maxlocal_allthread_tau))

      ! allocate(tau1indextmp(2,maxlocal_allthread_tau))
      ! allocate(tau2indextmp(2,maxlocal_allthread_tau))

      allocate (frcztau1d(9,maxlocal_allthread))
      allocate (frcytau1d(9,maxlocal_allthread))
      allocate (frcxtau1d(9,maxlocal_allthread))
      allocate (frcztau1p(9,maxlocal_allthread))
      allocate (frcytau1p(9,maxlocal_allthread))
      allocate (frcxtau1p(9,maxlocal_allthread))
      allocate (frcztau2d(9,maxlocal_allthread))
      allocate (frcytau2d(9,maxlocal_allthread))
      allocate (frcxtau2d(9,maxlocal_allthread))
      allocate (frcztau2p(9,maxlocal_allthread))
      allocate (frcytau2p(9,maxlocal_allthread))
      allocate (frcxtau2p(9,maxlocal_allthread))

      ! call LoadBal_ereal1d_3b_Perm_sentlist2_totfield_dEtensorOmp2( 
     & ! emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet) 
      !call LoadBalInner_thg4omp_mpi_72015torque2 (
     & !emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
        ! if(int(npole/numtasks_emreal).le.500) then
        !   call thg5_small_vir_nolistsend(
     &  !  emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
           

           call thg5_small_vir_nolistsend_Nn_noewtens(
     &      fieldnpolet,fieldpnpolet)
       ! call thg5_small_vir_nolistsend_noewtens(
     & !  emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
        call mpi_isend(ntpair_dEtmp,1,mpi_integer,numtasks-1,1,
     &      mpi_comm_world,reqs21,ierr)
        call mpi_isend(dEd1tmp,9*ntpair_dEtmp,mpi_real8,numtasks-1,2,
     &      mpi_comm_world,reqs22,ierr)
        call mpi_isend(dEd2tmp,9*ntpair_dEtmp,mpi_real8,numtasks-1,3,
     &      mpi_comm_world,reqs23,ierr)
        call mpi_isend(dEp1tmp,9*ntpair_dEtmp,mpi_real8,numtasks-1,4,
     &      mpi_comm_world,reqs24,ierr)
        call mpi_isend(dEp2tmp,9*ntpair_dEtmp,mpi_real8,numtasks-1,5,
     &      mpi_comm_world,reqs25,ierr)

        call mpi_isend(ntpair_tau1tmp,1,mpi_integer,numtasks-1,6,
     &      mpi_comm_world,reqs26,ierr)
        call mpi_isend(ntpair_tau2tmp,1,mpi_integer,numtasks-1,7,
     &      mpi_comm_world,reqs27,ierr)
        call mpi_isend(taud1tmp,9*ntpair_tau1tmp,mpi_real8,numtasks-1,8,
     &      mpi_comm_world,reqs28,ierr)
        call mpi_isend(taud2tmp,9*ntpair_tau2tmp,mpi_real8,numtasks-1,9,
     &      mpi_comm_world,reqs29,ierr)
       call mpi_isend(taup1tmp,9*ntpair_tau1tmp,mpi_real8,numtasks-1,10,
     &      mpi_comm_world,reqs30,ierr)
       call mpi_isend(taup2tmp,9*ntpair_tau2tmp,mpi_real8,numtasks-1,11,
     &      mpi_comm_world,reqs31,ierr)


        call mpi_isend(frcztau1d,9*ntpair_dEtmp,mpi_real8,numtasks-1,12,
     &      mpi_comm_world,reqs32,ierr)
        call mpi_isend(frcytau1d,9*ntpair_dEtmp,mpi_real8,numtasks-1,13,
     &      mpi_comm_world,reqs33,ierr)
        call mpi_isend(frcxtau1d,9*ntpair_dEtmp,mpi_real8,numtasks-1,14,
     &      mpi_comm_world,reqs34,ierr)
        call mpi_isend(frcztau1p,9*ntpair_dEtmp,mpi_real8,numtasks-1,15,
     &      mpi_comm_world,reqs35,ierr)
        call mpi_isend(frcytau1p,9*ntpair_dEtmp,mpi_real8,numtasks-1,16,
     &      mpi_comm_world,reqs36,ierr)
        call mpi_isend(frcxtau1p,9*ntpair_dEtmp,mpi_real8,numtasks-1,17,
     &      mpi_comm_world,reqs37,ierr)

        call mpi_isend(frcztau2d,9*ntpair_dEtmp,mpi_real8,numtasks-1,18,
     &      mpi_comm_world,reqs38,ierr)
        call mpi_isend(frcytau2d,9*ntpair_dEtmp,mpi_real8,numtasks-1,19,
     &      mpi_comm_world,reqs39,ierr)
        call mpi_isend(frcxtau2d,9*ntpair_dEtmp,mpi_real8,numtasks-1,20,
     &      mpi_comm_world,reqs40,ierr)
        call mpi_isend(frcztau2p,9*ntpair_dEtmp,mpi_real8,numtasks-1,21,
     &      mpi_comm_world,reqs41,ierr)
        call mpi_isend(frcytau2p,9*ntpair_dEtmp,mpi_real8,numtasks-1,22,
     &      mpi_comm_world,reqs42,ierr)
        call mpi_isend(frcxtau2p,9*ntpair_dEtmp,mpi_real8,numtasks-1,23,
     &      mpi_comm_world,reqs43,ierr)

        call mpi_isend(largest_nelst,1,mpi_integer,
     &   numtasks-1,24,mpi_comm_world,reqs44,ierr)
        call mpi_isend(elstgrad_tmp,sz*largest_nelst,mpi_integer,
     &   numtasks-1,25,mpi_comm_world,reqs45,ierr)
        call mpi_isend(elsttau1_index_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,26,mpi_comm_world,reqs46,ierr)
        call mpi_isend(elsttau1_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,27,mpi_comm_world,reqs47,ierr) 
        call mpi_isend(elsttau2_index_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,28,mpi_comm_world,reqs48,ierr)
        call mpi_isend(elsttau2_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,29,mpi_comm_world,reqs49,ierr)


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
        call mpi_wait(reqs21,stat,ierr)
        call mpi_wait(reqs22,stat,ierr)
        call mpi_wait(reqs23,stat,ierr)
        call mpi_wait(reqs24,stat,ierr)
        call mpi_wait(reqs25,stat,ierr)
        call mpi_wait(reqs26,stat,ierr)
        call mpi_wait(reqs27,stat,ierr)
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
        call mpi_wait(reqs45,stat,ierr)
        call mpi_wait(reqs46,stat,ierr)
        call mpi_wait(reqs47,stat,ierr)
        call mpi_wait(reqs48,stat,ierr)
        call mpi_wait(reqs49,stat,ierr)
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


         ntpair_dEoffset=0
         ntpair_tau1offset=0
         ntpair_tau2offset=0

         allocate(ntpair_dE_array(0:numtasks_emreal2-1))
         allocate(ntpair_tau1_array(0:numtasks_emreal2-1))
         allocate(ntpair_tau2_array(0:numtasks_emreal2-1))

         do k=2,numtasks_emreal2+1
           call mpi_recv(ntpair_dEtmp_rcv,1,mpi_integer,k,1,
     &     mpi_comm_world,stat,ierr)
           ntpair_dEoffset=ntpair_dEoffset+ntpair_dEtmp_rcv
           ntpair_dE_array(k-2)=ntpair_dEtmp_rcv

           call mpi_recv(ntpair_tau1tmp_rcv,1,mpi_integer,k,6,
     &     mpi_comm_world,stat,ierr)
           ntpair_tau1offset=ntpair_tau1offset+ntpair_tau1tmp_rcv
           ntpair_tau1_array(k-2)=ntpair_tau1tmp_rcv

           call mpi_recv(ntpair_tau2tmp_rcv,1,mpi_integer,k,7,
     &     mpi_comm_world,stat,ierr)
           ntpair_tau2offset=ntpair_tau2offset+ntpair_tau2tmp_rcv
           ntpair_tau2_array(k-2)=ntpair_tau2tmp_rcv
         end do

         ntpair_dE=ntpair_dEoffset
         ntpair_tau1=ntpair_tau1offset
         ntpair_tau2=ntpair_tau2offset

         if(allocated(elstgrad)) deallocate(elstgrad)
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau1)) deallocate(elsttau1)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)
         if(allocated(elsttau2)) deallocate(elsttau2)

         if(allocated(dEd1)) deallocate(dEd1)
         if(allocated(dEd2)) deallocate(dEd2)
         if(allocated(dEp1)) deallocate(dEp1)
         if(allocated(dEp2)) deallocate(dEp2)

         if(allocated(taud1)) deallocate(taud1)
         if(allocated(taud2)) deallocate(taud2)
         if(allocated(taup1)) deallocate(taup1)
         if(allocated(taup2)) deallocate(taup2)

         if(allocated(frcztau1dtot)) deallocate(frcztau1dtot)
         if(allocated(frcytau1dtot)) deallocate(frcytau1dtot)
         if(allocated(frcxtau1dtot)) deallocate(frcxtau1dtot)
         if(allocated(frcztau1ptot)) deallocate(frcztau1ptot)
         if(allocated(frcytau1ptot)) deallocate(frcytau1ptot)
         if(allocated(frcxtau1ptot)) deallocate(frcxtau1ptot)
         if(allocated(frcztau2dtot)) deallocate(frcztau2dtot)
         if(allocated(frcytau2dtot)) deallocate(frcytau2dtot)
         if(allocated(frcxtau2dtot)) deallocate(frcxtau2dtot)
         if(allocated(frcztau2ptot)) deallocate(frcztau2ptot)
         if(allocated(frcytau2ptot)) deallocate(frcytau2ptot)
         if(allocated(frcxtau2ptot)) deallocate(frcxtau2ptot)

         allocate(elstgrad(largest_nelst2,npole))
         allocate(elsttau1_index(4,largest_nelst2,npole))
         allocate(elsttau1(4,largest_nelst2,npole))
         allocate(elsttau2_index(4,largest_nelst2,npole))
         allocate(elsttau2(4,largest_nelst2,npole))

         allocate(dEd1(9,ntpair_dE))
         allocate(dEd2(9,ntpair_dE))
         allocate(dEp1(9,ntpair_dE))
         allocate(dEp2(9,ntpair_dE))

         allocate(taud1(9,ntpair_tau1))
         allocate(taup1(9,ntpair_tau1))
         allocate(taud2(9,ntpair_tau2))
         allocate(taup2(9,ntpair_tau2))

         allocate(frcztau1dtot(9,ntpair_dE))
         allocate(frcytau1dtot(9,ntpair_dE))
         allocate(frcxtau1dtot(9,ntpair_dE))
         allocate(frcztau1ptot(9,ntpair_dE))
         allocate(frcytau1ptot(9,ntpair_dE))
         allocate(frcxtau1ptot(9,ntpair_dE))
         allocate(frcztau2dtot(9,ntpair_dE))
         allocate(frcytau2dtot(9,ntpair_dE))
         allocate(frcxtau2dtot(9,ntpair_dE))
         allocate(frcztau2ptot(9,ntpair_dE))
         allocate(frcytau2ptot(9,ntpair_dE))
         allocate(frcxtau2ptot(9,ntpair_dE))

         ntpair_dEoffset=0
         ntpair_tau1offset=0
         ntpair_tau2offset=0

         do k=2,numtasks_emreal2+1
           m=ntpair_dE_array(k-2)
           allocate(dEd1tmp_rcv(9,m))
           allocate(dEd2tmp_rcv(9,m))
           allocate(dEp1tmp_rcv(9,m))
           allocate(dEp2tmp_rcv(9,m))

           allocate(frcztau1d_rcv(9,m))
           allocate(frcytau1d_rcv(9,m))
           allocate(frcxtau1d_rcv(9,m))
           allocate(frcztau1p_rcv(9,m))
           allocate(frcytau1p_rcv(9,m))
           allocate(frcxtau1p_rcv(9,m))

           allocate(frcztau2d_rcv(9,m))
           allocate(frcytau2d_rcv(9,m))
           allocate(frcxtau2d_rcv(9,m))
           allocate(frcztau2p_rcv(9,m))
           allocate(frcytau2p_rcv(9,m))
           allocate(frcxtau2p_rcv(9,m))
           
c         NEED TO CORRECT THE TAGS
           call mpi_recv(dEd1tmp_rcv,9*m,mpi_real8,k,2,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEd2tmp_rcv,9*m,mpi_real8,k,3,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEp1tmp_rcv,9*m,mpi_real8,k,4,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEp2tmp_rcv,9*m,mpi_real8,k,5,
     &     mpi_comm_world,stat,ierr)

           call mpi_recv(frcztau1d_rcv,9*m,mpi_real8,k,12,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau1d_rcv,9*m,mpi_real8,k,13,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau1d_rcv,9*m,mpi_real8,k,14,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcztau1p_rcv,9*m,mpi_real8,k,15,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau1p_rcv,9*m,mpi_real8,k,16,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau1p_rcv,9*m,mpi_real8,k,17,
     &     mpi_comm_world,stat,ierr)

           call mpi_recv(frcztau2d_rcv,9*m,mpi_real8,k,18,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau2d_rcv,9*m,mpi_real8,k,19,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau2d_rcv,9*m,mpi_real8,k,20,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcztau2p_rcv,9*m,mpi_real8,k,21,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau2p_rcv,9*m,mpi_real8,k,22,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau2p_rcv,9*m,mpi_real8,k,23,
     &     mpi_comm_world,stat,ierr)

           do i=1,m
              k2=ntpair_dEoffset+i
              do j=1,9
                dEd1(j,k2)=dEd1tmp_rcv(j,i)
                dEd2(j,k2)=dEd2tmp_rcv(j,i)
                dEp1(j,k2)=dEp1tmp_rcv(j,i)
                dEp2(j,k2)=dEp2tmp_rcv(j,i)

                frcztau1dtot(j,k2)=frcztau1d_rcv(j,i)
                frcytau1dtot(j,k2)=frcytau1d_rcv(j,i)
                frcxtau1dtot(j,k2)=frcxtau1d_rcv(j,i)
                frcztau1ptot(j,k2)=frcztau1p_rcv(j,i)
                frcytau1ptot(j,k2)=frcytau1p_rcv(j,i)
                frcxtau1ptot(j,k2)=frcxtau1p_rcv(j,i)

                frcztau2dtot(j,k2)=frcztau2d_rcv(j,i)
                frcytau2dtot(j,k2)=frcytau2d_rcv(j,i)
                frcxtau2dtot(j,k2)=frcxtau2d_rcv(j,i)
                frcztau2ptot(j,k2)=frcztau2p_rcv(j,i)
                frcytau2ptot(j,k2)=frcytau2p_rcv(j,i)
                frcxtau2ptot(j,k2)=frcxtau2p_rcv(j,i)
              end do
           end do

           deallocate(dEd1tmp_rcv)
           deallocate(dEd2tmp_rcv)
           deallocate(dEp1tmp_rcv)
           deallocate(dEp2tmp_rcv)

           deallocate(frcztau1d_rcv)
           deallocate(frcytau1d_rcv)
           deallocate(frcxtau1d_rcv)
           deallocate(frcztau1p_rcv)
           deallocate(frcytau1p_rcv)
           deallocate(frcxtau1p_rcv)
           deallocate(frcztau2d_rcv)
           deallocate(frcytau2d_rcv)
           deallocate(frcxtau2d_rcv)
           deallocate(frcztau2p_rcv)
           deallocate(frcytau2p_rcv)
           deallocate(frcxtau2p_rcv)


           m1=ntpair_tau1_array(k-2)
           m2=ntpair_tau2_array(k-2)
           allocate(taud1tmp_rcv(9,m1))
           allocate(taud2tmp_rcv(9,m2))
           allocate(taup1tmp_rcv(9,m1))
           allocate(taup2tmp_rcv(9,m2))

           call mpi_recv(taud1tmp_rcv,9*m1,mpi_real8,k,8,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(taud2tmp_rcv,9*m2,mpi_real8,k,9,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(taup1tmp_rcv,9*m1,mpi_real8,k,10,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(taup2tmp_rcv,9*m2,mpi_real8,k,11,
     &     mpi_comm_world,stat,ierr)

           do i=1,m1
              k2=ntpair_tau1offset+i
              do j=1,9
               taud1(j,k2)=taud1tmp_rcv(j,i)
               taup1(j,k2)=taup1tmp_rcv(j,i)
              end do
           end do

           do i=1,m2
              k2=ntpair_tau2offset+i
              do j=1,9
               taud2(j,k2)=taud2tmp_rcv(j,i)
               taup2(j,k2)=taup2tmp_rcv(j,i)
              end do
           end do

           deallocate(taud1tmp_rcv)
           deallocate(taud2tmp_rcv)
           deallocate(taup1tmp_rcv)
           deallocate(taup2tmp_rcv)

           call mpi_recv(largest_nelst,1,mpi_integer,k,24,
     &     mpi_comm_world,stat,ierr)

           sz=last_emreal2(k-2)-start_emreal2(k-2)+1

           allocate(elstgrad_rcv(largest_nelst,sz))
           allocate(elsttau1_index_rcv(4,largest_nelst,sz))
           allocate(elsttau1_rcv(4,largest_nelst,sz))
           allocate(elsttau2_index_rcv(4,largest_nelst,sz))
           allocate(elsttau2_rcv(4,largest_nelst,sz))

           call mpi_recv(elstgrad_rcv,sz*largest_nelst,
     &     mpi_integer,k,25,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau1_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,26,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau1_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,27,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau2_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,28,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau2_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,29,mpi_comm_world,stat,ierr)
 
           do i=start_emreal2(k-2),last_emreal2(k-2)
              iter=i-start_emreal2(k-2)+1
              do k2=1,nelst2(i)
                 elstgrad(k2,i)=elstgrad_rcv(k2,iter)+ntpair_dEoffset
                do k3=1,4
                  elsttau1_index(k3,k2,i)=elsttau1_index_rcv(k3,k2,iter)
                  elsttau2_index(k3,k2,i)=elsttau2_index_rcv(k3,k2,iter)

                  if(elsttau1_index_rcv(k3,k2,iter).ne.0) then
                    elsttau1(k3,k2,i)=elsttau1_rcv(k3,k2,iter)
     &              +ntpair_tau1offset
                  else
                    elsttau1(k3,k2,i)=elsttau1_rcv(k3,k2,iter) 
                  end if

                  if(elsttau2_index_rcv(k3,k2,iter).ne.0) then
                    elsttau2(k3,k2,i)=elsttau2_rcv(k3,k2,iter)
     &              +ntpair_tau2offset
                  else
                    elsttau2(k3,k2,i)=elsttau2_rcv(k3,k2,iter)
                  end if
                end do
              end do
           end do

           deallocate(elstgrad_rcv)
           deallocate(elsttau1_index_rcv)
           deallocate(elsttau1_rcv)
           deallocate(elsttau2_index_rcv)
           deallocate(elsttau2_rcv)

           ntpair_dEoffset=ntpair_dEoffset+m
           ntpair_tau1offset=ntpair_tau1offset+m1
           ntpair_tau2offset=ntpair_tau2offset+m2
         end do

      end if

      if(allocated(dEd1tmp)) deallocate(dEd1tmp)
      if(allocated(dEd2tmp)) deallocate(dEd2tmp)
      if(allocated(dEp1tmp)) deallocate(dEp1tmp)
      if(allocated(dEp2tmp)) deallocate(dEp2tmp)

      if(allocated(taud1tmp)) deallocate(taud1tmp)
      if(allocated(taud2tmp)) deallocate(taud2tmp)
      if(allocated(taup1tmp)) deallocate(taup1tmp)
      if(allocated(taup2tmp)) deallocate(taup2tmp)

      !if(allocated(tau1indextmp)) deallocate(tau1indextmp)
      !if(allocated(tau2indextmp)) deallocate(tau2indextmp)

      if(allocated(frcztau1d)) deallocate (frcztau1d)
      if(allocated(frcytau1d)) deallocate (frcytau1d)
      if(allocated(frcxtau1d)) deallocate (frcxtau1d)
      if(allocated(frcztau1p)) deallocate (frcztau1p)
      if(allocated(frcytau1p)) deallocate (frcytau1p)
      if(allocated(frcxtau1p)) deallocate (frcxtau1p)
      if(allocated(frcztau2d)) deallocate (frcztau2d)
      if(allocated(frcytau2d)) deallocate (frcytau2d)
      if(allocated(frcxtau2d)) deallocate (frcxtau2d)
      if(allocated(frcztau2p)) deallocate (frcztau2p)
      if(allocated(frcytau2p)) deallocate (frcytau2p)
      if(allocated(frcxtau2p)) deallocate (frcxtau2p)

      if(allocated(elstgrad_tmp)) deallocate(elstgrad_tmp)
      if(allocated(elsttau1_index_tmp)) deallocate(elsttau1_index_tmp)
      if(allocated(elsttau1_tmp)) deallocate(elsttau1_tmp)
      if(allocated(elsttau2_index_tmp)) deallocate(elsttau2_index_tmp)
      if(allocated(elsttau2_tmp)) deallocate(elsttau2_tmp)

         largest_nelst2=0

         do i=1,npole
           if(nelst2(i).gt.largest_nelst2) then
             largest_nelst2=nelst2(i)
           end if 
         end do

       call mpi_bcast(ntpair_dE,1,mpi_integer,numtasks-1,
     &   mpi_comm_world,ierr)
       call mpi_bcast(ntpair_tau1,1,mpi_integer,numtasks-1,
     &   mpi_comm_world,ierr)
       call mpi_bcast(ntpair_tau2,1,mpi_integer,numtasks-1,
     &   mpi_comm_world,ierr)

      if(taskid.ne.numtasks-1) then   
         if(allocated(elstgrad)) deallocate(elstgrad)
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau1)) deallocate(elsttau1)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)
         if(allocated(elsttau2)) deallocate(elsttau2)

         if(allocated(dEd1)) deallocate(dEd1)
         if(allocated(dEd2)) deallocate(dEd2)
         if(allocated(dEp1)) deallocate(dEp1)
         if(allocated(dEp2)) deallocate(dEp2)

         if(allocated(taud1)) deallocate(taud1)
         if(allocated(taud2)) deallocate(taud2)
         if(allocated(taup1)) deallocate(taup1)
         if(allocated(taup2)) deallocate(taup2)

         if(allocated(frcztau1dtot)) deallocate(frcztau1dtot)
         if(allocated(frcytau1dtot)) deallocate(frcytau1dtot)
         if(allocated(frcxtau1dtot)) deallocate(frcxtau1dtot)
         if(allocated(frcztau1ptot)) deallocate(frcztau1ptot)
         if(allocated(frcytau1ptot)) deallocate(frcytau1ptot)
         if(allocated(frcxtau1ptot)) deallocate(frcxtau1ptot)
         if(allocated(frcztau2dtot)) deallocate(frcztau2dtot)
         if(allocated(frcytau2dtot)) deallocate(frcytau2dtot)
         if(allocated(frcxtau2dtot)) deallocate(frcxtau2dtot)
         if(allocated(frcztau2ptot)) deallocate(frcztau2ptot)
         if(allocated(frcytau2ptot)) deallocate(frcytau2ptot)
         if(allocated(frcxtau2ptot)) deallocate(frcxtau2ptot)

         allocate(elstgrad(largest_nelst2,npole))
         allocate(elsttau1_index(4,largest_nelst2,npole))
         allocate(elsttau1(4,largest_nelst2,npole))
         allocate(elsttau2_index(4,largest_nelst2,npole))
         allocate(elsttau2(4,largest_nelst2,npole))

         allocate(dEd1(9,ntpair_dE))
         allocate(dEd2(9,ntpair_dE))
         allocate(dEp1(9,ntpair_dE))
         allocate(dEp2(9,ntpair_dE))

         allocate(taud1(9,ntpair_tau1))
         allocate(taup1(9,ntpair_tau1))
         allocate(taud2(9,ntpair_tau2))
         allocate(taup2(9,ntpair_tau2))

         allocate(frcztau1dtot(9,ntpair_dE))
         allocate(frcytau1dtot(9,ntpair_dE))
         allocate(frcxtau1dtot(9,ntpair_dE)) 
         allocate(frcztau1ptot(9,ntpair_dE)) 
         allocate(frcytau1ptot(9,ntpair_dE)) 
         allocate(frcxtau1ptot(9,ntpair_dE)) 
         allocate(frcztau2dtot(9,ntpair_dE)) 
         allocate(frcytau2dtot(9,ntpair_dE)) 
         allocate(frcxtau2dtot(9,ntpair_dE)) 
         allocate(frcztau2ptot(9,ntpair_dE)) 
         allocate(frcytau2ptot(9,ntpair_dE)) 
         allocate(frcxtau2ptot(9,ntpair_dE)) 
      end if


      call mpi_bcast(dEd1,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dEd2,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)      
      call mpi_bcast(dEp1,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dEp2,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(taud1,9*ntpair_tau1,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taud2,9*ntpair_tau2,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taup1,9*ntpair_tau1,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taup2,9*ntpair_tau2,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(frcztau1dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau1dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau1dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcztau1ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau1ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau1ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcztau2dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau2dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau2dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcztau2ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau2ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau2ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(elstgrad,largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau1_index,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau1,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau2_index,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau2,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)

      call mpi_bcast(fieldnpole,3*npole,mpi_real8,master1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(fieldpnpole,3*npole,mpi_real8,master1,
     &   mpi_comm_world,ierr)
      return
      end

