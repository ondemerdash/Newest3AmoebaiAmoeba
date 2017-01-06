      subroutine grad_covemrecip_vandr_reduce_totfield_dEtensor_vir6 
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

      !print*,"Grad5!"      

      allocate (demrealt(3,npole))

      allocate (devtmp(3,nvdw))

      f = electric / dielec
      evtmp=0.0d0
      emrealt=0.0d0
      do i=1,n
         do j=1,3
            demrealt(j,i)=0.0d0
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
c         largest_nelst=0
c           do i=1,npole
c              if(nelst2(i).gt.largest_nelst) then
c                largest_nelst=nelst2(i)
c              end if
c           end do
c
c
c      allocate (elsttau1_index(4,largest_nelst,npole))
c      allocate (elsttau2_index(4,largest_nelst,npole))
c
c      allocate(dEd1_3(9,largest_nelst,npole))
c      allocate(dEd2_3(9,largest_nelst,npole))
c      allocate(dEp1_3(9,largest_nelst,npole))
c      allocate(dEp2_3(9,largest_nelst,npole))
c
c      allocate(taud1_3(9,4,largest_nelst,npole))
c      allocate(taud2_3(9,4,largest_nelst,npole))
c      allocate(taup1_3(9,4,largest_nelst,npole))
c      allocate(taup2_3(9,4,largest_nelst,npole))
c
c
c      allocate (frcztau1dtot_3(9,largest_nelst,npole))
c      allocate (frcytau1dtot_3(9,largest_nelst,npole))
c      allocate (frcxtau1dtot_3(9,largest_nelst,npole))
c      allocate (frcztau1ptot_3(9,largest_nelst,npole))
c      allocate (frcytau1ptot_3(9,largest_nelst,npole))
c      allocate (frcxtau1ptot_3(9,largest_nelst,npole))
c      allocate (frcztau2dtot_3(9,largest_nelst,npole))
c      allocate (frcytau2dtot_3(9,largest_nelst,npole))
c      allocate (frcxtau2dtot_3(9,largest_nelst,npole))
c      allocate (frcztau2ptot_3(9,largest_nelst,npole))
c      allocate (frcytau2ptot_3(9,largest_nelst,npole))
c      allocate (frcxtau2ptot_3(9,largest_nelst,npole))

      if(taskid.eq.master) then
         largest_nelst=0
           do i=1,npole
              if(nelst2(i).gt.largest_nelst) then
                largest_nelst=nelst2(i)
              end if
           end do


      allocate (elsttau1_index(4,largest_nelst,npole))
      allocate (elsttau2_index(4,largest_nelst,npole))

      allocate(dEd1_3(9,largest_nelst,npole))
      allocate(dEd2_3(9,largest_nelst,npole))
      allocate(dEp1_3(9,largest_nelst,npole))
      allocate(dEp2_3(9,largest_nelst,npole))

      allocate(taud1_3(9,4,largest_nelst,npole))
      allocate(taud2_3(9,4,largest_nelst,npole))
      allocate(taup1_3(9,4,largest_nelst,npole))
      allocate(taup2_3(9,4,largest_nelst,npole))


      allocate (frcztau1dtot_3(9,largest_nelst,npole))
      allocate (frcytau1dtot_3(9,largest_nelst,npole))
      allocate (frcxtau1dtot_3(9,largest_nelst,npole))
      allocate (frcztau1ptot_3(9,largest_nelst,npole))
      allocate (frcytau1ptot_3(9,largest_nelst,npole))
      allocate (frcxtau1ptot_3(9,largest_nelst,npole))
      allocate (frcztau2dtot_3(9,largest_nelst,npole))
      allocate (frcytau2dtot_3(9,largest_nelst,npole))
      allocate (frcxtau2dtot_3(9,largest_nelst,npole))
      allocate (frcztau2ptot_3(9,largest_nelst,npole))
      allocate (frcytau2ptot_3(9,largest_nelst,npole))
      allocate (frcxtau2ptot_3(9,largest_nelst,npole))

         call gradient_covalent2
         call thg5_small_vir_nolistsend_Nn_noewtens_notenssends2

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

       call LoadBalInnerloop_ereal1d_3b_Perm_nolistsend(
     & emrealt,viremrealt,demrealt)

       !print*,"PermElec Real!"



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


      deallocate(demrealt)
      deallocate(devtmp)


         


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


      return
      end

