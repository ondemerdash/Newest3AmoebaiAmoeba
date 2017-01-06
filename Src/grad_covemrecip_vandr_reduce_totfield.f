      subroutine grad_covemrecip_vandr_reduce_totfield
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
      use polar
      use aprx
      use deriv3b
      use chgpot   
      use paremneigh
      use neigh
      use atoms
      use vdw
      use pme
      use chunks
      implicit none
      include 'mpif.h'
      real*8 emrealt,viremrealt(3,3)
      real*8, allocatable :: demrealt(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6
      real*8 t11,t21,t31,t41,t51,t61
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
      logical first
      save first
      data first  / .true. /

      allocate (demrealt(3,npole))
c      allocate (fieldnpolet(3,npole))
c      allocate (fieldpnpolet(3,npole))

      allocate (devtmp(3,nvdw))

      f = electric / dielec
      evtmp=0.0d0
      emrealt=0.0d0
      do i=1,n
         do j=1,3
            demrealt(j,i)=0.0d0
c            fieldnpolet(j,i)=0.0d0
c            fieldpnpolet(j,i) = 0.0d0
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
      !print*,"numtasks_emreal2=",numtasks_emreal2

      if(taskid.eq.master) then
         call gradient_covalent2

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

        if(embedtyp.ne.'O') then
        call emrecip1_3b_Perm2_totfield
        else
        !  t5=mpi_Wtime()
        call emrecip1_3b_Perm2
        !call emrecip1_3b_Perm2_noqpole_nodpole
        !  t6=mpi_Wtime()    
        !print*,"Computation PermElec PME taskid",taskid,"time=",t6-t5

        end if

       ! t5=mpi_Wtime()
      
      call mpi_isend(em,1,mpi_real8,master,1,mpi_comm_world,
     &      reqs1,ierr)
      call mpi_isend(dem,3*n,mpi_real8,master,2,mpi_comm_world,
     &      reqs2,ierr)
      call mpi_isend(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
     &      reqs3,ierr)
      !  t6=mpi_Wtime()
      !print*,"Isend PermElec PME taskid",taskid,"time=",t6-t5

c       call mpi_send(em,1,mpi_real8,master,1,mpi_comm_world,
c     &      ierr)
c      call mpi_send(dem,3*n,mpi_real8,master,2,mpi_comm_world,
c     &      ierr)
c      call mpi_send(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
c     &      ierr)
 
c      if(taskid.lt.numtasks_emreal2) then
      !else if(taskid.eq.master2) then
       !call dfield0b_totfield
      else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then

c           call thg5_small_vir_nolistsend(
c     &   emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
      !t1=mpi_Wtime()
      call LoadBalInnerloop_ereal1d_3b_Perm_nolistsend(
     & emrealt,viremrealt,demrealt)
      !t2=mpi_Wtime()
      !print*,"Computation PermElec Real taskid",taskid,"time=",t2-t1
      
       !   if((taskid.eq.master2).and.(embedtyp.ne.'O')) then
       !     call dfield0b_totfield
       !   end if
       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      else if((taskid.ge.numtasks_emreal+2)
     &         .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) then
      !t3=mpi_Wtime()

            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)
      !t4=mpi_Wtime()
      !print*,"Computation VdW taskid",taskid,"time=",t4-t3
      end if


       !call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
      !t11=mpi_Wtime()
       call mpi_ireduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs13,ierr)
       call mpi_ireduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs14,ierr)
       call mpi_ireduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs15,ierr)
      !t21=mpi_Wtime()
      !print*,"Ireduce PermElec Real taskid",taskid,"time=",t21-t11

c       print*,"PermElec Real After all reduce and send!"

       !call mpi_reduce(evtmp,ev,1,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(devtmp,dev,nvdw*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
       !call mpi_reduce(virevtmp,virev,3*3,mpi_real8,
     & !     mpi_sum,master,mpi_comm_world,ierr)
      !t31=mpi_Wtime()
       call mpi_ireduce(evtmp,ev,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs16,ierr)
       call mpi_ireduce(devtmp,dev,nvdw*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs17,ierr)
       call mpi_ireduce(virevtmp,virev,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs18,ierr)
      !t41=mpi_Wtime()
      !print*,"Ireduce VdW taskid",taskid,"time=",t41-t31

      deallocate(demrealt)
      deallocate(devtmp)
c      call mpi_barrier(mpi_comm_world,ierr)
      if(taskid.eq.master) then
c      call mpi_recv(em,1,mpi_real8,master1,1,mpi_comm_world,stat,ierr)
c      call mpi_recv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
c     &     stat,ierr)
c      call mpi_recv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
c     &     stat,ierr)
c      print*,"Before irecv"
      !t51=mpi_Wtime()
      call mpi_irecv(em,1,mpi_real8,master1,1,mpi_comm_world,
     &    reqs4,ierr)
      call mpi_irecv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
     &    reqs5,ierr)
      call mpi_irecv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
     &    reqs6,ierr)
      !t61=mpi_Wtime()
      !print*,"Irecv PME taskid",taskid,"time=",t61-t51

c       print*,"After irecv"
c
      end if

       !if(embedtyp.ne.'O') then
       !  call mpi_bcast(fieldnpole,3*n,mpi_real8,master2,
     & !     mpi_comm_world,ierr)
       !  call mpi_bcast(fieldpnpole,3*n,mpi_real8,master2,
     & !     mpi_comm_world,ierr)
       !end if
      return
      end

