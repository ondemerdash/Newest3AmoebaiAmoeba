      subroutine grad_vandr_ireduce_totfield_dEtensor_vir 
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
      use neigh
      use atoms
      use vdw
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
      
      sz=last_emreal2(taskid)-start_emreal2(taskid)+1
      
c      maxlocal_allthread=int(dble(sz)*dble(maxsize_elst(taskid)))
         maxlocal_allthread=0
         
           do i=start_emreal2(taskid),last_emreal2(taskid)
              !j=i-start_emreal2(taskid)+1
c             print*,"i",i,"j",j,"nelst_recv(i)",nelst_recv(j)
              !maxlocal_allthread=maxlocal_allthread+nelst_recv(j)     
              maxlocal_allthread=maxlocal_allthread+nelst(i)
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
           call thg5_small_vir_nolistsend(
     &   emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
       !print*,"PermElec Real!"

c        call thg6_save_small(
c     & emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

        ! else
        !    call thg5_vir(
     &  !  emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)
c        call thg6_save(
c     & emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

        ! end if

       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      else if((taskid.ge.numtasks_emreal)
     &         .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal))) then

c            call loadbal_ehal1c_half(
c     &        evtmp,devtmp,virevtmp)
            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)
       !print*,"VDW!"

      end if
            call mpi_barrier(mpi_comm_world,ierr)

         !    t1=mpi_Wtime()
      emrealt=emrealt+evtmp
      do i=1,n
         do j=1,3
            demrealt(j,i)=demrealt(j,i)+devtmp(j,i)
         end do
      end do
      do i=1,3
        do j=1,3
           viremrealt(j,i)=viremrealt(j,i)+virevtmp(j,i)
        end do
      end do

       call mpi_ireduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs(1),ierr)
       call mpi_ireduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs(2),ierr)
       call mpi_ireduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,reqs(3),ierr)

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

      end if

      
      deallocate(demrealt)
      deallocate(fieldnpolet)
      deallocate(fieldpnpolet)
      deallocate(devtmp)

c      call mpi_barrier(mpi_comm_world,ierr)

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

