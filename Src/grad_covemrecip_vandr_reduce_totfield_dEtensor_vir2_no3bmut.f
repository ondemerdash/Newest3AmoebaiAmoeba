      subroutine grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2no3b
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
         print*,"Finished gradient_covalent2 eub=",eub
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

c        t5=mpi_Wtime()
        call emrecip1_3b_Perm2_totfield
          print*,"After emrecip1_3b_Perm2_totfield",em
c        t6=mpi_Wtime()
c        print*,"PME time",t6-t5

      call mpi_isend(em,1,mpi_real8,master,1,mpi_comm_world,
     &      reqs1,ierr)
      call mpi_isend(dem,3*n,mpi_real8,master,2,mpi_comm_world,
     &      reqs2,ierr)
      call mpi_isend(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
     &      reqs3,ierr)

c       call mpi_send(em,1,mpi_real8,master,1,mpi_comm_world,
c     &      ierr)
c      call mpi_send(dem,3*n,mpi_real8,master,2,mpi_comm_world,
c     &      ierr)
c      call mpi_send(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
c     &      ierr)
       print*,"After mpi_isend of emrecip1 vars"
 
c      if(taskid.lt.numtasks_emreal2) then
      else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then

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
      
      sz=last_emreal2(taskid-2)-start_emreal2(taskid-2)+1
      
c      maxlocal_allthread=int(dble(sz)*dble(maxsize_elst(taskid)))
         maxlocal_allthread=0
         
           do i=start_emreal2(taskid-2),last_emreal2(taskid-2)
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
c        t1=mpi_Wtime()

           call ereal1d_Perm_savetensor(
     &   emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

c        t2=mpi_Wtime()
c        print*,"TensPermElec Real taskid",taskid,"time=",t2-t1

       !print*,"PermElec Real!"



       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      else if((taskid.ge.numtasks_emreal+2)
     &         .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) then

c            call loadbal_ehal1c_half(
c     &        evtmp,devtmp,virevtmp)
c        t3=mpi_Wtime()

            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)
c        t4=mpi_Wtime()
c        print*,"VDW taskid",taskid,"time=",t4-t3

       !print*,"VDW!"

      end if
c            call mpi_barrier(mpi_comm_world,ierr)

         !    t1=mpi_Wtime()
c      emrealt=emrealt+evtmp
c      do i=1,n
c         do j=1,3
c            demrealt(j,i)=demrealt(j,i)+devtmp(j,i)
c         end do
c      end do
c      do i=1,3
c        do j=1,3
c           viremrealt(j,i)=viremrealt(j,i)+virevtmp(j,i)
c        end do
c      end do

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

c      call mpi_barrier(mpi_comm_world,ierr)
c      print*,"longrangepoldir=",longrangepoldir

      if(taskid.eq.master) then
c      call mpi_recv(em,1,mpi_real8,master1,1,mpi_comm_world,stat,ierr)
c      call mpi_recv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
c     &     stat,ierr)
c      call mpi_recv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
c     &     stat,ierr)
      print*,"Before irecv"
      call mpi_irecv(em,1,mpi_real8,master1,1,mpi_comm_world,
     &    reqs4,ierr)
      call mpi_irecv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
     &    reqs5,ierr)
      call mpi_irecv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
     &    reqs6,ierr)
       print*,"After irecv"
c
      end if


      if(taskid.eq.master1) then
        ep3b=0.0d0
        
        if(uzepmedirpolz) then
           term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
           do i = 1, npole
             do j = 1, 3
c            fieldnpole(j,i)=fieldnpolet(j,i)
c            fieldpnpole(j,i)=fieldpnpolet(j,i)
c            demreal_tmp(j,i)=demrealt(j,i)
                fieldnpole(j,i) = fieldnpole(j,i) +fieldnpole_rcp(j,i)
                fieldpnpole(j,i) = fieldpnpole(j,i) +fieldnpole_rcp(j,i)
                fieldnpole(j,i) = fieldnpole(j,i) + term*rpole(j+1,i)
                fieldpnpole(j,i) = fieldpnpole(j,i) + term*rpole(j+1,i)
             end do
         
c              if(longrangepoldir) then
              uind(1,i)=polarity(i) * fieldnpole(1,i)
              uind(2,i)=polarity(i) * fieldnpole(2,i)
              uind(3,i)=polarity(i) * fieldnpole(3,i)
              uinp(1,i)=polarity(i) * fieldpnpole(1,i)
              uinp(2,i)=polarity(i) * fieldpnpole(2,i)
              uinp(3,i)=polarity(i) * fieldpnpole(3,i)
              ep3b=ep3b - 0.5d0*f*( uind(1,i)*fieldpnpole(1,i) +
     &         uind(2,i)*fieldpnpole(2,i) + uind(3,i)*fieldpnpole(3,i))
c              end if
           end do


        else
              term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
             do i=1,npole

c               if(longrangepoldir) then
    
                 uind(1,i)=polarity(i) * fieldnpole(1,i)
                 uind(2,i)=polarity(i) * fieldnpole(2,i)
                 uind(3,i)=polarity(i) * fieldnpole(3,i)

                 uinp(1,i)=polarity(i) * fieldpnpole(1,i)
                 uinp(2,i)=polarity(i) * fieldpnpole(2,i)
                 uinp(3,i)=polarity(i) * fieldpnpole(3,i)
                 ep3b=ep3b - 0.5d0*f*( uind(1,i)*fieldpnpole(1,i) +
     &         uind(2,i)*fieldpnpole(2,i) + uind(3,i)*fieldpnpole(3,i))
c               end if
             end do
        end if 
      end if
             !t2=mpi_Wtime()
       !print*,"TensPermElec Real taskid",taskid,"time=",t2-t1

c      if(.not.longrangepoldir) then
c         call mpi_bcast(fieldnpole,3*n,mpi_real8,master1,
c     &      mpi_comm_world,ierr)
c         call mpi_bcast(fieldpnpole,3*n,mpi_real8,master1,
c     &      mpi_comm_world,ierr) 
cccc         call mpi_bcast(thetai1,4*bsorder*n,mpi_real8,master1,
cccc     &       mpi_comm_world,ierr)
cccc         call mpi_bcast(thetai2,4*bsorder*n,mpi_real8,master1,
cccc     &       mpi_comm_world,ierr)
cccc         call mpi_bcast(thetai3,4*bsorder*n,mpi_real8,master1,
cccc     &       mpi_comm_world,ierr)
cccc         call mpi_bcast(qfac,nfft1*nfft2*nfft3,mpi_real8,master1,
cccc     &       mpi_comm_world,ierr)
cccc         call mpi_bcast(pmetable,n*nchunk,mpi_integer,master1,
cccc     &       mpi_comm_world,ierr)
cccc         call mpi_bcast(igrid,3*n,mpi_integer,master1,
cccc     &       mpi_comm_world,ierr)

c         allocate(dep3bt1(3,npole))
c         allocate(uindt(3,npole))
c         allocate(uinpt(3,npole))
c
c         do i=1,n
c           do j=1,3
c             dep3bt1(j,i)=0.0d0
c             uindt(j,i) = 0.0d0
c             uinpt(j,i) = 0.0d0
c           end do
c         end do
c
c
c         do i=1,3
c           do j=1,3
c            virep3bt1(j,i)=0.0d0
c           end do
c         end do
c         if(taskid.eq.master1) then    
c            do i=1,n
c              do j=1,3
c                uind(j,i) = 0.0d0
c                uinp(j,i) = 0.0d0
c              end do
c            end do
c         end if
c
c         if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then
c           call empole1c_3bPolar_totfield_dEtensor1bmut(
c     &          dep3bt1,virep3bt1,uindt,uinpt)
c         end if
c
c        call mpi_reduce(dep3bt1,dep3b1,npole*3,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c        call mpi_reduce(virep3bt1,virep3b1,3*3,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c
c        call mpi_reduce(uindt,uind,npole*3,mpi_real8,
c     &       mpi_sum,master1,mpi_comm_world,ierr)
c        call mpi_reduce(uinpt,uinp,npole*3,mpi_real8,
c     &       mpi_sum,master1,mpi_comm_world,ierr)
        
c        if(taskid.eq.master1) then
c          ep3b=0.0d0
c          do i=1,npole
c           ep3b=ep3b - 0.5d0*f*( uind(1,i)*fieldpnpole(1,i) +
c     &     uind(2,i)*fieldpnpole(2,i) + uind(3,i)*fieldpnpole(3,i))
c          end do 
c        end if
c        deallocate(dep3bt1)
c        deallocate(uindt)
c        deallocate(uinpt) 
c      end if
      return
      end

