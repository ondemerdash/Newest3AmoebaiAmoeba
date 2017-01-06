      subroutine ewtotfieldsmooth2_gradient_polar_load_balanced_dEtensor
      use mpole
      use mpidat
      use deriv3b
      use totfield
      use molcul
      use aprx
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do
      ntript=0
      if(taskid.lt.numtasks_polar) then
        if(approxmode.eq.'2BODYMODE') then
       !  call totfieldNpoleloadbalsmoothInner_1c_2bPolar(
     & ! ep3bt,virep3bt,dep3bt)
         call totfieldNpoleloadbalsmoothInner_1c_2bPolar2(
     &  ep3bt,virep3bt,dep3bt)
       !  call totfieldNpoleloadbalsmoothInner_1c_2bPolar2_511(
     &  !ep3bt,virep3bt,dep3bt)
        ! call totfieldNpoleloadbalsmoothInner_1c_2bPolar2_513(
     &  ! ep3bt,virep3bt,dep3bt)
        else
         call totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi(
     &   ep3bt,virep3bt,dep3bt)
        end if
      end if

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

      ! print*,"totfield_smooth_gradpolar After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

c      subroutine ewtotfieldsmallsmooth2_gradient_polar_1body_dEtensor(
c     & start,last,moli1rmndr)
      subroutine ewtotfieldsmallsmooth2_gradient_polar_1body_dEtensor
      use mpole
      use mpidat
      use deriv3b
      use totfield
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
c      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do


c      print*,"Before empole1c_3b_Polar_totfield_dEtensor1b_nouselist"
      call empole1c_3b_Polar_totfield_dEtensor1b_nouselist(
     &  dep3bt,virep3bt)
c                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c      print*,"After empole1c_3b_Polar_totfield_dEtensor1b_nouselist"
      deallocate(dep3bt)
      return
      end

      subroutine ewtotfieldsmooth2_gradient_polar_1body_dEtensor_mpi
      use mpole
      use mpidat
      use deriv3b
      use totfield
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
      real*8 t1,t2,t3,t4,t5,t6

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
c      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do
      
       !t1=mpi_Wtime()
      if(taskid.lt.numtasks_emreal2) then    
      !call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpi(
     & ! dep3bt,virep3bt)
      call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
     &  dep3bt,virep3bt)
      end if
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

        !t2=mpi_Wtime()
       !print*,"Tensor Empole1c_3b_Polar taskid",taskid,"time=",t2-t1
      deallocate(dep3bt)
      return
      end

      subroutine ewtotfieldsmooth2_gradient_polar_1body_dEtensor_mpioff
      use mpole
      use mpidat
      use deriv3b
      use totfield
      use aprx
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
      real*8 t1,t2,t3,t4,t5,t6

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
c      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do

c      if(longrangepoldir) then
         if(uzepmedirpolz) then
          if(taskid.eq.0) then
            call emrecip1_3b_Polar_usemlist_direct_virialonly
          else if(taskid.eq.1)then
           call emrecip1_3b_Polar_usemlist_direct
          else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then
         call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
     &     dep3bt,virep3bt)
          end if
         else
      
           if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then
         call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
     &     dep3bt,virep3bt)
           end if
         end if 
c      else
c         if(uzepmedirpolz) then
c          if(taskid.eq.0) then
c            call emrecip1_3b_Polar_usemlist_mutual_virialonly
c          else if(taskid.eq.1)then
c           call emrecip1_3b_Polar_usemlist_mutual
c          else if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then
c         call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
c     &     dep3bt,virep3bt)
c          end if
c         else
c           if((taskid.gt.1).and.(taskid.lt.numtasks_emreal2+2)) then
c         call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
c     &     dep3bt,virep3bt)
c           end if
c         end if
c      end if 

                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

        !t2=mpi_Wtime()
       !print*,"Tensor Empole1c_3b_Polar taskid",taskid,"time=",t2-t1
      deallocate(dep3bt)
      return
      end


      subroutine ewtotfieldsmooth2_gradient_polar_1body_dEtensor_mpiall
      use mpole
      use mpidat
      use deriv3b
      use totfield
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
      real*8 t1,t2,t3,t4,t5,t6

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
c      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do

       !t1=mpi_Wtime()
      if(taskid.lt.numtasks_emreal2) then
      !call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpi(
     & ! dep3bt,virep3bt)
      call empole1c_3b_Polar_totfield_dEtensor1b_nouselist_mpivir(
     &  dep3bt,virep3bt)
      end if
                  call mpi_allreduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,mpi_comm_world,ierr)
                  call mpi_allreduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,mpi_comm_world,ierr)

        !t2=mpi_Wtime()
       !print*,"Tensor Empole1c_3b_Polar taskid",taskid,"time=",t2-t1
      deallocate(dep3bt)
      return
      end


