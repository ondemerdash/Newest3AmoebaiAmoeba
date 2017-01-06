      subroutine ewtotfieldsmallsmooth2_gradient_polar_serial
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
      !if(taskid.lt.numtasks_polar) then
      start_polar(taskid)=1
      last_polar(taskid)=nmol 
        if(approxmode.eq.'2BODYMODE') then
       !  call totfieldNpoleloadbalsmoothInner_1c_2bPolar(
     & ! ep3bt,virep3bt,dep3bt)
         call totfieldNpoleloadbalsmoothInner_1c_2bPolar2(
     &  ep3bt,virep3bt,dep3bt)
       !  call totfieldNpoleloadbalsmoothInner_1c_2bPolar2_511(
     &  !ep3bt,virep3bt,dep3bt)
        ! call totfieldNpoleloadbalsmoothInner_1c_2bPolar2_513(
     &  !ep3bt,virep3bt,dep3bt)
        else
         call totfieldNpoleloadbalsmoothInner_1c_3bPolar(
     &  ep3bt,virep3bt,dep3bt)
        end if
      !end if

              !    call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &        !  mpi_sum,master,mpi_comm_world,ierr)
        ep3b=ep3bt

             !     call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &       !   mpi_sum,master,mpi_comm_world,ierr)
          do i=1,npole
              dep3b(1,i)=dep3bt(1,i)
              dep3b(2,i)=dep3bt(2,i)
              dep3b(3,i)=dep3bt(3,i)
          end do
          !        call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &    !      mpi_sum,master,mpi_comm_world,ierr)
          do i=1,3
             do j=1,3
               virep3b(j,i)=virep3bt(j,i)
             end do
          end do
              !    call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &        !  mpi_sum,master,mpi_comm_world,ierr)

      ! print*,"totfield_smooth_gradpolar After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

