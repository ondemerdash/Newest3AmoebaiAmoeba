      subroutine clustsmooth2_gradient_polar_load_balanced(start,
     & last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      use neigh2clust
      use aprx
      use deriv1bmat
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
      
      if(do2waterclustlist.and.(taskid.lt.numtasks_polar)) then
         if(approxmode.eq.'2BODYMODE') then
         call ClustLoadBalNew2SmoothInnerloop1_2b_wskinOrig(
     &       ep3bt,virep3bt,dep3bt,ntript)
         else
         call ClustLoadBalNew2SmoothInnerloop1_2cut_wskinOrig(
     &       ep3bt,virep3bt,dep3bt,ntript)
         end if
      else if((.not.do2waterclustlist).and.
     &         (taskid.lt.numtasks_polar1)) then
         if(approxmode.eq.'1BODYMODE') then
          call ClustNoListSmoothInnerloop1_1b_Orig
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr)
         else
          call ClustNoListSmoothInnerloop1_2cut_wskinOrig
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr)
         end if
        ! call ClustNoListVanilla2body(
       !if(taskid.lt.2) then
       !call EwClustNoListSmoothInnerloop1_2cut_wskinOrig(
c     &     ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr)
       !end if
      end if

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"In clust After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end


