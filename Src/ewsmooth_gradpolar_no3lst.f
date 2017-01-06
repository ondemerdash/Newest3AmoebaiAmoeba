
      subroutine ewsmallsmooth2_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
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
c        call LoadBalSmoothInnerloop1_ew(ep3bt,virep3bt,
c     &  dep3bt)
      end if
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

c      print*,"Ewloadbal gradient_polar After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

      subroutine ewNOsmallsmooth2_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
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
c        call LoadBalNOSmoothInnerloop1_ew(ep3bt,virep3bt,
c     &  dep3bt)
      end if
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

c      print*,"Ewloadbal gradient_polar After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

