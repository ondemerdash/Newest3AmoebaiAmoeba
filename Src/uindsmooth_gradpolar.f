      subroutine uindsmallsmooth2_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
      real*8, allocatable :: uindt(:,:)

      allocate (dep3bt(3,npole))
      allocate (uindt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
         uindt(1,i)=0.0d0
         uindt(2,i)=0.0d0
         uindt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do
      ntript=0

      if(taskid.lt.numtasks_polar) then
        call uindLoadBalNew2SmoothInnerloop1_2cut_wskinOrig(
     &  ep3bt,virep3bt,dep3bt,ntript,uindt) 

      end if
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)
                call mpi_reduce(uindt,uind3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In loadbal gradient_polar After mpi_reduce


      deallocate(dep3bt)
      deallocate(uindt)
      return
      end

      subroutine uindsmooth2_gradient_polar(start,last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
      real*8, allocatable :: dep3bt(:,:)
      real*8, allocatable :: uindt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
c      integer master,ierr
      allocate (dep3bt(3,npole))
      allocate (uindt(3,npole))
c      print*,"Beginning of gradient_polar"
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
         uindt(1,i)=0.0d0
         uindt(2,i)=0.0d0
         uindt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do

c        call All2bAll3bNew2SmoothInnerloop1_empole1c(start,last,
c     &            ep3bt,virep3bt,dep3bt,moli1rmndr)

      call uindInnerloop1_Nocut(start,last,
     &  ep3bt,virep3bt,dep3bt,uindt,moli1rmndr)

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                call mpi_reduce(uindt,uind3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In loadbal gradient_polar After mpi_reduce


      deallocate(dep3bt)
      deallocate(uindt)
      return
      end

