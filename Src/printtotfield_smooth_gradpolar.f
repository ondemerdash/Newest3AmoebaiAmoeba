      subroutine PrintNpoletotfield(
     & start,last,moli1rmndr)
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

         call PrinttotfieldNpoleLoadBalNew2SmoothInnerloop1_w1body(
     &   start,last,ep3bt,virep3bt,dep3bt,moli1rmndr) 


                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

       print*,"PrintNpoletotfield After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

      subroutine Print3Npoletotfield(
     & start,last,moli1rmndr)
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

         call Printtotfield3NpoleLoadBalNew2SmoothInnerloop1_w1body(
     &   start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)


                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

       print*,"PrintNpoletotfield After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end


      subroutine Printtotfield(
     & start,last,moli1rmndr)
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

       call totfieldLoadBalNew2SmoothInnerloop1_w1body_rtrndeptemp2
     &  (start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

       print*,"Printtotfield After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end


