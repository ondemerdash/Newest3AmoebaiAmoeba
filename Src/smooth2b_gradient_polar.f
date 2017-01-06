      subroutine smooth2b_gradient_polar_load_balanced
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
      
c        call LoadBalCOMSmooth2bInnerloop1_2cut_wskin(ep3bt,virep3bt,
c     &  dep3bt)
      if(taskid.lt.numtasks_polar) then
       call LoadBalNew2SmoothInnerloop1_2cut_wskin2b(
     &  ep3bt,virep3bt,dep3bt,ntript)
      end if

            call mpi_reduce(ep3bt,ep3bmut,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
            call mpi_reduce(dep3bt,dep3bmut,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
            call mpi_reduce(virep3bt,virep3bmut,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

      deallocate(dep3bt)
      return
      end

      subroutine smallsmooth2b_gradient_polar_load_balanced
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
c      print*,"numtasks_polar=",numtasks_polar
      if(taskid.lt.numtasks_polar) then
        call LoadBalCOMSmooth2bInnerloop1_2cut_wskin(ep3bt,virep3bt,
     &  dep3bt)
      end if

       call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)

      deallocate(dep3bt)
      return
      end

      subroutine smooth2b_gradient_polar(start,last,moli1)
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
      real*8, allocatable :: dep3bt_tot(:,:)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1,i,j,ierr,ntript
c      integer master,ierr
      allocate (dep3bt(3,npole))
      allocate (dep3bt_tot(3,npole))
c      print*,"Beginning of gradient_polar"
                ep3bt_tot=0.0d0
c                ep3bt_tot2=0.0d0
c                ep3bt_tot3=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt_tot(j,i) = 0.0d0
c                      dep3bt_tot2(j,i) = 0.0d0
c                      dep3bt_tot3(j,i) = 0.0d0
                   end do
                end do
                do i=1,3
                   do j=1,3
                      virep3bt_tot(j,i)=0.0d0
c                      virep3bt_tot2(i,j)=0.0d0
c                      virep3bt_tot3(i,j)=0.0d0
                   end do
                end do

        ! call COMSmooth2bInnerloop1_2cut_wskin(start,
     &   !           last,ep3bt,virep3bt,dep3bt,ntript)

c              print*,"After call to Innerloop1_2cut"
                 do i=1,3
                   do j=1,3
                 virep3bt_tot(j,i)=virep3bt_tot(j,i)+virep3bt(j,i)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot=ep3bt_tot+ep3bt

       if(moli1.ne.0) then
      !call COMSmooth2bInnerloop1_2cut_wskin_single(
     & !             moli1,ep3bt,virep3bt,dep3bt,ntript)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot(j,i)=virep3bt_tot(j,i)+virep3bt(j,i)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot=ep3bt_tot+ep3bt
       end if
c       print*,"In gradient_polar",ep3bt_tot
c       print*,"In gradient_polar Before mpi_reduce ",ep3bt_tot
                  call mpi_reduce(ep3bt_tot,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(dep3bt_tot)
      deallocate(dep3bt)
      return
      end

