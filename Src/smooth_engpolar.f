
      subroutine smallsmooth2_eng_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript

      ep3bt=0.0d0

      if(taskid.lt.numtasks_polar) then
        call LoadBalNew2SmoothInnerloop0_2cut_wskinOrig(
     &  ep3bt)
      end if
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

      return
      end
