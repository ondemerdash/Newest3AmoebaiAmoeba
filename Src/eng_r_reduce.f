
      subroutine eng_emreal_reduce_load_balanced
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      implicit none
      include 'mpif.h'
      real*8 emrealt,viremrealt(3,3)
      real*8, allocatable :: demrealt(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6

      emrealt=0.0d0

      if(taskid.lt.numtasks_emreal2) then 
      call LoadBalInnerloop_ereal0d_3b_Perm_sentlist2(
     & emrealt)
 
      end if 

       call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)

      return
      end
