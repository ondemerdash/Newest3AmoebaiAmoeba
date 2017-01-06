      subroutine gradient_emreal_reduce(startem,lastem,atomind)
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

      allocate (demrealt(3,npole))

      emrealt=0.0d0
      do i=1,npole
         do j=1,3
            demrealt(j,i)=0.0d0
         end do 
      end do  
      do i=1,3
        do j=1,3
           viremrealt(j,i)=0.0d0
        end do  
      end do 

c       call Innerloop_ereal1d_3b_Perm_fulllist2(atomind,startem,
c     &      lastem,emrealt,viremrealt,demrealt)

       call Innerloop_ereal1d_3b_Perm_sentlist2(atomind,startem,
     &      lastem,emrealt,viremrealt,demrealt)

       call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(demrealt)
      return
      end

      subroutine gradient_emreal_reduce_load_balanced
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

      allocate (demrealt(3,npole))

      emrealt=0.0d0
      do i=1,npole
         do j=1,3
            demrealt(j,i)=0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           viremrealt(j,i)=0.0d0
        end do
      end do
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"

      if(taskid.lt.numtasks_emreal2) then 
c             t3=mpi_Wtime()
      call LoadBalInnerloop_ereal1d_3b_Perm_sentlist2(
     & emrealt,viremrealt,demrealt)
c             t4=mpi_Wtime()
c       print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
c       print*,"emrealt=",emrealt,"taskid=",taskid
      end if 
c            call mpi_barrier(mpi_comm_world,ierr)

c             t1=mpi_Wtime()
       call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)

c             t2=mpi_Wtime()
c       print*,"PermElec Real time MPIReduceOnly,taskid",t2-t1,taskid

c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(demrealt)
      return
      end

c      subroutine gradient_emreal_reduce_load_balanced2
c      use mpole
c      use mpidat
c      use deriv
c      use energi
c      use virial
c      implicit none
c      include 'mpif.h'
c      real*8 emrealt,viremrealt(3,3)
c      real*8, allocatable :: demrealt(:,:)
c      integer startem,lastem,atomind,i,j,ierr
c
c      allocate (demrealt(3,npole))
c
c      emrealt=0.0d0
c      do i=1,npole
c         do j=1,3
c            demrealt(j,i)=0.0d0
c         end do
c      end do
c      do i=1,3
c        do j=1,3
c           viremrealt(j,i)=0.0d0
c        end do
c      end do
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"
c      call LoadBalInnerloop_ereal1d_3b_Perm_fulllist2(
c     & emrealt,viremrealt,demrealt)
c
c       call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
c     &      mpi_sum,master,mpi_comm_world,ierr)
c       call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
c     &      mpi_sum,master,mpi_comm_world,ierr)
c       call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
c     &      mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

c      deallocate(demrealt)
c      return
c      end
 
