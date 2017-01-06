
      subroutine gradient_vdw_reduce_load_balanced
      use mpidat
      use deriv
      use energi
      use virial
      use vdw
      implicit none
      include 'mpif.h'
      real*8 evtmp,virevtmp(3,3)
      real*8, allocatable :: devtmp(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6

      allocate (devtmp(3,nvdw))
      !print*,"nvdw=",nvdw
      evtmp=0.0d0
      do i=1,nvdw
         do j=1,3
            devtmp(j,i)=0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           virevtmp(j,i)=0.0d0
        end do
      end do
c             t1=mpi_Wtime()

      if(taskid.lt.numtasks_vdw2) then 
c      call loadbal_ehal1c(
c     & evtmp,devtmp,virevtmp)
            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)

      end if 
c             t2=mpi_Wtime()
       !print*,"evtmp=",evtmp 
c        print*,"Vdw time taskid",t2-t1,taskid
       call mpi_reduce(evtmp,ev,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(devtmp,dev,nvdw*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(virevtmp,virev,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)

         !    t2=mpi_Wtime()
       !print*,"evtmp=",evtmp 
        !print*,"Vdw time taskid",t2-t1,taskid


      deallocate(devtmp)
      return
      end

      subroutine gradient_vdw_reduce_load_balanced_allreduce
      use mpidat
      use deriv
      use energi
      use virial
      use vdw
      implicit none
      include 'mpif.h'
      real*8 evtmp,virevtmp(3,3)
      real*8, allocatable :: devtmp(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6

      allocate (devtmp(3,nvdw))
      !print*,"nvdw=",nvdw
      evtmp=0.0d0
      do i=1,nvdw
         do j=1,3
            devtmp(j,i)=0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           virevtmp(j,i)=0.0d0
        end do
      end do
c             t1=mpi_Wtime()

      if(taskid.lt.numtasks_vdw2) then
c      call loadbal_ehal1c(
c     & evtmp,devtmp,virevtmp)
            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)

      end if
c             t2=mpi_Wtime()
       !print*,"evtmp=",evtmp 
c        print*,"Vdw time taskid",t2-t1,taskid
       call mpi_allreduce(evtmp,ev,1,mpi_real8,
     &      mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(devtmp,dev,nvdw*3,mpi_real8,
     &      mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(virevtmp,virev,3*3,mpi_real8,
     &      mpi_sum,mpi_comm_world,ierr)

         !    t2=mpi_Wtime()
       !print*,"evtmp=",evtmp 
        !print*,"Vdw time taskid",t2-t1,taskid


      deallocate(devtmp)
      return
      end

