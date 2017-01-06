
      subroutine gradient_vdw_reduce_load_balanced_master2_send
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
      integer stat(MPI_STATUS_SIZE)
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
             t1=mpi_Wtime()

      if(taskid.lt.numtasks_vdw2) then 
      call loadbal_ehal1c(
     & evtmp,devtmp,virevtmp)
      end if 

       call mpi_reduce(evtmp,ev,1,mpi_real8,
     &      mpi_sum,2,mpi_comm_world,ierr)
       call mpi_reduce(devtmp,dev,nvdw*3,mpi_real8,
     &      mpi_sum,2,mpi_comm_world,ierr)
       call mpi_reduce(virevtmp,virev,3*3,mpi_real8,
     &      mpi_sum,2,mpi_comm_world,ierr)

             t2=mpi_Wtime()
       !print*,"evtmp=",evtmp 
       print*,"Vdw time taskid",t2-t1,taskid


      deallocate(devtmp)

      if(taskid.eq.2) then
      call mpi_send(ev,1,mpi_real8,master,1,mpi_comm_world,ierr)
      call mpi_send(dev,3*nvdw,mpi_real8,master,2,mpi_comm_world,ierr)
      call mpi_send(virev,3*3,mpi_real8,master,3,mpi_comm_world,ierr) 
      end if
      if(taskid.eq.master) then
      call mpi_recv(ev,1,mpi_real8,2,1,mpi_comm_world,stat,ierr)
      call mpi_recv(dev,3*nvdw,mpi_real8,2,2,mpi_comm_world,stat,ierr)
      call mpi_recv(virev,3*3,mpi_real8,2,3,mpi_comm_world,stat,ierr)
      end if
      return
      end
