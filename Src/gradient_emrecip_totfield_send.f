      subroutine gradient_emrecip_totfield_send
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use virial
      use totfield
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,ierr
      logical first
      integer stat(MPI_STATUS_SIZE)
      save first
      data first  / .true. /

      if(taskid.eq.master1) then     
      em = 0.0d0

      if (first) then
       first = .false.
       if (.not. allocated(dem))  allocate (dem(3,n))
      end if

      do i = 1, n
         do j=1,3
            dem(j,i) =0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            viremrecip(j,i) =0.0d0
         end do
      end do

        call emrecip1_3b_Perm2_totfield 
c      print*,"From gradient_emrecip em=",em
      !print*,"Before Emrecip Sends"
      call mpi_send(em,1,mpi_real8,master,1,mpi_comm_world,ierr)
      call mpi_send(dem,3*n,mpi_real8,master,2,mpi_comm_world,ierr)
      call mpi_send(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
     & ierr)
      end if

      if(taskid.eq.master) then
      call mpi_recv(em,1,mpi_real8,master1,1,mpi_comm_world,stat,ierr)
      !     print*,"PME em received",em
      call mpi_recv(dem,3*n,mpi_real8,master1,2,mpi_comm_world,
     &     stat,ierr)
      call mpi_recv(viremrecip,3*3,mpi_real8,master1,3,mpi_comm_world,
     &     stat,ierr)

      end if       
      ! print*,"After Emrecip Sends" 
      return
      end

