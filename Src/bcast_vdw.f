      subroutine bcast_vdw
      use mpidat
      use sizes
      use atoms
      use vdw
      use vdwpot
      use limits
      use potent 
      use mutant
      implicit none
      include 'mpif.h'
      integer ierr
          call mpi_bcast(nvdw,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ivdw,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(jvdw,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ired,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(kred,n,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(radmin,maxclass*maxclass,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(epsilon,maxclass*maxclass,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(radmin4,maxclass*maxclass,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(epsilon4,maxclass*maxclass,mpi_real8,master,
     &   mpi_comm_world,ierr)


      call mpi_bcast(vdwcut,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
c      call mpi_bcast(use_vlist,1,mpi_logical,master,
c     &   mpi_comm_world,ierr)
c      call mpi_bcast(dovlst,1,mpi_logical,master,
c     &   mpi_comm_world,ierr)
      call mpi_bcast(use_vdw,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
c      call mpi_bcast(vbuf2,1,mpi_real8,master,
c     &   mpi_comm_world,ierr)
c      call mpi_bcast(vbufx,1,mpi_real8,master,
c     &   mpi_comm_world,ierr)
      call mpi_bcast(vdwtaper,1,mpi_real8,master,
     &   mpi_comm_world,ierr)

      call mpi_bcast(ghal,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dhal,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(v2scale,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(v3scale,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(v4scale,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(v5scale,1,mpi_real8,master,
     &   mpi_comm_world,ierr)

      call mpi_bcast(mut,n,mpi_logical,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(lambda,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(vlambda,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(elambda,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(scexp,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(scalpha,1,mpi_real8,master,
     &   mpi_comm_world,ierr)


      return
      end
