c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns particle mesh Ewald parameters and options
c     for a periodic system
c
c
      subroutine kewaldreg3b
      use sizes
      use limits 
      use math
      use boxes
      use ewald
      use ewreg3bpolz
      implicit none

c
c
c     estimate an optimal value for the Ewald coefficient
c
      call ewaldcof (aewald3b,ewaldcut3b)

      jmax3b=int(ewaldcut3b*aewald3b*xbox*aewald3b/pi)
      kmax3b=int(ewaldcut3b*aewald3b*ybox*aewald3b/pi)
      lmax3b=int(ewaldcut3b*aewald3b*zbox*aewald3b/pi)

      print*,"jmax3bcalc kmax3bcalc lmax3bcalc",jmax3b,kmax3b,lmax3b
c      jmax3b=min(maxvec,jmax3b)
c      kmax3b=min(maxvec,kmax3b)
c      lmax3b=min(maxvec,lmax3b)
      jmax3b=12
      kmax3b=12
      lmax3b=12      
      print*,"jmax3bused kmax3bused lmax3bused",jmax3b,kmax3b,lmax3b
      return
      end

      subroutine bcast_ewald3b
      use mpidat
      use ewald
      use limits
      use ewreg3bpolz
      implicit none
      include 'mpif.h'
      integer ierr

      call mpi_bcast(ewaldcut3b,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(aewald3b,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(jmax3b,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(kmax3b,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(lmax3b,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      return
      end

