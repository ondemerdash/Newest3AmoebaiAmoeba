      subroutine gradient_emrecip_totfield_mpi
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
c      use rigid
c      use vdwpot
      use virial
      use totfield
c      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

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

        call emrecip1_3b_Perm2_totfield_mpi
c      print*,"From gradient_emrecip em=",em

      return
      end

      subroutine splitloop(lstart,lend,loop)
      use mpidat
      implicit none
      integer, intent(out) :: lstart
      integer, intent(out) :: lend
      integer, intent(in) :: loop
      integer :: avg

c      avg = loop/nprocs
c      lstart = avg*rank + 1
c      lend = avg*(rank+1) 
c      if(rank .eq. nprocs-1) lend = loop

      avg = loop/numtasks
      lstart=avg*taskid+1
      lend=avg*(taskid+1)
      if(taskid .eq. numtasks-1) lend=loop
      end 


      subroutine gradient_emrecip_totfield_mpi_gridmpole
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
c      use rigid
c      use vdwpot
      use virial
      use totfield
c      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

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

        call emrecip1_3b_Perm2_totfield_mpi_gridmpole
c      print*,"From gradient_emrecip em=",em

      return
      end


      subroutine gradient_emrecip_totfield_mpi_fphimpole
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
c      use rigid
c      use vdwpot
      use virial
      use totfield
c      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

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

        call emrecip1_3b_Perm2_totfield_mpi_fphimpole
c      print*,"From gradient_emrecip em=",em

      return
      end


