      subroutine bcast_pme
      use sizes
      use mpidat
      use chunks
      use pme
      use boxes
      implicit none
      include 'mpif.h'
      integer ierr

      call mpi_bcast(nchunk,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nchk1,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nchk2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nchk3,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(ngrd1,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(ngrd2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(ngrd3,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nlpts,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nrpts,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(grdoff,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(lvec,3*3,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(recip,3*3,mpi_real8,master,
     &   mpi_comm_world,ierr)   
      call mpi_bcast(nfft1,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nfft2,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(nfft3,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(bsorder,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(bsmod1,maxfft,mpi_real8,master,
     &   mpi_comm_world,ierr)   
      call mpi_bcast(bsmod2,maxfft,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(bsmod3,maxfft,mpi_real8,master,
     &   mpi_comm_world,ierr)
      return
      end

      subroutine bsend_pme
      use sizes
      use atoms
      use mpidat
      use pme
      use chunks
      implicit none
      include 'mpif.h'
      integer ierr

      call mpi_bsend(thetai1,4*bsorder*n,mpi_real8,master2,1,
     &     mpi_comm_world,ierr)
      call mpi_bsend(thetai2,4*bsorder*n,mpi_real8,master2,2,
     &     mpi_comm_world,ierr)
      call mpi_bsend(thetai3,4*bsorder*n,mpi_real8,master2,3,
     &     mpi_comm_world,ierr)
      call mpi_bsend(qgrid,2*nfft1*nfft2*nfft3,mpi_real8,master2,4,
     &     mpi_comm_world,ierr)
      call mpi_bsend(qfac,nfft1*nfft2*nfft3,mpi_real8,master2,5,
     &     mpi_comm_world,ierr)
      call mpi_bsend(pmetable,n*nchunk,mpi_real8,master2,6,
     &     mpi_comm_world,ierr)

      return
      end

      subroutine recv_pme
      use sizes
      use atoms
      use mpidat
      use pme
      use chunks
      implicit none
      include 'mpif.h'
      integer ierr
      integer stat(MPI_STATUS_SIZE)

      call mpi_recv(thetai1,4*bsorder*n,mpi_real8,master,1,
     &     mpi_comm_world,stat,ierr)
      call mpi_recv(thetai2,4*bsorder*n,mpi_real8,master,2,
     &     mpi_comm_world,stat,ierr)
      call mpi_recv(thetai3,4*bsorder*n,mpi_real8,master,3,
     &     mpi_comm_world,stat,ierr)
      call mpi_recv(qgrid,2*nfft1*nfft2*nfft3,mpi_real8,master,4,
     &     mpi_comm_world,stat,ierr)
      call mpi_recv(qfac,nfft1*nfft2*nfft3,mpi_real8,master,5,
     &     mpi_comm_world,stat,ierr)
      call mpi_recv(pmetable,n*nchunk,mpi_real8,master,6,
     &     mpi_comm_world,stat,ierr)
      return
      end

      subroutine bcast_ewald
      use mpidat
      use sizes
      use atoms
      use couple
      use neigh3b
      use molcul
      use usage
      use boxes
      use group
      use mpole
      use polar
      use polgrp
      use bound
      use ewald
      use limits
      implicit none
      include 'mpif.h'
      integer ierr

      call mpi_bcast(ewaldcut,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(aewald,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(mpolecut,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(use_ewald,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
      return
      end

      subroutine bcast_ewaldreg
      use mpidat
      use sizes
      use atoms
      use couple
      use neigh3b
      use molcul
      use usage
      use boxes
      use group
      use mpole
      use polar
      use polgrp
      use bound
      use ewald
      use limits
      use parewreg
      use ewaldneigh
      implicit none
      include 'mpif.h'
      integer ierr

      call mpi_bcast(ewaldcut,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(aewald,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(jmax,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(kmax,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(lmax,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(remainder_bcast_alloc,1,mpi_logical,master,
     &   mpi_comm_world,ierr)

      return
      end


      subroutine bcast_ewald_mlistvar
      use mpidat
      use sizes
      use atoms
      use couple
      use neigh3b
      use molcul
      use usage
      use boxes
      use group
      use mpole
      use polar
      use polgrp
      use bound
      use ewald
      use limits
      use neigh
      implicit none
      include 'mpif.h'
      integer ierr

      call mpi_bcast(ewaldcut,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(aewald,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(use_mlist,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(domlst,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(use_ewald,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(lbuffer,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(lbuf2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(mbuf2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(mbufx,1,mpi_real8,master,
     &   mpi_comm_world,ierr)

      return
      end

      subroutine allocPermElec
      use sizes
      use atoms
      use deriv
      use energi
      use virial
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      emreal=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(demreal)) allocate (demreal(3,n))
      end if
        do i=1,3
          do j=1,3
            viremreal(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            demreal(j,i) = 0.0d0
         end do
      end do

      return
      end

      subroutine allocPermElec1
      use sizes
      use atoms
      use deriv
      use energi
      use virial
      implicit none
      logical first3
      save first3
      data first3  / .true. /
      integer i,j

      emreal_tmp=0.0d0

      if (first3) then
         first3 = .false.
         if (.not. allocated(demreal_tmp)) allocate (demreal_tmp(3,n))
      end if
        do i=1,3
          do j=1,3
            viremreal_tmp(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            demreal_tmp(j,i) = 0.0d0
         end do
      end do

      return
      end

      subroutine allocVdW1
      use sizes
      use atoms
      use deriv
      use energi
      use virial
      implicit none
      logical first3
      save first3
      data first3  / .true. /
      integer i,j

      ev=0.0d0

      if (first3) then
         first3 = .false.
         if (.not. allocated(dev)) allocate (dev(3,n))
      end if
        do i=1,3
          do j=1,3
            virev(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dev(j,i) = 0.0d0
         end do
      end do

      return
      end


      subroutine allocPermElec1NoZero
      use sizes
      use atoms
      use deriv
      use energi
      use virial
      implicit none
      logical first3
      save first3
      data first3  / .true. /
      integer i,j

c      emreal_tmp=0.0d0

      if (first3) then
         first3 = .false.
         if (.not. allocated(demreal_tmp)) allocate (demreal_tmp(3,n))
      end if
c        do i=1,3
c          do j=1,3
c            viremreal_tmp(j,i)=0.0d0
c         end do
c        end do
c      do i = 1, n
c         do j=1,3
c            demreal_tmp(j,i) = 0.0d0
c         end do
c      end do

      return
      end

      subroutine allocPermElec1_mlist
      use sizes
      use atoms
      use deriv
      use energi
      use virial
      implicit none
      logical first3
      save first3
      data first3  / .true. /
      integer i,j

      emreal_tmp=0.0d0

      if (first3) then
         first3 = .false.
         if (.not. allocated(demreal_tmp)) allocate (demreal_tmp(3,n))
      end if
        do i=1,3
          do j=1,3
            viremreal_tmp(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            demreal_tmp(j,i) = 0.0d0
         end do
      end do
      call mlist
      return
      end


      subroutine gradient_covalent_vdw
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      implicit none
      include 'mpif.h'
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
c      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time
      real*8 t1,t2,t3,t4,t5,t6

      save first
      data first  / .true. /

      !call system_clock(t1,clock_rate)

      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0
      ev = 0.0d0

      if (first) then
         first = .false.
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
      end if


      do i = 1, n
         do j=1,3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            vir(j,i) =0.0d0
         end do
      end do

      !call ebond1
      !call eangle1
      !call estrbnd1
      !call eurey1
      !call eopbend1
      !call etors1
      !call epitors1
      !call etortor1
c
c     call the local geometry energy and gradient routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_angtor)  call eangtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy and gradient routines
c
      !call system_clock(t1,clock_rate)
        t1=mpi_Wtime()

      call ehal1
      print*,"After vdw from gradient_covalent_vdw vdw=",ev
      !call system_clock(t2,clock_rate)
      t2=mpi_Wtime()
       !tot_time =  (t2-t1)/real(clock_rate)
       !tot_time = t2 - t1
      print*, "MPI WTime Timing  Just vdw", t2-t1

      return
      end

      subroutine gradient_covalent_respa
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0

      if (first) then
         first = .false.
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
      end if


      do i = 1, n
         do j=1,3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            vir(j,i) =0.0d0
         end do
      end do
      call ebond1
      call eangle1
      call estrbnd1
      call eurey1
      call eopbend1
      call etors1
      call epitors1
      call etortor1
c      print*,"After vdw from gradient_covalent_vdw vdw=",ev
      return
      end


      subroutine gradient_emreal
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,ierr
c      real*8 energy,cutoff
c      real*8 derivs(3,*)
c      logical first
c      save first
c      data first  / .true. /
c      real*8, allocatable :: demreal_tmp(:,:)
c      real *8 emreal_tmp,viremreal_tmp(3,3)
c
c      emreal_tmp=0.0d0
c
c      allocate (demreal_tmp(3,n))
c
c      do i = 1, n
c         do j=1,3
c            demreal_tmp(j,i)=0.0d0
c         end do
c      end do
c
c      do i = 1, 3
c         do j=1,3
c            viremreal_tmp(j,i)=0.0d0
c         end do
c      end do
c        call ereal1d_3b_Perm2(emreal_tmp,demreal_tmp,viremreal_tmp)
c        print*,"Before ereal1d_3b_Perm2 in gradient_emreal",emreal_tmp
        call ereal1d_3b_Perm2

c        call mpi_reduce(emreal_tmp,emreal,1,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c        call mpi_reduce(demreal_tmp,demreal,n*3,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c        call mpi_reduce(viremreal_tmp,viremreal,3*3,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c      print*,"From gradient_emreal before sends emreal_tmp=",emreal_tmp

        call mpi_bsend(emreal_tmp,1,mpi_real8,master,4,
     &       mpi_comm_world,ierr)
        call mpi_bsend(demreal_tmp,3*n,mpi_real8,master,5,
     &       mpi_comm_world,ierr)
        call mpi_bsend(viremreal_tmp,3*3,mpi_real8,master,6,
     &       mpi_comm_world,ierr)
c      print*,"From gradient_emreal emreal_tmp=",emreal_tmp
c      print*,"From gradient_emreal after sends emreal_tmp=",emreal_tmp
c      deallocate(demreal_tmp)
      return
      end

      subroutine gradient_emrecip
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
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

      print*,"B4 gradient_emrecip em=",em
      
        call emrecip1_3b_Perm2 
c      call emrecip1_3b_Perm2_noqpole
c      call emrecip1_3b_Perm2_noqpole_nodpole
      print*,"From gradient_emrecip em=",em

      return
      end


      subroutine gradient_emrecip_send
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
c      use deriv3b
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,ierr
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

        call emrecip1_3b_Perm2
      print*,"From gradient_emrecip em=",em

        call mpi_send(em,1,mpi_real8,master,1,
     &       mpi_comm_world,ierr)

        call mpi_send(dem,3*n,mpi_real8,master,2,
     &      mpi_comm_world,ierr)

        call mpi_send(viremrecip,3*3,mpi_real8,master,3,
     &    mpi_comm_world,ierr)

      print*,"AFTER SENDS gradient_emrecip em=",em

      return
      end




c      subroutine gradient_erecip
c      use sizes
c      use atoms
c      use bound
c      use deriv
c      use energi
c      use inter
c      use iounit
c      use limits
c      use potent
c      use rigid
c      use vdwpot
c      use virial
c      use deriv3b
c      implicit none
c      integer i,j
c      real*8 energy,cutoff
c      logical first
c      save first
c      data first  / .true. /
c
c      em = 0.0d0
c
c      if (first) then
c       first = .false.
c       if (.not. allocated(dem))  allocate (dem(3,n))
c      end if
c
c      do i = 1, n
c         do j=1,3
c            dem(j,i) =0.0d0
c         end do
c      end do
c      do i = 1, 3
c         do j=1,3
c            viremrecip(j,i) =0.0d0
c         end do
c      end do
c
c        call erecip1_3b_Perm
c      call parallel_erecip1_3b_Perm
c      call parallel2_erecip1_3b_Perm
c      call parallel3_erecip1_3b_Perm
c      print*,"From gradient_emrecip em=",em

c      return
c      end



      subroutine allocPermElec_forbcast
      use sizes
      use atoms
      use deriv
      use energi
      use virial
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      emreal_tmp=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(demreal_tmp)) allocate (demreal_tmp(3,n))
      end if
        do i=1,3
          do j=1,3
            viremreal_tmp(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            demreal_tmp(j,i) = 0.0d0
         end do
      end do

      return
      end


      subroutine gradient_emreal_bcast
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,ierr
c      real*8 energy,cutoff
c      real*8 derivs(3,*)
c      logical first
c      save first
c      data first  / .true. /
c      real*8, allocatable :: demreal_tmp(:,:)
c      real *8 emreal_tmp,viremreal_tmp(3,3)
c

      emreal_tmp=0.0d0
c
c      allocate (demreal_tmp(3,n))
c
      do i = 1, n
         do j=1,3
            demreal_tmp(j,i)=0.0d0
         end do
      end do
c
      do i = 1, 3
         do j=1,3
            viremreal_tmp(j,i)=0.0d0
         end do
      end do
c        call ereal1d_3b_Perm2(emreal_tmp,demreal_tmp,viremreal_tmp)
c        print*,"Before ereal1d_3b_Perm2 in gradient_emreal",emreal_tmp
        call ereal1d_3b_Perm2

c        call mpi_reduce(emreal_tmp,emreal,1,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c        call mpi_reduce(demreal_tmp,demreal,n*3,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c        call mpi_reduce(viremreal_tmp,viremreal,3*3,mpi_real8,
c     &       mpi_sum,master,mpi_comm_world,ierr)
c      print*,"From gradient_emreal before sends emreal_tmp=",emreal_tmp

c        call mpi_bcast(emreal_tmp,1,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c        call mpi_bcast(demreal_tmp,3*n,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c        call mpi_bcast(viremreal_tmp,3*3,mpi_real8,master1,
c     &       mpi_comm_world,ierr)
c      print*,"After ereal1d_3b_Perm2  emreal_tmp=",emreal_tmp
c      print*,"From gradient_emreal after sends emreal_tmp=",emreal_tmp
c      deallocate(demreal_tmp)
      return
      end

      subroutine gradient_vdw
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

      ev = 0.0d0

      if (first) then
         first = .false.
         if (.not. allocated(dev))  allocate (dev(3,n))
      end if


      do i = 1, n
         do j=1,3
            dev(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            vir(j,i) =0.0d0
         end do
      end do

c
c     call the van der Waals energy component routines
c
c      print*,"Before vdw from gradient_covalent_vdw=",ev
      call ehal1
c      print*,"After vdw from gradient_covalent_vdw vdw=",ev
      return
      end


      subroutine gradient_covalent_vdwTest
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      eat = 0.0d0
      ett = 0.0d0
      ev = 0.0d0

      if (first) then
         first = .false.
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deaa))  allocate (deaa(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(deopd))  allocate (deopd(3,n))
         if (.not. allocated(deid))  allocate (deid(3,n))
         if (.not. allocated(deit))  allocate (deit(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(debt))  allocate (debt(3,n))
         if (.not. allocated(deat))  allocate (deat(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
      end if


      do i = 1, n
         do j=1,3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            deat(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            vir(j,i) =0.0d0
         end do
      end do

c
c     call the van der Waals energy component routines
c
c      print*,"Before vdw from gradient_covalent_vdw=",ev
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_angtor)  call eangtor1
      if (use_tortor)  call etortor1
c      call ehal1
      call ehal1c_noswitch
c      print*,"After vdw from gradient_covalent_vdw vdw=",ev
      return
      end


      subroutine gradient_emrecip_send2
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
c      use deriv3b
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,ierr
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      integer stat(MPI_STATUS_SIZE)
      logical first
      save first
      data first  / .true. /

      if(taskid.eq.master1) then
        call emrecip1_3b_Perm2
      print*,"From gradient_emrecip em=",em

        call mpi_send(em,1,mpi_real8,master,1,
     &       emreal_comm,ierr)

        call mpi_send(dem,3*n,mpi_real8,master,2,
     &      emreal_comm,ierr)

        call mpi_send(viremrecip,3*3,mpi_real8,master,3,
     &    emreal_comm,ierr)

      print*,"AFTER SENDS gradient_emrecip em=",em
      end if

            if(taskid.eq.master) then
                 call mpi_recv(em,1,mpi_real8,master1,1,
     &          emreal_comm,stat,ierr)

                 call mpi_recv(dem,3*n,mpi_real8,master1,2,
     &          emreal_comm,stat,ierr)

                 call mpi_recv(viremrecip,3*3,mpi_real8,master1,3,
     &          emreal_comm,stat,ierr)
            end if

      return
      end


