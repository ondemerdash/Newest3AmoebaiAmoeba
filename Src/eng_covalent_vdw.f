      subroutine eng_covalent_vdw
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

      if (use_bond)  call ebond
      if (use_angle)  call eangle
      if (use_strbnd)  call estrbnd
      if (use_urey)  call eurey
      if (use_angang)  call eangang
      if (use_opbend)  call eopbend
      if (use_opdist)  call eopdist
      if (use_improp)  call eimprop
      if (use_imptor)  call eimptor
      if (use_tors)  call etors
      if (use_pitors)  call epitors
      if (use_strtor)  call estrtor
      if (use_angtor)  call eangtor
      if (use_tortor)  call etortor
c
c     call the van der Waals energy and gradient routines
c
      !call system_clock(t1,clock_rate)
c        t1=mpi_Wtime()

      call ehal
c      print*,"After vdw from gradient_covalent_vdw vdw=",ev
      !call system_clock(t2,clock_rate)
c      t2=mpi_Wtime()
       !tot_time =  (t2-t1)/real(clock_rate)
       !tot_time = t2 - t1
c      print*, "MPI WTime Timing  Just vdw", t2-t1

      return
      end

