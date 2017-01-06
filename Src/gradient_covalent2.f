      subroutine gradient_covalent2
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use iounit
      use limits
      use potent
      use virial
      implicit none
      integer i,j
      real*8 energy,cutoff
c      real*8 derivs(3,*)
      logical first
      integer t1,t2,t3,t4
      integer clock_rate
      real tot_time

      save first
      data first  / .true. /

c      call system_clock(t1,clock_rate)

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
            virev(j,i)=0.0d0
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
c      print*,"use_angang=",use_angang
c      print*,"use_opbend=",use_opbend
c      print*,"use_opdist=",use_opdist
c      print*,"use_improp=",use_improp
c      print*,"use_imptor=",use_imptor
c      print*,"use_strtor=",use_strtor
      if (use_bond)  call ebond1
c      print*,"After ebond1, after all the allocs",eb
      if (use_angle)  call eangle1
c      print*,"After eangle1, after all the allocs",ea
      if (use_strbnd)  call estrbnd1
c      print*,"After estrbnd1, after all the allocs",eba
      if (use_urey)  call eurey1
c      print*,"After eurey1, after all the allocs",eub
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

c      print*,"Form gradient_covalent eub=",eub
c      call system_clock(t2,clock_rate)
c      tot_time =  (t2-t1)/real(clock_rate)
c      print*, "Timing  Just covalent", tot_time

      return
      end

