      subroutine totfieldsmallsmooth2_gradient_polar_load_balanced(
     & start,last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      use totfield
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do
      ntript=0



c       call totfieldLoadBalNew2SmoothInnerloop1_rtrndeptemp2
c     &  (start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
       if(embedtyp.eq.'I') then
        print*,"embedtyp=",embedtyp
        call totfieldLoadBalNew2SmoothInnerloop1_w1body_rtrndeptemp2
     &  (start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
       else if(embedtyp.eq.'II') then
        print*,"embedtyp=",embedtyp
         call  totfieldNpoleLoadBalNew2SmoothInnerloop1_2cut_w1body(
     &   start,last,ep3bt,virep3bt,dep3bt,moli1rmndr) 
       else if(embedtyp.eq.'III') then
        print*,"embedtyp=",embedtyp
          call totfield3NpoleLoadBalNew2SmoothInnerloop1_2cut_w1body(
     &   start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
       end if
c       call All2bAll3bNew2SmoothInnerloop1_empole1c_totfield(
c     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)

c          call totfield3NpoleLoadBalNew2SmoothInnerloop1_2cut_w1body(
c     &   start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
c         call  totfieldNpoleLoadBalNew2SmoothInnerloop1_2cut_w1body(
c     &   start,last,ep3bt,virep3bt,dep3bt,moli1rmndr) 
c      end if


c      end if
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"totfield_smooth_gradpolar After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end


      !subroutine ewtotfieldsmallsmooth2_gradient_polar_load_balanced(
     & !start,last,moli1rmndr)
      subroutine ewtotfieldsmallsmooth2_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
      use totfield
      use molcul
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do
      ntript=0
      !start =1
      !last = nmol
      !moli1rmndr = 0
       !call All2bAll3bNew2SmoothInnerloop1_empole1c_totfield(
     & ! start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
      !if(taskid.eq.master) then
      !  call totfieldNpolesmoothInner_1c_3bPolar(
     &!  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
      !end if
      if(taskid.lt.numtasks_polar) then
        call totfieldNpoleloadbalsmoothInner_1c_3bPolar(
     &  ep3bt,virep3bt,dep3bt)
      end if

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

       !print*,"totfield_smooth_gradpolar After mpi_reduce",taskid,ep3bt


      deallocate(dep3bt)
      return
      end


      subroutine stotfieldsmallsmooth2_gradient_polar_load_balanced(
     & start,last,moli1rmndr) 
      use mpole
      use mpidat
      use deriv3b
      use totfield
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript

      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do
      ntript=0

       if(embedtyp.eq.'I') then

        print*,"embedtyp=",embedtyp      
        if(taskid.eq.0) then
        !call s2btotfieldLoadBalNew2SmoothInnerloop1_rtrndeptemp2
        !call s1btotfieldLoadBalNew2SmoothInnerloop1_rtrndeptemp2
         call stotfieldLoadBalNew2SmoothInnerloop1w1body_rtrndeptemp2 
     &  (1,4,ep3bt,virep3bt,dep3bt,0)
        else if(taskid.eq.1) then
        !call s2btotfieldLoadBalNew2SmoothInnerloop1_rtrndeptemp2
        ! call s1btotfieldLoadBalNew2SmoothInnerloop1_rtrndeptemp2
         call stotfieldLoadBalNew2SmoothInnerloop1w1body_rtrndeptemp2 
     &  (5,8,ep3bt,virep3bt,dep3bt,0)
        end if

       else if(embedtyp.eq.'II') then

        print*,"embedtyp=",embedtyp
        if(taskid.eq.0) then
        !call s2btotfieldNpoleLoadBalNew2SmoothInnerloop1_w1body
        !call s1btotfieldNpoleLoadBalNew2SmoothInnerloop1_w1body
        call stotfieldNpoleLoadBalNew2SmoothInnerloop1_2cut_w1body 
     &  (1,4,ep3bt,virep3bt,dep3bt,0)
        else if(taskid.eq.1) then
        !call s2btotfieldNpoleLoadBalNew2SmoothInnerloop1_w1body
        !call s1btotfieldNpoleLoadBalNew2SmoothInnerloop1_w1body
        call stotfieldNpoleLoadBalNew2SmoothInnerloop1_2cut_w1body
     &  (5,8,ep3bt,virep3bt,dep3bt,0)
        end if

       else if(embedtyp.eq.'III') then

        print*,"embedtyp=",embedtyp
        if(taskid.eq.0) then
         !call stotfield3NpoleLoadBalNew2SmoothInnerloop1_2cut_w1body
        !call s2btotfield3NpoleLoadBalNew2SmoothInnerloop1_w1body 
        !call s1btotfield3NpoleLoadBalNew2SmoothInnerloop1_w1body
        call stotfield3NpoleLoadBalNew2SmoothInnerloop1_2cut_w1body
     &  (1,4,ep3bt,virep3bt,dep3bt,0)
        else if(taskid.eq.1) then
         !call stotfield3NpoleLoadBalNew2SmoothInnerloop1_2cut_w1body
        ! call s2btotfield3NpoleLoadBalNew2SmoothInnerloop1_w1body
        !call s1btotfield3NpoleLoadBalNew2SmoothInnerloop1_w1body
        call stotfield3NpoleLoadBalNew2SmoothInnerloop1_2cut_w1body
     &  (5,8,ep3bt,virep3bt,dep3bt,0)
        end if

       end if
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"In loadbal gradient_polar After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

