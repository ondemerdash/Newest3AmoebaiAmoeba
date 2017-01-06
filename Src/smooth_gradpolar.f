      subroutine smooth_gradient_polar(start,last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
      real*8, allocatable :: dep3bt_tot(:,:)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
c      integer master,ierr
      allocate (dep3bt(3,npole))
      allocate (dep3bt_tot(3,npole))
c      print*,"Beginning of gradient_polar"
                ep3bt_tot=0.0d0
c                ep3bt_tot2=0.0d0
c                ep3bt_tot3=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt_tot(j,i) = 0.0d0
c                      dep3bt_tot2(j,i) = 0.0d0
c                      dep3bt_tot3(j,i) = 0.0d0
                   end do
                end do
                do i=1,3
                   do j=1,3
                      virep3bt_tot(j,i)=0.0d0
c                      virep3bt_tot2(i,j)=0.0d0
c                      virep3bt_tot3(i,j)=0.0d0
                   end do
                end do

        ! call COMSmooth2bInnerloop1_2cut_wskin(start,
     &   !           last,ep3bt,virep3bt,dep3bt,ntript)

c              print*,"After call to Innerloop1_2cut"
                 do i=1,3
                   do j=1,3
                 virep3bt_tot(j,i)=virep3bt_tot(j,i)+virep3bt(j,i)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot=ep3bt_tot+ep3bt

       if(moli1rmndr.ne.0) then
      !call COMSmooth2bInnerloop1_2cut_wskin_single(
     & !             moli1rmndr,ep3bt,virep3bt,dep3bt,ntript)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot(j,i)=virep3bt_tot(j,i)+virep3bt(j,i)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot=ep3bt_tot+ep3bt
       end if

c       print*,"In gradient_polar Before mpi_reduce ",ep3bt_tot
                  call mpi_reduce(ep3bt_tot,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(dep3bt_tot)
      deallocate(dep3bt)
      return
      end

      subroutine print_smooth2_gradient_polar(start,last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
      real*8, allocatable :: dep3bt_tot(:,:)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
c      integer master,ierr
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

       !call PrintAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(start,last,
     & !           ep3bt,virep3bt,dep3bt,ntript,moli1rmndr)

        if(taskid.eq.0) then
        !call sAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(1,4,
     &  !          ep3bt,virep3bt,dep3bt,0)
       !  call sAll2bNew2SmoothInnerloop1_2cut_wskin(1,4,
     &  !          ep3bt,virep3bt,dep3bt,0)

        else if(taskid.eq.1) then
        ! call sAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(5,8,
     &   !         ep3bt,virep3bt,dep3bt,0)
       !  call sAll2bNew2SmoothInnerloop1_2cut_wskin(5,8,
     & !           ep3bt,virep3bt,dep3bt,0)
        end if
        !call All2bAll3bNew2SmoothInnerloop1_2cut_wskin(start,last,
     &  !          ep3bt,virep3bt,dep3bt,moli1rmndr)


c        call All2bNew2SmoothInnerloop1_2cut_wskin(start,last,
c     &            ep3bt,virep3bt,dep3bt,moli1rmndr)


c       print*,"In gradient_polar Before mpi_reduce ",ep3bt_tot
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot


      deallocate(dep3bt)
      return
      end


      subroutine smooth2_gradient_polar(start,last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
      real*8, allocatable :: dep3bt(:,:)
      real*8, allocatable :: uindt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
c      integer master,ierr
      allocate (dep3bt(3,npole))
c      allocate (uindt(3,npole))
c      print*,"Beginning of gradient_polar"
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
c      ntript=0

c         call NewL3ongTaprR123Smooth3b2b2Innerloop1_2cut_wskin(start,
c     &              last,ep3bt,virep3bt,dep3bt,ntript,moli1rmndr)
c          call NoOMPL3ongTaprR123Smooth3b2b2Innerloop1_2cut_wskin(start,
c     &              last,ep3bt,virep3bt,dep3bt,ntript,moli1rmndr)
c         call NewSmoothInnerloop1_2cut_wskin(start,last,ep3bt,
c     &            virep3bt,dep3bt,ntript,moli1rmndr)

c         call New2SmoothInnerloop1_2cut_wskin(start,last,ep3bt,
c     &            virep3bt,dep3bt,ntript,moli1rmndr)

c         call AllatomNew2SmoothInnerloop1_2cut_wskin(start,last,ep3bt,
c     &            virep3bt,dep3bt,ntript,moli1rmndr)

c        call UindAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(start,last,
c     &            ep3bt,virep3bt,dep3bt,ntript,moli1rmndr,uindt)

        call All2bAll3bNew2SmoothInnerloop1_2cut_wskin(start,last,
     &            ep3bt,virep3bt,dep3bt,moli1rmndr)

        !call All2bAll3bNew2SmoothInnerloop1_empole1c(start,last,
     &  !          ep3bt,virep3bt,dep3bt,moli1rmndr)

        !if(taskid.eq.0) then
        !  call sAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(1,4,
     &  !          ep3bt,virep3bt,dep3bt,0)
         !call sAll2bNew2SmoothInnerloop1_2cut_wskin(1,4,
     &   !         ep3bt,virep3bt,dep3bt,0) 
        !else if(taskid.eq.1) then
         !  call sAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(5,8,
     &   !         ep3bt,virep3bt,dep3bt,0)
         ! call sAll2bNew2SmoothInnerloop1_2cut_wskin(5,8,
     &   !         ep3bt,virep3bt,dep3bt,0)
        !end if


        !call All2bNew2SmoothInnerloop1_2cut_wskin(start,last,
     &  !          ep3bt,virep3bt,dep3bt,moli1rmndr)
c        call All2bAll3bAll4bNew2SmoothInnerloop1_2cut_wskin(start,last,
c     &            ep3bt,virep3bt,dep3bt,moli1rmndr)

c       print*,"In gradient_polar Before mpi_reduce ",ep3bt_tot
              !    call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &        !  mpi_sum,master,mpi_comm_world,ierr)
              !    call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &        !  mpi_sum,master,mpi_comm_world,ierr)
              !    call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &        !  mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(uindt,uind3b,npole*3,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
                  
              !    call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &        !  mpi_sum,master,mpi_comm_world,ierr)

                  call mpi_reduce(ep3bt,ep3bmut,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3bmut,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3bmut,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

c      deallocate(uindt)
      deallocate(dep3bt)
      return
      end

      subroutine smooth2_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
      use aprx
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

      if(taskid.lt.numtasks_polar) then
        if(longrangepoldir) then
        call LoadBalNew2SmoothInnerloop1_2cut_wskin(ep3bt,virep3bt,
     &  dep3bt,ntript)
        else
        call LoadBalNew2SmoothInnerloop1_2cut_wskin1bmut(ep3bt,
     &   virep3bt,dep3bt,ntript)
        end if
c        call b2test(ep3bt,virep3bt,
c     &  dep3bt,ntript)
      end if
                  call mpi_reduce(ep3bt,ep3bmut,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3bmut,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3bmut,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"In loadbal gradient_polar After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

      subroutine smallsmooth2_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
      use aprx
      use mpidat3
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript
      real*8 t1,t2,t3,t4,t5,t6

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

      if(taskid.lt.numtasks_polar) then
c        call AnintLoadBalNew2SmoothInnerloop1_2cut_wskin(ep3bt,virep3bt,
c     &  dep3bt,ntript)
        !t1=mpi_Wtime()
        if(approxmode.eq.'2BODYMODE') then
        call LoadBalNew2SmoothInnerloop1_2cut_wskinOrig2b(ep3bt,
     &  virep3bt,dep3bt,ntript)
        else
        call LoadBalNew2SmoothInnerloop1_2cut_wskinOrig(ep3bt,virep3bt,
     &  dep3bt,ntript)
        end if
        !t2=mpi_Wtime()
       !print*,"Computation 2b/3b nnPolz taskid",taskid,"time=",t2-t1
      end if

       ! t3=mpi_Wtime()
c                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

                  call mpi_ireduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,reqs21,ierr)
                  call mpi_ireduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,reqs22,ierr)
                  call mpi_ireduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,reqs23,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

       ! t4=mpi_Wtime()
       !print*,"Ireduce 2b/3b nnPolz taskid",taskid,"time=",t4-t3

c       print*,"In loadbal gradient_polar After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

      subroutine smallsmooth2body_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
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

      if(taskid.lt.numtasks_polar) then
        !call LoadBalSmoothInnerloop2body(ep3bt,virep3bt,
     &  ! dep3bt)
      end if
                  call mpi_reduce(ep3bt,ep2b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep2b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep2b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(ntript,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr)

c       print*,"In loadbal gradient_polar After mpi_reduce ",taskid,ep3bt


      deallocate(dep3bt)
      return
      end

      subroutine smallsmooth3body_gradient_polar_load_balanced
      use mpole
      use mpidat
      use deriv3b
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

      if(taskid.lt.numtasks_polar3) then
        !call LoadBalSmoothInnerloop3body(ep3bt,virep3bt,
     &  !dep3bt,ntript)
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

