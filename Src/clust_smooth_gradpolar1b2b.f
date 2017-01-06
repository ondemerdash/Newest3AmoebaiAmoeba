      subroutine clustsmooth2_gradient_polar_load_balanced1b2b(start,
     & last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      use neigh2clust
      use aprx
      use deriv1bmat
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript,k
      real*8, allocatable :: dep1btmat(:,:,:)
      real*8, allocatable :: ep1btmat(:)
      real*8, allocatable :: virep1btmat(:,:,:)

      allocate (dep3bt(3,npole))
 
      allocate(ep1btmat(clustcount))
      allocate(dep1btmat(3,3*maxsizeclust,clustcount))
      allocate(virep1btmat(3,3,clustcount))

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

      do i=1,clustcount
         ep1btmat(i)=0.0d0
         do j=1,3*maxsizeclust
            do k=1,3
              dep1btmat(k,j,i)=0.0d0
            end do
         end do
         do j=1,3
            do k=1,3
              virep1btmat(k,j,i)=0.0d0
            end do
         end do  
      end do
      
      ntript=0

c      print*,"numtasks_polar1=",numtasks_polar1
      
      if(taskid.lt.numtasks_polar1) then
        if(ompouterloop3b) then
            call ClustNoListSmoothInnerloop1_1b_Orig2
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat)
        else
          if(mlistclust) then
            call ClustNoListSmoothInnerloop1_1b_Orig2NoOmplist
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat)
c          print*,"Done w ClustNoList1_1b_Orig2NoOmplist",ep3bt,taskid
          else
            call ClustNoListSmoothInnerloop1_1b_Orig2NoOmp
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat)
          end if
        end if
      end if

c      print*,"In clust1b2b After 1b taskid=",taskid,"ep3bt=",ep3bt

       call mpi_allreduce(ep1btmat,ep1bmat,clustcount,mpi_real8,
     &       mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(dep1btmat,dep1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(virep1btmat,virep1bmat,3*3*clustcount,
     &      mpi_real8,mpi_sum,mpi_comm_world,ierr)

      deallocate(ep1btmat)
      deallocate(dep1btmat)
      deallocate(virep1btmat)

      if(do2waterclustlist.and.(taskid.lt.numtasks_polar)) then
        if(ompouterloop3b) then
         call ClustLoadBalNew2SmoothInnerloop1_2b_wskinOrigPrecalc1b(
     &       ep3bt,virep3bt,dep3bt,ntript)
        else
          if(mlistclust) then
           call ClustLoadBalNew2SmoothInnerloop1_2b_Precalc1bNoOmplist(
     &       ep3bt,virep3bt,dep3bt,ntript)
          else
           call ClustLoadBalNew2SmoothInnerloop1_2b_Precalc1bNoOmp(
     &       ep3bt,virep3bt,dep3bt,ntript)
          end if
        end if
      end if

c       print*,"In clust1b2b After 2b taskid=",taskid,"ep3bt=",ep3bt

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

      deallocate(dep3bt)

      return
      end


