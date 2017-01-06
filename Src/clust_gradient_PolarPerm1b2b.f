      subroutine clust_gradient_PolarPerm1b2b(start,
     & last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      use neigh2clust
      use aprx
      use deriv1bmat
      use energi
      use deriv
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1rmndr,i,j,ierr,ntript,k
      real*8, allocatable :: dep1btmat(:,:,:)
      real*8, allocatable :: ep1btmat(:)
      real*8, allocatable :: virep1btmat(:,:,:)
      real*8 t1,t2,t3,t4,t5,t6,t7,t8
      real*8 em3bt,virem3bt(3,3)
      real*8, allocatable :: dem3bt(:,:)
      real*8, allocatable :: dem1btmat(:,:,:)
      real*8, allocatable :: em1btmat(:)
      real*8, allocatable :: virem1btmat(:,:,:)


      allocate (dep3bt(3,npole))
 
      allocate(ep1btmat(clustcount))
      allocate(dep1btmat(3,3*maxsizeclust,clustcount))
      allocate(virep1btmat(3,3,clustcount))

      allocate (dem3bt(3,npole))

      allocate(em1btmat(clustcount))
      allocate(dem1btmat(3,3*maxsizeclust,clustcount))
c      allocate(virem1btmat(3,3,clustcount))

      do i=1,npole
         dep3bt(1,i)=0.0d0  
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0   
         dem3bt(1,i)=0.0d0
         dem3bt(2,i)=0.0d0
         dem3bt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
      em3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
c            virem3bt(j,i)=0.0d0
         end do 
      end do 

      do i=1,clustcount
         ep1btmat(i)=0.0d0
         em1btmat(i)=0.0d0
         do j=1,3*maxsizeclust
            do k=1,3
              dep1btmat(k,j,i)=0.0d0
              dem1btmat(k,j,i)=0.0d0
            end do
         end do
         do j=1,3
            do k=1,3
              virep1btmat(k,j,i)=0.0d0
c              virem1btmat(k,j,i)=0.0d0
            end do
         end do  
      end do
      
      ntript=0

c      print*,"numtasks_polar1=",numtasks_polar1
      
      if(taskid.lt.numtasks_polar1) then
      !if(taskid.eq.7) then
          if(mlistclust) then
c            t1=mpi_Wtime()
            call ClustEmpole1_1b_NoOmplist
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat,em3bt,dem3bt,em1btmat,dem1btmat
     &  )
c            t2=mpi_Wtime()
c          print*,"Done w ClustEmpole1_1b_NoOmplist",ep3bt,taskid
c           print*,"taskid=",taskid,"1body Computation Cost=",t2-t1
          else
            call ClustEmpole1_1b_NoOmp
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat,em3bt,dem3bt,em1btmat,dem1btmat
     &  )
c          print*,"Done w ClustEmpole1_1b_NoOmp",ep3bt,taskid
          end if
      end if

c      print*,"In clust1b2b After 1b taskid=",taskid,"ep3bt=",ep3bt

c       call mpi_barrier(mpi_comm_world,ierr)

c            t3=mpi_Wtime()
       call mpi_allreduce(ep1btmat,ep1bmat,clustcount,mpi_real8,
     &       mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(dep1btmat,dep1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(virep1btmat,virep1bmat,3*3*clustcount,
     &      mpi_real8,mpi_sum,mpi_comm_world,ierr)
c
       call mpi_allreduce(em1btmat,em1bmat,clustcount,mpi_real8,
     &       mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(dem1btmat,dem1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,mpi_comm_world,ierr)

      deallocate(ep1btmat)
      deallocate(dep1btmat)
      deallocate(virep1btmat)
      deallocate(em1btmat)
      deallocate(dem1btmat)
c      deallocate(virem1btmat)
c            t4=mpi_Wtime()
c            print*,"taskid=",taskid,"Allreduce Cost=",t4-t3


      if(do2waterclustlist.and.(taskid.lt.numtasks_polar)) then
          if(mlistclust) then
c           t5=mpi_Wtime()
           call ClustEmpole1_2b_NoOmplist(
     &       ep3bt,virep3bt,dep3bt,ntript,em3bt,dem3bt)
c           t6=mpi_Wtime()
         !   print*,"taskid=",taskid,"2body Computation Cost=",t6-t5
c        print*,"Done w ClustEmpole1_2b_NoOmplist",ep3bt,taskid
          else
           call ClustEmpole1_2b_NoOmp(
     &       ep3bt,virep3bt,dep3bt,ntript,em3bt,dem3bt)
c        print*,"Done w ClustEmpole1_2b_NoOmp",ep3bt,taskid
          end if
      end if

c       print*,"In clust1b2b After 2b taskid=",taskid,"ep3bt=",ep3bt
c       call mpi_barrier(mpi_comm_world,ierr)

c           t7=mpi_Wtime()
                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

                  call mpi_reduce(em3bt,em,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dem3bt,dem,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
c           t8=mpi_Wtime()
c            print*,"taskid=",taskid,"Reduce Cost=",t8-t7


      deallocate(dep3bt)
      deallocate(dem3bt)
      return
      end


