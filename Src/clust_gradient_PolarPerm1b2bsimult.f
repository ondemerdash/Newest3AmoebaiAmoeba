      subroutine clust_gradient_PolarPerm1b2bsimult(start,
     & last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      use neigh2clust
      use aprx
      use deriv1bmat
      use energi
      use deriv
      use deriv2bmat
      use molcul
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
      real*8, allocatable :: dep2btmat(:,:,:,:)
      real*8, allocatable :: ep2btmat(:,:)
      real*8, allocatable :: virep2btmat(:,:,:,:)
      real*8, allocatable :: dem2btmat(:,:,:,:)
      real*8, allocatable :: em2btmat(:,:)
      integer j1,clust1,npole3b,np1,l1,clust2,np2,l2,moli1
      integer, allocatable :: pnum(:) 

c      allocate (dep3bt(3,npole))
 
      allocate(ep1btmat(clustcount))
      allocate(dep1btmat(3,3*maxsizeclust,clustcount))
      allocate(virep1btmat(3,3,clustcount))

c      allocate (dem3bt(3,npole))

      allocate(em1btmat(clustcount))
      allocate(dem1btmat(3,3*maxsizeclust,clustcount))


      allocate(ep2btmat(clustcount,clustcount))
      allocate(dep2btmat(3,2*3*maxsizeclust,clustcount,clustcount))
      allocate(virep2btmat(3,3,clustcount,clustcount))
      allocate(em2btmat(clustcount,clustcount))
      allocate(dem2btmat(3,2*3*maxsizeclust,clustcount,clustcount))

c      do i=1,npole
c         dep3bt(1,i)=0.0d0  
c         dep3bt(2,i)=0.0d0
c         dep3bt(3,i)=0.0d0   
c         dem3bt(1,i)=0.0d0
c         dem3bt(2,i)=0.0d0
c         dem3bt(3,i)=0.0d0
c      end do
c      ep3bt=0.0d0
c      em3bt=0.0d0
c      do i=1,3
c         do j=1,3
c            virep3bt(j,i)=0.0d0
c         end do 
c      end do 

      do i=1,clustcount
         do j1=1,clustcount
            ep2btmat(j1,i)=0.0d0
            em2btmat(j1,i)=0.0d0
            do j=1,2*3*maxsizeclust
               do k=1,3
                dep2btmat(k,j,j1,i)=0.0d0
                dem2btmat(k,j,j1,i)=0.0d0
               end do
            end do
            do j=1,3
              do k=1,3
               virep2btmat(k,j,j1,i)=0.0d0
              end do
            end do
         end do
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
            end do
         end do  
      end do

      
      ntript=0

c      print*,"numtasks_polar1=",numtasks_polar1
      
      if(taskid.lt.numtasks_polar1) then
            call ClustEmpole1_1b_NoOmplistsimult
     &  (start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat,em1btmat,dem1btmat
     &  )
      else if(do2waterclustlist.and.(taskid.ge.numtasks_polar1)
     &     .and.(taskid.lt.(numtasks_polar+numtasks_polar1))) then
           call ClustEmpole1_2b_NoOmplistsimult(
     &       ep2btmat,virep2btmat,dep2btmat,em2btmat,dem2btmat)
      end if



       call mpi_reduce(ep1btmat,ep1bmat,clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(dep1btmat,dep1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(virep1btmat,virep1bmat,3*3*clustcount,
     &      mpi_real8,mpi_sum,master,mpi_comm_world,ierr)
c
       call mpi_reduce(em1btmat,em1bmat,clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(dem1btmat,dem1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)

       call mpi_reduce(ep2btmat,ep2bmat,clustcount*clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(dep2btmat,dep2bmat,
     &      3*2*3*maxsizeclust*clustcount*clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(virep2btmat,virep2bmat,3*3*clustcount*clustcount,
     &      mpi_real8,mpi_sum,master,mpi_comm_world,ierr)
c
       call mpi_reduce(em2btmat,em2bmat,clustcount*clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(dem2btmat,dem2bmat,
     &      3*2*3*maxsizeclust*clustcount*clustcount,mpi_real8,
     &       mpi_sum,master,mpi_comm_world,ierr)


      deallocate(ep1btmat)
      deallocate(dep1btmat)
      deallocate(virep1btmat)
      deallocate(em1btmat)
      deallocate(dem1btmat)

      deallocate(ep2btmat)
      deallocate(dep2btmat)
      deallocate(virep2btmat)
      deallocate(em2btmat)
      deallocate(dem2btmat)

      if(taskid.eq.master) then
        allocate(pnum(2*3*maxsizeclust))
        do clust1=1,clustcount

           npole3b=0
           np1=0
           do i=1,sizeclust(clust1)
             moli1=clust(i,clust1)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           npole3b=np1

           ep3b = ep3b + ep1bmat(clust1)
           em = em + em1bmat(clust1)
           do l1 = 1, npole3b
             i = pnum(l1)
             do j = 1, 3
               dep3b(j,i)=  dep3b(j,i) + dep1bmat(j,l1,clust1)
               dem(j,i) = dem(j,i) + dem1bmat(j,l1,clust1)
             end do
           end do
           do i=1,3
              do j=1,3
                virep3b(j,i)=virep3b(j,i) + virep1bmat(j,i,clust1)
              end do
           end do

           print*,"clst1=",clust1,"ep1b=",ep1bmat(clust1)
           do clust2=clust1+1,clustcount
c              np2=0
c              do i=1,sizeclust(clust2)
c               moli1=clust(i,clust2)
c               pnum(np1+3*(i-1)+1)=imol(1,moli1)
c               pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
c               pnum(np1+3*(i-1)+3)=imol(2,moli1)
c               np2=np2+3
c              end do
c              npole3b=np1+np2
c
              ep3b = ep3b - ep1bmat(clust1)
              em = em - em1bmat(clust1)   
c              do l1=1,np1
c                do j=1,3
c                i=pnum(l1)
c                dep3b(j,i) = dep3b(j,i) - dep1bmat(j,l1,clust1)
c                dem(j,i) = dem(j,i) - dem1bmat(j,l1,clust1) 
c                end do
c              end do      
c
c              do i=1,3
c                do j=1,3
c                 virep3b(j,i)=virep3b(j,i)-virep1bmat(j,i,clust1)
c                end do
c              end do
c
c              npole3b=np2
c
              ep3b = ep3b - ep1bmat(clust2)
              em = em - em1bmat(clust2) 
c              do l1 =1, npole3b
c                 l2=l1+np1
c                 i=pnum(l2)
c                 do j=1,3
c                  dep3b(j,i) = dep3b(j,i)-dep1bmat(j,l1,clust2)
c                  dem(j,i) = dem(j,i)-dem1bmat(j,l1,clust2)
c                 end do
c              end do
c              do i=1,3
c                 do j=1,3
c                  virep3b(j,i)=virep3b(j,i)-virep1bmat(j,i,clust2)
c                 end do
c              end do
c
c              npole3b=np1+np2
c            
              ep3b = ep3b + ep2bmat(clust2,clust1)
      print*,"cst1=",clust1,"cst2=",clust2,"ep2b",ep2bmat(clust2,clust1)
c              em =  em + em2bmat(clust2,clust1)
c              do l1 = 1, npole3b
c                 i = pnum(l1)
c                 do j=1,3
c                  dep3b(j,i) = dep3b(j,i) + dep2bmat(j,l1,clust2,clust1)
c                  dem(j,i) = dem(j,i) + dem2bmat(j,l1,clust2,clust1)
c                 end do
c              end do
c              do i=1,3
c                do j=1,3
c                  virep3b(j,i)=virep3b(j,i) 
c     &               + virep2bmat(j,i,clust2,clust1)
c                end do
c              end do
c
           end do
        end do
        deallocate(pnum)
      end if

      return
      end
