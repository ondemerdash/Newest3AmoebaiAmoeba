      subroutine clust_gradient_PolarPerm1b2bsimultmod(start,
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
      use neigh3b
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
      integer iter1,kouter,k3,cnt,k1
      real*8 xr1,yr1,zr1,r1_2,r1,Rcut2b,Rcut2b2


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
        Rcut2b=cut2b_input
        Rcut2b2=Rcut2b*Rcut2b

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
        end do 

      allocate (dep3bt(3,npole))
      allocate (dem3bt(3,npole))
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
         end do 
      end do 

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,
!$OMP& nmollst,mollst,
!$OMP& Rcut2b,molnew,
!$OMP& Rcut2b2,sizeclust,clust,clust_cm,ep1bmat,dep1bmat,virep1bmat,
!$OMP& start_mollst_polar,last_mollst_polar,num_mollst_chunk,
!$OMP& numtasks_polar,em1bmat,dem1bmat,
!$OMP& ep2bmat,dep2bmat,virep2bmat,em2bmat,dem2bmat,em3bt,dem3bt)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,em3bt,dem3bt)
!$OMP& schedule(guided)

        do iter1=0,numtasks_polar-1
           do kouter=start_polar(iter1),last_polar(iter1)
              clust1=molnew(kouter)

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

              do k3=1,num_mollst_chunk(clust1,iter1)
               do k1=start_mollst_polar(k3,clust1,iter1),
     &                     last_mollst_polar(k3,clust1,iter1)
                     clust2=mollst(k1,clust1)
                 xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
                 yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
                 zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
                 call image(xr1,yr1,zr1)
                 r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
c                 r1=sqrt(r1_2)
c                 cnt=cnt+1

                if(r1_2.le.Rcut2b2) then

                  np2=0
                  do i=1,sizeclust(clust2)
                    moli1=clust(i,clust2)
                    pnum(np1+3*(i-1)+1)=imol(1,moli1)
                    pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
                    pnum(np1+3*(i-1)+3)=imol(2,moli1)
                    np2=np2+3
                  end do

                  npole3b=np1+np2
c                  ep3b = ep3b - ep1bmat(clust1)
c                  em = em - em1bmat(clust1)
                  ep3bt = ep3bt - ep1bmat(clust1)
                  em3bt = em3bt - em1bmat(clust1)

                  do l1 = 1, np1
                    i = pnum(l1)
                    do j = 1, 3
c                     dep3b(j,i) = dep3b(j,i)-dep1bmat(j,l1,clust1)
c                     dem(j,i) = dem(j,i)-dem1bmat(j,l1,clust1)

                     dep3bt(j,i) = dep3bt(j,i)-dep1bmat(j,l1,clust1)
                     dem3bt(j,i) = dem3bt(j,i)-dem1bmat(j,l1,clust1)
                    end do
                  end do

                  do i=1,3
                    do j=1,3
c                       virep3b(j,i)=virep3b(j,i)-virep1bmat(j,i,clust1)
                      virep3bt(j,i)=virep3bt(j,i)-virep1bmat(j,i,clust1)
                    end do
                  end do

                  npole3b=np2

c                  ep3b = ep3b - ep1bmat(clust2)
c                  em = em - em1bmat(clust2)
                  ep3bt = ep3bt - ep1bmat(clust2)
                  em3bt = em3bt - em1bmat(clust2)

                  do l1 =1, npole3b
                     l2=l1+np1
                     i=pnum(l2)
                     do j=1,3
c                      dep3b(j,i) = dep3b(j,i)-dep1bmat(j,l1,clust2)
c                      dem(j,i) = dem(j,i)-dem1bmat(j,l1,clust2)

                      dep3bt(j,i) = dep3bt(j,i)-dep1bmat(j,l1,clust2)
                      dem3bt(j,i) = dem3bt(j,i)-dem1bmat(j,l1,clust2)
                     end do
                  end do
                  do i=1,3
                     do j=1,3
c                       virep3b(j,i)=virep3b(j,i)-virep1bmat(j,i,clust2)
                      virep3bt(j,i)=virep3bt(j,i)-virep1bmat(j,i,clust2)
                     end do
                  end do

                  npole3b=np1+np2
c                  ep3b=ep3b+ep2bmat(clust2,clust1)
c                  em=em+em2bmat(clust2,clust1)

                  ep3bt=ep3bt+ep2bmat(clust2,clust1)
                  em3bt=em3bt+em2bmat(clust2,clust1)

                  do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
c                  dep3b(j,i) = dep3b(j,i) + dep2bmat(j,l1,clust2,clust1)
c                  dem(j,i) = dem(j,i) + dem2bmat(j,l1,clust2,clust1)

                  dep3bt(j,i)=dep3bt(j,i) + dep2bmat(j,l1,clust2,clust1)
                  dem3bt(j,i)=dem3bt(j,i) + dem2bmat(j,l1,clust2,clust1)
                    end do
                  end do
                  do i=1,3
                     do j=1,3
c                        virep3b(j,i)=virep3b(j,i)
c     &                    +virep2bmat(j,i,clust2,clust1)
                        virep3bt(j,i)=virep3bt(j,i)
     &                    +virep2bmat(j,i,clust2,clust1)

                     end do
                  end do

                end if
               end do
              end do

           end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        ep3b=ep3b+ep3bt
        em=em+em3bt
        do i=1,3
           do j=1,3
              virep3b(j,i)=virep3b(j,i)+virep3bt(j,i)
           end do
        end do       
        do i=1,npole
           do j=1,3
              dep3b(j,i)=dep3b(j,i)+dep3bt(j,i)
              dem(j,i)=dem(j,i)+dem3bt(j,i)
           end do
        end do

        deallocate (dep3bt)
        deallocate (dem3bt)
        deallocate(pnum)
      end if

      return
      end
