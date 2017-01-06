      subroutine clust_gradient_PolarPerm1b2b3bsimult_reducematcommsmall
     & (start,
     & last,moli1rmndr)
      use mpole
      use mpidat
      use deriv3b
      use neigh2clust
      use aprx
      use deriv1bmat
      use energi
      use deriv
      use mpidat2
      use vdw
      use atoms
      use virial
      use mpidat3
      use molcul
      use neigh3b
      use deriv2bmat
      use cobar
      implicit none
      include 'mpif.h'
      integer stat(MPI_STATUS_SIZE)
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
      real*8, allocatable :: devtmp(:,:)
      real*8 evtmp,virevtmp(3,3)
      real*8, allocatable :: dep2bt(:,:)
      real*8, allocatable :: dem2bt(:,:)
      integer j1,clust1,npole3b,np1,l1,clust2,np2,l2,moli1
      integer, allocatable :: pnum(:)
      integer iter1,kouter,k3,cnt,k1
      real*8 xr1,yr1,zr1,r1_2,r1,Rcut2b,Rcut2b2
      real*8 ep2bt,em2bt,virep2bt(3,3),em2b
      real*8, allocatable :: dem2b(:,:)
      real*8, allocatable :: dep1bt(:,:)
      real*8, allocatable :: dem1bt(:,:)
      real*8 ep1bt,em1bt,virep1bt(3,3),em1b
      real*8, allocatable :: dep2btmat(:,:,:)
      real*8, allocatable :: ep2btmat(:)
      real*8, allocatable :: virep2btmat(:,:,:)
      integer clust3,ind,np3
      real*8 xr2,yr2,zr2,r2_2,r2,xr3,yr3,zr3,r3_2,r3,r
      
      if(taskid.eq.master) then
        ep2b=0.0d0
        em2b=0.0d0
        ep1b=0.0d0
        !em1b=0.0d0
         allocate (dep2b(3,n))
         allocate (dem2b(3,n))
         allocate (dep1b(3,n))
        ! allocate (dem1b(3,n))

        do i=1,3
          do j=1,3
            virep2b(j,i)=0.0d0
            virep1b(j,i)=0.0d0
          end do
        end do
        do i = 1, n
         do j=1,3
            dep2b(j,i) = 0.0d0
            dem2b(j,i) = 0.0d0
            dep1b(j,i) = 0.0d0
          !  dem1b(j,i) = 0.0d0
         end do
        end do
      end if
     

      if(taskid.lt.numtasks_polar1) then 
      allocate (dep1bt(3,npole))
      allocate(ep1btmat(clustcount))
      allocate(dep1btmat(3,3*maxsizeclust,clustcount))
      allocate(virep1btmat(3,3,clustcount))
      allocate (dem1bt(3,npole))
      allocate(em1btmat(clustcount))
      allocate(dem1btmat(3,3*maxsizeclust,clustcount))
c      allocate(virem1btmat(3,3,clustcount))
     
      do i=1,npole
         dep1bt(1,i)=0.0d0  
         dep1bt(2,i)=0.0d0
         dep1bt(3,i)=0.0d0   
         dem1bt(1,i)=0.0d0
         dem1bt(2,i)=0.0d0
         dem1bt(3,i)=0.0d0
      end do
      ep1bt=0.0d0
      em1bt=0.0d0
      do i=1,3
         do j=1,3
            virep1bt(j,i)=0.0d0
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
      end if

      allocate (dem2bt(3,npole))
      allocate (dep2bt(3,npole))
      do i=1,npole
         dep2bt(1,i)=0.0d0
         dep2bt(2,i)=0.0d0
         dep2bt(3,i)=0.0d0
         dem2bt(1,i)=0.0d0
         dem2bt(2,i)=0.0d0
         dem2bt(3,i)=0.0d0
      end do
      ep2bt=0.0d0
      em2bt=0.0d0
      do i=1,3
         do j=1,3
            virep2bt(j,i)=0.0d0
c            virem3bt(j,i)=0.0d0
         end do
      end do

      allocate(ep2btmat(num2bsave))
      allocate(dep2btmat(3,2*3*maxsizeclust,num2bsave))
      allocate(virep2btmat(3,3,num2bsave))
      do i=1,num2bsave
         ep2btmat(i)=0.0d0
         do j=1,2*3*maxsizeclust
            do k=1,3
              dep2btmat(k,j,i)=0.0d0
            end do
         end do
         do j=1,3
            do k=1,3
              virep2btmat(k,j,i)=0.0d0
            end do
         end do
      end do

      
c      allocate (dem3bt(3,npole))
      allocate (dep3bt(3,npole))
      do i=1,npole
         dep3bt(1,i)=0.0d0
         dep3bt(2,i)=0.0d0
         dep3bt(3,i)=0.0d0
c         dem3bt(1,i)=0.0d0
c         dem3bt(2,i)=0.0d0
c         dem3bt(3,i)=0.0d0
      end do
      ep3bt=0.0d0
c      em3bt=0.0d0
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
c            virem3bt(j,i)=0.0d0
         end do
      end do

      

c      print*,"numtasks_polar1=",numtasks_polar1
      
      if(taskid.lt.numtasks_polar1) then
      !if(taskid.eq.7) then          
       !   if(mlistclust) then
c            t1=mpi_Wtime()
            call ClustEmpole1_1b_NoOmplist
     &  (ep1bt,virep1bt,dep1bt,ntript,start,last,moli1rmndr,ep1btmat,
     &  dep1btmat,virep1btmat,em1bt,dem1bt,em1btmat,dem1btmat
     &  )
c            t2=mpi_Wtime()
c          print*,"Done w ClustEmpole1_1b_NoOmplist",ep3bt,taskid
c           print*,"taskid=",taskid,"1body Computation Cost=",t2-t1
       !   else
       !     call ClustEmpole1_1b_NoOmp
     & ! (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,ep1btmat,
     & ! dep1btmat,virep1btmat,em3bt,dem3bt,em1btmat,dem1btmat
     & ! )
       !   print*,"Done w ClustEmpole1_1b_NoOmplist",ep1bt,taskid
       !   end if
      else if( (taskid.ge.numtasks_polar1)
     &   .and.(taskid.lt.(numtasks_polar+numtasks_polar1))) then
          if(mlistclust) then
c           t5=mpi_Wtime()
           call ClustEmpole1_2b_NoOmplistRtrn2BmatNoSubtr1bSimult(
     &       ep2bt,virep2bt,dep2bt,ntript,em2bt,dem2bt,ep2btmat,
     &       dep2btmat,virep2btmat)
c           t6=mpi_Wtime()
         !   print*,"taskid=",taskid,"2body Computation Cost=",t6-t5
c        print*,"Done w ClustEmpole1_2b_NoOmplist",ep3bt,taskid
          else
c           call ClustEmpole1_2b_NoOmp(
c     &       ep3bt,virep3bt,dep3bt,ntript,em3bt,dem3bt)
c        print*,"Done w ClustEmpole1_2b_NoOmp",ep2bt,taskid
          end if
       ! print*,"Done w ClustEmpole1_2b_NoOmplist",ep2bt,taskid
        
     &  
      else if( (taskid.ge.(numtasks_polar+numtasks_polar1))
     &    .and. (taskid.lt.(numtasks_polar+numtasks_polar1
     &         +numtasks_polar3b))) then 
           call ClustEmpole1_3b_NoOmplistRtrn3BNoSubtr1b2bSimult(
     &       ep3bt,virep3bt,dep3bt,ntript)
       ! print*,"Done w ClustEmpole1_3b_NoOmplist",ep3bt,taskid
      end if


c      print*,"In clust1b2b After 1b taskid=",taskid,"ep3bt=",ep3bt

c       call mpi_barrier(mpi_comm_world,ierr)

c            t3=mpi_Wtime()
       if(taskid.lt.numtasks_polar1) then
       call mpi_reduce(ep1btmat,ep1bmat,clustcount,mpi_real8,
     &       mpi_sum,master,emreal_comm,reqs19,ierr)
       call mpi_reduce(dep1btmat,dep1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,master,emreal_comm,reqs20,ierr)
       call mpi_reduce(virep1btmat,virep1bmat,3*3*clustcount,
     &      mpi_real8,mpi_sum,master,emreal_comm,reqs21,ierr)
c
       call mpi_reduce(em1btmat,em1bmat,clustcount,mpi_real8,
     &       mpi_sum,master,emreal_comm,reqs22,ierr)
       call mpi_reduce(dem1btmat,dem1bmat,
     &      3*3*maxsizeclust*clustcount,mpi_real8,
     &       mpi_sum,master,emreal_comm,reqs23,ierr)

                  call mpi_reduce(ep1bt,ep1b,1,mpi_real8,
     &          mpi_sum,master,emreal_comm,reqs24,ierr)
                  call mpi_reduce(dep1bt,dep1b,npole*3,mpi_real8,
     &          mpi_sum,master,emreal_comm,reqs25,ierr)
                  call mpi_reduce(virep1bt,virep1b,3*3,mpi_real8,
     &          mpi_sum,master,emreal_comm,reqs26,ierr)

                  call mpi_reduce(em1bt,em,1,mpi_real8,
     &          mpi_sum,master,emreal_comm,reqs27,ierr)
                  call mpi_reduce(dem1bt,dem,npole*3,mpi_real8,
     &          mpi_sum,master,emreal_comm,reqs28,ierr)
c           t8=mpi_Wtime()
c            print*,"taskid=",taskid,"Reduce Cost=",t8-t7
       end if

      if(taskid.lt.numtasks_polar1) then
      deallocate(ep1btmat)
      deallocate(dep1btmat)
      deallocate(virep1btmat)
      deallocate(em1btmat)
      deallocate(dem1btmat)
      deallocate(dep1bt)
      deallocate(dem1bt)
      end if

                  call mpi_reduce(ep2bt,ep2b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep2bt,dep2b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep2bt,virep2b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

                  call mpi_reduce(em2bt,em2b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dem2bt,dem2b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

                  call mpi_reduce(ep2btmat,ep2bmatmod,num2bsave,
     &            mpi_real8,mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep2btmat,dep2bmatmod,
     &                3*2*3*maxsizeclust*num2bsave,mpi_real8,
     &                mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep2btmat,virep2bmatmod,
     &                3*3*num2bsave,mpi_real8,
     &                mpi_sum,master,mpi_comm_world,ierr)

                  call mpi_reduce(ep3bt,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

c           t8=mpi_Wtime()
c            print*,"taskid=",taskid,"Reduce Cost=",t8-t7


      deallocate(dep2bt)
      deallocate(dem2bt)
      deallocate(ep2btmat)
      deallocate(dep2btmat)
      deallocate(virep2btmat)
      deallocate(dep3bt)

      if(taskid.eq.master) then
        Rcut2b=cut2b_input
        Rcut2b2=Rcut2b*Rcut2b
        allocate(pnum(3*3*maxsizeclust))
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

c                  npole3b=np1+np2
c                  ep3b = ep3b - ep1bmat(clust1)
c                  em = em - em1bmat(clust1)
                    ep2b = ep2b - ep1bmat(clust1)
                    em2b = em2b - em1bmat(clust1)

                    do l1 = 1, np1
                     i = pnum(l1)
                     do j = 1, 3
c                     dep3b(j,i) = dep3b(j,i)-dep1bmat(j,l1,clust1)
c                     dem(j,i) = dem(j,i)-dem1bmat(j,l1,clust1)

                      dep2b(j,i) = dep2b(j,i)-dep1bmat(j,l1,clust1)
                      dem2b(j,i) = dem2b(j,i)-dem1bmat(j,l1,clust1)
                     end do
                    end do

                    do i=1,3
                     do j=1,3
c                       virep3b(j,i)=virep3b(j,i)-virep1bmat(j,i,clust1)
                      virep2b(j,i)=virep2b(j,i)-virep1bmat(j,i,clust1)
                     end do
                    end do
                    npole3b=np2

c                  ep3b = ep3b - ep1bmat(clust2)
c                  em = em - em1bmat(clust2)
                    ep2b = ep2b - ep1bmat(clust2)
                    em2b = em2b - em1bmat(clust2)

                     do l1 =1, npole3b
                      l2=l1+np1
                      i=pnum(l2)
                      do j=1,3
c                      dep3b(j,i) = dep3b(j,i)-dep1bmat(j,l1,clust2)
c                      dem(j,i) = dem(j,i)-dem1bmat(j,l1,clust2)

                       dep2b(j,i) = dep2b(j,i)-dep1bmat(j,l1,clust2)
                       dem2b(j,i) = dem2b(j,i)-dem1bmat(j,l1,clust2)
                      end do
                     end do
                     do i=1,3
                      do j=1,3
c                       virep3b(j,i)=virep3b(j,i)-virep1bmat(j,i,clust2)
                       virep2b(j,i)=virep2b(j,i)-virep1bmat(j,i,clust2)
                      end do
                     end do
               
                  end if
               end do
              end do

           end do
        end do

        !ep3b=ep3b+ep2b
        em=em+em2b
        !do i=1,3
        !   do j=1,3
        !      virep3b(j,i)=virep3b(j,i)+virep2b(j,i)
        !   end do
        !end do
        do i=1,npole
           do j=1,3
        !      dep3b(j,i)=dep3b(j,i)+dep2b(j,i)
              dem(j,i)=dem(j,i)+dem2b(j,i)
           end do
        end do
        !deallocate(dep2b)
        deallocate(dem2b)




      do iter1=0,numtasks_polar3b-1
       do kouter=start_polar3b(iter1),last_polar3b(iter1)
              clust1=molnew3b(kouter)

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

        do k3=1,num_mollst_chunk3b(clust1,iter1)
         do k1=start_mollst_polar3b(k3,clust1,iter1),
     &                  last_mollst_polar3b(k3,clust1,iter1),2
           clust2=mollst3(k1,clust1)
           clust3=mollst3(k1+1,clust1)

           xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
           yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
           zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
           call image(xr1,yr1,zr1)
           r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
           r1=sqrt(r1_2)

           xr2=clust_cm(1,clust3)-clust_cm(1,clust1)
           yr2=clust_cm(2,clust3)-clust_cm(2,clust1)
           zr2=clust_cm(3,clust3)-clust_cm(3,clust1)
           call image(xr2,yr2,zr2)
           r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
           r2=sqrt(r2_2)

           xr3=clust_cm(1,clust3)-clust_cm(1,clust2)
           yr3=clust_cm(2,clust3)-clust_cm(2,clust2)
           zr3=clust_cm(3,clust3)-clust_cm(3,clust2)
           call image(xr3,yr3,zr3)
           r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
           r3=sqrt(r3_2)

            if((r1.lt.r3).and.(r2.lt.r3) ) then
                 r=r1+r2
              ! print*,"c1=",clust1,"c2=",clust2,"c3=",clust3,"Cobar=",r
            else if ((r1.lt.r2).and.(r3.lt.r2)) then
                 r=r1+r3
              ! print*,"c1=",clust1,"c2=",clust2,"c3=",clust3,"Cobar=",r
            else if ((r2.lt.r1).and.(r3.lt.r1)) then
                 r=r2+r3
              ! print*,"c1=",clust1,"c2=",clust2,"c3=",clust3,"Cobar=",r
            end if
           cnt=cnt+1

          if(r.le.Cobarcut3b) then
          np2=0
          do i=1,sizeclust(clust2)
             moli1=clust(i,clust2)
             pnum(np1+3*(i-1)+1)=imol(1,moli1)
             pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
             pnum(np1+3*(i-1)+3)=imol(2,moli1)
             np2=np2+3
          end do

          np3=0
          do i=1,sizeclust(clust3)
             moli1=clust(i,clust3)
             pnum(np1+np2+3*(i-1)+1)=imol(1,moli1)
             pnum(np1+np2+3*(i-1)+2)=imol(1,moli1)+1
             pnum(np1+np2+3*(i-1)+3)=imol(2,moli1)
             np3=np3+3
          end do

          npole3b=np1+np2

          ind=ind2b_counter(clust2,clust1)
          ep3b =  ep3b - ep2bmatmod(ind)          
          do l1=1,npole3b
             i=pnum(l1)   
             do j=1,3
                dep3b(j,i)=dep3b(j,i)-dep2bmatmod(j,l1,ind)
             end do
          end do
          do i=1,3
            do j=1,3
              virep3b(j,i)=virep3b(j,i)-virep2bmatmod(j,i,ind)
            end do
          end do 

          ind=ind2b_counter(clust3,clust1)
          ep3b = ep3b - ep2bmatmod(ind)
          do l1=1,np1
             i=pnum(l1)
             do j=1,3
               dep3b(j,i)=dep3b(j,i)-dep2bmatmod(j,l1,ind)
             end do             
          end do
          do l1=1,np3
             i=pnum(np1+np2+l1)
             l2=np1+l1
             do j=1,3
                dep3b(j,i)=dep3b(j,i)-dep2bmatmod(j,l2,ind)
             end do  
          end do
          do i=1,3
            do j=1,3
              virep3b(j,i)=virep3b(j,i)-virep2bmatmod(j,i,ind)
            end do
          end do

          npole3b=np2+np3
          ind=ind2b_counter(clust3,clust2)
          ep3b = ep3b - ep2bmatmod(ind)
          do l1=1,npole3b
             i=pnum(np1+l1)
             do j=1,3
                dep3b(j,i)=dep3b(j,i)-dep2bmatmod(j,l1,ind)
             end do
          end do
          do i=1,3
             do j=1,3
                virep3b(j,i)=virep3b(j,i)-virep2bmatmod(j,i,ind)
             end do
          end do


          ep3b = ep3b + ep1bmat(clust1)
          do l1 =1,np1
             i=pnum(l1)
             do j=1,3
                dep3b(j,i)=dep3b(j,i)+dep1bmat(j,l1,clust1)
             end do
          end do
          do i=1,3
             do j=1,3
                virep3b(j,i)=virep3b(j,i)+virep1bmat(j,i,clust1)
             end do
          end do

          ep3b = ep3b + ep1bmat(clust2)
          do l1=1,np2
             i=pnum(np1+l1)
             do j=1,3
                dep3b(j,i)=dep3b(j,i)+dep1bmat(j,l1,clust2)
             end do
          end do
          do i=1,3
             do j=1,3
                virep3b(j,i)=virep3b(j,i)+virep1bmat(j,i,clust2)
             end do
          end do

          ep3b = ep3b + ep1bmat(clust3)
          do l1=1,np3
             i=pnum(np1+np2+l1)
             do j=1,3
                dep3b(j,i)=dep3b(j,i)+dep1bmat(j,l1,clust3)
             end do
          end do
          do i=1,3
             do j=1,3
                virep3b(j,i)=virep3b(j,i)+virep1bmat(j,i,clust3)
             end do
          end do

          end if
         end do
        end do
       end do
      end do

      ep3b=ep1b+ep2b+ep3b
      do i=1,npole
        do j=1,3
           dep3b(j,i)=dep1b(j,i)+dep2b(j,i)+dep3b(j,i) 
        end do
      end do
      do i=1,3
        do j=1,3
          virep3b(j,i)=virep1b(j,i)+virep2b(j,i)+virep3b(j,i)
        end do
      end do
      deallocate(dep2b)
      deallocate(dep1b)

      end if
      return
      end

