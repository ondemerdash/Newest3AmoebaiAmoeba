      subroutine polar_load_balance_clust_splitlist_offset3b_simple
      use mpidat
      use neigh3b
      use molcul
      use neigh2clust
      implicit none
      include 'mpif.h'
      integer i,j,nmollst_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load,i1,counter,k1
      logical first
      integer numtasks_offset
      save first
      data first  / .true. /
            
      nmollst_tot=0
      nmol_nonzero=0

      
      if (first) then
         if (.not. allocated(molnew3b))
     &          allocate (molnew3b(clustcount))
      end if      

      max_nmollst3b=0
      do i=1,clustcount
        if(nmollst3(i).ne.0) then
          if(nmollst3(i).gt.max_nmollst3b) then 
              max_nmollst3b=nmollst3(i)   
          end if
c          do j=1,nmollst3(i),2
c             k=mollst3(j,i)
c             k2=mollst3(j+1,i)
c             nmollst_tot = nmollst_tot+sizeclust(k2)+sizeclust(k)
c     &                     +sizeclust(i)
c          end do

          nmol_nonzero=nmol_nonzero+1
          molnew3b(nmol_nonzero)=i
      !    print*,"cluster_index, nmollst3(i)",i,nmollst3(i)
        end if
      end do

        numtasks_offset=cnt3bwork

      if (first) then
         first = .false.
        if (.not. allocated(start_polar3b))
     &          allocate (start_polar3b(0:numtasks-1))

        if (.not. allocated(start_mollst_polar3b))
     &          allocate (start_mollst_polar3b(max_nmollst3b,clustcount,
     &          0:numtasks-1))

        if (.not. allocated(last_polar3b))
     &          allocate (last_polar3b(0:numtasks-1))

        if (.not. allocated(last_mollst_polar3b))
     &          allocate (last_mollst_polar3b(max_nmollst3b,clustcount,
     &          0:numtasks-1))
        if (.not. allocated(num_mollst_chunk3b))
     &       allocate (num_mollst_chunk3b(clustcount,0:numtasks-1))
      end if



      
      taskcount=0
      counter=0 
      do i1=1,nmol_nonzero
         i=molnew3b(i1)
           do j=1,nmollst3(i),2
             ! k=mollst3(j,i)
             ! k2=mollst3(j+1,i)
              !if(j.eq.1) then
              start_polar3b(taskcount)=i1
              last_polar3b(taskcount)=i1
              start_mollst_polar3b(1,i,taskcount)=j  
              last_mollst_polar3b(1,i,taskcount)=j+1
              num_mollst_chunk3b(i,taskcount)=1
              taskcount=taskcount+1
           end do
      end do           
               

      !do i=0,taskcount-1
      !   do k1=start_polar3b(i), last_polar3b(i)
      !      j=molnew3b(k1)
      ! print*,"3btask=",i,"clust=",j,"numchunk=",num_mollst_chunk3b(j,i)
      !   end do
      !end do

      !do i=0,taskcount-1
      !print*,"3btask molnwstart,molnewend",i,molnew3b(start_polar3b(i)),
     &!           molnew3b(last_polar3b(i))
      !   do k1=start_polar3b(i), last_polar3b(i)
      !      j=molnew3b(k1)
      !      do k=1,num_mollst_chunk3b(j,i)
      !        print*,"3btskid clust start_mollst last_mollst",i,j,
     &!          start_mollst_polar3b(k,j,i),last_mollst_polar3b(k,j,i)
      !      end do 
      !   end do
      !end do

      numtasks_polar3b=taskcount
      print*,"Listloadbal: numtasks_polar3b=",numtasks_polar3b
      return
      end

      subroutine bcast_polar_load_balance_clust_splitlist3b
      use mpidat
      use molcul
      use neigh2clust
      implicit none
      include 'mpif.h'
      integer ierr
      logical first
      save first
      data first  / .true. /

      call mpi_bcast(max_nmollst3b,1,mpi_integer,master,
     &   mpi_comm_world,ierr)


      if (first.and.(taskid.ne.master)) then
         first = .false.
c         if (.not. allocated(start_mollst_polar))
c     &          allocate (start_mollst_polar(0:numtasks-1))
c         if (.not. allocated(last_mollst_polar))
c     &          allocate (last_mollst_polar(0:numtasks-1))
         if (.not. allocated(molnew3b))
     &          allocate (molnew3b(clustcount))

         if (.not. allocated(start_polar3b))
     &          allocate (start_polar3b(0:numtasks-1))

         if (.not. allocated(start_mollst_polar3b))
     &          allocate (start_mollst_polar3b(max_nmollst3b,clustcount,
     &          0:numtasks-1))

         if (.not. allocated(last_polar3b))
     &          allocate (last_polar3b(0:numtasks-1))

         if (.not. allocated(last_mollst_polar3b))
     &          allocate (last_mollst_polar3b(max_nmollst3b,clustcount,
     &          0:numtasks-1))
         if (.not. allocated(num_mollst_chunk3b))
     &          allocate (num_mollst_chunk3b(clustcount,0:numtasks-1))

      end if
c         print*,"numtasks_polar in bcast",numtasks_polar
         call mpi_bcast(numtasks_polar3b,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(start_polar3b,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(last_polar3b,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(molnew3b,clustcount,mpi_integer,master,
     &   mpi_comm_world,ierr)

c       call mpi_bcast(start_mollst_polar,numtasks,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c       call mpi_bcast(last_mollst_polar,numtasks,mpi_integer,master,
c     &   mpi_comm_world,ierr)


      call mpi_bcast(start_mollst_polar3b,
     &    numtasks*clustcount*max_nmollst3b
     &   ,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(last_mollst_polar3b,
     &    numtasks*clustcount*max_nmollst3b
     &   ,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(num_mollst_chunk3b,numtasks*clustcount
     &   ,mpi_integer,master,mpi_comm_world,ierr)

      !print*,"Successful Load Bal Bcasts! Have a nice day!"
      return
      end

