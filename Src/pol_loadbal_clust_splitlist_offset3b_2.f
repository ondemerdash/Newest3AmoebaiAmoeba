      subroutine polar_load_balance_clust_splitlist_offset3b_2
      use mpidat
      use neigh3b
      use molcul
      use neigh2clust
      implicit none
      include 'mpif.h'
      integer i,j,nmollst_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load,i1,counter,k1
      logical first
      integer numtasks_offset,k2
      save first
      data first  / .true. /
            
      nmollst_tot=0
      nmol_nonzero=0

      
      if (first) then
         if (.not. allocated(molnew3b))
     &          allocate (molnew3b(clustcount))
      end if      

      print*,"cnt3bwork=",cnt3bwork

      max_nmollst3b=0
      do i=1,clustcount
        if(nmollst3(i).ne.0) then
          if(nmollst3(i).gt.max_nmollst3b) then 
              max_nmollst3b=nmollst3(i)   
          end if
          do j=1,nmollst3(i),2
             k=mollst3(j,i)
             k2=mollst3(j+1,i)
             nmollst_tot = nmollst_tot+sizeclust(k2)+sizeclust(k)
     &                     +sizeclust(i)
          end do

          nmol_nonzero=nmol_nonzero+1
          molnew3b(nmol_nonzero)=i
          print*,"cluster_index, nmollst3(i)",i,nmollst3(i)
        end if
      end do

c        numtasks_offset=cnt3bwork
        numtasks_offset=numtasks/2

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


c      do i=1,nmol_nonzero
c       print*,"molnew",molnew(i),"nmollst3(molnew(i))",nmollst3(molnew(i))
c      end do

c      print*,"nmol_nonzero",nmol_nonzero

      avg_load=nmollst_tot/numtasks_offset
      print*,"avg_load=",avg_load

      if(mod(nmollst_tot,numtasks_offset).gt.0) then
        avg_load=avg_load+1
        print*,"avg_load rmndr",avg_load
      end if


c      countmol=1
c      taskcount=0
c      do while (countmol.le.nmol_nonzero)
c         load=0
c         do while (load.lt.avg_load)
c            if(load.eq.0) then
c              start_polar(taskcount)=countmol
c            end if        
c
c            k=molnew(countmol) 
c            load=load+nmollst3(k)          
c            countmol=countmol+1
c
c            if(countmol.eq.nmol_nonzero+1) then
c               goto 31 
c            end if
c         end do
c   31             continue
c         print*,"taskcount load",taskcount,load
c         last_polar(taskcount)=countmol-1
c         taskcount=taskcount+1
c      end do

      
      load=0
      taskcount=0
      counter=0
c      do i=1,clustcount
      do i1=1,nmol_nonzero
         i=molnew3b(i1)
         counter=0
           do j=1,nmollst3(i),2
              k=mollst3(j,i)  
              k2=mollst3(j+1,i)            
              if(load.eq.0) then
                start_polar3b(taskcount)=i1
                counter=counter+1
                start_mollst_polar3b(counter,i,taskcount)=j
                num_mollst_chunk3b(i,taskcount)=counter
              else if((load.ne.0).and.(j.eq.1)) then
                counter=counter+1
                start_mollst_polar3b(counter,i,taskcount)=j
                num_mollst_chunk3b(i,taskcount)=counter
              end if

              load=load+sizeclust(k2)+sizeclust(k)+sizeclust(i)
              if(load.ge.avg_load) then
                last_polar3b(taskcount)=i1
                last_mollst_polar3b(counter,i,taskcount)=j+1
                taskcount=taskcount+1
                counter=0
              end if 
              if((load.lt.avg_load).and.(j.eq.nmollst3(i)-1)) then
                last_mollst_polar3b(counter,i,taskcount)=j+1
              end if

c              if( (load.ge.avg_load) .or. 
c     &           ((load.lt.avg_load).and.(j.eq.nmollst3(i))) ) then 
c                 counter=counter+1
                 if(load.ge.avg_load) then
                   load=0
                 end if
c              end if
           end do
      end do
c      print*,"taskcount=",taskcount
c      if(taskcount.gt.numtasks) then
c        print*,"Tilt! Taskcount exceeds total number of tasks!"
c      end if
      do i=0,taskcount-1
         do k1=start_polar3b(i), last_polar3b(i)
            j=molnew3b(k1)
       print*,"3btask=",i,"clust=",j,"numchunk=",num_mollst_chunk3b(j,i)
         end do
      end do

      do i=0,taskcount-1
c      print*,"taskid start_polar last_polar",i,start_polar(i),
c     &   last_polar(i)
      print*,"3btask molnwstart,molnewend",i,molnew3b(start_polar3b(i)),
     &           molnew3b(last_polar3b(i))
         do k1=start_polar3b(i), last_polar3b(i)
            j=molnew3b(k1)
            do k=1,num_mollst_chunk3b(j,i)
              print*,"3btskid clust start_mollst last_mollst",i,j,
     &          start_mollst_polar3b(k,j,i),last_mollst_polar3b(k,j,i)
            end do 
         end do
      end do

      numtasks_polar3b=taskcount
      print*,"Listloadbal: numtasks_polar3b=",numtasks_polar3b
      return
      end

c      subroutine bcast_polar_load_balance_clust_splitlist3b
c      use mpidat
c      use molcul
c      use neigh2clust
c      implicit none
c      include 'mpif.h'
c      integer ierr
c      logical first
c      save first
c      data first  / .true. /
c
c      call mpi_bcast(max_nmollst3b,1,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c
c
c      if (first.and.(taskid.ne.master)) then
c         first = .false.
cc         if (.not. allocated(start_mollst_polar))
cc     &          allocate (start_mollst_polar(0:numtasks-1))
cc         if (.not. allocated(last_mollst_polar))
cc     &          allocate (last_mollst_polar(0:numtasks-1))
c         if (.not. allocated(molnew3b))
c     &          allocate (molnew3b(clustcount))
c
c         if (.not. allocated(start_polar3b))
c     &          allocate (start_polar3b(0:numtasks-1))
c
c         if (.not. allocated(start_mollst_polar3b))
c     &          allocate (start_mollst_polar3b(max_nmollst3b,clustcount,
c     &          0:numtasks-1))
c
c         if (.not. allocated(last_polar3b))
c     &          allocate (last_polar3b(0:numtasks-1))
c
c         if (.not. allocated(last_mollst_polar3b))
c     &          allocate (last_mollst_polar3b(max_nmollst3b,clustcount,
c     &          0:numtasks-1))
c         if (.not. allocated(num_mollst_chunk3b))
c     &          allocate (num_mollst_chunk3b(clustcount,0:numtasks-1))
c
c      end if
cc         print*,"numtasks_polar in bcast",numtasks_polar
c         call mpi_bcast(numtasks_polar3b,1,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(start_polar3b,numtasks,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(last_polar3b,numtasks,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(molnew3b,clustcount,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c
cc       call mpi_bcast(start_mollst_polar,numtasks,mpi_integer,master,
cc     &   mpi_comm_world,ierr)
cc       call mpi_bcast(last_mollst_polar,numtasks,mpi_integer,master,
cc     &   mpi_comm_world,ierr)
c
c
c      call mpi_bcast(start_mollst_polar3b,
c     &    numtasks*clustcount*max_nmollst3b
c     &   ,mpi_integer,master,mpi_comm_world,ierr)
c      call mpi_bcast(last_mollst_polar3b,
c     &    numtasks*clustcount*max_nmollst3b
c     &   ,mpi_integer,master,mpi_comm_world,ierr)
c      call mpi_bcast(num_mollst_chunk3b,numtasks*clustcount
c     &   ,mpi_integer,master,mpi_comm_world,ierr)
c
c      !print*,"Successful Load Bal Bcasts! Have a nice day!"
c      return
c      end

