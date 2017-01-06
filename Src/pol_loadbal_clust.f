      subroutine polar_load_balance_clust
      use mpidat
      use neigh3b
      use molcul
      use neigh2clust
      implicit none
      include 'mpif.h'
      integer i,j,nmollst_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load
c      integer, allocatable :: start_polar(:)
c      integer, allocatable :: last_polar(:)
c      integer, allocatable :: molnew(:)      
      logical first
      save first
      data first  / .true. /
            
      nmollst_tot=0
      nmol_nonzero=0

      if (first) then
         first = .false.
         if (.not. allocated(start_polar)) 
     &          allocate (start_polar(0:numtasks-1))
         if (.not. allocated(last_polar)) 
     &          allocate (last_polar(0:numtasks-1))
         if (.not. allocated(molnew))
     &          allocate (molnew(clustcount))
      end if      

      do i=1,clustcount
        if(nmollst(i).ne.0) then
          nmollst_tot = nmollst_tot + nmollst(i)
          nmol_nonzero=nmol_nonzero+1
          molnew(nmol_nonzero)=i
          print*,"cluster_index, nmollst(i)",i,nmollst(i)
        end if
      end do

      do i=1,nmol_nonzero
       print*,"molnew",molnew(i),"nmollst(molnew(i))",nmollst(molnew(i))
      end do

c      print*,"nmol_nonzero",nmol_nonzero

      avg_load=nmollst_tot/numtasks
c      print*,"avg_load=",avg_load

      if(mod(nmollst_tot,numtasks).gt.0) then
        avg_load=avg_load+1
        print*,"avg_load rmndr",avg_load
      end if


      countmol=1
      taskcount=0

      do while (countmol.le.nmol_nonzero)
         load=0
         do while (load.lt.avg_load)
            if(load.eq.0) then
              start_polar(taskcount)=countmol
            end if        

            k=molnew(countmol) 
            load=load+nmollst(k)          
            countmol=countmol+1

            if(countmol.eq.nmol_nonzero+1) then
               goto 31 
            end if
         end do
   31             continue
c         print*,"taskcount load",taskcount,load
         last_polar(taskcount)=countmol-1
         taskcount=taskcount+1
      end do

c      print*,"taskcount=",taskcount
c      if(taskcount.gt.numtasks) then
c        print*,"Tilt! Taskcount exceeds total number of tasks!"
c      end if

      do i=0,taskcount-1
         print*,"taskid nmol_nonzero start_polar",i,start_polar(i)
         print*,"taskid nmol start",i,molnew(start_polar(i))
         print*,"taskid nmol_nonzero last_polar",i,last_polar(i)
         print*,"taskid nmol end",i,molnew(last_polar(i))
      end do

      numtasks_polar=taskcount
      print*,"Listloadbal: numtasks_polar=",numtasks_polar
      return
      end

      subroutine bcast_polar_load_balance_clust
      use mpidat
      use molcul
      use neigh2clust
      implicit none
      include 'mpif.h'
      integer ierr
      logical first
      save first
      data first  / .true. /

      if (first.and.(taskid.ne.master)) then
         first = .false.
         if (.not. allocated(start_polar))
     &          allocate (start_polar(0:numtasks-1))
         if (.not. allocated(last_polar))
     &          allocate (last_polar(0:numtasks-1))
         if (.not. allocated(molnew))
     &          allocate (molnew(clustcount))
      end if
c         print*,"numtasks_polar in bcast",numtasks_polar
         call mpi_bcast(numtasks_polar,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(start_polar,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(last_polar,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(molnew,clustcount,mpi_integer,master,
     &   mpi_comm_world,ierr)
      !print*,"Successful Load Bal Bcasts! Have a nice day!"
      return
      end   
