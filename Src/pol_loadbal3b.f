      subroutine polar_load_balance3b_sizelist
      use mpidat
      use neigh3b
      use molcul
      implicit none
      include 'mpif.h'
      integer i,j,nmollst3_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load,itercount,maxsize
      logical first
      save first
      data first  / .true. /
            
      nmollst3_tot=0
c      nmol_nonzero=0

      if (first) then
         first = .false.
         if (.not. allocated(start_polar3)) 
     &          allocate (start_polar3(0:numtasks-1))
         if (.not. allocated(last_polar3)) 
     &          allocate (last_polar3(0:numtasks-1))
c         if (.not. allocated(molnew3))
c     &          allocate (molnew3(nmol))
         if (.not. allocated(maxsize_mollst3))
     &         allocate (maxsize_mollst3(0:numtasks-1))

      end if      

      do i=1,nmol
          nmollst3_tot = nmollst3_tot + nmollst3(i)
      end do


      avg_load=nmollst3_tot/numtasks
c      print*,"avg_load=",avg_load

      if(mod(nmollst3_tot,numtasks).gt.0) then
        avg_load=avg_load+1
c        print*,"avg_load rmndr",avg_load
      end if


      countmol=1
      taskcount=0

      do while (countmol.lt.nmol)
         load=0
         itercount=0
         maxsize=0
         do while (load.lt.avg_load)
            if(itercount.eq.0) then
               start_polar3(taskcount)=countmol
               maxsize=nmollst3(countmol)
            else
               if(nmollst3(countmol).gt.maxsize) then
                  maxsize=nmollst3(countmol)
               end if
            end if        

           ! k=molnew3(countmol) 
            load=load+nmollst3(countmol)          
            countmol=countmol+1
            itercount=itercount+1

            if(countmol.eq.nmol+1) then
               goto 31 
            end if
         end do
   31             continue
c         print*,"taskcount load3body",taskcount,load
         last_polar3(taskcount)=countmol-1
         maxsize_mollst3(taskcount)=maxsize
         taskcount=taskcount+1
      end do


      numtasks_polar3=taskcount
c      print*,"numtasks_polar3=",numtasks_polar3
      return
      end

      subroutine bcast_polar_load_balance3b
      use mpidat
      use molcul
      implicit none
      include 'mpif.h'
      integer ierr
      logical first
      save first
      data first  / .true. /

      if (first.and.(taskid.ne.master)) then
         first = .false.
         if (.not. allocated(start_polar3))
     &          allocate (start_polar3(0:numtasks-1))
         if (.not. allocated(last_polar3))
     &          allocate (last_polar3(0:numtasks-1))
c         if (.not. allocated(molnew3))
c     &          allocate (molnew3(nmol))
         if (.not. allocated(maxsize_mollst3))
     &         allocate (maxsize_mollst3(0:numtasks-1))
      end if
c         print*,"numtasks_polar in bcast",numtasks_polar
         call mpi_bcast(numtasks_polar3,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(start_polar3,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(last_polar3,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(maxsize_mollst3,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
c      print*,"Successful 3BODY Load Bal Bcasts! Have a nice day!"
      return
      end  

      subroutine  sendmollst3_load_balanced_sizelist
      use sizes
      use mpidat
      use neigh3b
      use molcul
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt
      integer nmollst3_tot,nmol_nonzero,avg_load
      integer countpole,taskcount,offset_cum
      integer load,itercount,offset_polar3
      logical first,first2
      save first2
      data first2  / .true. /

      if(first2) then
         first2 =  .false.
         offset_polar3=last_polar3(taskid)-start_polar3(taskid)+1
         allocate (mollst3_recv(maxsize_mollst3(taskid),offset_polar3))
         allocate (nmollst3_recv(offset_polar3))
      else
         offset_polar3=last_polar3(taskid)-start_polar3(taskid)+1
          if(allocated(mollst3_recv)) deallocate(mollst3_recv)
          if(allocated(nmollst3_recv)) deallocate(nmollst3_recv)
          allocate (mollst3_recv(maxsize_mollst3(taskid),offset_polar3))
          allocate (nmollst3_recv(offset_polar3))
c           print*,"After recv allocs",taskid,offset_emreal3
      end if


         if(taskid.eq.master) then
c           print*,"numtasks_emreal",numtasks_emreal
           offset_cum=0
           do i=0,numtasks_polar3-1
             tag1=2*i+1
             tag2=2*i+2
             offset=last_polar3(i)-start_polar3(i)+1
c             print*,"BEFORE SEND i offset=",i,offset
c         print*,"i offset_emreal2=",i,offset_emreal2(i)
            if(i.eq.master) then
              cnt=0
c              do j=i*offset+1,(i+1)*offset
              do j=offset_cum+1,offset_cum+offset
                 cnt=cnt+1
                 nmollst3_recv(cnt)=nmollst3(j)
                 do k=1,maxsize_mollst3(i)
                    mollst3_recv(k,cnt)=mollst3(k,j)
                 end do
              end do
c              print*,"3BODYNOSEND:) taskid offset=",i,offset
            else

           call mpi_send(mollst3(1:maxsize_mollst3(i)
     &                        ,offset_cum+1:offset_cum+offset),
     &      maxsize_mollst3(i)*offset,mpi_integer,i,
     &        tag1,mpi_comm_world,ierr)
             call mpi_send(nmollst3(offset_cum+1:offset_cum+offset),
     &      offset,mpi_integer,i,tag2,mpi_comm_world,ierr)

c             print*,"AFTER 3BODYSEND taskid offset=",i,offset
            end if
            offset_cum=offset_cum+offset
           end do
         end if

         if(taskid.ne.master) then
           call mpi_recv(mollst3_recv,
     &                 maxsize_mollst3(taskid)*offset_polar3,
     &       mpi_integer,master,2*taskid+1,mpi_comm_world,stat,ierr)
           call mpi_recv(nmollst3_recv,offset_polar3,
     &       mpi_integer,master,2*taskid+2,mpi_comm_world,stat,ierr)
c            print*,"After 3BODYRECV"
         end if

      return
      end
      
