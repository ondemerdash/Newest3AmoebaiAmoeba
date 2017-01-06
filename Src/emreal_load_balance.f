      subroutine emreal_load_balance 
      use sizes
      use mpidat 
      use neigh
      use paremneigh
      use mpole
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt
      integer nelst_tot,npole_nonzero,avg_load
      integer countpole,taskcount
      integer load,itercount,offset_emreal2
      logical first,first2
      save first2
      data first2  / .true. /

             nelst_tot=0
             npole_nonzero=0
             print*,"numtasks_emreal=",numtasks_emreal
             do i=1,npole
               nelst_tot=nelst_tot+nelst(i) 
             end do 
        
             avg_load = nelst_tot/numtasks_emreal
             print*,"PERM ELEC avg_load",avg_load

             if(mod(nelst_tot,numtasks_emreal).gt.0) then
               avg_load=avg_load+1
              print*,"PERM ELEC avg_load rmndr",avg_load
             end if

             countpole=1
             taskcount=0

             do while (countpole.lt.npole)
                load=0
                itercount=0
                do while (load.lt.avg_load)
                   if(itercount.eq.0) then
                     start_emreal2(taskcount)=countpole
                   end if
                   load=load+nelst(countpole)
                   countpole=countpole+1
                   itercount=itercount+1

                   if(countpole.eq.npole+1) then
                     goto 31
                   end if
                end do
   31             continue
               print*,"PERM ELEC taskcount load",taskcount,load
               last_emreal2(taskcount)=countpole-1
c            offset_emreal2(taskcount)=last_emreal2(taskcount)
c     &        - start_emreal2(taskcount)+1
               taskcount=taskcount+1
             end do 
             print*,"PERM ELEC taskcount",taskcount
             do i=0,taskcount-1
              print*,"i start_emreal2",i,start_emreal2(i)
              print*,"i last_emreal2",i,last_emreal2(i)
             end do
      return
      end


      subroutine sendmlist_load_balanced3
      use sizes
      use mpidat
      use neigh
      use paremneigh
      use mpole
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt
      integer nelst_tot,npole_nonzero,avg_load
      integer countpole,taskcount,offset_cum
      integer load,itercount,offset_emreal2
      logical first,first2
      save first2
      data first2  / .true. /

c      print*,"taskid start_emreal2(taskid)",taskid,start_emreal2(taskid)
c      print*,"taskid last_emreal2(taskid)",taskid,last_emreal2(taskid)

         if(first2) then
           first2= .false.
           offset_emreal3=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           allocate (elst_recv(maxelst,offset_emreal3))
           allocate (nelst_recv(offset_emreal3))
           print*,"After first recv allocs",taskid,offset_emreal3
         else
           offset_emreal3=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           if(allocated(elst_recv)) deallocate(elst_recv)
           if(allocated(nelst_recv)) deallocate(nelst_recv)
           allocate (elst_recv(maxelst,offset_emreal3))
           allocate (nelst_recv(offset_emreal3))
           print*,"After recv allocs",taskid,offset_emreal3
         end if 

c         call mpi_barrier(mpi_comm_world,ierr)

         if(taskid.eq.master1) then
           print*,"numtasks_emreal",numtasks_emreal
           offset_cum=0
           do i=0,numtasks_emreal-1
             tag1=2*i+1
             tag2=2*i+2
             offset=last_emreal2(i)-start_emreal2(i)+1
             print*,"BEFORE SEND i offset=",i,offset
c         print*,"i offset_emreal2=",i,offset_emreal2(i)
            if(i.eq.master1) then
              cnt=0
c              do j=i*offset+1,(i+1)*offset
              do j=offset_cum+1,offset_cum+offset
                 cnt=cnt+1
                 nelst_recv(cnt)=nelst(j)
                 do k=1,maxelst
                    elst_recv(k,cnt)=elst(k,j)
                 end do 
              end do 
              print*,"NOSEND:) taskid offset=",i,offset
            else
c             call mpi_send(elst(1:maxelst,i*offset+1:(i+1)*offset),
c     &      maxelst*offset,mpi_integer,i,tag1,emreal_comm,ierr)
c             call mpi_send(nelst(i*offset+1:(i+1)*offset),
c     &      offset,mpi_integer,i,tag2,emreal_comm,ierr)

           call mpi_send(elst(1:maxelst,offset_cum+1:offset_cum+offset),
     &      maxelst*offset,mpi_integer,i,tag1,emreal_comm,ierr)
             call mpi_send(nelst(offset_cum+1:offset_cum+offset),
     &      offset,mpi_integer,i,tag2,emreal_comm,ierr)

             print*,"AFTER SEND taskid offset=",i,offset
            end if
            offset_cum=offset_cum+offset
           end do
         end if

         if(taskid.ne.master1) then
            
           call mpi_recv(elst_recv,maxelst*offset_emreal3,
     &       mpi_integer,master1,2*taskid+1,emreal_comm,stat,ierr) 
           call mpi_recv(nelst_recv,offset_emreal3,
     &       mpi_integer,master1,2*taskid+2,emreal_comm,stat,ierr)
           print*,"AFTER RECV taskid offset=",taskid,offset_emreal3
         end if

      return
      end
