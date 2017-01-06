      subroutine sendmlist_load_balanced3_sizelist_serial
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
         !if(taskid.lt.numtasks_emreal2) then
         if(first2) then
           first2= .false.
           offset_emreal3=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           allocate (elst_recv(maxsize_elst(taskid),offset_emreal3))
           allocate (nelst_recv(offset_emreal3))
c           print*,"After first recv allocs",taskid,offset_emreal3
         else
           offset_emreal3=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           if(allocated(elst_recv)) deallocate(elst_recv)
           if(allocated(nelst_recv)) deallocate(nelst_recv)
           allocate (elst_recv(maxsize_elst(taskid),offset_emreal3))
           allocate (nelst_recv(offset_emreal3))
c           print*,"After recv allocs",taskid,offset_emreal3
         end if
         !end if

         !if(taskid.eq.master1) then
c           print*,"numtasks_emreal",numtasks_emreal
           offset_cum=0
           do i=0,numtasks_emreal2-1
             tag1=2*i+1
             tag2=2*i+2
             offset=last_emreal2(i)-start_emreal2(i)+1
c             print*,"BEFORE SEND i offset=",i,offset
c         print*,"i offset_emreal2=",i,offset_emreal2(i)
            !if(i.eq.master1) then
            if(i.eq.master) then
              cnt=0
              do j=offset_cum+1,offset_cum+offset
                 cnt=cnt+1
                 nelst_recv(cnt)=nelst(j)
                 do k=1,maxsize_elst(i)
                    elst_recv(k,cnt)=elst(k,j)
                 end do
              end do
c              print*,"NOSEND:) taskid offset=",i,offset
            else

           !call mpi_send(elst(1:maxsize_elst(i)
     &     !                   ,offset_cum+1:offset_cum+offset),
     &     ! maxsize_elst(i)*offset,mpi_integer,i,tag1,emreal_comm,ierr)
           !  call mpi_send(nelst(offset_cum+1:offset_cum+offset),
     &     ! offset,mpi_integer,i,tag2,emreal_comm,ierr)

c             print*,"AFTER SEND taskid offset=",i,offset
            end if
            offset_cum=offset_cum+offset
           end do
         !end if

         !if((taskid.ne.master1).and.(taskid.lt.numtasks_emreal2)) then

         !  call mpi_recv(elst_recv,maxsize_elst(taskid)*offset_emreal3,
     &   !    mpi_integer,master1,2*taskid+1,emreal_comm,stat,ierr)
         !  call mpi_recv(nelst_recv,offset_emreal3,
     &   !    mpi_integer,master1,2*taskid+2,emreal_comm,stat,ierr)
c           print*,"AFTER RECV taskid offset=",taskid,offset_emreal3
         !end if

         return
         end
