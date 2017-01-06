      subroutine sendmlist
      use sizes
      use mpidat 
      use neigh
      use paremneigh
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt

      offset=offset_emreal
      if(taskid.lt.numtasks_emreal) then
        if (.not. allocated(elst_recv)) 
     &      allocate (elst_recv(maxelst,offset_emreal))
        if (.not. allocated(nelst_recv)) 
     &      allocate (nelst_recv(offset_emreal))
        if(taskid.lt.remainder_emreal) then
          if (.not. allocated(elst_recv_single)) 
     &      allocate(elst_recv_single(maxelst,1))
c          allocate(nelst_recv_single(1))
        end if
      end if

      tot=numtasks_emreal*offset

      if(taskid.eq.master1) then
        do i=0,numtasks_emreal-1
          tag1=2*i+1
          tag2=2*i+2
         if(i.eq.master1) then
           cnt=0
           do j=i*offset+1,(i+1)*offset
              cnt=cnt+1
              nelst_recv(cnt)=nelst(j)
              do k=1,maxelst
                 elst_recv(k,cnt)=elst(k,j)
              end do 
           end do 
         else
          call mpi_send(elst(1:maxelst,i*offset+1:(i+1)*offset),
     &    maxelst*offset,mpi_integer,i,tag1,emreal_comm,ierr)
          call mpi_send(nelst(i*offset+1:(i+1)*offset),
     &    offset,mpi_integer,i,tag2,emreal_comm,ierr)
         end if
        end do
       if(remainder_emreal.gt.0) then
        do i=0,remainder_emreal-1
          tag1=2*numtasks_emreal+2*i+1
          tag2=2*numtasks_emreal+2*i+2 
          if(i.eq.master1) then
            nelst_recv_single=nelst(tot+i+1)
            do k=1,maxelst
              elst_recv_single(k,1)=elst(k,tot+i+1)
            end do
          else
           call mpi_send(elst(1:maxelst,tot+i+1:tot+i+1),
     &     maxelst,mpi_integer,i,tag1,emreal_comm,ierr)
           call mpi_send(nelst(tot+i+1:tot+i+1),
     &     1,mpi_integer,i,tag2,emreal_comm,ierr)
          end if        
        end do
       end if
      end if

      if(taskid.lt.numtasks_emreal.and.taskid.ne.master1) then
        call mpi_recv(elst_recv,maxelst*offset,mpi_integer,master1,
     &       2*taskid+1,emreal_comm,stat,ierr) 
        call mpi_recv(nelst_recv,offset,mpi_integer,master1,
     &       2*taskid+2,emreal_comm,stat,ierr)
        if(taskid.lt.remainder_emreal) then
          call mpi_recv(elst_recv_single,maxelst,mpi_integer,master1,
     &      2*numtasks_emreal+2*taskid+1,emreal_comm,stat,ierr)
          call mpi_recv(nelst_recv_single,1,mpi_integer,master1,
     &      2*numtasks_emreal+2*taskid+2,emreal_comm,stat,ierr)
        end if
        print*,"After recvs, offsets",taskid,offset
      end if
      
      return
      end
