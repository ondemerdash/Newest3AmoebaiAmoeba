
      subroutine vdw_load_balance_sizelist_half 
      use sizes
      use mpidat 
      use neigh
      use parvdwneigh
      use vdw
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt
      integer nvlst_tot,nvdw_nonzero,avg_load
      integer countvdw,taskcount
      integer load,itercount,offset_vdw2
      integer maxsize
      logical first,first2       
      save first2
      data first2  / .true. /

             nvlst_tot=0
             nvdw_nonzero=0
c             print*,"numtasks_emreal=",numtasks_emreal
             do i=1,nvdw
               nvlst_tot=nvlst_tot+nvlst(i) 
             end do 
        
             avg_load = nvlst_tot/numtasks_vdw
c             print*,"VDW avg_load",avg_load

             if(mod(nvlst_tot,numtasks_vdw).gt.0) then
               avg_load=avg_load+1
c              print*,"VDW avg_load rmndr",avg_load
             end if

             countvdw=1
             taskcount=0

             do while (countvdw.lt.nvdw)
                load=0
                itercount=0
                maxsize=0
                do while (load.lt.avg_load)
                   if(itercount.eq.0) then
                     start_vdw2(taskcount)=countvdw
                     maxsize=nvlst(countvdw)
                   else
                     if(nvlst(countvdw).gt.maxsize) then
                       maxsize=nvlst(countvdw)
                     end if
                   end if
                   load=load+nvlst(countvdw)
                   
                   countvdw=countvdw+1
                   itercount=itercount+1

                   if(countvdw.eq.nvdw+1) then
                     goto 31
                   end if
                end do
   31             continue
c               print*,"VDW taskcount load",taskcount,load
               last_vdw2(taskcount)=countvdw-1
               maxsize_vlst(taskcount)=maxsize
c            offset_emreal2(taskcount)=last_emreal2(taskcount)
c     &        - start_emreal2(taskcount)+1
               taskcount=taskcount+1
             end do 
c             print*,"VDW Total taskcount",taskcount

             numtasks_vdw2=taskcount
c             do i=0,taskcount-1
c              print*,"taskid start_vdw2",i,start_vdw2(i)
c              print*,"taskid last_vdw2",i,last_vdw2(i)
c             end do
      return
      end


      subroutine sendvlist_load_balanced3_sizelist_half
      use sizes
      use mpidat
      use neigh
      use parvdwneigh
      use vdw
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt
      integer nvlst_tot,nvdw_nonzero,avg_load
      integer countvdw,taskcount,offset_cum
      integer load,itercount,offset_vdw2
      logical first,first2
      real*8 t1,t2,t3,t4
      save first2
      data first2  / .true. /

c         if(taskid.lt.numtasks_vdw2) then
         if((taskid.ge.numtasks_emreal)
     &   .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal))) then
           offset_vdw3=last_vdw2(taskid-numtasks_emreal)
     &        - start_vdw2(taskid-numtasks_emreal)+1
           if(allocated(vlst_recv)) deallocate(vlst_recv)
           if(allocated(nvlst_recv)) deallocate(nvlst_recv)
           allocate (vlst_recv(maxsize_vlst(taskid-numtasks_emreal)
     &             ,offset_vdw3))
           allocate (nvlst_recv(offset_vdw3))
         end if
c         call mpi_barrier(mpi_comm_world,ierr)

         if(taskid.eq.master) then
c           print*,"numtasks_emreal",numtasks_emreal
             t3=mpi_Wtime()
           offset_cum=0
           do i=0,numtasks_vdw2-1
             tag1=2*(i+numtasks_emreal)+1
             tag2=2*(i+numtasks_emreal)+2
             offset=last_vdw2(i)-start_vdw2(i)+1
c            if(i.eq.master) then
c              cnt=0
c              do j=offset_cum+1,offset_cum+offset
c                 cnt=cnt+1
c                 nvlst_recv(cnt)=nvlst(j)
c                 do k=1,maxsize_vlst(i)
c                    vlst_recv(k,cnt)=vlst(k,j)
c                 end do 
c              end do 
c            else

           call mpi_send(vlst(1:maxsize_vlst(i)
     &                        ,offset_cum+1:offset_cum+offset),
     &      maxsize_vlst(i)*offset,mpi_integer,i+numtasks_emreal,tag1,
     &                  mpi_comm_world,ierr)
             call mpi_send(nvlst(offset_cum+1:offset_cum+offset),
     &      offset,mpi_integer,i+numtasks_emreal,tag2,
     &                  mpi_comm_world,ierr)

c             print*,"AFTER VDW SEND taskid offset=",i,offset
c            end if
            offset_cum=offset_cum+offset
           end do
             t4=mpi_Wtime()
        print*,"VDWSEND Timing taskid",taskid,"time=",t4-t3

         end if

      if((taskid.ge.numtasks_emreal)
     &    .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal)))then
             t1=mpi_Wtime()
            
      call mpi_recv(vlst_recv,
     &       maxsize_vlst(taskid-numtasks_emreal)*offset_vdw3,
     &       mpi_integer,master,2*taskid+1,mpi_comm_world,stat,ierr)
      call mpi_recv(nvlst_recv,offset_vdw3,
     &       mpi_integer,master,2*taskid+2,mpi_comm_world,stat,ierr)
c           print*,"AFTER VDW RECV taskid offset=",taskid,offset_vdw3
             t2=mpi_Wtime()
        print*,"VDWRECV Timing taskid",taskid,"time=",t2-t1

      end if

      return
      end
