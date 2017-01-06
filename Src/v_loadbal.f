      subroutine vdw_sizelist_only
      use mpidat
      use neigh
      implicit none
      integer i,j
     
c      print*,"BEGIN emreal_sizelist_only"
c      do i=0,numtasks_emreal-1
      do i=0,numtasks_vdw2-1
         do j=start_vdw2(i),last_vdw2(i)
            if(nvlst(j).gt.maxsize_vlst(i)) then
              maxsize_vlst(i)=nvlst(j)
            end if            
         end do
      end do
c      print*,"END emreal_sizelist_only"
      return
      end

      subroutine parallel_vdw_sizelist_only
      use mpidat
      use neigh
      implicit none
      integer i,j
      print*,"BEGIN parallel_vdw_sizelist_only"
!$OMP PARALLEL DO default(private) shared(numtasks_vdw2,
!$OMP& start_vdw2,last_vdw2,nvlst,maxsize_vlst)
      do i=0,numtasks_vdw2-1
         do j=start_vdw2(i),last_vdw2(i)
            if(nvlst(j).gt.maxsize_vlst(i)) then
              maxsize_vlst(i)=nvlst(j)
            end if
         end do
      end do
!$OMP END PARALLEL DO
      print*,"END parallel_vdw_sizelist_only"
      return
      end


      subroutine vdw_load_balance_sizelist 
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
        
             avg_load = nvlst_tot/numtasks
c             print*,"VDW avg_load",avg_load

             if(mod(nvlst_tot,numtasks).gt.0) then
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
c             print*,"VDW taskcount",taskcount

             numtasks_vdw2=taskcount
c             do i=0,taskcount-1
c              print*,"taskid start_vdw2",i,start_vdw2(i)
c              print*,"taskid last_vdw2",i,last_vdw2(i)
c             end do
      return
      end


      subroutine sendvlist_load_balanced3_sizelist
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
      save first2
      data first2  / .true. /

         if(taskid.lt.numtasks_vdw2) then
         if(first2) then
           first2= .false.
           offset_vdw3=last_vdw2(taskid)
     &        - start_vdw2(taskid)+1
           allocate (vlst_recv(maxsize_vlst(taskid),offset_vdw3))
           allocate (nvlst_recv(offset_vdw3))
           !print*,"After first recv allocs",taskid,offset_emreal3
         else
           offset_vdw3=last_vdw2(taskid)
     &        - start_vdw2(taskid)+1
           if(allocated(vlst_recv)) deallocate(vlst_recv)
           if(allocated(nvlst_recv)) deallocate(nvlst_recv)
           allocate (vlst_recv(maxsize_vlst(taskid),offset_vdw3))
           allocate (nvlst_recv(offset_vdw3))
         end if 
         end if
c         call mpi_barrier(mpi_comm_world,ierr)

         if(taskid.eq.master) then
c           print*,"numtasks_emreal",numtasks_emreal
           offset_cum=0
           do i=0,numtasks_vdw2-1
             tag1=2*i+1
             tag2=2*i+2
             offset=last_vdw2(i)-start_vdw2(i)+1
            if(i.eq.master) then
              cnt=0
              do j=offset_cum+1,offset_cum+offset
                 cnt=cnt+1
                 nvlst_recv(cnt)=nvlst(j)
                 do k=1,maxsize_vlst(i)
                    vlst_recv(k,cnt)=vlst(k,j)
                 end do 
              end do 
            else

           call mpi_send(vlst(1:maxsize_vlst(i)
     &                        ,offset_cum+1:offset_cum+offset),
     &      maxsize_vlst(i)*offset,mpi_integer,i,tag1,
     &                  mpi_comm_world,ierr)
             call mpi_send(nvlst(offset_cum+1:offset_cum+offset),
     &      offset,mpi_integer,i,tag2,mpi_comm_world,ierr)

c             print*,"AFTER VDW SEND taskid offset=",i,offset
            end if
            offset_cum=offset_cum+offset
           end do
         end if

         if((taskid.ne.master).and.(taskid.lt.numtasks_vdw2)) then
            
           call mpi_recv(vlst_recv,maxsize_vlst(taskid)*offset_vdw3,
     &       mpi_integer,master,2*taskid+1,mpi_comm_world,stat,ierr)
           call mpi_recv(nvlst_recv,offset_vdw3,
     &       mpi_integer,master,2*taskid+2,mpi_comm_world,stat,ierr)
c           print*,"AFTER VDW RECV taskid offset=",taskid,offset_vdw3
         end if

      return
      end
