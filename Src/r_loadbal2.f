      subroutine emreal_load_balance_sizelist2 
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
      integer maxsize
      logical first,first2       
      save first2
      data first2  / .true. /

             nelst_tot=0
             npole_nonzero=0
c             print*,"numtasks_emreal=",numtasks_emreal
             do i=1,npole
               nelst_tot=nelst_tot+nelst(i) 
             end do 
        
             avg_load = nelst_tot/numtasks_emreal
c             print*,"PERM ELEC avg_load",avg_load

             if(mod(nelst_tot,numtasks_emreal).gt.0) then
               avg_load=avg_load+1
c              print*,"PERM ELEC avg_load rmndr",avg_load
             end if

             countpole=1
             taskcount=0
             maxsize_elst2=0
             do while (countpole.lt.npole)
                load=0
                itercount=0
                maxsize=0
                do while (load.lt.avg_load)
                   if(itercount.eq.0) then
                     start_emreal2(taskcount)=countpole
                     maxsize=nelst(countpole)
                   else
                     if(nelst(countpole).gt.maxsize) then
                       maxsize=nelst(countpole)
                     end if
                   end if
                   load=load+nelst(countpole)
                   
                   countpole=countpole+1
                   itercount=itercount+1

                   if(countpole.eq.npole+1) then
                     goto 31
                   end if
                end do
   31             continue
c               print*,"PERM ELEC taskcount load",taskcount,load
               last_emreal2(taskcount)=countpole-1
               !maxsize_elst(taskcount)=maxsize
               if(maxsize.gt.maxsize_elst2) then
                 maxsize_elst2 = maxsize
c                print*,"maxsize=",maxsize,"maxsize_elst2=",maxsize_elst2
               end if 
               taskcount=taskcount+1
             end do 
c             print*,"PERM ELEC taskcount",taskcount
c             print*,"Final maxsize_elst2=",maxsize_elst2
             numtasks_emreal2=taskcount
             do i=0,taskcount-1
c              print*,"taskid start_emreal2",i,start_emreal2(i)
c              print*,"taskid last_emreal2",i,last_emreal2(i)
             end do
      return
      end

      subroutine emreal_sizelist_only2
      use mpidat
      use neigh
      implicit none
      integer i,j

      do i=0,numtasks_emreal-1
         do j=start_emreal2(i),last_emreal2(i)
            if(nelst(j).gt.maxsize_elst2) then
              maxsize_elst2=nelst(j)
            end if
         end do
      end do

      return
      end

