      subroutine emreal_sizelist_only
      use mpidat
      use neigh
      implicit none
      integer i,j
     
c      print*,"BEGIN emreal_sizelist_only"
c      do i=0,numtasks_emreal-1
      do i=0,numtasks_emreal2-1
         do j=start_emreal2(i),last_emreal2(i)
            if(nelst(j).gt.maxsize_elst(i)) then
              maxsize_elst(i)=nelst(j)
            end if            
         end do
      end do
c      print*,"END emreal_sizelist_only"
      return
      end

      subroutine parallel_emreal_sizelist_only
      use mpidat
      use neigh
      implicit none
      integer i,j
      print*,"BEGIN parallel_emreal_sizelist_only"
!$OMP PARALLEL DO default(private) shared(numtasks_emreal,
!$OMP& start_emreal2,last_emreal2,nelst,maxsize_elst)
      do i=0,numtasks_emreal2-1
         do j=start_emreal2(i),last_emreal2(i)
            if(nelst(j).gt.maxsize_elst(i)) then
              maxsize_elst(i)=nelst(j)
            end if
         end do
      end do
!$OMP END PARALLEL DO
      print*,"END parallel_emreal_sizelist_only"
      return
      end


      subroutine emreal_load_balance_sizelist 
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
             !print*,"numtasks_emreal=",numtasks_emreal
             do i=1,npole
               nelst_tot=nelst_tot+nelst(i) 
             end do 
        
             avg_load = nelst_tot/numtasks_emreal
             !print*,"PERM ELEC avg_load",avg_load

             if(mod(nelst_tot,numtasks_emreal).gt.0) then
               avg_load=avg_load+1
            !  print*,"PERM ELEC avg_load rmndr",avg_load
             end if

             countpole=1
             taskcount=0

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
              ! print*,"PERM ELEC taskcount load",taskcount,load
               last_emreal2(taskcount)=countpole-1
               maxsize_elst(taskcount)=maxsize
c            offset_emreal2(taskcount)=last_emreal2(taskcount)
c     &        - start_emreal2(taskcount)+1
               taskcount=taskcount+1
             end do 
          !   print*,"PERM ELEC taskcount",taskcount

             numtasks_emreal2=taskcount
           !  do i=0,taskcount-1
           !   print*,"taskid start_emreal2",i,start_emreal2(i)
           !   print*,"taskid last_emreal2",i,last_emreal2(i)
           !   print*,"taskid maxsize_elst",i,maxsize_elst(i)
           !  end do
c             do i=1,npole
c                print*,"i",i,"nelst(i)",nelst(i)
c             end do
      return
      end


      subroutine sendmlist_load_balanced3_sizelist
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
      real*8 t1,t2,t3,t4
      save first2
      data first2  / .true. /

c      print*,"taskid start_emreal2(taskid)",taskid,start_emreal2(taskid)
c      print*,"taskid last_emreal2(taskid)",taskid,last_emreal2(taskid)
       print*,"In sendmlist routine numtasks_emreal2=",numtasks_emreal2

         if(taskid.lt.numtasks_emreal2) then
         !if(first2) then
         !  first2= .false.
         !  offset_emreal3=last_emreal2(taskid)
     &   !     - start_emreal2(taskid)+1
         !  allocate (elst_recv(maxsize_elst(taskid),offset_emreal3))
         !  allocate (nelst_recv(offset_emreal3))
           !print*,"After first recv allocs",taskid,offset_emreal3
         !else
           offset_emreal3=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           if(allocated(elst_recv)) deallocate(elst_recv)
           if(allocated(nelst_recv)) deallocate(nelst_recv)
           allocate (elst_recv(maxsize_elst(taskid),offset_emreal3))
           allocate (nelst_recv(offset_emreal3))
           !print*,"After recv allocs",taskid,offset_emreal3
         !end if 
         end if
c         call mpi_barrier(mpi_comm_world,ierr)

         if(taskid.eq.master1) then
          ! print*,"numtasks_emreal2",numtasks_emreal2
           offset_cum=0
          !   t3=mpi_Wtime()
           do i=0,numtasks_emreal2-1
             tag1=2*i+1
             tag2=2*i+2
             offset=last_emreal2(i)-start_emreal2(i)+1
c             print*,"BEFORE SEND i offset=",i,offset
c         print*,"i offset_emreal2=",i,offset_emreal2(i)
            if(i.eq.master1) then
              cnt=0
              do j=offset_cum+1,offset_cum+offset
                 cnt=cnt+1
                 nelst_recv(cnt)=nelst(j)
c              print*,"j",j,"cnt=",cnt,"nelst_recv(cnt)",nelst_recv(cnt)
                 do k=1,maxsize_elst(i)
                    elst_recv(k,cnt)=elst(k,j)
                 end do 

c                 do k=1,nelst(j)
c                    elst_recv(k,cnt)=elst(k,j)
c                 end do 
c                 do k=nelst(j)+1,maxsize_elst(i)
c                    elst_recv(k,cnt)=0
c                 end do
              end do 
c              print*,"MPOLE NOSEND:) taskid offset=",i,offset
            else

           call mpi_send(elst(1:maxsize_elst(i)
     &                        ,offset_cum+1:offset_cum+offset),
     &      maxsize_elst(i)*offset,mpi_integer,i,tag1,
     &                  mpi_comm_world,ierr)
             call mpi_send(nelst(offset_cum+1:offset_cum+offset),
     &      offset,mpi_integer,i,tag2,mpi_comm_world,ierr)

c             print*,"AFTER MPOLE SEND taskid offset=",i,offset
            end if
            offset_cum=offset_cum+offset
           end do
         !    t4=mpi_Wtime()
        !print*,"MPOLESEND Timing taskid",taskid,"time=",t4-t3
         end if

c         call mpi_barrier(mpi_comm_world,ierr)

         if((taskid.ne.master1).and.(taskid.lt.numtasks_emreal2)) then
         !    t1=mpi_Wtime()
           call mpi_recv(elst_recv,maxsize_elst(taskid)*offset_emreal3,
     &       mpi_integer,master1,2*taskid+1,mpi_comm_world,stat,ierr)
           call mpi_recv(nelst_recv,offset_emreal3,
     &       mpi_integer,master1,2*taskid+2,mpi_comm_world,stat,ierr)
c           print*,"AFTER MPOLE RECV taskid offset",taskid,offset_emreal3
c           do i=start_emreal2(taskid),last_emreal2(taskid)
c             print*,"i",i,"nelst_recv(i)",nelst_recv(i)
c           end do
c         else if(taskid.eq.master1) then
c           do i=start_emreal2(taskid),last_emreal2(taskid)
c              j=i-start_emreal2(taskid)+1
c             print*,"i",i,"j",j,"nelst_recv(i)",nelst_recv(j)
c           end do
        !     t2=mpi_Wtime()
       ! print*,"MPOLERECV Timing taskid",taskid,"time=",t2-t1
         end if

      return
      end

      subroutine sendmlist_load_balanced3_sizelist_master
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
      real*8 t1,t2,t3,t4
      save first2
      data first2  / .true. /

c      print*,"taskid start_emreal2(taskid)",taskid,start_emreal2(taskid)
c      print*,"taskid last_emreal2(taskid)",taskid,last_emreal2(taskid)

c         if(first2) then
c           first2= .false.
c           offset_emreal3=last_emreal2(taskid)
c     &        - start_emreal2(taskid)+1
c           allocate (elst_recv(maxsize_elst(taskid),offset_emreal3))
c           allocate (nelst_recv(offset_emreal3))
c           print*,"After first recv allocs",taskid,offset_emreal3
c         else
           offset_emreal3=last_emreal2(taskid)
     &        - start_emreal2(taskid)+1
           if(allocated(elst_recv)) deallocate(elst_recv)
           if(allocated(nelst_recv)) deallocate(nelst_recv)
           allocate (elst_recv(maxsize_elst(taskid),offset_emreal3))
           allocate (nelst_recv(offset_emreal3))
c           print*,"After recv allocs",taskid,offset_emreal3
c         end if

c         call mpi_barrier(mpi_comm_world,ierr)

         if(taskid.eq.master) then
c           print*,"numtasks_emreal",numtasks_emreal
           offset_cum=0
             t3=mpi_Wtime()
           do i=0,numtasks_emreal-1
             tag1=2*i+1
             tag2=2*i+2
             offset=last_emreal2(i)-start_emreal2(i)+1
c             print*,"BEFORE SEND i offset=",i,offset
c         print*,"i offset_emreal2=",i,offset_emreal2(i)
            if(i.eq.master) then
              cnt=0
c              do j=i*offset+1,(i+1)*offset
              do j=offset_cum+1,offset_cum+offset
                 cnt=cnt+1
                 nelst_recv(cnt)=nelst(j)
                 do k=1,maxsize_elst(i)
                    elst_recv(k,cnt)=elst(k,j)
                 end do
              end do
c              print*,"NOSEND:) taskid offset=",i,offset
            else

           call mpi_send(elst(1:maxsize_elst(i)
     &                        ,offset_cum+1:offset_cum+offset),
     &      maxsize_elst(i)*offset,mpi_integer,i,tag1,
     &               mpi_comm_world,ierr)
             call mpi_send(nelst(offset_cum+1:offset_cum+offset),
     &      offset,mpi_integer,i,tag2,mpi_comm_world,ierr)

c             print*,"AFTER SEND taskid offset=",i,offset
            end if
            offset_cum=offset_cum+offset
           end do
             t4=mpi_Wtime()
        print*,"MPOLESEND Timing taskid",taskid,"time=",t4-t3
         end if

         if(taskid.ne.master) then
           t1=mpi_Wtime()
           call mpi_recv(elst_recv,maxsize_elst(taskid)*offset_emreal3,
     &       mpi_integer,master,2*taskid+1,mpi_comm_world,stat,ierr)
           call mpi_recv(nelst_recv,offset_emreal3,
     &       mpi_integer,master,2*taskid+2,mpi_comm_world,stat,ierr)
c           print*,"AFTER RECV taskid offset=",taskid,offset_emreal3
             t2=mpi_Wtime()
        print*,"MPOLERECV Timing taskid",taskid,"time=",t2-t1
             
         end if

      return
      end

      subroutine emreal_load_balance_sizelist_mol
      use sizes
      use mpidat
      use neigh
      use paremneigh
      use mpole
      use molcul
      implicit none
      include 'mpif.h'
      integer ierr,i,offset,cnti,j,k,tot,tag1,tag2
      integer stat(MPI_STATUS_SIZE),cnt
      integer nelst_tot,npole_nonzero,avg_load
      integer countpole,taskcount
      integer load,itercount,offset_emreal2
      integer maxsize
      logical first,first2
      integer lenmol,countmol
      save first2
      data first2  / .true. /

             nelst_tot=0
             npole_nonzero=0
             !print*,"numtasks_emreal=",numtasks_emreal
             do i=1,npole
               nelst_tot=nelst_tot+nelst(i)
             end do

             avg_load = nelst_tot/numtasks_emreal
             !print*,"PERM ELEC avg_load",avg_load

             if(mod(nelst_tot,numtasks_emreal).gt.0) then
               avg_load=avg_load+1
            !  print*,"PERM ELEC avg_load rmndr",avg_load
             end if

             countpole=1
             countmol=1
             taskcount=0

c             do while (countpole.lt.npole)
             do while (countmol.lt.nmol)
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

                   
                   lenmol=imol(2,countmol)-imol(1,countmol)+1

                   do k=1,lenmol
                    load=load+nelst(countpole)
                    countpole = countpole + 1
                   end do

                   countmol=countmol+1
                   itercount=itercount+1

                   if(countpole.eq.npole+1) then
                     goto 31
                   end if
                end do
   31             continue
              ! print*,"PERM ELEC taskcount load",taskcount,load
               last_emreal2(taskcount)=countpole-1
               maxsize_elst(taskcount)=maxsize
c            offset_emreal2(taskcount)=last_emreal2(taskcount)
c     &        - start_emreal2(taskcount)+1
               taskcount=taskcount+1
             end do
           print*,"PERM ELECMol taskcount",taskcount

             numtasks_emreal2=taskcount
             do i=0,taskcount-1
            print*,"taskid mol start_emreal2",i,start_emreal2(i)
            print*,"taskid mol last_emreal2",i,last_emreal2(i)
            ! print*,"taskid maxsize_elst",i,maxsize_elst(i)
c             end do
c             do i=1,npole
c                print*,"i",i,"nelst(i)",nelst(i)
             end do
      return
      end

