      subroutine polar_load_balance_clust_splitlist_offset
      use mpidat
      use neigh3b
      use molcul
      use neigh2clust
      use aprx
      implicit none
      include 'mpif.h'
      integer i,j,nmollst_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load,i1,counter,k1
      logical first
      integer numtasks_offset
      save first
      data first  / .true. /
            
      nmollst_tot=0
      nmol_nonzero=0

      
      if (first) then
         if (.not. allocated(molnew))
     &          allocate (molnew(clustcount))
      end if      

      max_nmollst=0
      do i=1,clustcount
        if(nmollst(i).ne.0) then
          if(nmollst(i).gt.max_nmollst) then 
              max_nmollst=nmollst(i)   
          end if
          do j=1,nmollst(i)
             k=mollst(j,i)
             nmollst_tot = nmollst_tot+sizeclust(k)+sizeclust(i)
          end do

          nmol_nonzero=nmol_nonzero+1
          molnew(nmol_nonzero)=i
        !  print*,"cluster_index, nmollst(i)",i,nmollst(i)
        end if
      end do

      if(approxmode.eq.'3BODYMODE') then
        numtasks_offset=numtasks-clustcount-cnt3bwork
      else 
        numtasks_offset=numtasks-clustcount
      end if

      if (first) then
         first = .false.
        if (.not. allocated(start_polar))
     &          allocate (start_polar(0:numtasks-1))

        if (.not. allocated(start_mollst_polar))
     &          allocate (start_mollst_polar(max_nmollst,clustcount,
     &          0:numtasks-1))

        if (.not. allocated(last_polar))
     &          allocate (last_polar(0:numtasks-1))

        if (.not. allocated(last_mollst_polar))
     &          allocate (last_mollst_polar(max_nmollst,clustcount,
     &          0:numtasks-1))
        if (.not. allocated(num_mollst_chunk))
     &       allocate (num_mollst_chunk(clustcount,0:numtasks-1))
      end if


c      do i=1,nmol_nonzero
c       print*,"molnew",molnew(i),"nmollst(molnew(i))",nmollst(molnew(i))
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
c            load=load+nmollst(k)          
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
         i=molnew(i1)
         counter=0
           do j=1,nmollst(i)
              k=mollst(j,i)              
              if(load.eq.0) then
                start_polar(taskcount)=i1
                counter=counter+1
                start_mollst_polar(counter,i,taskcount)=j
                num_mollst_chunk(i,taskcount)=counter
              else if((load.ne.0).and.(j.eq.1)) then
                counter=counter+1
                start_mollst_polar(counter,i,taskcount)=j
                num_mollst_chunk(i,taskcount)=counter
              end if

              load=load+sizeclust(k)+sizeclust(i)
              if(load.ge.avg_load) then
                last_polar(taskcount)=i1
                last_mollst_polar(counter,i,taskcount)=j
                taskcount=taskcount+1
                counter=0
              end if 
              if((load.lt.avg_load).and.(j.eq.nmollst(i))) then
                last_mollst_polar(counter,i,taskcount)=j
              end if

c              if( (load.ge.avg_load) .or. 
c     &           ((load.lt.avg_load).and.(j.eq.nmollst(i))) ) then 
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
      !do i=0,taskcount-1
      !   do k1=start_polar(i), last_polar(i)
      !      j=molnew(k1)
      !   print*,"task=",i,"clust=",j,"numchunk=",num_mollst_chunk(j,i)
      !   end do
      !end do

      !do i=0,taskcount-1
      !print*,"taskid molnewstart,molnewend",i,molnew(start_polar(i)),
     &!           molnew(last_polar(i))
      !   do k1=start_polar(i), last_polar(i)
      !      j=molnew(k1)
      !      do k=1,num_mollst_chunk(j,i)
      !         print*,"tskid clust start_mollst last_mollst",i,j,
     &!          start_mollst_polar(k,j,i),last_mollst_polar(k,j,i)
      !      end do 
      !   end do
      !end do

      numtasks_polar=taskcount
      print*,"Listloadbal: numtasks_polar=",numtasks_polar
      return
      end

