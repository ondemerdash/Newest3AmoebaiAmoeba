      subroutine polar_load_balance6
      use mpidat
      use neigh3b
      use molcul
      implicit none
      include 'mpif.h'
      integer i,j,nmollst_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load, new_load, count
      integer subtasks, add, part
      real scal
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
     &          allocate (molnew(nmol))
      end if      

      do i=1,nmol
        if(nmollst(i).ne.0) then
          nmollst_tot = nmollst_tot + nmollst(i)
          nmol_nonzero=nmol_nonzero+1
          molnew(nmol_nonzero)=i
c          print*,"Nonzero nmollst",i,nmollst(i)
        end if
      end do

c      print*,"nmol_nonzero",nmol_nonzero

      avg_load=nmollst_tot/numtasks
c      print*,"avg_load=",avg_load

      if(mod(nmollst_tot,numtasks).gt.0) then
        avg_load=avg_load+1
c        print*,"avg_load rmndr",avg_load
      end if


      countmol=1
      taskcount=0
      subtasks = 3.0/4.0 * numtasks +1
      count = subtasks +1
      scal = 0
      new_load = 0 
      part = 0 
      add = -2
      
      do while (countmol.lt.nmol_nonzero)
         load=0
         if(taskcount .lt. subtasks)then
            scal = 1.0 - (1.0/(taskcount+2)) -  
     &           (taskcount-2.0)/(subtasks*subtasks)
         else 
            scal=0
            if (taskcount .lt.(numtasks - numtasks/8.0)) then
               part = 2
            else
               part = 4
            end if
            do i=1,part 
               scal = scal + 1.0/count + add/(subtasks*subtasks*1.0)    
               count= count - 1
               add= add + 1
            end do
            scal = scal + 1.0 
         end if
        
         new_load = avg_load*scal
c         print*, "scal", scal 
c         print*, "new_load for task",new_load, taskcount
         
         do while (load.lt.new_load)
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
                  
 31         continue
c         print*,"taskcount load",taskcount,load
         last_polar(taskcount)=countmol-1
         taskcount=taskcount+1
         
      end do

c      print*,"taskcount=",taskcount
c      if(taskcount.gt.numtasks) then
c        print*,"Tilt! Taskcount exceeds total number of tasks!"
c      end if

c      do i=0,taskcount-1
c         print*,"i nmol_nonzero start_polar",start_polar(i)
c         print*,"i nmol start",molnew(start_polar(i))
c         print*,"i nmol_nonzero last_polar",last_polar(i)
c         print*,"i nmol end",molnew(last_polar(i))
c      end do

      numtasks_polar=taskcount
      return
      end
