      subroutine ewpolar_load_balance7
      use mpidat
      use neigh3b
      use molcul
      implicit none
      include 'mpif.h'
      integer i,j,nmollst_tot,load,countmol,taskcount
      integer nmol_nonzero,k,avg_load
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
c         if (.not. allocated(molnew))
c     &          allocate (molnew(nmol))
      end if      

      do i=1,nmol
c        if(nmollst(i).ne.0) then
          do j=1,nmollst(i)
             k=mollst(j,i)
            nmollst_tot = nmollst_tot + nmollst(k)
          end do
c            nmol_nonzero=nmol_nonzero+1
c            molnew(nmol_nonzero)=i
         ! print*,"molecule_index, nmollst(i)",i,nmollst(i)
c        end if
      end do

      !do i=1,nmol_nonzero
      ! print*,"molnew",molnew(i),"nmollst(molnew(i))",nmollst(molnew(i))
      !end do

c      print*,"nmol_nonzero",nmol_nonzero

      avg_load=nmollst_tot/numtasks
c      print*,"avg_load=",avg_load

      if(mod(nmollst_tot,numtasks).gt.0) then
        avg_load=avg_load+1
c        print*,"avg_load rmndr",avg_load
      end if


      countmol=1
      taskcount=0

c      do while (countmol.lt.nmol_nonzero)
      do while (countmol.lt.nmol)
         load=0
         do while (load.lt.avg_load)
            if(load.eq.0) then
              start_polar(taskcount)=countmol
            end if        

c            i=molnew(countmol) 
c            do j=1,nmollst(i)
            do j=1,nmollst(countmol)
c               k=mollst(j,i)
               k=mollst(j,countmol)
               load=load+nmollst(k)
            end do 
            countmol=countmol+1

c            if(countmol.eq.nmol_nonzero+1) then
            if(countmol.eq.nmol+1) then
               goto 31 
            end if
         end do
   31             continue
c         print*,"taskcount load",taskcount,load
         last_polar(taskcount)=countmol-1
         taskcount=taskcount+1
      end do

c      print*,"taskcount=",taskcount
c      if(taskcount.gt.numtasks) then
c        print*,"Tilt! Taskcount exceeds total number of tasks!"
c      end if

      do i=0,taskcount-1
         print*,"taskid start_polar Edwinload",i,start_polar(i)
      !   print*,"taskid nmol start",i,molnew(start_polar(i))
         print*,"taskid last_polar Edwinload",i,last_polar(i)
      !   print*,"taskid nmol end",i,molnew(last_polar(i))
      end do

      numtasks_polar=taskcount
      return
      end

