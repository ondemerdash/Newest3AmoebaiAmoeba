c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns particle mesh Ewald parameters and options
c     for a periodic system
c
c
      subroutine kewaldreg
      use limits 
      use math
      use boxes
      use ewald
c      use ewreg
      use parewreg
      implicit none

c
c
c     estimate an optimal value for the Ewald coefficient
c
      call ewaldcof (aewald,ewaldcut)

      jmax=int(ewaldcut*aewald*xbox*aewald/pi)
      kmax=int(ewaldcut*aewald*ybox*aewald/pi)
      lmax=int(ewaldcut*aewald*zbox*aewald/pi)

      print*,"jmax calc kmax calc lmax calc",jmax,kmax,lmax
      jmax=min(maxvec,jmax)
      kmax=min(maxvec,kmax)
      lmax=min(maxvec,lmax)

c
c     set the system extent for nonperiodic Ewald summation
c
      return
      end

      subroutine ewaldreglist
      use ewaldneigh
      use parewreg
      use mpidat
      implicit none
      integer i,j,k,l,ksq,lmin,kmin
      integer totiter,taskcnt,taskcnt2,taskcnt_overflo
      integer offseterecip,remainder,numchunk

c      kcut=105
    
      print*,"kcut=",kcut      
      print*,"jmax kmax lmax in ewaldreglist",jmax,kmax,lmax
      if (.not.allocated(njlst)) allocate(njlst(0:jmax))
      if (.not.allocated(jlst)) 
     &  allocate(jlst( (2*kmax+1)*(2*lmax+1)*2,0:jmax))

      if (.not.allocated(starterecip)) 
     &    allocate(starterecip(0:numtasks-1))
      if (.not.allocated(enderecip)) allocate(enderecip(0:numtasks-1))
      if (.not.allocated(jiter)) allocate(jiter(0:numtasks-1))
c      print*,"After allocs"
c      if (.not.allocated(doremainder)) 
c     &       allocate(doremainder(0:numtasks-1))
      if (.not.allocated(starterecip2))
     &    allocate(starterecip2(0:numtasks-1))
      if (.not.allocated(enderecip2)) allocate(enderecip2(0:numtasks-1))
      if (.not.allocated(jiter2)) allocate(jiter2(0:numtasks-1))
      do j=0,jmax
         njlst(j)=0
      end do

      remainder_bcast_alloc=.false. 
c      do i=0,numtasks-1
c         doremainder(i)=.false.
c      end do

      j=0
      k=0
      do l=1,lmax
         ksq=j*j + k*k + l*l
         if (ksq .le. kcut) then
            njlst(j)=njlst(j)+1
            jlst(njlst(j),j)=k
            njlst(j)=njlst(j)+1
            jlst(njlst(j),j)=l
         end if
      end do
     
      lmin=-lmax
   
      do k=1,kmax
         do l=lmin,lmax
          ksq=j*j + k*k + l*l
           if (ksq .le. kcut) then
            njlst(j)=njlst(j)+1
            jlst(njlst(j),j)=k
            njlst(j)=njlst(j)+1
            jlst(njlst(j),j)=l
           end if
         end do 
      end do

      kmin=-kmax

      do j=1,jmax
        do k=kmin,kmax
          do l=lmin,lmax
             ksq=j*j + k*k + l*l
             if (ksq .le. kcut) then
               njlst(j)=njlst(j)+1
               jlst(njlst(j),j)=k
               njlst(j)=njlst(j)+1
               jlst(njlst(j),j)=l
             end if
          end do  
        end do   
      end do
      print*,"Absolute max size for jlist",(2*kmax+1)*(2*lmax+1)*2
      totiter=0
      do j=0,jmax
         totiter=totiter+njlst(j)  
        print*,"j Numiters Njlst(i)",j,njlst(j)/2,njlst(j)
      end do 

      offseterecip=int(totiter/numtasks)
c      print*,"Badoffseterecip",offseterecip

      if( mod(offseterecip,2).ne.0) then
      print*,"Badoffseterecip",offseterecip
        offseterecip=offseterecip+1
      end if
      print*,"Goodoffseterecip",offseterecip
 
      taskcnt=0
      taskcnt_overflo=0      
      taskcnt2=0
      do i=0,jmax
         if(njlst(i).ne.0) then
           taskcnt2=0
              if(njlst(i).gt.offseterecip) then
                numchunk=int(njlst(i)/offseterecip)
                remainder=njlst(i)-numchunk*offseterecip
         print*,"j njlst nchunk remainder",i,njlst(i),numchunk,remainder
                do j=1,numchunk
                  if(taskcnt.le.numtasks-1) then
                   starterecip(taskcnt)=taskcnt2*offseterecip+1
c            enderecip(taskcnt)=taskcnt2*offseterecip+1+offseterecip-1
                   enderecip(taskcnt)=taskcnt2*offseterecip+offseterecip
                   jiter(taskcnt)=i
                   taskcnt=taskcnt+1
                   taskcnt2=taskcnt2+1                 
                  else
                   remainder_bcast_alloc=.true.
                   starterecip2(taskcnt_overflo)=taskcnt2*offseterecip+1
                   enderecip2(taskcnt_overflo)=taskcnt2*offseterecip
     &                                         +offseterecip
                   jiter2(taskcnt_overflo)=i
                   !doremainder(taskcnt_overflo)=.true.
                   taskcnt_overflo=taskcnt_overflo+1
                   taskcnt2=taskcnt2+1
                   taskcnt=taskcnt+1
                  end if
                end do 
                if(taskcnt.le.numtasks) then
                 enderecip(taskcnt-1)=enderecip(taskcnt-1)+remainder 
                else
                 enderecip2(taskcnt_overflo-1)=
     &            enderecip2(taskcnt_overflo-1)+remainder
                end if
              else
                numchunk=1
                remainder=0
         print*,"j njlst nchunk remainder",i,njlst(i),numchunk,remainder
                 if(taskcnt.le.numtasks-1) then
                   starterecip(taskcnt)=taskcnt2*offseterecip+1
                   enderecip(taskcnt)=njlst(i)
                   jiter(taskcnt)=i
                   taskcnt=taskcnt+1
                   taskcnt2=taskcnt2+1
                 else
                   remainder_bcast_alloc=.true.
                   starterecip2(taskcnt_overflo)=taskcnt2*offseterecip+1
                   enderecip2(taskcnt_overflo)=njlst(i)
                   jiter2(taskcnt_overflo)=i
                   !doremainder(taskcnt_overflo)=.true.
                   taskcnt=taskcnt+1
                   taskcnt_overflo=taskcnt_overflo+1
                   taskcnt2=taskcnt2+1
                 end if 
              end if 
         end if
      end do 
      
      taskcnterecip=taskcnt
      if(taskcnt.gt.numtasks) then
        taskcnterecip=numtasks
      end if

      taskcnterecip_overflo=taskcnt_overflo
      print*,"Task count total",taskcnt
      print*,"Taskcount erecip",taskcnterecip
      print*,"Taskcount erecip overflo",taskcnterecip_overflo
      print*,"remainder_bcast_alloc",remainder_bcast_alloc
      if(taskcnterecip.gt.numtasks) then
        print*,"Tilt! Num tasks for Kewald exceeds num MPI tasks!"
         call fatal
      end if

      return
      end

      subroutine bcast_ewaldreglist
      use ewaldneigh
      use parewreg
      use mpidat
      implicit none
      include 'mpif.h'
      integer ierr

c      call mpi_bcast(njlst,jmax+1,mpi_integer,master,
c     &   mpi_comm_world,ierr)
      call mpi_bcast(jlst,(2*kmax+1)*(2*lmax+1)*2*(jmax+1),mpi_integer,
     &   master,mpi_comm_world,ierr)
      call mpi_bcast(taskcnterecip,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taskcnterecip_overflo,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(starterecip,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(enderecip,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
      call mpi_bcast(jiter,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
c      print*,"remainder_bcast_alloc bcast_ewaldreglist",
c     &  remainder_bcast_alloc
      if(remainder_bcast_alloc) then
        call mpi_bcast(starterecip2,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
        call mpi_bcast(enderecip2,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
        call mpi_bcast(jiter2,numtasks,mpi_integer,master,
     &   mpi_comm_world,ierr)
      end if
      return
      end

