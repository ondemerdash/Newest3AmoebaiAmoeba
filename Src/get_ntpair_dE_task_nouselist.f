      subroutine get_ntpair_dE_task_nouselist
      use mpidat
      use dEtensor
      use neigh
      use molcul
      use mpole
      implicit none
      integer offset,remainder,pairoffset_tot,pairoffset
      integer start,last,i,j,k,pnum(3),totelst

      offset=int(ntpair_dE/numtasks)
      remainder = mod(ntpair_dE,numtasks)
      pairoffset_tot=0

      do i=0,numtasks-1
         ntpair_start_rmndr(i)=0
      end do

      do i=0,numtasks-1
         start = i*offset+1
         last  = start+offset-1
         ntpair_start(i)=start
         ntpair_last(i)=last
      end do

      do i=0,remainder-1
         start=numtasks*offset+i+1         
         ntpair_start_rmndr(i)=start
      end do

c      print*,"ntpair_dE=",ntpair_dE
c      do i=0,numtasks-1
c        print*,"i=",i,"ntpair_start(i)=",ntpair_start(i)
c        print*,"i=",i,"ntpair_last(i)=",ntpair_last(i)
c      end do

c      do i=0,remainder-1
c        print*,"i=",i,"ntpair_start_rmndr(i)=",ntpair_start_rmndr(i)
c        print*,"i=",i,"ntpair_last_rmndr(i)=",ntpair_last_rmndr(i)
c      end do 

      return
      end
