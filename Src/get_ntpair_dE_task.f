      subroutine get_ntpair_dE_task
      use mpidat
      use dEtensor
      use neigh
      use molcul
      use mpole
      implicit none
      integer offset,remainder,pairoffset_tot,pairoffset
      integer start,last,i,j,k,pnum(3),totelst
      offset=int(nmol/numtasks)
      remainder = mod(nmol,numtasks)

      pairoffset_tot=0

      do i=0,numtasks-1
         start = i*offset+1
         last  = start+offset-1
         ntpair_start(i)=pairoffset_tot+1
         pairoffset=0
         do j=start,last
            pnum(1)=imol(1,j)
            pnum(2)=imol(1,j)+1     
            pnum(3)=imol(2,j)
            do k=1,3              
              pairoffset_tot=pairoffset_tot+nelst(pnum(k)) 
              pairoffset=pairoffset+nelst(pnum(k))  
            end do
         end do
         ntpair_last(i)=ntpair_start(i)+pairoffset-1
      end do

      do i=0,remainder-1
         start=numtasks*offset+i+1         
         ntpair_start_rmndr(i)=pairoffset_tot+1
          pnum(1)=imol(1,start)
          pnum(2)=imol(1,start)+1
          pnum(3)=imol(2,start)
          pairoffset=0 
           do k=1,3
             pairoffset_tot=pairoffset_tot+nelst(pnum(k))
             pairoffset=pairoffset+nelst(pnum(k))
           end do         
         ntpair_last_rmndr(i)=ntpair_start_rmndr(i)+pairoffset-1  
      end do
      totelst=0
      do i=1,npole
         totelst=totelst+nelst(i)
      end do
      print*,"totelst",totelst
      print*,"ntpair_dE=",ntpair_dE
      do i=0,numtasks-1
        print*,"i=",i,"ntpair_start(i)=",ntpair_start(i)
        print*,"i=",i,"ntpair_last(i)=",ntpair_last(i)
      end do

      do i=0,remainder-1
        print*,"i=",i,"ntpair_start_rmndr(i)=",ntpair_start_rmndr(i)
        print*,"i=",i,"ntpair_last_rmndr(i)=",ntpair_last_rmndr(i)
      end do 

      return
      end
