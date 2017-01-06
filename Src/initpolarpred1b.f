      subroutine initpolarpred1b (start,last,moli1rmndr)
      use sizes
      use molcul
      use neigh2clust
      use limits
      use neigh
      use neigh3b
      use uprior
      use uprior1b
      implicit none
      integer num1blists,start,last,moli1rmndr,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,counter,limit
      integer, allocatable :: pnum(:)
      logical first
      integer j,k
      save first
      data first  / .true. /

c      if(first) then
c         first=.false.
         if(moli1rmndr.ne.0) then
           num1blists=last-start+2      
         else 
           num1blists=last-start+1
         end if

         if(.not.allocated(nualt1b)) allocate(nualt1b(num1blists))
         if(.not.allocated(udalt1b)) 
     &       allocate (udalt1b(maxualt,3,3*maxsizeclust,num1blists))
         if(.not.allocated(upalt1b)) 
     &       allocate (upalt1b(maxualt,3,3*maxsizeclust,num1blists))
         if(.not.allocated(bpred1b)) 
     &       allocate (bpred1b(maxualt,num1blists))
         if(.not.allocated(bpredp1b)) 
     &       allocate (bpredp1b(maxualt,num1blists))


c
c     set the Gear predictor binomial coefficients
c
         gear(1) = 6.0d0
         gear(2) = -15.0d0
         gear(3) = 20.0d0
         gear(4) = -15.0d0
         gear(5) = 6.0d0
         gear(6) = -1.0d0
         gear(7) = 0.0d0
c
c     set always stable predictor-corrector (ASPC) coefficients
c
         aspc(1) = 22.0d0 / 7.0d0
         aspc(2) = -55.0d0 / 14.0d0
         aspc(3) = 55.0d0 / 21.0d0
         aspc(4) = -22.0d0 / 21.0d0
         aspc(5) = 5.0d0 / 21.0d0
         aspc(6) = -1.0d0 / 42.0d0
         aspc(7) = 0.0d0
c      end if

      counter=0
      do clust1=start,last
         counter=counter+1
         nualt1b(counter)=0
         do i=1,3*maxsizeclust
            do j=1,3
               do k=1, maxualt
                 udalt1b(k,j,i,counter)=0.0d0
                 upalt1b(k,j,i,counter)=0.0d0
               end do
            end do
         end do           
      end do

      if(moli1rmndr.ne.0) then
         counter=counter+1
         nualt1b(counter)=0
         do i=1,3*maxsizeclust
            do j=1,3
               do k=1, maxualt
                 udalt1b(k,j,i,counter)=0.0d0
                 upalt1b(k,j,i,counter)=0.0d0
               end do
            end do
         end do
      end if    

      return
      end 
