

      subroutine clust1blists (start,last,moli1rmndr)
      use sizes
      use molcul
      use neigh2clust
      use limits
      use neigh
      use neigh3b
      implicit none
      integer num1blists,start,last,moli1rmndr,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,counter,limit
      integer, allocatable :: pnum(:)
      logical first
      save first
      data first  / .true. /

      if(first) then
         first=.false.
         if(moli1rmndr.ne.0) then
           num1blists=last-start+2      
         else 
           num1blists=last-start+1
         end if

         maxnpole3b=0

         do clust1=start,last
            npole3b=0
            do i=1,sizeclust(clust1)
               npole3b=npole3b+3
            end do
            if(npole3b.gt.maxnpole3b) then
              maxnpole3b=npole3b
            end if
         end do

         lbuffer = 2.0d0
         pbuffer = 2.0d0
         mbuf2clust = (mpolecut+lbuffer)**2
         mbufxclust = (mpolecut+2.0d0*lbuffer)**2
         ubuf2 = (usolvcut+pbuffer)**2
         pbuf2 = (0.5d0*pbuffer)**2
         ubufx = (usolvcut+2.0d0*pbuffer)**2

         maxulst2 = 500
         limit = int(0.5d0*sqrt(ubuf2)**3) + 50
         maxulst2 = min(limit,maxulst2)


         if(.not.allocated(elst1b)) 
     &       allocate(elst1b(maxelst,maxnpole3b,num1blists))
         if(.not.allocated(nelst1b))
     &        allocate(nelst1b(maxnpole3b,num1blists))
         if(.not.allocated(domlstclust1b))
     &        allocate(domlstclust1b(num1blists))    

         if(listbcast) then
         do i=1,num1blists
            domlstclust1b(i)=.true.
         end do     
         end if

         if(.not.allocated(xmold1b))
     &      allocate(xmold1b(maxnpole3b,num1blists))
         if(.not.allocated(ymold1b))
     &      allocate(ymold1b(maxnpole3b,num1blists))
         if(.not.allocated(zmold1b))
     &      allocate(zmold1b(maxnpole3b,num1blists))


         if(.not.allocated(ulst1b))
     &       allocate(ulst1b(maxulst2,maxnpole3b,num1blists))
         if(.not.allocated(nulst1b))
     &        allocate(nulst1b(maxnpole3b,num1blists))
         if(.not.allocated(doulstclust1b))
     &        allocate(doulstclust1b(num1blists))

         if(listbcast) then
         do i=1,num1blists
            doulstclust1b(i)=.true.
         end do
         end if

         if(.not.allocated(xuold1b))
     &      allocate(xuold1b(maxnpole3b,num1blists))
         if(.not.allocated(yuold1b))
     &      allocate(yuold1b(maxnpole3b,num1blists))
         if(.not.allocated(zuold1b))
     &      allocate(zuold1b(maxnpole3b,num1blists))

      end if

      allocate(pnum(3*maxsizeclust))

      counter=0
      do clust1=start,last
         counter=counter+1

           npole3b=0
           np1=0
           do i=1,sizeclust(clust1)
             moli1=clust(i,clust1)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           npole3b=np1
           
           call mlist_1bclust(npole3b,pnum,counter)
           call ulist_1bclust(npole3b,pnum,counter)
      end do

      if(moli1rmndr.ne.0) then
         counter=counter+1

          npole3b=0
           np1=0
           do i=1,sizeclust(moli1rmndr)
             moli1=clust(i,moli1rmndr)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           npole3b=np1

           call mlist_1bclust(npole3b,pnum,counter)
           call ulist_1bclust(npole3b,pnum,counter)
      end if    

      deallocate(pnum)
      return
      end 
