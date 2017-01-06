      subroutine clust3blists_offset
      use sizes
      use molcul
      use neigh2clust
      use mpidat
      use limits
      use neigh
      use neigh3b
      implicit none
      integer num3blists,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,np2,clust2
      integer, allocatable :: pnum(:)
      integer kouter,k3,k1,counter,limit
      real*8 xr1,yr1,zr1,r1_2,Rcut3b,Rcut3b2
      logical first3b
      integer taskid_offset,clust3,np3
      save first3b
      data first3b  / .true. /


c      Rcut3b=cut3b_input
c      Rcut3b2=Rcut3b*Rcut3b

c      if(first3b) then
         first3b=.false.

         !print*,"maxsizeclust=",maxsizeclust
         num3blists=0         
         maxnpole3b=0

         taskid_offset=taskid-clustcount-numtasks_polar
         
         do kouter=start_polar3b(taskid_offset),
     &                last_polar3b(taskid_offset)
            clust1=molnew3b(kouter)

            np1=0
            do i=1,sizeclust(clust1)
               np1=np1+3
            end do
            
            do k3=1,num_mollst_chunk3b(clust1,taskid_offset)

               do k1=start_mollst_polar3b(k3,clust1,taskid_offset),
     &                   last_mollst_polar3b(k3,clust1,taskid_offset),2
                  clust2=mollst3(k1,clust1)
                  clust3=mollst3(k1+1,clust1)


                    num3blists=num3blists+1
                    np2=0
                    do i=1,sizeclust(clust2)
                       np2=np2+3
                    end do

                    np3=0
                    do i=1,sizeclust(clust3)
                       np3=np3+3
                    end do

                    npole3b=np1+np2+np3
                    if(npole3b.gt.maxnpole3b) then
                     maxnpole3b=npole3b
                    end if
               end do
            end do
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

         if(listbcast) then
           if(allocated(elst3b)) deallocate(elst3b)
           if(allocated(nelst3b)) deallocate(nelst3b)
           if(allocated(domlstclust3b)) deallocate(domlstclust3b)
           if(allocated(xmold3b)) deallocate(xmold3b)
           if(allocated(ymold3b)) deallocate(ymold3b)
           if(allocated(zmold3b)) deallocate(zmold3b)
           if(allocated(ulst3b)) deallocate(ulst3b)
           if(allocated(nulst3b)) deallocate(nulst3b)
           if(allocated(doulstclust3b)) deallocate(doulstclust3b)
           if(allocated(xuold3b)) deallocate(xuold3b)
           if(allocated(yuold3b)) deallocate(yuold3b)
           if(allocated(zuold3b)) deallocate(zuold3b)
         end if


         if(.not.allocated(elst3b)) 
     &       allocate(elst3b(maxelst,maxnpole3b,num3blists))
         if(.not.allocated(nelst3b))
     &        allocate(nelst3b(maxnpole3b,num3blists))
         if(.not.allocated(domlstclust3b))
     &        allocate(domlstclust3b(num3blists))    

         if(listbcast) then
         do i=1,num3blists
            domlstclust3b(i)=.true.
         end do     
         end if

         if(.not.allocated(xmold3b))
     &      allocate(xmold3b(maxnpole3b,num3blists))
         if(.not.allocated(ymold3b))
     &      allocate(ymold3b(maxnpole3b,num3blists))
         if(.not.allocated(zmold3b))
     &      allocate(zmold3b(maxnpole3b,num3blists))


         if(.not.allocated(ulst3b))
     &       allocate(ulst3b(maxulst2,maxnpole3b,num3blists))
         if(.not.allocated(nulst3b))
     &        allocate(nulst3b(maxnpole3b,num3blists))
         if(.not.allocated(doulstclust3b))
     &        allocate(doulstclust3b(num3blists))

         if(listbcast) then
         do i=1,num3blists
            doulstclust3b(i)=.true.
         end do
         end if

         if(.not.allocated(xuold3b))
     &      allocate(xuold3b(maxnpole3b,num3blists))
         if(.not.allocated(yuold3b))
     &      allocate(yuold3b(maxnpole3b,num3blists))
         if(.not.allocated(zuold3b))
     &      allocate(zuold3b(maxnpole3b,num3blists))

c      end if
      !print*,"task=",taskid,"maxnpole3b=",maxnpole3b
      !print*,"task=",taskid,"num3blists=",num3blists
      allocate(pnum(maxnpole3b))

      counter=0

      do kouter=start_polar3b(taskid_offset),last_polar3b(taskid_offset)
         clust1=molnew3b(kouter)

           np1=0
           do i=1,sizeclust(clust1)
             moli1=clust(i,clust1)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           
           do k3=1,num_mollst_chunk3b(clust1,taskid_offset)
               do k1=start_mollst_polar3b(k3,clust1,taskid_offset),
     &                    last_mollst_polar3b(k3,clust1,taskid_offset),2
                  clust2=mollst3(k1,clust1) 
                  clust3=mollst3(k1+1,clust1)

                  counter=counter+1
                  np2=0
                  do i=1,sizeclust(clust2)
                    moli1=clust(i,clust2)
                    pnum(np1+3*(i-1)+1)=imol(1,moli1)
                    pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
                    pnum(np1+3*(i-1)+3)=imol(2,moli1)
                    np2=np2+3
                  end do

                  np3=0
                  do i=1,sizeclust(clust3)
                    moli1=clust(i,clust3)
                    pnum(np1+np2+3*(i-1)+1)=imol(1,moli1)
                    pnum(np1+np2+3*(i-1)+2)=imol(1,moli1)+1
                    pnum(np1+np2+3*(i-1)+3)=imol(2,moli1)
                    np3=np3+3
                  end do

                  npole3b=np1+np2+np3
                  call mlist_3bclust(npole3b,pnum,counter)
                  call ulist_3bclust(npole3b,pnum,counter)

               end do
           end do
      end do

      deallocate(pnum)
      return
      end
