      subroutine clust2blists
      use sizes
      use molcul
      use neigh2clust
      use mpidat
      use limits
      use neigh
      use neigh3b
      implicit none
      integer num2blists,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,np2,clust2
      integer, allocatable :: pnum(:)
      integer kouter,k3,k1,counter,limit
      real*8 xr1,yr1,zr1,r1_2,Rcut2b,Rcut2b2
      logical first2b
      save first2b
      data first2b  / .true. /


      Rcut2b=cut2b_input
      Rcut2b2=Rcut2b*Rcut2b

c      if(first2b) then
         first2b=.false.

         !print*,"maxsizeclust=",maxsizeclust
         num2blists=0         
         maxnpole3b=0

         do kouter=start_polar(taskid),last_polar(taskid)
            clust1=molnew(kouter)

            np1=0
            do i=1,sizeclust(clust1)
               np1=np1+3
            end do
            
            do k3=1,num_mollst_chunk(clust1,taskid)
c                    num2blists=num2blists+1
               do k1=start_mollst_polar(k3,clust1,taskid),
     &                         last_mollst_polar(k3,clust1,taskid)
                  clust2=mollst(k1,clust1)
c                 xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
c                 yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
c                 zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
c                 call image(xr1,yr1,zr1)
c                 r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1

c                 if(r1_2.le.Rcut2b2) then

                    num2blists=num2blists+1
                    np2=0
                    do i=1,sizeclust(clust2)
                       np2=np2+3
                    end do

                    npole3b=np1+np2
                    if(npole3b.gt.maxnpole3b) then
                     maxnpole3b=npole3b
                    end if
c                 end if
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

         if(.not.allocated(elst2b)) 
     &       allocate(elst2b(maxelst,maxnpole3b,num2blists))
         if(.not.allocated(nelst2b))
     &        allocate(nelst2b(maxnpole3b,num2blists))
         if(.not.allocated(domlstclust2b))
     &        allocate(domlstclust2b(num2blists))    
         do i=1,num2blists
            domlstclust2b(i)=.true.
         end do     

         if(.not.allocated(xmold2b))
     &      allocate(xmold2b(maxnpole3b,num2blists))
         if(.not.allocated(ymold2b))
     &      allocate(ymold2b(maxnpole3b,num2blists))
         if(.not.allocated(zmold2b))
     &      allocate(zmold2b(maxnpole3b,num2blists))


         if(.not.allocated(ulst2b))
     &       allocate(ulst2b(maxulst2,maxnpole3b,num2blists))
         if(.not.allocated(nulst2b))
     &        allocate(nulst2b(maxnpole3b,num2blists))
         if(.not.allocated(doulstclust2b))
     &        allocate(doulstclust2b(num2blists))
         do i=1,num2blists
            doulstclust2b(i)=.true.
         end do

         if(.not.allocated(xuold2b))
     &      allocate(xuold2b(maxnpole3b,num2blists))
         if(.not.allocated(yuold2b))
     &      allocate(yuold2b(maxnpole3b,num2blists))
         if(.not.allocated(zuold2b))
     &      allocate(zuold2b(maxnpole3b,num2blists))

c      end if
      !print*,"task=",taskid,"maxnpole3b=",maxnpole3b
      !print*,"task=",taskid,"num2blists=",num2blists
      allocate(pnum(maxnpole3b))

      counter=0

      do kouter=start_polar(taskid),last_polar(taskid)
         clust1=molnew(kouter)

           np1=0
           do i=1,sizeclust(clust1)
             moli1=clust(i,clust1)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           
           do k3=1,num_mollst_chunk(clust1,taskid)
               do k1=start_mollst_polar(k3,clust1,taskid),
     &                         last_mollst_polar(k3,clust1,taskid)
                  clust2=mollst(k1,clust1) 

c                  xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
c                  yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
c                  zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
c                  call image(xr1,yr1,zr1)
c                  r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1

c                 if(r1_2.le.Rcut2b2) then
                  counter=counter+1
                  np2=0
                  do i=1,sizeclust(clust2)
                    moli1=clust(i,clust2)
                    pnum(np1+3*(i-1)+1)=imol(1,moli1)
                    pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
                    pnum(np1+3*(i-1)+3)=imol(2,moli1)
                    np2=np2+3
                  end do
                  npole3b=np1+np2 
                  call mlist_2bclust(npole3b,pnum,counter)
                  call ulist_2bclust(npole3b,pnum,counter)
c                 end if
               end do
           end do
      end do

      deallocate(pnum)
      return
      end
