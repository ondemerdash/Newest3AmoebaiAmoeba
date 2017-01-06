      subroutine initpolarpred2b
      use sizes
      use molcul
      use neigh2clust
      use mpidat
      use limits
      use neigh
      use neigh3b
      use uprior
      use uprior2b
      implicit none
      integer num2blists,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,np2,clust2
      integer, allocatable :: pnum(:)
      integer kouter,k3,k1,counter,limit
      real*8 xr1,yr1,zr1,r1_2,Rcut2b,Rcut2b2
      logical first2b
      integer j,k
      save first2b
      data first2b  / .true. /


      Rcut2b=cut2b_input
      Rcut2b2=Rcut2b*Rcut2b

c      if(first2b) then
         first2b=.false.

c         print*,"maxsizeclust=",maxsizeclust
         num2blists=0         
         maxnpole3b=0

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
         
         do kouter=start_polar(taskid),last_polar(taskid)
            clust1=molnew(kouter)

            np1=0
            do i=1,sizeclust(clust1)
               np1=np1+3
            end do
            
            do k3=1,num_mollst_chunk(clust1,taskid)
               do k1=start_mollst_polar(k3,clust1,taskid),
     &                        last_mollst_polar(k3,clust1,taskid)
                  clust2=mollst(k1,clust1)

                    num2blists=num2blists+1
                    np2=0
                    do i=1,sizeclust(clust2)
                       np2=np2+3
                    end do

                    npole3b=np1+np2
                    if(npole3b.gt.maxnpole3b) then
                     maxnpole3b=npole3b
                    end if
               end do
            end do
         end do

         if(allocated(udalt2b)) deallocate(udalt2b)
         if(allocated(upalt2b)) deallocate(upalt2b)
         if(allocated(nualt2b)) deallocate(nualt2b)
         if(allocated(bpred2b)) deallocate(bpred2b)
         if(allocated(bpredp2b)) deallocate(bpredp2b)

         if(.not.allocated(udalt2b)) 
     &       allocate (udalt2b(maxualt,3,maxnpole3b,num2blists))
         if(.not.allocated(upalt2b)) 
     &       allocate (upalt2b(maxualt,3,maxnpole3b,num2blists))
         if(.not.allocated(nualt2b)) allocate (nualt2b(num2blists))

         if(.not.allocated(bpred2b)) 
     &       allocate (bpred2b(maxualt,num2blists))
         if(.not.allocated(bpredp2b))
     &       allocate (bpredp2b(maxualt,num2blists))
c      end if
c      print*,"task=",taskid,"maxnpole3b=",maxnpole3b
c      print*,"task=",taskid,"num2blists=",num2blists

      counter=0

      do kouter=start_polar(taskid),last_polar(taskid)
            clust1=molnew(kouter)
           do k3=1,num_mollst_chunk(clust1,taskid)
               do k1=start_mollst_polar(k3,clust1,taskid),
     &                       last_mollst_polar(k3,clust1,taskid)
                   counter=counter+1
                   if(counter.gt.num2blists) then
                    print*,"Tilt! tsk=",taskid
                   end if
                   nualt2b(counter)=0
                   do i=1,maxnpole3b
                    do j=1,3
                     do k=1, maxualt
                       udalt2b(k,j,i,counter)=0.0d0
                       upalt2b(k,j,i,counter)=0.0d0
                     end do
                    end do
                   end do
               end do
           end do
      end do

      return
      end
