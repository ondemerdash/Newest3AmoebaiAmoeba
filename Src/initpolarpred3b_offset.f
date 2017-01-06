      subroutine initpolarpred3b_offset 
      use sizes
      use molcul
      use neigh2clust
      use mpidat
      use limits
      use neigh
      use neigh3b
      use uprior
      use uprior3b
      implicit none
      integer num3blists,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,np2,clust2
      integer, allocatable :: pnum(:)
      integer kouter,k3,k1,counter,limit
      real*8 xr1,yr1,zr1,r1_2,Rcut3b,Rcut3b2
      logical first3b
      integer taskid_offset
      integer j,k,clust3,np3
      save first3b
      data first3b  / .true. /


      !Rcut3b=cut3b_input
      !Rcut3b2=Rcut3b*Rcut3b

c      if(first3b) then
         first3b=.false.

         !print*,"maxsizeclust=",maxsizeclust
         num3blists=0         
         maxnpole3b=0

         taskid_offset=taskid-clustcount-numtasks_polar
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
         
         do kouter=start_polar3b(taskid_offset),
     &                last_polar3b(taskid_offset)
            clust1=molnew3b(kouter)

            np1=0
            do i=1,sizeclust(clust1)
               np1=np1+3
            end do
            
            do k3=1,num_mollst_chunk3b(clust1,taskid_offset)
               do k1=start_mollst_polar3b(k3,clust1,taskid_offset),
     &                    last_mollst_polar3b(k3,clust1,taskid_offset),2
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

         if(listbcast) then
           if(allocated(udalt3b)) deallocate(udalt3b)
           if(allocated(upalt3b)) deallocate(upalt3b)
           if(allocated(nualt3b)) deallocate(nualt3b)
           if(allocated(bpred3b)) deallocate(bpred3b)
           if(allocated(bpredp3b)) deallocate(bpredp3b)
         end if

         if(.not.allocated(udalt3b)) 
     &       allocate (udalt3b(maxualt,3,maxnpole3b,num3blists))
         if(.not.allocated(upalt3b)) 
     &       allocate (upalt3b(maxualt,3,maxnpole3b,num3blists))
         if(.not.allocated(nualt3b)) allocate (nualt3b(num3blists))

         if(.not.allocated(bpred3b)) 
     &       allocate (bpred3b(maxualt,num3blists))
         if(.not.allocated(bpredp3b))
     &       allocate (bpredp3b(maxualt,num3blists))

c      end if
      !print*,"task=",taskid,"maxnpole3b=",maxnpole3b
      !print*,"task=",taskid,"num3blists=",num3blists

      counter=0

      do kouter=start_polar3b(taskid_offset),last_polar3b(taskid_offset)
            clust1=molnew3b(kouter)
           do k3=1,num_mollst_chunk3b(clust1,taskid_offset)
               do k1=start_mollst_polar3b(k3,clust1,taskid_offset),
     &                   last_mollst_polar3b(k3,clust1,taskid_offset),2
                   counter=counter+1
           !        print*,"task=",taskid,"counter=",counter 
                   nualt3b(counter)=0
                   do i=1,maxnpole3b
                    do j=1,3
                     do k=1, maxualt
                       udalt3b(k,j,i,counter)=0.0d0
                       upalt3b(k,j,i,counter)=0.0d0
                     end do
                    end do
                   end do
               end do
           end do
      end do

      return
      end
