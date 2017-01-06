      subroutine clust2blatticeboxes 
      use sizes
      use molcul
      use neigh2clust
      use mpidat
      use limits
      use neigh
      use neigh3b
      use boxes2bclust
      use cell2bclust
      use atoms
      implicit none
      integer num2blists,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,np2,clust2
      integer, allocatable :: pnum(:)
      integer kouter,k3,k1,counter,limit
      real*8 xr1,yr1,zr1,r1_2,Rcut2b,Rcut2b2
      logical first2b
      real*8 xmax,ymax,zmax,xdiff,ydiff,zdiff
      integer l1,l3,k
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

         if(.not.allocated(xbox2b))
     &       allocate(xbox2b(num2blists))
         if(.not.allocated(ybox2b))
     &        allocate(ybox2b(num2blists))
         if(.not.allocated(zbox2b))
     &        allocate(zbox2b(num2blists))
         if(.not.allocated(xcell2b))
     &       allocate(xcell2b(num2blists))
         if(.not.allocated(ycell2b))
     &        allocate(ycell2b(num2blists))
         if(.not.allocated(zcell2b))
     &        allocate(zcell2b(num2blists))
         if(.not.allocated(xbox2_2b))
     &       allocate(xbox2_2b(num2blists))
         if(.not.allocated(ybox2_2b))
     &        allocate(ybox2_2b(num2blists))
         if(.not.allocated(zbox2_2b))
     &        allocate(zbox2_2b(num2blists))
         if(.not.allocated(xcell2_2b))
     &       allocate(xcell2_2b(num2blists))
         if(.not.allocated(ycell2_2b))
     &        allocate(ycell2_2b(num2blists))
         if(.not.allocated(zcell2_2b))
     &        allocate(zcell2_2b(num2blists))
         if(.not.allocated(ncell2b))
     &       allocate(ncell2b(num2blists))
         if(.not.allocated(icell2b))
     &       allocate(icell2b(3,maxcell,num2blists))
         if(.not.allocated(volbox2b))
     &      allocate(volbox2b(num2blists))
         if(.not.allocated(lvec2b))
     &      allocate(lvec2b(3,3,num2blists))
         if(.not.allocated(recip2b))
     &      allocate(recip2b(3,3,num2blists))


c      end if
      !print*,"task=",taskid,"maxnpole3b=",maxnpole3b
      print*,"Boxtask=",taskid,"num2blists=",num2blists
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
                  xmax=0.0d0
                  ymax=0.0d0
                  zmax=0.0d0
                  do l1=1,npole3b-1
                     i=pnum(l1)
                     do l3=l1+1,npole3b
                       k=pnum(l3)
                       xdiff=abs(x(k)-x(i))
                       ydiff=abs(y(k)-y(i))
                       zdiff=abs(z(k)-z(i))
                       if(xdiff.gt.xmax) then
                         xmax=xdiff
                       end if
                       if(ydiff.gt.ymax) then
                         ymax=ydiff
                       end if
                       if(zdiff.gt.zmax) then
                         zmax=zdiff
                       end if
                     end do
                  end do
                  xbox2b(counter)=62.0d0!xmax+20.0d0
                  ybox2b(counter)=62.0d0!ymax+20.0d0
                  zbox2b(counter)=62.0d0!zmax+20.0d0

                  call lattice2bclust(counter)

c                 end if
               end do
           end do
      end do

      deallocate(pnum)
      return
      end
