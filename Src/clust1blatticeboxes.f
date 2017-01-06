

      subroutine clust1blatticeboxes (start,last,moli1rmndr)
      use sizes
      use molcul
      use neigh2clust
      use limits
      use neigh
      use neigh3b
      use boxes1bclust
      use cell1bclust
      use atoms
      use mpidat
      implicit none
      integer num1blists,start,last,moli1rmndr,clust1,maxnpole3b
      integer i,npole3b,np1,moli1,counter,limit
      integer, allocatable :: pnum(:)
      logical first
      real*8 xmax,ymax,zmax,xdiff,ydiff,zdiff
      integer l1,l3,k
      save first
      data first  / .true. /

      if(first) then
         first=.false.
         if(moli1rmndr.ne.0) then
           num1blists=last-start+2      
         else 
           num1blists=last-start+1
         end if

         !maxnpole3b=0

         !do clust1=start,last
         !   npole3b=0
         !   do i=1,sizeclust(clust1)
         !      npole3b=npole3b+3
         !   end do
         !   if(npole3b.gt.maxnpole3b) then
         !     maxnpole3b=npole3b
         !   end if
         !end do


         if(.not.allocated(xbox1b)) 
     &       allocate(xbox1b(num1blists))
         if(.not.allocated(ybox1b))
     &        allocate(ybox1b(num1blists))
         if(.not.allocated(zbox1b))
     &        allocate(zbox1b(num1blists))    

         if(.not.allocated(xcell1b))
     &       allocate(xcell1b(num1blists))
         if(.not.allocated(ycell1b))
     &        allocate(ycell1b(num1blists))
         if(.not.allocated(zcell1b))
     &        allocate(zcell1b(num1blists))

c         if(.not.allocated(alpha1b))
c     &      allocate(alpha1b(num1blists))
c         if(.not.allocated(beta1b))
c     &      allocate(beta1b(num1blists))
c         if(.not.allocated(gamma1b))
c     &      allocate(gamma1b(num1blists))
         if(.not.allocated(xbox2_1b))
     &       allocate(xbox2_1b(num1blists))
         if(.not.allocated(ybox2_1b))
     &        allocate(ybox2_1b(num1blists))
         if(.not.allocated(zbox2_1b))
     &        allocate(zbox2_1b(num1blists))
         
         if(.not.allocated(xcell2_1b))
     &       allocate(xcell2_1b(num1blists))
         if(.not.allocated(ycell2_1b))
     &        allocate(ycell2_1b(num1blists))
         if(.not.allocated(zcell2_1b))
     &        allocate(zcell2_1b(num1blists))

         if(.not.allocated(ncell1b))
     &       allocate(ncell1b(num1blists))
         if(.not.allocated(icell1b))
     &       allocate(icell1b(3,maxcell,num1blists))
c	         if(.not.allocated(box34_1b))
c     &      allocate(box34_1b(num1blists))
         if(.not.allocated(volbox1b))
     &      allocate(volbox1b(num1blists))
c         if(.not.allocated(beta_sin1b))
c     &      allocate(beta_sin1b(num1blists))
c         if(.not.allocated(beta_cos1b))
c     &      allocate(beta_cos1b(num1blists))

c         if(.not.allocated(gamma_sin1b))
c     &      allocate(gamma_sin1b(num1blists))
c         if(.not.allocated(gamma_cos1b))
c     &      allocate(gamma_cos1b(num1blists))
c         if(.not.allocated(beta_term1b))
c     &     allocate(beta_term1b(num1blists))
c         if(.not.allocated(gamma_term1b))
c     &     allocate(gamma_term1b(num1blists))
         if(.not.allocated(lvec1b))
     &      allocate(lvec1b(3,3,num1blists))
         if(.not.allocated(recip1b))
     &      allocate(recip1b(3,3,num1blists))
c         if(.not.allocated(orthogonal1b))
c     &      allocate(orthogonal1b(num1blists))
c         if(.not.allocated(monoclinic1b))
c     &     allocate(monoclinic1b(num1blists))
c         if(.not.allocated(triclinic1b))
c     &     allocate(triclinic1b(num1blists))
c         if(.not.allocated(octahedron1b))
c     &     allocate(octahedron1b(num1blists))
c         if(.not.allocated(spacegrp1b))
c     &     allocate(spacegrp1b(num1blists))

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
           xbox1b(counter)=62.0d0!xmax+20.0d0
           ybox1b(counter)=62.0d0!ymax+20.0d0
           zbox1b(counter)=62.0d0!zmax+20.0d0

           !call unitcell1bclust(counter,xmax,ymax,zmax)
           call lattice1bclust(counter)
           print*,"task=",taskid,"xbox1b_2=",xbox2_1b(counter)
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
           xbox1b(counter)=62.0d0!xmax
           ybox1b(counter)=62.0d0!ymax
           zbox1b(counter)=62.0d0!zmax

           call lattice1bclust(counter)
      end if    
      print*,"Boxtask=",taskid,"num1blists=",num1blists
      
      deallocate(pnum)
      return
      end 
