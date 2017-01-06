c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readxyz  --  input of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readxyz" gets a set of Cartesian coordinates from
c     an external disk file
c
c
      subroutine readclust (iclust)
      use sizes
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use iounit
      use titles
      use neigh2clust
      use molcul
      use mpidat
      implicit none
      integer i,j,k,m
      integer iclust,nmax
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      integer, allocatable :: list(:)
      real*8 xlen,ylen,zlen
      real*8 aang,bang,gang
      logical exist,opened
      logical quit,reorder
      logical clash
      character*120 xyzfile
      character*120 record
      character*120 string
      integer molclust,moli1,lenmol,l1,k1,j1
      real*8 xcm,ycm,zcm,x1,y1,z1
c
c
c     initialize the total number of atoms in the system
c
c      n = 0
c
c     open the input file if it has not already been done
c
      inquire (unit=iclust,opened=opened)
c      if (.not. opened) then
c         xyzfile = filename(1:leng)//'.xyz'
c         call version (xyzfile,'old')
c         inquire (file=xyzfile,exist=exist)
c         if (exist) then
c            open (unit=ixyz,file=xyzfile,status='old')
c            rewind (unit=ixyz)
c         else
c            write (iout,10)
c   10       format (/,' READXYZ  --  Unable to Find the Cartesian',
c     &                 ' Coordinates File')
c            call fatal
c         end if
c      end if
c
c     read first line and return if already at end of file
c
c      quit = .false.
c      abort = .true.
c      size = 0
c      do while (size .eq. 0)
c         read (iclust,20,err=80,end=80)  record
c   20    format (a120)
c         size = trimtext (record)
c      end do
c      abort = .false.
c      quit = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
c      next = 1
c      call gettext (record,string,next)
         read (iclust,20,err=80,end=80)  record
   20    format (a120)
      read (record,*,err=80,end=80) clustcount,maxsizeclust
      !print*,"Clustcount read from Kmeans file=",clustcount
      !print*,"Maxsixe of Clusters read from Kmeans file=",maxsizeclust

      if(.not.allocated(clust)) allocate(clust(maxsizeclust,clustcount))
      if(.not.allocated(sizeclust)) allocate(sizeclust(clustcount))
      if(.not.allocated(clust_cm)) allocate(clust_cm(3,clustcount))
      if(.not.allocated(distmax)) allocate(distmax(clustcount))

      
c
c     extract the title and determine its length
c
c      string = record(next:120)
c      first = nexttext (string)
c      last = trimtext (string)
c      if (last .eq. 0) then
c         title = ' '
c         ltitle = 0
c      else
c         title = string(first:last)
c         ltitle = trimtext (title)
c      end if



c
c     check for too few or too many total atoms in the file
c
c
c     initialize coordinates and connectivities for each atom
c


c     initialize clust_cm and sizeclust

      do i=1,clustcount
         sizeclust(i)=0
         do j=1,3
            clust_cm(j,i)=0.0d0
         end do
      end do


c     Read in cluster assignment for each molecule
      do i = 1, nmol
         read (iclust,50,err=80,end=80)  record
   50       format (a120)
         read (record,*,err=80,end=80) molclust
c       print*,"Reading in cluster index for each mol",i,molclust
         sizeclust(molclust)=sizeclust(molclust)+1
         clust(sizeclust(molclust),molclust)=i
      end do

      do i=1,clustcount
         xcm=0.0d0
         ycm=0.0d0
         zcm=0.0d0
         do j=1,sizeclust(i)
            moli1=clust(j,i)
            lenmol=imol(2,moli1)-imol(1,moli1)+1

            do l1 = 1, lenmol
               k1=l1-1
               j1=imol(1,moli1)+k1
               if(name(j1).eq.'O') then
                 x1 = x(j1)
                 y1 = y(j1)
                 z1 = z(j1)
               end if
            end do
            
            xcm=xcm+x1
            ycm=ycm+y1
            zcm=zcm+z1
         end do
         clust_cm(1,i)=xcm/sizeclust(i)
         clust_cm(2,i)=ycm/sizeclust(i)
         clust_cm(3,i)=zcm/sizeclust(i)
      end do      


      if(taskid.eq.0) then
        do i=1,clustcount
          print*,"Clustind=",i,"sizeclust(i)=",sizeclust(i)
        end do
      end if

   80 continue

c
c     read the coordinates and connectivities for each atom
c
c
c     an error occurred in reading the coordinate file
c
c
c     for each atom, count and sort its attached atoms
c
c
c     perform dynamic allocation of some local arrays
c
c
c     check for scrambled atom order and attempt to renumber
c
c
c     perform deallocation of some local arrays
c
c
c     check for atom pairs with identical coordinates
c
c
c     make sure that all connectivities are bidirectional
c
      return
      end
