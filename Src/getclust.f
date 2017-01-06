c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getclust 
c     ##                                                             ##
c     #################################################################
c
c
c     "getclust" reads in a list of kmeans clusters, 1 for each molecule 
c
c
      subroutine getclust
      use inform
      use iounit
      use output
      implicit none
      integer iclust
      integer freeunit
      logical exist
      character*120 clustfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (clustfile,exist)
      if (exist) then
         call basefileclust (clustfile)
         call suffix (clustfile,'txt2','old')
         inquire (file=clustfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Clust File Name :  ',$)
         read (input,20)  clustfile
   20    format (a120)
         call basefileclust (clustfile)
         call suffix (clustfile,'txt2','old')
         inquire (file=clustfile,exist=exist)
      end do
c
c     first open and then read the Cartesian coordinates file
c
c      coordtype = 'CARTESIAN'
      iclust = freeunit ()
      open (unit=iclust,file=clustfile,status='old')
      rewind (unit=iclust)
      call readclust (iclust)
      close (unit=iclust)
c
c     quit if the Cartesian coordinates file contains no atoms
c
      if (abort) then
         write (iout,30)
   30    format (/,' Error With Cluster Input file',
     &              ' Tilt!')
         call fatal
      end if
      return
      end
