c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine replica  --  periodicity via cell replication  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "replica" decides between images and replicates for generation
c     of periodic boundary conditions, and sets the cell replicate
c     list if the replicates method is to be used
c
c
      subroutine replicabox2bclust (cutoff,cnt)
      use sizes
      use bound
      use boxes
      use cell
      use boxes2bclust
      use cell2bclust
      use inform
      use iounit
      implicit none
      integer i,j,k
      integer nx,ny,nz
      real*8 cutoff,maximage
      real*8 xlimit,ylimit,zlimit
      integer cnt
c
c
c     only necessary if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     find the maximum sphere radius inscribed in periodic box
c
      if (orthogonal) then
         xlimit = xbox2_2b(cnt)
         ylimit = ybox2_2b(cnt)
         zlimit = zbox2_2b(cnt)
C   WE HAVE NO NEED FOR OTHER TYPES OF PERIODIC BOXES FOR THE CLUSTER 
C   SUBSYSTEMS RIGHT NOW.
c      else if (monoclinic) then
c         xlimit = xbox2 * beta_sin
c         ylimit = ybox2
c         zlimit = zbox2 * beta_sin
c      else if (triclinic) then
c         xlimit = xbox2 * beta_sin * gamma_sin
c         ylimit = ybox2 * gamma_sin
c         zlimit = zbox2 * beta_sin
c      else if (octahedron) then
c         xlimit = (sqrt(3.0d0)/4.0d0) * xbox
c         ylimit = xlimit
c         zlimit = xlimit
      end if
      maximage = min(xlimit,ylimit,zlimit)
c
c     use replicate method to handle cutoffs too large for images
c
      if (cutoff .le. maximage) then
         use_replica = .false.
      else
         use_replica = .true.
      end if
c
c     truncated octahedron cannot use the replicates method
c
      if (octahedron .and. use_replica) then
         write (iout,10)
   10    format (/,' REPLICA  --  Truncated Octahedron',
     &              ' cannot be Replicated')
         call fatal
      end if
c
c     find the number of replicates needed based on cutoff
c
      nx = int(cutoff/xlimit)
      ny = int(cutoff/ylimit)
      nz = int(cutoff/zlimit)
      if (cutoff .gt. dble(nx)*xlimit)  nx = nx + 1
      if (cutoff .gt. dble(ny)*ylimit)  ny = ny + 1
      if (cutoff .gt. dble(nz)*zlimit)  nz = nz + 1
      if (nx .lt. 1)  nx = 1
      if (ny .lt. 1)  ny = 1
      if (nz .lt. 1)  nz = 1
c
c     set the replicated cell length and the half width
c
      xcell2b(cnt) = dble(nx) * xbox2b(cnt)
      ycell2b(cnt) = dble(ny) * ybox2b(cnt)
      zcell2b(cnt) = dble(nz) * zbox2b(cnt)
      xcell2_2b(cnt) = 0.5d0 * xcell2b(cnt)
      ycell2_2b(cnt) = 0.5d0 * ycell2b(cnt)
      zcell2_2b(cnt) = 0.5d0 * zcell2b(cnt)
c
c     check the total number of replicated unit cells
c
      ncell2b(cnt) = nx*ny*nz - 1
      if (ncell2b(cnt) .gt. maxcell) then
         write (iout,20)
   20    format (/,' REPLICA  --  Increase MAXCELL or Decrease',
     &              ' the Interaction Cutoffs')
         call fatal
      end if
c
c     assign indices to the required cell replicates
c
      ncell2b(cnt) = 0
      do k = 0, nz-1
         do j = 0, ny-1
            do i = 0, nx-1
               if (k.ne.0 .or. j.ne.0 .or. i.ne.0) then
                  ncell2b(cnt) = ncell2b(cnt) + 1
                  icell2b(1,ncell,cnt) = i
                  icell2b(2,ncell,cnt) = j
                  icell2b(3,ncell,cnt) = k
               end if
            end do
         end do
      end do
c
c     print a message indicating the number of replicates used
c
      if (debug .and. ncell2b(cnt).ne.0) then
         write (iout,30)  nx,ny,nz
   30    format (/,' REPLICA  --  Period Boundary via',i3,' x',
     &              i3,' x',i3,' Set of Cell Replicates')
      end if
      return
      end
