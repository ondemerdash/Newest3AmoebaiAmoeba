c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole_omp
      use sizes
      use mpole
      implicit none
      integer i
      real*8 a(3,3)
c
c
c     rotate the atomic multipoles at each site in turn
c
!$OMP PARALLEL DO default(shared) private(i,a)
      do i = 1, npole
         call rotmat (i,a)
         call rotsite (i,a)
      end do
!$OMP END PARALLEL DO
      return
      end
c
c
