c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module cell  --  replicated cell periodic boundaries  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ncell    total number of cell replicates for periodic boundaries
c     icell    offset along axes for each replicate periodic cell
c     xcell    length of the a-axis of the complete replicated cell
c     ycell    length of the b-axis of the complete replicated cell
c     zcell    length of the c-axis of the complete replicated cell
c     xcell2   half the length of the a-axis of the replicated cell
c     ycell2   half the length of the b-axis of the replicated cell
c     zcell2   half the length of the c-axis of the replicated cell
c
c
      module cell2bclust
      use sizes
      implicit none
      integer, allocatable :: ncell2b(:)
c      integer icell(3,maxcell)
      integer, allocatable :: icell2b(:,:,:)
      real*8, allocatable ::  xcell2b(:)
      real*8, allocatable ::  ycell2b(:)
      real*8, allocatable ::  zcell2b(:)
      real*8, allocatable ::  xcell2_2b(:)
      real*8, allocatable ::  ycell2_2b(:)
      real*8, allocatable ::  zcell2_2b(:)
      save
      end
