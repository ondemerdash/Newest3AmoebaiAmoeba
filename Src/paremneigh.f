c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module neigh  --  pairwise neighbor list indices & storage  ##
c     ##                                                              ##
c     ##################################################################
c
c
c
c
      module paremneigh
      implicit none
      integer, allocatable :: elst_recv(:,:)
      integer, allocatable :: nelst_recv(:)
      integer, allocatable :: elst_recv_single(:,:)
      integer nelst_recv_single
      save
      end
