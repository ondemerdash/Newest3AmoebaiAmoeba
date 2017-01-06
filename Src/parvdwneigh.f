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
      module parvdwneigh
      implicit none
      integer, allocatable :: vlst_recv(:,:)
      integer, allocatable :: nvlst_recv(:)
      integer, allocatable :: vlst_recv_single(:,:)
      integer nvlst_recv_single
      save
      end
