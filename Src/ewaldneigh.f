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
      module ewaldneigh
      implicit none
      integer, allocatable :: njlst(:)
      integer, allocatable :: jlst(:,:)
      integer taskcnterecip,taskcnterecip_overflo
      integer, allocatable :: starterecip(:)
      integer, allocatable :: enderecip(:)
      integer, allocatable :: jiter(:)
      integer, allocatable :: starterecip2(:)
      integer, allocatable :: enderecip2(:)
      integer, allocatable :: jiter2(:)
      logical, allocatable :: doremainder(:)
      integer, allocatable :: jlst_recv(:)
      integer, allocatable :: jlst2_recv(:)
      integer jiter_recv,jiter2_recv
      integer starterecip_recv,enderecip_recv
      integer starterecip2_recv,enderecip2_recv
      logical remainder_bcast_alloc
      save
      end
