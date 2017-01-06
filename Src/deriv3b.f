c
c
      module deriv3b
      implicit none
      real*8, allocatable :: dep3b(:,:)
      real*8, allocatable :: dep2b(:,:)
      real*8 ep2b,virep2b(3,3),virep3b(3,3),ep3b
      integer ntriples
      real*8, allocatable :: uind3b(:,:)
      real*8, allocatable :: dep3b_recip(:,:)
      real*8 ep3b_recip,virep3b_recip(3,3)
      real*8, allocatable :: dep3bmut(:,:)
      real*8 ep3bmut,virep3bmut(3,3)
      real*8 virep3b_recip_pt1(3,3)
      real*8, allocatable :: dep3b1(:,:)
      real*8 virep3b1(3,3)
      real*8 ep1b,virep1b(3,3)
      real*8, allocatable :: dep1b(:,:)
      real*8, allocatable :: uind1b(:,:)
      real*8, allocatable :: uind2b(:,:)
      save 
      end
