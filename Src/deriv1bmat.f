c
c
      module deriv1bmat
      implicit none
      real*8, allocatable :: ep1bmat(:)
      real*8, allocatable :: dep1bmat(:,:,:)
      real*8, allocatable ::  virep1bmat(:,:,:)
      real*8, allocatable :: em1bmat(:)
      real*8, allocatable :: dem1bmat(:,:,:)
      real*8, allocatable ::  virem1bmat(:,:,:)
      integer numtasks_polar1
      save 
      end
