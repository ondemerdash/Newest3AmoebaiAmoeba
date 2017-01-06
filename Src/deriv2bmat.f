c
c
      module deriv2bmat
      implicit none
      real*8, allocatable :: dep2bmat(:,:,:,:)
      real*8, allocatable :: ep2bmat(:,:)
      real*8, allocatable :: virep2bmat(:,:,:,:)
      real*8, allocatable :: dem2bmat(:,:,:,:)
      real*8, allocatable :: em2bmat(:,:)      
      real*8, allocatable :: dep2bmatmod(:,:,:)
      real*8, allocatable :: ep2bmatmod(:)
      real*8, allocatable :: virep2bmatmod(:,:,:)
      save 
      end
