      subroutine init1bmatPolarPerm
      use deriv1bmat
      use neigh2clust
      implicit none
      integer i,j,k

      do i=1,clustcount
         ep1bmat(i)=0.0d0
         em1bmat(i)=0.0d0
         do j=1,3*maxsizeclust
            do k=1,3
              dep1bmat(k,j,i)=0.0d0
              dem1bmat(k,j,i)=0.0d0
            end do
         end do
         do j=1,3
            do k=1,3
              virep1bmat(k,j,i)=0.0d0
c              virem1bmat(k,j,i)=0.0d0
            end do
         end do
      end do

      return
      end 
