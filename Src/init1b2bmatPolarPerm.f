      subroutine init1b2bmatPolarPerm
      use deriv1bmat
      use deriv2bmat
      use neigh2clust
      implicit none
      integer i,j,k,j1

      do i=1,clustcount
         do j1=1,clustcount
            ep2bmat(j1,i)=0.0d0
            em2bmat(j1,i)=0.0d0
            do j=1,2*3*maxsizeclust
               do k=1,3
                dep2bmat(k,j,j1,i)=0.0d0
                dem2bmat(k,j,j1,i)=0.0d0
               end do
            end do
            do j=1,3
              do k=1,3
               virep2bmat(k,j,j1,i)=0.0d0
              end do
            end do
         end do
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
