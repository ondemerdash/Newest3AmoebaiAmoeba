      subroutine init1bmat_uind
      use uind1bmatrix
      use neigh2clust
      implicit none
      integer i,j,k

      do i=1,clustcount
         do j=1,3*maxsizeclust
            do k=1,3
              uind1bmat(k,j,i)=0.0d0
            end do
         end do
      end do

      return
      end 
