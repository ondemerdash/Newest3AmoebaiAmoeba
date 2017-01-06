      subroutine init2bmatmodPolar
      use deriv1bmat
      use deriv2bmat
      use neigh2clust
      implicit none
      integer i,j,k,j1

      do i=1,num2bsave
         ep2bmatmod(i)=0.0d0
         do j=1,2*3*maxsizeclust
            do k=1,3
              dep2bmatmod(k,j,i)=0.0d0
            end do
         end do
         do j=1,3
            do k=1,3
              virep2bmatmod(k,j,i)=0.0d0
            end do
         end do
      end do

      return
      end 
