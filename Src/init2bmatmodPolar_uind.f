      subroutine init2bmatmodPolar_uind
      use uind2bmat
      use neigh2clust
      implicit none
      integer i,j,k,j1

      do i=1,num2bsave
         do j=1,2*3*maxsizeclust
            do k=1,3
              uind2bmatmod(k,j,i)=0.0d0
            end do
         end do
      end do

      return
      end 
