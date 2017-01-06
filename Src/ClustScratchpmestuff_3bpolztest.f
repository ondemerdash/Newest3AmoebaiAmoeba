c
c
      subroutine cmp_to_fmp3bclust1b (npole3b,pnum,cmp,fmp,cnt)
      use sizes
      use mpole
      implicit none
      integer i,j,k,npole3b,pnum(*),l1
      real*8 ctf(10,10)
      real*8 cmp(10,*)
      real*8 fmp(10,*)
      integer cnt
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac1bclust (ctf,cnt)
c
c     apply the transformation to get the fractional multipoles
c
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
c         fmp(1,i) = ctf(1,1) * cmp(1,i)
         fmp(1,l1) = ctf(1,1) * cmp(1,l1)
         do j = 2, 4
c            fmp(j,i) = 0.0d0
            fmp(j,l1) = 0.0d0
            do k = 2, 4
c               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
               fmp(j,l1) = fmp(j,l1) + ctf(j,k)*cmp(k,l1)
            end do
         end do
         do j = 5, 10
c            fmp(j,i) = 0.0d0
            fmp(j,l1) = 0.0d0
            do k = 5, 10
c               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
               fmp(j,l1) = fmp(j,l1) + ctf(j,k)*cmp(k,l1)
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_to_cphi  --  transformation of potential  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_to_cphi" transforms the reciprocal space potential from
c     fractional to Cartesian coordinates
c
c
      subroutine fphi_to_cphi3bclust1b (npole3b,pnum,fphi,cphi,cnt)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ftc(10,10)
      real*8 cphi(10,*)
      real*8 fphi(20,*)
      integer npole3b,pnum(*),l1,cnt
c
c
c     find the matrix to convert fractional to Cartesian
c
      call frac_to_cart1bclust (ftc,cnt)
c
c     apply the transformation to get the Cartesian potential
c
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
c         cphi(1,i) = ftc(1,1) * fphi(1,i)
         cphi(1,l1) = ftc(1,1) * fphi(1,l1)
         do j = 2, 4
c            cphi(j,i) = 0.0d0
            cphi(j,l1) = 0.0d0
            do k = 2, 4
c               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               cphi(j,l1) = cphi(j,l1) + ftc(j,k)*fphi(k,l1)
            end do
         end do
         do j = 5, 10
c            cphi(j,i) = 0.0d0
            cphi(j,l1) = 0.0d0
            do k = 5, 10
c               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               cphi(j,l1) = cphi(j,l1) + ftc(j,k)*fphi(k,l1)
            end do
         end do
      end do
      return
      end
c
c
c
c
c
c
c
      subroutine cmp_to_fmp3bclust2b (npole3b,pnum,cmp,fmp,cnt)
      use sizes
      use mpole
      implicit none
      integer i,j,k,npole3b,pnum(*),l1
      real*8 ctf(10,10)
      real*8 cmp(10,*)
      real*8 fmp(10,*)
      integer cnt
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac2bclust (ctf,cnt)
c
c     apply the transformation to get the fractional multipoles
c
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
c         fmp(1,i) = ctf(1,1) * cmp(1,i)
         fmp(1,l1) = ctf(1,1) * cmp(1,l1)
         do j = 2, 4
c            fmp(j,i) = 0.0d0
            fmp(j,l1) = 0.0d0
            do k = 2, 4
c               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
               fmp(j,l1) = fmp(j,l1) + ctf(j,k)*cmp(k,l1)
            end do
         end do
         do j = 5, 10
c            fmp(j,i) = 0.0d0
            fmp(j,l1) = 0.0d0
            do k = 5, 10
c               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
               fmp(j,l1) = fmp(j,l1) + ctf(j,k)*cmp(k,l1)
            end do
         end do
      end do
      return
      end


c
c     "fphi_to_cphi" transforms the reciprocal space potential from
c     fractional to Cartesian coordinates
c
c
      subroutine fphi_to_cphi3bclust2b (npole3b,pnum,fphi,cphi,cnt)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ftc(10,10)
      real*8 cphi(10,*)
      real*8 fphi(20,*)
      integer npole3b,pnum(*),l1,cnt
c
c
c     find the matrix to convert fractional to Cartesian
c
      call frac_to_cart2bclust (ftc,cnt)
c
c     apply the transformation to get the Cartesian potential
c
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
c         cphi(1,i) = ftc(1,1) * fphi(1,i)
         cphi(1,l1) = ftc(1,1) * fphi(1,l1)
         do j = 2, 4
c            cphi(j,i) = 0.0d0
            cphi(j,l1) = 0.0d0
            do k = 2, 4
c               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               cphi(j,l1) = cphi(j,l1) + ftc(j,k)*fphi(k,l1)
            end do
         end do
         do j = 5, 10
c            cphi(j,i) = 0.0d0
            cphi(j,l1) = 0.0d0
            do k = 5, 10
c               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               cphi(j,l1) = cphi(j,l1) + ftc(j,k)*fphi(k,l1)
            end do
         end do
      end do
      return
      end

