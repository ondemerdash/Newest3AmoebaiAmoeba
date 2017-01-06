c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspred3b(npole3b,cnt)
      use sizes
      use mpole
      use uprior
      use uprior3b
      implicit none
      integer i,j,k,m
      real*8 coeff,udk,upk
      real*8 amax,apmax
      real*8 b(maxualt)
      real*8 bp(maxualt)
      real*8 a(maxualt*(maxualt+1)/2)
      real*8 ap(maxualt*(maxualt+1)/2)
      real*8 c(maxualt,maxualt)
      real*8 cp(maxualt,maxualt)
      integer npole3b,cnt
c
c
c     set the Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
         do i = 1, nualt3b(cnt)
            coeff = gear(i)
            bpred3b(i,cnt) = coeff
            bpredp3b(i,cnt) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred .eq. 'ASPC') then
         do i = 1, nualt3b(cnt)
            coeff = aspc(i)
            bpred3b(i,cnt) = coeff
            bpredp3b(i,cnt) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         do k = 1, nualt3b(cnt)
            b(k) = 0.0d0
            bp(k) = 0.0d0
            do m = k, nualt3b(cnt)
               c(k,m) = 0.0d0
               cp(k,m) = 0.0d0
            end do
         end do
         do i = 1, npole3b
            do j = 1, 3
               do k = 1, nualt3b(cnt)
                  udk = udalt3b(k,j,i,cnt)
                  upk = upalt3b(k,j,i,cnt)
                  do m = k, nualt3b(cnt)
                     c(k,m) = c(k,m) + udk*udalt3b(m,j,i,cnt)
                     cp(k,m) = cp(k,m) + upk*upalt3b(m,j,i,cnt)
                  end do
               end do
            end do
         end do
         i = 0
         do k = 2, nualt3b(cnt)
            b(k-1) = c(1,k)
            bp(k-1) = cp(1,k)
            do m = k, nualt3b(cnt)
               i = i + 1
               a(i) = c(k,m)
               ap(i) = cp(k,m)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt3b(cnt) - 1
         amax = 0.0d0
         apmax = 0.0d0
         do i = 1, k*(k+1)/2
            amax = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
         if (amax .ne. 0.0d0)  call cholesky (k,a,b)
         if (apmax .ne. 0.0d0)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt3b(cnt)-1
            bpred3b(k,cnt) = b(k)
            bpredp3b(k,cnt) = bp(k)
         end do
         bpred3b(nualt3b(cnt),cnt) = 0.0d0
         bpredp3b(nualt3b(cnt),cnt) = 0.0d0
      end if
      return
      end
c
