c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0a  --  conjugate gradient dipole solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0a" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c
c
      subroutine induce0a_nobranch
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 udsum,upsum
      real*8 a,ap,b,bp
      real*8 sum,sump
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
c
c     get the electrostatic field due to permanent multipoles
c
       print*,"poleps=",poleps
c      if (use_ewald) then
c         call dfield0c (field,fieldp)
c      else if (use_mlist) then
c         call dfield0b (field,fieldp)
c      else
         call dfield0a (field,fieldp)
c      end if

c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            call ulspred
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
c
c     get the electrostatic field due to induced dipoles
c
c         if (use_ewald) then
c            call ufield0c (field,fieldp)
c         else if (use_mlist) then
c            call ufield0b (field,fieldp)
c         else
            call ufield0a (field,fieldp)
c         end if
c
c     set initial conjugate gradient residual and conjugate vector
c
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                       + fieldp(j,i)
            end do
         end do
         mode = 'BUILD'
c         if (use_mlist) then
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c         else
            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
            mode = 'APPLY'
            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c         end if
         do i = 1, npole
            do j = 1, 3
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  vecp(j,i) = uinp(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
               end do
            end do
c            if (use_ewald) then
c               call ufield0c (field,fieldp)
c            else if (use_mlist) then
c               call ufield0b (field,fieldp)
c            else
               call ufield0a (field,fieldp)
c            end if
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
               end do
            end do
c            if (use_mlist) then
c               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            else
               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c            end if
            b = 0.0d0
            bp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (conj)
         deallocate (conjp)
         deallocate (vec)
         deallocate (vecp)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
      return
      end
c
c
