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
c      subroutine induce0a_3b(npole3b,pnum,uind,uinp,qgrid,qfac)
      subroutine induce0a_3b_totfield(npole3b,pnum,uind,uinp,qgrid3b,
     &                 planf3b,planb3b,iprime3b,ffttable3b) 
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use mpole
      use polar, only: polarity, thole, pdamp
      use polpot
      use potent
      use units
      use uprior
c      use pme, only: nfft1, nfft2, nfft3
      use pme
      use chunks
      use fft
      use totfield
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
      integer npole3b,pnum(*),l1,l3
      real*8 uind(3,npole3b),uinp(3,npole3b)
      real*8, allocatable :: minv(:)
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      real*8 qfac3b(nfft1,nfft2,nfft3)
      real*8 thetai1_3b(4,bsorder,n)
      real*8 thetai2_3b(4,bsorder,n)
      real*8 thetai3_3b(4,bsorder,n)
      integer pmetable3b(n,nchunk)
      integer igrid3b(3,n)
      integer*8 planf3b
      integer*8 planb3b
      integer iprime3b(maxprime,3)
      real*8 ffttable3b(maxtable,3)
c
c
c     zero out the induced dipoles at each site
c
c      do i = 1, npole
      do l1 = 1, npole3b
         do j = 1, 3
c            uind(j,i) = 0.0d0
c            uinp(j,i) = 0.0d0
            uind(j,l1) = 0.0d0
            uinp(j,l1) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
c      allocate (field(3,npole))
c      allocate (fieldp(3,npole))
c      allocate (udir(3,npole))
c      allocate (udirp(3,npole))

      allocate (field(3,npole3b))
      allocate (fieldp(3,npole3b))
      allocate (udir(3,npole3b))
      allocate (udirp(3,npole3b))
      allocate (minv(6*maxulst*npole3b))

         allocate (poli(npole3b))
         allocate (rsd(3,npole3b))
         allocate (rsdp(3,npole3b))
         allocate (zrsd(3,npole3b))
         allocate (zrsdp(3,npole3b))
         allocate (conj(3,npole3b))
         allocate (conjp(3,npole3b))
         allocate (vec(3,npole3b))
         allocate (vecp(3,npole3b))

      do l1=1,npole3b
         poli(l1)=0.0d0
         do j=1,3
            field(j,l1)=0.0d0
            fieldp(j,l1)=0.0d0
            uind(j,l1)=0.0d0
            uinp(j,l1)=0.0d0
            udir(j,l1)=0.0d0
            udirp(j,l1)=0.0d0
            rsd(j,l1)=0.0d0
            rsdp(j,l1)=0.0d0
            zrsd(j,l1)=0.0d0
            zrsdp(j,l1)=0.0d0
            conj(j,l1)=0.0d0
            conjp(j,l1)=0.0d0
            vec(j,l1)=0.0d0
            vecp(j,l1)=0.0d0
         end do
      end do
      
c
c     get the electrostatic field due to permanent multipoles
c
c      print*,"use_ewald in induce",use_ewald
c      print*,"use_mlist in induce",use_mlist
c
c     set induced dipoles to polarizability times direct field
c
         do l1 =1,npole3b
            i=pnum(l1)
            do j = 1, 3
               udir(j,l1) = polarity(i) * fieldnpole(j,i)
               udirp(j,l1) = polarity(i) * fieldpnpole(j,i)
               uind(j,l1) = udir(j,l1)
               uinp(j,l1) = udirp(j,l1)
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
c         print*,"poltyp=",poltyp
c
c     estimated induced dipoles from polynomial predictor
c


c
c     perform dynamic allocation of some local arrays
c

c
c     get the electrostatic field due to induced dipoles
c
         call ufield0c_3b_totfield (npole3b,pnum,field,fieldp,uind,uinp,
     &      qgrid3b,
     &      planf3b,planb3b,iprime3b,ffttable3b)

c      subroutine ufield0c_3b_totfield(npole3b,pnum,field,fieldp,uind,
c     &      uinp,qgrid3b,
c     &      planf3b,planb3b,iprime3b,ffttable3b)

c
c     set initial conjugate gradient residual and conjugate vector
c

c         do i = 1, npole
         do l1= 1,npole3b
            i=pnum(l1)
            poli(l1) = max(polmin,polarity(i))
            do j = 1, 3
c               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
c     &                       + field(j,i)
c               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
c     &                       + fieldp(j,i)
               rsd(j,l1) = (udir(j,l1)-uind(j,l1))/poli(l1)
     &                       + field(j,l1)
               rsdp(j,l1) = (udirp(j,l1)-uinp(j,l1))/poli(l1)
     &                       + fieldp(j,l1)

            end do
         end do
         mode = 'BUILD'
c         if (use_mlist) then
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c         else
            call uscale0a_3b (npole3b,pnum,mode,rsd,rsdp,zrsd,zrsdp,
     &              minv)
            mode = 'APPLY'
            call uscale0a_3b (npole3b,pnum,mode,rsd,rsdp,zrsd,zrsdp,
     &              minv)
c         end if
c         do i = 1, npole
         do l1 =1,npole3b
            i=pnum(l1)
            do j = 1, 3
c               conj(j,i) = zrsd(j,i)
c               conjp(j,i) = zrsdp(j,i)
               conj(j,l1) = zrsd(j,l1)
               conjp(j,l1) = zrsdp(j,l1)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
c         print*,"POLAR-EPS=",poleps

         do while (.not. done)
            iter = iter + 1
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  vec(j,i) = uind(j,i)
c                  vecp(j,i) = uinp(j,i)
c                  uind(j,i) = conj(j,i)
c                  uinp(j,i) = conjp(j,i)

                  vec(j,l1) = uind(j,l1)
                  vecp(j,l1) = uinp(j,l1)
                  uind(j,l1) = conj(j,l1)
                  uinp(j,l1) = conjp(j,l1)
               end do
            end do
            if (use_ewald) then
         call ufield0c_3b_totfield (npole3b,pnum,field,fieldp,uind,uinp,
     &      qgrid3b,
     &      planf3b,planb3b,iprime3b,ffttable3b) 
c            else if (use_mlist) then
c               call ufield0b (field,fieldp)
            else
               call ufield0a (npole3b,pnum,field,fieldp)
            end if
c            do i = 1, npole
            do l1 =1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  uind(j,i) = vec(j,i)
c                  uinp(j,i) = vecp(j,i)
c                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
c                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  uind(j,l1) = vec(j,l1)
                  uinp(j,l1) = vecp(j,l1)
                  vec(j,l1) = conj(j,l1)/poli(l1) - field(j,l1)
                  vecp(j,l1) = conjp(j,l1)/poli(l1) - fieldp(j,l1)

               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  a = a + conj(j,i)*vec(j,i)
c                  ap = ap + conjp(j,i)*vecp(j,i)
c                  sum = sum + rsd(j,i)*zrsd(j,i)
c                  sump = sump + rsdp(j,i)*zrsdp(j,i)

                  a = a + conj(j,l1)*vec(j,l1)
                  ap = ap + conjp(j,l1)*vecp(j,l1)
                  sum = sum + rsd(j,l1)*zrsd(j,l1)
                  sump = sump + rsdp(j,l1)*zrsdp(j,l1)

               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1)
               do j = 1, 3
c                  uind(j,i) = uind(j,i) + a*conj(j,i)
c                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
c                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
c                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                  uind(j,l1) = uind(j,l1) + a*conj(j,l1)
                  uinp(j,l1) = uinp(j,l1) + ap*conjp(j,l1)
                  rsd(j,l1) = rsd(j,l1) - a*vec(j,l1)
                  rsdp(j,l1) = rsdp(j,l1) - ap*vecp(j,l1)

               end do
            end do
c            if (use_mlist) then
c               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            else
               call uscale0a_3b (npole3b,pnum,mode,rsd,rsdp,zrsd,zrsdp,
     &                          minv)
c            end if
            b = 0.0d0
            bp = 0.0d0
c            do i = 1, npole
            do l1 =1, npole3b
               i=pnum(l1)
               do j = 1, 3
c                  b = b + rsd(j,i)*zrsd(j,i)
c                  bp = bp + rsdp(j,i)*zrsdp(j,i)
                  b = b + rsd(j,l1)*zrsd(j,l1)
                  bp = bp + rsdp(j,l1)*zrsdp(j,l1)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
c            do i = 1, npole
            do l1 = 1,npole3b
               i=pnum(l1) 
               do j = 1, 3
c                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
c                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
c                  epsd = epsd + rsd(j,i)*rsd(j,i)
c                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  conj(j,l1) = zrsd(j,l1) + b*conj(j,l1)
                  conjp(j,l1) = zrsdp(j,l1) + bp*conjp(j,l1)
                  epsd = epsd + rsd(j,l1)*rsd(j,l1)
                  epsp = epsp + rsdp(j,l1)*rsdp(j,l1)

               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp)
c            eps = debye * sqrt(eps/dble(npolar))
            eps = debye * sqrt(eps/dble(npole3b))

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
         !deallocate (poli)
         !deallocate (rsd)
         !deallocate (rsdp)
         !deallocate (zrsd)
         !deallocate (zrsdp)
         !deallocate (conj)
         !deallocate (conjp)
         !deallocate (vec)
         !deallocate (vecp)
         !deallocate (minv)
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
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (conj)
         deallocate (conjp)
         deallocate (vec)
         deallocate (vecp)
         deallocate (minv)

      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
      return
      end
c
c
c
c
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0c  --  mutual induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0c" computes the mutual electrostatic field due to
c     induced dipole moments via Ewald summation
c
c
c      subroutine ufield0c_3b (npole3b,pnum,field,fieldp,uind,uinp,
c     &      thetai1_3b,thetai2_3b,thetai3_3b,qgrid3b,qfac3b,pmetable3b,
c     &      igrid3b,planf3b,planb3b,iprime3b,ffttable3b)
      subroutine ufield0c_3b_totfield(npole3b,pnum,field,fieldp,uind,
     &      uinp,qgrid3b,
     &      planf3b,planb3b,iprime3b,ffttable3b)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar, only: polarity, thole, pdamp
c      use pme, only: nfft1,nfft2,nfft3
      use pme
      use chunks
      use fft
      implicit none
      integer i,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
c      real*8 field(3,*)
c      real*8 fieldp(3,*)
      integer npole3b,pnum(*),l1,l3
      real*8 uind(3,npole3b),uinp(3,npole3b)
      real*8 field(3,npole3b),fieldp(3,npole3b)
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
c      real*8 qfac3b(nfft1,nfft2,nfft3)
c      real*8 thetai1_3b(4,bsorder,n)
c      real*8 thetai2_3b(4,bsorder,n)
c      real*8 thetai3_3b(4,bsorder,n)
c      integer pmetable3b(n,nchunk)
c      integer igrid3b(3,n)
      integer*8 planf3b
      integer*8 planb3b
      integer iprime3b(maxprime,3)
      real*8 ffttable3b(maxtable,3)

c
c
c     zero out the electrostatic field at each site
c
c      do i = 1, npole
      do l1=1,npole3b
         do j = 1, 3
            field(j,l1) = 0.0d0
            fieldp(j,l1) = 0.0d0
         end do
      end do
c
c     get the reciprocal space part of the electrostatic field
c

c      call umutual1_3b (npole3b,pnum,field,fieldp,uind,uinp,
c     &          thetai1_3b,thetai2_3b,thetai3_3b,qgrid3b,qfac3b,
c     &         pmetable3b,igrid3b,planf3b,planb3b,iprime3b,ffttable3b)

      call umutual1_3b_totfield (npole3b,pnum,field,fieldp,uind,uinp,
     &       qgrid3b,
     &      planf3b,planb3b,iprime3b,ffttable3b)

c
c     get the real space portion of the electrostatic field
c
c      if (use_mlist) then
c         call umutual2b (field,fieldp)
c      else
         call umutual2a_3b (npole3b,pnum,field,fieldp,uind,uinp)
c      end if
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
         do j = 1, 3
c            field(j,i) = field(j,i) + term*uind(j,i)
c            fieldp(j,i) = fieldp(j,i) + term*uinp(j,i)
            field(j,l1) = field(j,l1) + term*uind(j,l1)
            fieldp(j,l1) = fieldp(j,l1) + term*uinp(j,l1)

         end do
      end do
c
c     compute the cell dipole boundary correction to the field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
            ucellp(i) = 0.0d0
         end do
         do i = 1, npole
            do j = 1, 3
               ucell(j) = ucell(j) + uind(j,i)
               ucellp(j) = ucellp(j) + uinp(j,i)
            end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
            end do
         end do
      end if
      return
      end
c
c
c
c     perform deallocation of some local arrays
c
c
c
c
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine umutual1  --  Ewald recip mutual induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
c      subroutine umutual1_3b (npole3b,pnum,field,fieldp,uind,uinp,
c     & qgrid,qfac)
c      subroutine umutual1_3b (npole3b,pnum,field,fieldp,uind,uinp,
c     &    thetai1_3b,thetai2_3b,thetai3_3b,qgrid3b,qfac3b,
c     &    pmetable3b,igrid3b,planf3b,planb3b,iprime3b,ffttable3b)
      subroutine umutual1_3b_totfield (npole3b,pnum,field,fieldp,uind,
     &      uinp,qgrid3b,
     &      planf3b,planb3b,iprime3b,ffttable3b)
      use sizes
      use boxes
      use ewald
      use math
      use mpole
c      use pme, only: nfft1,nfft2,nfft3,maxorder,bsorder,igrid,thetai1,
c     & thetai2,thetai3
      use pme
      use atoms
      use chunks
      use polar, only: polarity, thole, pdamp
      use fft
      implicit none
      integer i,j,k
      real*8 term
      real*8 a(3,3)
c      real*8 field(3,*)
c      real*8 fieldp(3,*)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
      integer npole3b,pnum(*),l1
      real*8 uind(3,npole3b),uinp(3,npole3b)
      real*8 field(3,npole3b),fieldp(3,npole3b)
c      real*8 qgrid(2,nfft1,nfft2,nfft3)
c      real*8 qfac(nfft1,nfft2,nfft3)
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
c      real*8 qfac3b(nfft1,nfft2,nfft3)
c      real*8 thetai1_3b(4,bsorder,n)
c      real*8 thetai2_3b(4,bsorder,n)
c      real*8 thetai3_3b(4,bsorder,n)
c      integer pmetable3b(n,nchunk)
c      integer igrid3b(3,n)
      integer*8 planf3b
      integer*8 planb3b
      integer iprime3b(maxprime,3)
      real*8 ffttable3b(maxtable,3)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c

      allocate (fuind(3,npole3b))
      allocate (fuinp(3,npole3b))
      allocate (fdip_phi1(10,npole3b))
      allocate (fdip_phi2(10,npole3b))
      allocate (fdip_sum_phi(20,npole3b))
      allocate (dipfield1(3,npole3b))
      allocate (dipfield2(3,npole3b))

c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do l1=1,npole3b
         i=pnum(l1)
         do k = 1, 3

            fuind(k,l1) = a(k,1)*uind(1,l1) + a(k,2)*uind(2,l1)
     &                      + a(k,3)*uind(3,l1)
            fuinp(k,l1) = a(k,1)*uinp(1,l1) + a(k,2)*uinp(2,l1)
     &                      + a(k,3)*uinp(3,l1)

         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
c      call grid_uind3b (npole3b,pnum,fuind,fuinp,qgrid3b,
c     &       thetai1_3b,thetai2_3b,thetai3_3b,pmetable3b,igrid3b)

      call grid_uind3b_513 (npole3b,pnum,fuind,fuinp,qgrid3b)
      call fftfront3b(planf3b,qgrid3b,iprime3b,ffttable3b)
c
c     complete the transformation of the PME grid
c
      qfac(1,1,1) = 0.0d0
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               !term = qfac3b(i,j,k)
               term = qfac(i,j,k)
               qgrid3b(1,i,j,k) = term * qgrid3b(1,i,j,k)
               qgrid3b(2,i,j,k) = term * qgrid3b(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback3b(planb3b,qgrid3b,iprime3b,ffttable3b)
c      LEFT OFF HERE
c      call fphi_uind3b (npole3b,pnum,fdip_phi1,fdip_phi2,fdip_sum_phi,
c     & qgrid3b,igrid3b,thetai1_3b,thetai2_3b,thetai3_3b)

      call fphi_uind3b_513 (npole3b,pnum,fdip_phi1,fdip_phi2,
     & fdip_sum_phi,qgrid3b)
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
         do k = 1, 3

            dipfield1(k,l1) = a(k,1)*fdip_phi1(2,l1)
     &                          + a(k,2)*fdip_phi1(3,l1)
     &                          + a(k,3)*fdip_phi1(4,l1)
            dipfield2(k,l1) = a(k,1)*fdip_phi2(2,l1)
     &                          + a(k,2)*fdip_phi2(3,l1)
     &                          + a(k,3)*fdip_phi2(4,l1)

         end do
      end do
c
c     increment the field at each multipole site
c
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
         do k = 1, 3
c            field(k,i) = field(k,i) - dipfield1(k,i)
c            fieldp(k,i) = fieldp(k,i) - dipfield2(k,i)
            field(k,l1) = field(k,l1) - dipfield1(k,l1)
            fieldp(k,l1) = fieldp(k,l1) - dipfield2(k,l1)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end
c
c
