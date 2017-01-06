c     ##  subroutine ptest  --  find pressure via finite-difference  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ptest" compares the virial-based value of dE/dV to an estimate
c     from finite-difference volume changes; also finds the isotropic
c     pressure via finite-differences
c
c     original version written by John D. Chodera, University of
c     California, Berkeley, December 2010
c
c
      subroutine ptest_parallel
      use sizes
      use atoms
      use bath
      use bound
      use boxes
      use iounit
      use units
      use virial
      implicit none
      integer i
      real*8 energy,third
      real*8 delta,step,scale
      real*8 vold,xboxold
      real*8 yboxold,zboxold
      real*8 epos,eneg
      real*8 dedv_vir,dedv_fd
      real*8 pres_vir,pres_fd
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c
c
c     set relative volume change for finite-differences
c
      if (.not. use_bounds)  return
      delta = 0.000001d0
      step = volbox * delta
      print*,"delta=",delta      
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     store original box dimensions and coordinate values
c
      xboxold = xbox
      yboxold = ybox
      zboxold = zbox
      vold = volbox
      do i = 1, n
         xold(i) = x(i)
         yold(i) = y(i)
         zold(i) = z(i)
      end do
c
c     get scale factor to reflect a negative volume change
c
      volbox = vold + step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
c
c     set new box dimensions and coordinate values
c
      xbox = xboxold * scale
      ybox = yboxold * scale
      zbox = zboxold * scale
      call lattice
      do i = 1, n
         x(i) = xold(i) * scale
         y(i) = yold(i) * scale
         z(i) = zold(i) * scale
      end do
c
c     compute potential energy for negative volume change
c
c      eneg = energy ()
c
c     get scale factor to reflect a positive volume change
c
c      volbox = vold + step
c      third = 1.0d0 / 3.0d0
c      scale = (volbox/vold)**third
c
c     set new box dimensions and coordinate values
c
c      xbox = xboxold * scale
c      ybox = yboxold * scale
c      zbox = zboxold * scale
c      call lattice
c      do i = 1, n
c         x(i) = xold(i) * scale
c         y(i) = yold(i) * scale
c         z(i) = zold(i) * scale
c      end do
c
c     compute potential energy for positive volume change
c
c      epos = energy ()
c
c     restore original box dimensions and coordinate values
c
c      xbox = xboxold
c      ybox = yboxold
c      zbox = zboxold
c      call lattice
c      do i = 1, n
c         x(i) = xold(i)
c         y(i) = yold(i)
c         z(i) = zold(i)
c      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
c
c     get virial and finite difference values of dE/dV
c
c      dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
c      dedv_fd = (epos-eneg) / (2.0d0*delta*volbox)
c      write (iout,10)  dedv_vir
c   10 format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
c      write (iout,20)  dedv_fd
c   20 format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
c
c     compute analytical and finite-difference isotropic pressure
c
c      pres_vir = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_vir)
c      pres_fd = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_fd)
c      if (kelvin .eq. 0.0d0) then
c         write (iout,30)  pres_vir
c         write (iout,40)  pres_fd
c   30    format (/,' Pressure (Analytical, 0 K) :',5x,f15.3,
c     &              ' Atmospheres')
c   40    format (' Pressure (Numerical, 0 K) :',6x,f15.3,
c     &              ' Atmospheres')
c      else
c         write (iout,50)  nint(kelvin),pres_vir
c         write (iout,60)  nint(kelvin),pres_fd
c   50    format (/,' Pressure (Analytical,',i4,' K) :',3x,f15.3,
c     &              ' Atmospheres')
c   60    format (' Pressure (Numerical,',i4,' K) :',4x,f15.3,
c     &              ' Atmospheres')
c      end if
      return
      end
