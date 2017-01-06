c
c     Edited moments routine to only calculate dipole moment
c
      subroutine moments_onlyDipole_fixedcharge
      use sizes
      use atomid
      use atoms
      use bound
      use charge
      use dipole
      use limits
      use moment
      use mpole
      use polar
      use potent
      use rigid
      use solute
      use units
      use usage
      use deriv3b
      implicit none
      integer i,j,k
      real*8 weigh,qave
      real*8 xc,yc,zc
      real*8 xi,yi,zi,ri
      real*8 xmid,ymid,zmid
      real*8 xbnd,ybnd,zbnd
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)      
      real*8 a(3,3),b(3,3)
      real*8, allocatable :: totdipole(:,:)
c
c
c     zero out total charge, dipole and quadrupole components
c
      netchg = 0.0d0
      netdpl = 0.0d0
      xdpl = 0.0d0
      ydpl = 0.0d0
      zdpl = 0.0d0
c
c     maintain periodic boundaries and neighbor lists
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xcm(n))
      allocate (ycm(n))
      allocate (zcm(n))
      allocate (totdipole(3,npole))
c
c     find the center of mass of the set of active atoms
c
      weigh = 0.0d0
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      do i = 1, n
c         if (use(i)) then
            weigh = weigh + mass(i)
            xmid = xmid + x(i)*mass(i)
            ymid = ymid + y(i)*mass(i)
            zmid = zmid + z(i)*mass(i)
c         end if
      end do
      if (weigh .ne. 0.0d0) then
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
      end if
      do i = 1, n
         xcm(i) = x(i) - xmid
         ycm(i) = y(i) - ymid
         zcm(i) = z(i) - zmid
      end do
c
c     set the multipole moment components due to partial charges
c
c
c     set the multipole moment components due to bond dipoles
c
c
c     find atomic multipoles and induced dipoles in global frame
c
c
c     set the multipole moment components due to atomic multipoles
c
      do i = 1, npole
         k = ipole(i)
         totdipole(1,i) = rpole(2,i) + xcm(k)*rpole(1,i)
         totdipole(2,i) = rpole(3,i) + ycm(k)*rpole(1,i)
         totdipole(3,i) = rpole(4,i) + zcm(k)*rpole(1,i)
         print*,"fixedpole_x=",totdipole(1,i)
         print*,"fixedpole_y=",totdipole(2,i)
         print*,"fixedpole_z=",totdipole(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xcm)
      deallocate (ycm)
      deallocate (zcm)
c
c     convert the quadrupole from traced to traceless form
c
c
c     add the traceless atomic quadrupoles to total quadrupole
c
c
c     convert dipole to Debyes and quadrupole to Buckinghams
c

c
c     get dipole magnitude and diagonalize quadrupole tensor
c
      deallocate (totdipole)
      return
      end
