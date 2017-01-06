c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine unitcell  --  get periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "unitcell" gets the periodic boundary box size and related
c     values from an external keyword file
c
c
      subroutine unitcell1bclust(cnt,xmax,ymax,zmax)
      use sizes
      use bound
      use boxes1bclust
      use iounit
      use keys
      implicit none
      integer i,next
      real*8 boxmax
      logical nosymm
      character*20 keyword
      character*120 record
      character*120 string
      real*8 xmax,ymax,zmax
      integer cnt
c
c
c     set the default values for periodic boundary conditions
c
C     NOTE:  COMMENTED BELOW, SINCE ALREADY DEFINED PREVIOUSLY

c      use_bounds = .false.
c      use_replica = .false.

c
c     set the default values for the unitcell variables
c
C     NOTE:  COMMENTED BELOW, SINCE ALREADY DEFINED PREVIOUSLY

c      orthogonal = .false.
c      monoclinic = .false.
c      triclinic = .false.
c      octahedron = .false.
c      spacegrp = '          '
c      nosymm = .false.

c
c     get keywords containing crystal lattice dimensions
c     NOTE:  COMMENTED BELOW BECAUSE THEY HAVE EITHER ALREADY BEEN DEFINED OR 
C            ARE SPECIFIC TO THE SUBSYSTEM AT HAND.
c
c      do i = 1, nkey
c         next = 1
c         record = keyline(i)
c         call gettext (record,keyword,next)
c         call upcase (keyword)
c         string = record(next:120)
c         if (keyword(1:7) .eq. 'X-AXIS ') then
c            if (xbox .eq. 0.0d0)  read (string,*,err=10,end=10)  xbox
c         else if (keyword(1:7) .eq. 'Y-AXIS ') then
c            if (ybox .eq. 0.0d0)  read (string,*,err=10,end=10)  ybox
c         else if (keyword(1:7) .eq. 'Z-AXIS ') then
c            if (zbox .eq. 0.0d0)  read (string,*,err=10,end=10)  zbox
c         else if (keyword(1:7) .eq. 'A-AXIS ') then
c            if (xbox .eq. 0.0d0)  read (string,*,err=10,end=10)  xbox
c         else if (keyword(1:7) .eq. 'B-AXIS ') then
c            if (ybox .eq. 0.0d0)  read (string,*,err=10,end=10)  ybox
c         else if (keyword(1:7) .eq. 'C-AXIS ') then
c            if (zbox .eq. 0.0d0)  read (string,*,err=10,end=10)  zbox
c         else if (keyword(1:6) .eq. 'ALPHA ') then
c            if (alpha .eq. 0.0d0)  read (string,*,err=10,end=10)  alpha
c         else if (keyword(1:5) .eq. 'BETA ') then
c            if (beta .eq. 0.0d0)  read (string,*,err=10,end=10)  beta
c         else if (keyword(1:6) .eq. 'GAMMA ') then
c            if (gamma .eq. 0.0d0)  read (string,*,err=10,end=10)  gamma
c         else if (keyword(1:11) .eq. 'OCTAHEDRON ') then
c            octahedron = .true.
c         else if (keyword(1:11) .eq. 'SPACEGROUP ') then
c            call getword (record,spacegrp,next)
c         else if (keyword(1:12) .eq. 'NO-SYMMETRY ') then
c            nosymm = .true.
c         end if
c   10    continue
c      end do
c
c     use periodic boundary conditions if a cell was defined
c
      boxmax = max(xmax,ymax,zmax)
      xbox1b(cnt)=xmax
      ybox1b(cnt)=ymax
      zbox1b(cnt)=zmax

c      if (boxmax .ne. 0.0d0)  use_bounds = .true.
c
c     set unspecified periodic boundary box lengths and angles
c
c      if (use_bounds) then
c         if (xbox .eq. 0.0d0)  xbox = boxmax
c         if (ybox .eq. 0.0d0)  ybox = xbox
c         if (zbox .eq. 0.0d0)  zbox = xbox
c         if (alpha .eq. 0.0d0)  alpha = 90.0d0
c         if (beta .eq. 0.0d0)  beta = 90.0d0
c         if (gamma .eq. 0.0d0)  gamma = 90.0d0
c
c     determine the general periodic boundary lattice type
c
c         if (nosymm) then
c            triclinic = .true.
c         else if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
c     &               .and. gamma.eq.90.0d0) then
c            orthogonal = .true.
c         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
c            monoclinic = .true.
c         else
c            triclinic = .true.
c         end if
c      end if
c
c     check for proper use of truncated octahedron boundary
c
c      if (octahedron) then
c         if (xbox.eq.ybox .and. xbox.eq.zbox .and. orthogonal) then
c            orthogonal = .false.
c            monoclinic = .false.
c            triclinic = .false.
c         else
c            write (iout,20)
c   20       format (/,' UNITCELL  --  Truncated Octahedron',
c     &                 ' Incompatible with Defined Cell')
c            call fatal
c         end if
c      end if
      return
      end
