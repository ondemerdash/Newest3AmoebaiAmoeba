c     ##  subroutine cutoffs  --  set distance and Hessian cutoffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cutoffs" initializes and stores spherical energy cutoff
c     distance windows, Hessian element and Ewald sum cutoffs,
c     and the pairwise neighbor generation method
c
c
      subroutine cutoffs_vdw_nolist
      use sizes
      use atoms
      use bound
      use hescut
      use keys
      use limits
      use neigh
      use polpot
      use tarray
      implicit none
      integer i,next
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
      if (use_bounds) then
         vdwcut = 9.0d0
         chgcut = 9.0d0
         dplcut = 9.0d0
c         mpolecut = 9.0d0
         mpolecut = big

      else
         vdwcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
      end if
      ewaldcut = 7.0d0
      usolvcut = 4.5d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffers
c
      vdwtaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
      hesscut = 0.0d0
      lbuffer = 2.0d0
      pbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      truncate = .false.
      use_lights = .false.
      use_list = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      use_ulist = .false.
      dovlst = .true.
      doclst = .true.
      domlst = .true.
      doulst = .true.
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:6) .eq. 'EWALD ') then
            use_ewald = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
c
c     get cutoff for preconditioner of dipole solver
c
         else if (keyword(1:14) .eq. 'USOLVE-CUTOFF ') then
            read (string,*,err=10,end=10)  usolvcut
c
c     get values for the tapering style and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
            use_list = .true.
            use_vlist = .true.
            use_clist = .true.
            use_mlist = .true.
            use_ulist = .true.
         else if (keyword(1:9) .eq. 'VDW-LIST ') then
            use_list = .true.
            use_vlist = .true.
         else if (keyword(1:9) .eq. 'CHG-LIST ') then
            use_list = .true.
            use_clist = .true.
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.
c            use_ulist = .true.
c
c     get cutoff for the magnitude of Hessian elements
c
         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
c
c     get buffer width for use with pairwise neighbor lists
c
         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10)  lbuffer
         else if (keyword(1:14) .eq. 'USOLVE-BUFFER ') then
            read (string,*,err=10,end=10)  pbuffer
         end if
   10    continue
      end do
c
c     preconditioner list only needed for mutual polarization
c
      if (poltyp .ne. 'MUTUAL')  use_ulist = .false.
      if (use_list)  usolvcut = usolvcut - pbuffer
c
c     apply any Ewald cutoff to charge and multipole terms
c
c      if (use_ewald) then
c         chgcut = ewaldcut
c         mpolecut = ewaldcut
c      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
      end if
c
c     set buffer region limits for pairwise neighbor lists
c
      lbuf2 = (0.5d0*lbuffer)**2
      pbuf2 = (0.5d0*pbuffer)**2
      vbuf2 = (vdwcut+lbuffer)**2
      cbuf2 = (chgcut+lbuffer)**2
c      mbuf2 = (mpolecut+lbuffer)**2
      mbuf2 = (ewaldcut+lbuffer)**2
      ubuf2 = (usolvcut+pbuffer)**2
      vbufx = (vdwcut+2.0d0*lbuffer)**2
      cbufx = (chgcut+2.0d0*lbuffer)**2
c      mbufx = (mpolecut+2.0d0*lbuffer)**2
      mbufx = (ewaldcut+2.0d0*lbuffer)**2
      ubufx = (usolvcut+2.0d0*pbuffer)**2
c
c     perform dynamic allocation of some global arrays
c
c      if (use_vlist) then
c         if (.not.allocated(nvlst))  allocate (nvlst(n))
c         if (.not.allocated(vlst))  allocate (vlst(maxvlst,n))
c         if (.not.allocated(xvold))  allocate (xvold(n))
c         if (.not.allocated(yvold))  allocate (yvold(n))
c         if (.not.allocated(zvold))  allocate (zvold(n))
c      end if
c      if (use_clist .or. use_mlist) then
c         if (.not.allocated(nelst))  allocate (nelst(n))
c         if (.not.allocated(elst))  allocate (elst(maxelst,n))
c      end if
c      if (use_clist) then
c         if (.not.allocated(xcold))  allocate (xcold(n))
c         if (.not.allocated(ycold))  allocate (ycold(n))
c         if (.not.allocated(zcold))  allocate (zcold(n))
c      end if
c      if (use_mlist) then
c         if (.not.allocated(xmold))  allocate (xmold(n))
c         if (.not.allocated(ymold))  allocate (ymold(n))
c         if (.not.allocated(zmold))  allocate (zmold(n))
c         if (poltyp .eq. 'MUTUAL') then
c            if (.not.allocated(tindex))  allocate (tindex(2,n*maxelst))
c            if (.not.allocated(tdipdip))
c     &         allocate (tdipdip(6,n*maxelst))
c         end if
c      end if
c      if (use_ulist) then
c         if (.not.allocated(nulst))  allocate (nulst(n))
c         if (.not.allocated(ulst))  allocate (ulst(maxulst,n))
c         if (.not.allocated(xuold))  allocate (xuold(n))
c         if (.not.allocated(yuold))  allocate (yuold(n))
c         if (.not.allocated(zuold))  allocate (zuold(n))
c      end if
      return
      end

c
c
      subroutine cutoffs_vdw_parallel
      use sizes
      use atoms
      use bound
      use hescut
      use keys
      use limits
      use neigh
      use polpot
      use tarray
      implicit none
      integer i,next
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
c      if (use_bounds) then
c         vdwcut = 9.0d0
c         chgcut = 9.0d0
c         dplcut = 9.0d0
c         mpolecut = 9.0d0
c         mpolecut = big
c
c      else
c         vdwcut = big
c         chgcut = big
c         dplcut = big
c         mpolecut = big
c      end if
c      ewaldcut = 7.0d0
c      usolvcut = 4.5d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffers
c
c      vdwtaper = 0.90d0
c      chgtaper = 0.65d0
c      dpltaper = 0.75d0
c      mpoletaper = 0.65d0
c      hesscut = 0.0d0
c      lbuffer = 2.0d0
c      pbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
c      use_ewald = .false.
c      truncate = .false.
c      use_lights = .false.
c      use_list = .false.
c      use_vlist = .false.
c      use_clist = .false.
c      use_mlist = .false.
c      use_ulist = .false.
c      dovlst = .true.
c      doclst = .true.
c      domlst = .true.
c      doulst = .true.
c
c     search the keywords for various cutoff parameters
c
c      do i = 1, nkey
c         next = 1
c         record = keyline(i)
c         call gettext (record,keyword,next)
c         call upcase (keyword)
c         string = record(next:120)
c
c     get values related to use of Ewald summation
c
c         if (keyword(1:6) .eq. 'EWALD ') then
c            use_ewald = .true.
c         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
c            read (string,*,err=10,end=10)  ewaldcut
c
c     get cutoff for preconditioner of dipole solver
c
c         else if (keyword(1:14) .eq. 'USOLVE-CUTOFF ') then
c            read (string,*,err=10,end=10)  usolvcut
c
c     get values for the tapering style and neighbor method
c
c         else if (keyword(1:9) .eq. 'TRUNCATE ') then
c            truncate = .true.
c         else if (keyword(1:7) .eq. 'LIGHTS ') then
c            use_lights = .true.
c         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
c            use_list = .true.
c            use_vlist = .true.
c            use_clist = .true.
c            use_mlist = .true.
c            use_ulist = .true.
c         else if (keyword(1:9) .eq. 'VDW-LIST ') then
c            use_list = .true.
c            use_vlist = .true.
c         else if (keyword(1:9) .eq. 'CHG-LIST ') then
c            use_list = .true.
c            use_clist = .true.
c         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
c            use_list = .true.
c            use_mlist = .true.
c            use_ulist = .true.
c
c     get cutoff for the magnitude of Hessian elements
c
c         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
c            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
c         else if (keyword(1:7) .eq. 'CUTOFF ') then
c            read (string,*,err=10,end=10)  value
c            vdwcut = value
c            chgcut = value
c            dplcut = value
c            mpolecut = value
c            ewaldcut = value
c         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
c            read (string,*,err=10,end=10)  vdwcut
c         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
c            read (string,*,err=10,end=10)  chgcut
c         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
c            read (string,*,err=10,end=10)  dplcut
c         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
c            read (string,*,err=10,end=10)  mpolecut
c
c     get distance for initialization of energy switching
c
c         else if (keyword(1:6) .eq. 'TAPER ') then
c            read (string,*,err=10,end=10)  value
c            vdwtaper = value
c            chgtaper = value
c            dpltaper = value
c            mpoletaper = value
c         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
c            read (string,*,err=10,end=10)  vdwtaper
c         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
c            read (string,*,err=10,end=10)  chgtaper
c         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
c            read (string,*,err=10,end=10)  dpltaper
c         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
c            read (string,*,err=10,end=10)  mpoletaper
c
c     get buffer width for use with pairwise neighbor lists
c
c         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
c            read (string,*,err=10,end=10)  lbuffer
c         else if (keyword(1:14) .eq. 'USOLVE-BUFFER ') then
c            read (string,*,err=10,end=10)  pbuffer
c         end if
c   10    continue
c      end do
c
c     preconditioner list only needed for mutual polarization
c
c      if (poltyp .ne. 'MUTUAL')  use_ulist = .false.
c      if (use_list)  usolvcut = usolvcut - pbuffer
c
c     apply any Ewald cutoff to charge and multipole terms
c
c      if (use_ewald) then
c         chgcut = ewaldcut
c         mpolecut = ewaldcut
c      end if
c
c     convert any tapering percentages to absolute distances
c
c      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
c      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
c      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
c      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
c
c     apply truncation cutoffs if they were requested
c
c      if (truncate) then
c         vdwtaper = big
c         chgtaper = big
c         dpltaper = big
c         mpoletaper = big
c      end if
c
c     set buffer region limits for pairwise neighbor lists
c
c      lbuf2 = (0.5d0*lbuffer)**2
c      pbuf2 = (0.5d0*pbuffer)**2
c      vbuf2 = (vdwcut+lbuffer)**2
c      cbuf2 = (chgcut+lbuffer)**2
c      mbuf2 = (mpolecut+lbuffer)**2
c      mbuf2 = (ewaldcut+lbuffer)**2
c      ubuf2 = (usolvcut+pbuffer)**2
c      vbufx = (vdwcut+2.0d0*lbuffer)**2
c      cbufx = (chgcut+2.0d0*lbuffer)**2
c      mbufx = (mpolecut+2.0d0*lbuffer)**2
c      mbufx = (ewaldcut+2.0d0*lbuffer)**2
c      ubufx = (usolvcut+2.0d0*pbuffer)**2
c
c     perform dynamic allocation of some global arrays
c
      if (use_vlist) then
         if (.not.allocated(nvlst))  allocate (nvlst(n))
         if (.not.allocated(vlst))  allocate (vlst(maxvlst,n))
         if (.not.allocated(xvold))  allocate (xvold(n))
         if (.not.allocated(yvold))  allocate (yvold(n))
         if (.not.allocated(zvold))  allocate (zvold(n))
      end if
c      if (use_clist .or. use_mlist) then
c         if (.not.allocated(nelst))  allocate (nelst(n))
c         if (.not.allocated(elst))  allocate (elst(maxelst,n))
c      end if
c      if (use_clist) then
c         if (.not.allocated(xcold))  allocate (xcold(n))
c         if (.not.allocated(ycold))  allocate (ycold(n))
c         if (.not.allocated(zcold))  allocate (zcold(n))
c      end if
c      if (use_mlist) then
c         if (.not.allocated(xmold))  allocate (xmold(n))
c         if (.not.allocated(ymold))  allocate (ymold(n))
c         if (.not.allocated(zmold))  allocate (zmold(n))
c         if (poltyp .eq. 'MUTUAL') then
c            if (.not.allocated(tindex))  allocate (tindex(2,n*maxelst))
c            if (.not.allocated(tdipdip))
c     &         allocate (tdipdip(6,n*maxelst))
c         end if
c      end if
c      if (use_ulist) then
c         if (.not.allocated(nulst))  allocate (nulst(n))
c         if (.not.allocated(ulst))  allocate (ulst(maxulst,n))
c         if (.not.allocated(xuold))  allocate (xuold(n))
c         if (.not.allocated(yuold))  allocate (yuold(n))
c         if (.not.allocated(zuold))  allocate (zuold(n))
c      end if
      return
      end

