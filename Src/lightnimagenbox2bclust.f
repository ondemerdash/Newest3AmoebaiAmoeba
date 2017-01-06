c     ##                                                             ##
c     ##  subroutine lightn  --  method of lights for list building  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "lightn" computes the set of nearest neighbor interactions
c     using the method of lights algorithm
c
c     note this is a special version for neighbor list generation
c     which includes each pair in both directions, (A,B) and (B,A)
c
c     literature reference:
c
c     F. Sullivan, R. D. Mountain and J. O'Connell, "Molecular
c     Dynamics on Vector Computers", Journal of Computational
c     Physics, 61, 138-153 (1985)
c
c
      subroutine lightnbox2bclust (cutoff,nsite,xsort,ysort,zsort,cnt)
      use sizes
      use bound
      use boxes
      use cell
      use boxes2bclust
      use cell2bclust
      use iounit
      use light
      implicit none
      integer i,j,k
      integer nsite
      integer extent
      real*8 cutoff,box
      real*8 xcut,ycut,zcut
      real*8 xsort(*)
      real*8 ysort(*)
      real*8 zsort(*)
      real*8, allocatable :: xfrac(:)
      real*8, allocatable :: yfrac(:)
      real*8, allocatable :: zfrac(:)
      integer cnt
c
c
c     truncated octahedron periodicity is not handled at present
c
      if (use_bounds) then
         if (octahedron) then
            write (iout,10)
   10       format (/,' LIGHTS  --  Truncated Octahedron not',
     &                 ' Supported by Method of Lights')
            call fatal
         end if
      end if
c
c     set the light width based on input distance cutoff
c
      xcut = cutoff
      ycut = cutoff
      zcut = cutoff
      if (use_bounds) then
         if (monoclinic) then
            zcut = zcut / beta_sin
            xcut = xcut + zcut*abs(beta_cos)
         else if (triclinic) then
            zcut = zcut / gamma_term
            ycut = (ycut + zcut*abs(beta_term)) / gamma_sin
            xcut = xcut + ycut*abs(gamma_cos) + zcut*abs(beta_cos)
         end if
         xcut = min(xcut,xcell2_2b(cnt))
         ycut = min(ycut,ycell2_2b(cnt))
         zcut = min(zcut,zcell2_2b(cnt))
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xfrac(nsite))
      allocate (yfrac(nsite))
      allocate (zfrac(nsite))
c
c     find fractional coordinates for the unitcell atoms
c
      if (use_bounds) then
         if (orthogonal) then
            do i = 1, nsite
               zfrac(i) = zsort(i)
               yfrac(i) = ysort(i)
               xfrac(i) = xsort(i)
            end do
         else if (monoclinic) then
            do i = 1, nsite
               zfrac(i) = zsort(i) / beta_sin
               yfrac(i) = ysort(i)
               xfrac(i) = xsort(i) - zfrac(i)*beta_cos
            end do
         else if (triclinic) then
            do i = 1, nsite
               zfrac(i) = zsort(i) / gamma_term
               yfrac(i) = (ysort(i) - zfrac(i)*beta_term) / gamma_sin
               xfrac(i) = xsort(i) - yfrac(i)*gamma_cos
     &                       - zfrac(i)*beta_cos
            end do
         end if
      end if
c
c     use images to move coordinates into periodic cell
c
      if (use_bounds) then
         do i = 1, nsite
            xsort(i) = xfrac(i)
            ysort(i) = yfrac(i)
            zsort(i) = zfrac(i)
            do while (abs(xsort(i)) .gt. xcell2_2b(cnt))
               xsort(i) = xsort(i) - sign(xcell2b(cnt),xsort(i))
            end do
            do while (abs(ysort(i)) .gt. ycell2_2b(cnt))
               ysort(i) = ysort(i) - sign(ycell2b(cnt),ysort(i))
            end do
            do while (abs(zsort(i)) .gt. zcell2_2b(cnt))
               zsort(i) = zsort(i) - sign(zcell2b(cnt),zsort(i))
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xfrac)
      deallocate (yfrac)
      deallocate (zfrac)
c
c     perform dynamic allocation of some global arrays
c
      nlight = nsite
      extent = 0
      if (allocated(rgx))  extent = size(rgx)
      if (extent .lt. nlight) then
         if (allocated(kbx))  deallocate (kbx)
         if (allocated(kby))  deallocate (kby)
         if (allocated(kbz))  deallocate (kbz)
         if (allocated(kex))  deallocate (kex)
         if (allocated(key))  deallocate (key)
         if (allocated(kez))  deallocate (kez)
         if (allocated(locx))  deallocate (locx)
         if (allocated(locy))  deallocate (locy)
         if (allocated(locz))  deallocate (locz)
         if (allocated(rgx))  deallocate (rgx)
         if (allocated(rgy))  deallocate (rgy)
         if (allocated(rgz))  deallocate (rgz)
         allocate (kbx(nsite))
         allocate (kby(nsite))
         allocate (kbz(nsite))
         allocate (kex(nsite))
         allocate (key(nsite))
         allocate (kez(nsite))
         allocate (locx(nlight))
         allocate (locy(nlight))
         allocate (locz(nlight))
         allocate (rgx(nlight))
         allocate (rgy(nlight))
         allocate (rgz(nlight))
      end if
c
c     sort the coordinate components into ascending order
c
      call sort2 (nlight,xsort,locx)
      call sort2 (nlight,ysort,locy)
      call sort2 (nlight,zsort,locz)
c
c     index the position of each atom in the sorted coordinates
c
      do i = 1, nlight
         rgx(locx(i)) = i
         rgy(locy(i)) = i
         rgz(locz(i)) = i
      end do
c
c     find the negative x-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locx(i)
         do while (xsort(i)-xsort(j)+box .le. xcut)
            if (j .eq. 1) then
               if (use_bounds) then
                  j = nlight + 1
                  box = xcell2b(cnt)
               end if
            end if
            j = j - 1
            if (j .lt. 1)  goto 20
         end do
   20    continue
         j = j + 1
         if (j .gt. nlight) then
            j = 1
            box = 0.0d0
         end if
         kbx(k) = j
      end do
c
c     find the positive x-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locx(i)
         do while (xsort(j)-xsort(i)+box .lt. xcut)
            if (j .eq. nlight) then
               if (use_bounds) then
                  j = 0
                  box = xcell2b(cnt)
               end if
            end if
            j = j + 1
            if (j .gt. nlight)  goto 30
         end do
   30    continue
         j = j - 1
         if (j .lt. 1) then
            j = nlight
            box = 0.0d0
         end if
         kex(k) = j
      end do
c
c     find the negative y-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locy(i)
         do while (ysort(i)-ysort(j)+box .le. ycut)
            if (j .eq. 1) then
               if (use_bounds) then
                  j = nlight + 1
                  box = ycell2b(cnt)
               end if
            end if
            j = j - 1
            if (j .lt. 1)  goto 40
         end do
   40    continue
         j = j + 1
         if (j .gt. nlight) then
            j = 1
            box = 0.0d0
         end if
         kby(k) = j
      end do
c
c     find the positive y-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locy(i)
         do while (ysort(j)-ysort(i)+box .lt. ycut)
            if (j .eq. nlight) then
               if (use_bounds) then
                  j = 0
                  box = ycell2b(cnt)
               end if
            end if
            j = j + 1
            if (j .gt. nlight)  goto 50
         end do
   50    continue
         j = j - 1
         if (j .lt. 1) then
            j = nlight
            box = 0.0d0
         end if
         key(k) = j
      end do
c
c     find the negative z-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locz(i)
         do while (zsort(i)-zsort(j)+box .le. zcut)
            if (j .eq. 1) then
               if (use_bounds) then
                  j = nlight + 1
                  box = zcell2b(cnt)
               end if
            end if
            j = j - 1
            if (j .lt. 1)  goto 60
         end do
   60    continue
         j = j + 1
         if (j .gt. nlight) then
            j = 1
            box = 0.0d0
         end if
         kbz(k) = j
      end do
c
c     find the positive z-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locz(i)
         do while (zsort(j)-zsort(i)+box .lt. zcut)
            if (j .eq. nlight) then
               if (use_bounds) then
                  j = 0
                  box = zcell2b(cnt)
               end if
            end if
            j = j + 1
            if (j .gt. nlight)  goto 70
         end do
   70    continue
         j = j - 1
         if (j .lt. 1) then
            j = nlight
            box = 0.0d0
         end if
         kez(k) = j
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine imagen  --  neighbor minimum image distance  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "imagen" takes the components of pairwise distance between
c     two points and converts to the components of the minimum
c     image distance
c
c     note this is a fast version for neighbor list generation
c     which only returns the correct component magnitudes
c
c
      subroutine imagen2b (xr,yr,zr,cnt)
      use boxes
      use boxes2bclust
      implicit none
      real*8 xr,yr,zr
      integer cnt
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = abs(xr)
         yr = abs(yr)
         zr = abs(zr)
         if (xr .gt. xbox2_2b(cnt))  xr = xr - xbox2b(cnt)
         if (yr .gt. ybox2_2b(cnt))  yr = yr - ybox2b(cnt)
         if (zr .gt. zbox2_2b(cnt))  zr = zr - zbox2b(cnt)
C  WE'RE NOT USING ANY OTHER KIND OF PERIODIC BOX FOR THE CLUSTER 
C  SUBSYSTEMS AT THIS TIME.
c
c     for monoclinic lattice, convert "xr" and "zr" specially
c
c      else if (monoclinic) then
c         zr = zr / beta_sin
c         yr = abs(yr)
c         xr = xr - zr*beta_cos
c         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
c         if (yr .gt. ybox2)  yr = yr - ybox
c         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
c         xr = xr + zr*beta_cos
c         zr = zr * beta_sin
c
c     for triclinic lattice, use general conversion equations
c
c      else if (triclinic) then
c         zr = zr / gamma_term
c         yr = (yr - zr*beta_term) / gamma_sin
c         xr = xr - yr*gamma_cos - zr*beta_cos
c         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
c         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
c         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
c         xr = xr + yr*gamma_cos + zr*beta_cos
c         yr = yr*gamma_sin + zr*beta_term
c         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
c      else if (octahedron) then
c         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
c         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
c         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
c         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
c            xr = xr - sign(xbox2,xr)
c            yr = yr - sign(ybox2,yr)
c            zr = zr - sign(zbox2,zr)
c         end if
      end if
      return
      end
