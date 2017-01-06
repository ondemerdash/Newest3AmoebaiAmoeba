c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  compute the minimum image distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
      subroutine image2b (xr,yr,zr,cnt)
      use sizes
      use boxes
      use cell
      use boxes2bclust
      use cell2bclust
      implicit none
      real*8 xr,yr,zr
      integer cnt
c
c
c     for orthogonal lattice, find the desired image directly
c
c      if (orthogonal) then
         do while (abs(xr) .gt. xcell2_2b(cnt))
            xr = xr - sign(xcell2b(cnt),xr)
         end do
         do while (abs(yr) .gt. ycell2_2b(cnt))
            yr = yr - sign(ycell2b(cnt),yr)
         end do
         do while (abs(zr) .gt. zcell2_2b(cnt))
            zr = zr - sign(zcell2b(cnt),zr)
         end do
c  WE'RE NOT CONSIDERING ANY OTHER TYPES OF PERIODIC BOXES FOR THE 
C  CLUSTER SUBSYSTEMS RIGHT NOW.
c         print*,"xcell ycell zcell",xcell,ycell,zcell
c         print *,"xcell2 ycell2 zcell2",xcell2,ycell2,zcell2
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
c      else if (monoclinic) then
c         zr = zr / beta_sin
c         xr = xr - zr*beta_cos
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
c         xr = xr + zr*beta_cos
c         zr = zr * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
c      else if (triclinic) then
c         zr = zr / gamma_term
c         yr = (yr - zr*beta_term) / gamma_sin
c         xr = xr - yr*gamma_cos - zr*beta_cos
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
c         xr = xr + yr*gamma_cos + zr*beta_cos
c         yr = yr*gamma_sin + zr*beta_term
c         zr = zr * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
c      else if (octahedron) then
c         do while (abs(xr) .gt. xbox2)
c            xr = xr - sign(xbox,xr)
c         end do
c         do while (abs(yr) .gt. ybox2)
c            yr = yr - sign(ybox,yr)
c         end do
c         do while (abs(zr) .gt. zbox2)
c            zr = zr - sign(zbox,zr)
c         end do
c         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
c            xr = xr - sign(xbox2,xr)
c            yr = yr - sign(ybox2,yr)
c            zr = zr - sign(zbox2,zr)
c         end if
c      end if
      return
      end
c
