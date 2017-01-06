c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine lattice  --  setup periodic boundary conditions  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "lattice" stores the periodic box dimensions and sets angle
c     values to be used in computing fractional coordinates
c
c
      subroutine lattice2bclust(cnt)
      use sizes
      use boxes
      use boxes2bclust
      use cell
      use cell2bclust
      use inform
      use iounit
      use math
      implicit none
      real*8 alpha_cos
      real*8 ar1,ar2,ar3
      real*8 br1,br2,br3
      real*8 cr1,cr2,cr3
      integer cnt
c
c
c     compute and store the half box length values
c
      xbox2_2b(cnt) = 0.5d0 * xbox2b(cnt)
      ybox2_2b(cnt) = 0.5d0 * ybox2b(cnt)
      zbox2_2b(cnt) = 0.5d0 * zbox2b(cnt)
C     WE'RE NOT DOING ANY OCTAHEDRAL CELLS FOR THE SUBSYSTEMS RIGHT NOW.
c      if (octahedron)  box34 = 0.75d0 * xbox
c
c     set replicated cell dimensions equal to the unitcell
c
      xcell2b(cnt) = xbox2b(cnt)
      ycell2b(cnt) = ybox2b(cnt)
      zcell2b(cnt) = zbox2b(cnt)
      xcell2_2b(cnt) = xbox2_2b(cnt)
      ycell2_2b(cnt) = ybox2_2b(cnt)
      zcell2_2b(cnt) = zbox2_2b(cnt)
c
c     get values needed for fractional coordinate computations
c
C    NOTE: WE HAVE ALREADY DEFINED THE PARAMETERS BELOW, AND DON'T NEED TO
C          DEFINE SPECIAL VALUES OF THESE FOR THE CLUSTER SUBSYSTEMS AT PRESENT
c      if (orthogonal .or. octahedron) then
c         alpha_cos = 0.0d0
c         beta_sin = 1.0d0
c         beta_cos = 0.0d0
c         gamma_sin = 1.0d0
c         gamma_cos = 0.0d0
c         beta_term = 0.0d0
c         gamma_term = 1.0d0
c      else if (monoclinic) then
c         alpha_cos = 0.0d0
c         beta_sin = sin(beta/radian)
c         beta_cos = cos(beta/radian)
c         gamma_sin = 1.0d0
c         gamma_cos = 0.0d0
c         beta_term = 0.0d0
c         gamma_term = beta_sin
c      else if (triclinic) then
c         alpha_cos = cos(alpha/radian)
c         beta_sin = sin(beta/radian)
c         beta_cos = cos(beta/radian)
c         gamma_sin = sin(gamma/radian)
c         gamma_cos = cos(gamma/radian)
c         beta_term = (alpha_cos - beta_cos*gamma_cos) / gamma_sin
c         gamma_term = sqrt(beta_sin**2 - beta_term**2)
c      end if
c
c     determine the volume of the parent periodic box
c
      volbox2b(cnt) = 0.0d0
C     NOTE: WE ARE ONLY DOING ORTHOGONAL BOUNDARY CONDITIONS AT PRESENT.
      
c      if (orthogonal .or. octahedron) then
         volbox2b(cnt) = xbox2b(cnt) * ybox2b(cnt) * zbox2b(cnt)
c      else if (monoclinic) then
c         volbox = beta_sin * xbox * ybox * zbox
c      else if (triclinic) then
c         volbox = (gamma_sin*gamma_term) * xbox * ybox * zbox
c      end if
c
c     compute and store real space lattice vectors as rows
c
      ar1 = xbox2b(cnt)
      ar2 = 0.0d0
      ar3 = 0.0d0
      br1 = ybox2b(cnt) * gamma_cos
      br2 = ybox2b(cnt) * gamma_sin
      br3 = 0.0d0
      cr1 = zbox2b(cnt) * beta_cos
      cr2 = zbox2b(cnt) * beta_term
      cr3 = zbox2b(cnt) * gamma_term
      lvec2b(1,1,cnt) = ar1
      lvec2b(1,2,cnt) = ar2
      lvec2b(1,3,cnt) = ar3
      lvec2b(2,1,cnt) = br1
      lvec2b(2,2,cnt) = br2
      lvec2b(2,3,cnt) = br3
      lvec2b(3,1,cnt) = cr1
      lvec2b(3,2,cnt) = cr2
      lvec2b(3,3,cnt) = cr3
c
c     compute and store reciprocal lattice vectors as columns
c
      if (volbox2b(cnt) .ne. 0.0d0) then
         recip2b(1,1,cnt) = (br2*cr3 - cr2*br3) / volbox2b(cnt)
         recip2b(2,1,cnt) = (br3*cr1 - cr3*br1) / volbox2b(cnt)
         recip2b(3,1,cnt) = (br1*cr2 - cr1*br2) / volbox2b(cnt)
         recip2b(1,2,cnt) = (cr2*ar3 - ar2*cr3) / volbox2b(cnt)
         recip2b(2,2,cnt) = (cr3*ar1 - ar3*cr1) / volbox2b(cnt)
         recip2b(3,2,cnt) = (cr1*ar2 - ar1*cr2) / volbox2b(cnt)
         recip2b(1,3,cnt) = (ar2*br3 - br2*ar3) / volbox2b(cnt)
         recip2b(2,3,cnt) = (ar3*br1 - br3*ar1) / volbox2b(cnt)
         recip2b(3,3,cnt) = (ar1*br2 - br1*ar2) / volbox2b(cnt)
      end if
c
c     volume of truncated octahedron is half of cubic parent
c
c      if (octahedron)  volbox = 0.5d0 * volbox
      return
      end
