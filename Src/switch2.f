c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine switch  --  get switching function coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "switch" sets the coeffcients used by the fifth and seventh
c     order polynomial switching functions for spherical cutoffs
c
c
      subroutine switch2 (mode)
      use sizes
      use limits2
      use nonpol
      use shunt
      implicit none
      real*8 denom,term
      real*8 off3,off4,off5
      real*8 off6,off7
      real*8 cut3,cut4,cut5
      real*8 cut6,cut7
      character*6 mode
c
c
c     get the switching window for the current potential type
c
c         print*,"mpolecut in switch",mpolecut
      if (mode(1:6) .eq. 'EWALD2') then
         off = ewaldcutshort
         cut = ewaldcutshort
      end if
c
c     test for replicate periodic boundaries at this cutoff
c
      call replica (off)
c
c
c     store the powers of the switching window cutoffs
c
      off2 = off * off
c
c
      return
      end
