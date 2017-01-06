c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module limits  --  interaction taper & cutoff distances  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     vdwcut      cutoff distance for van der Waals interactions
c     chgcut      cutoff distance for charge-charge interactions
c     dplcut      cutoff distance for dipole-dipole interactions
c     mpolecut    cutoff distance for atomic multipole interactions
c     vdwtaper    distance at which van der Waals switching begins
c     chgtaper    distance at which charge-charge switching begins
c     dpltaper    distance at which dipole-dipole switching begins
c     mpoletaper  distance at which atomic multipole switching begins
c     ewaldcut    cutoff distance for direct space Ewald summation
c     usolvcut    cutoff distance for dipole solver preconditioner
c     use_ewald   logical flag governing use of Ewald summation
c     use_lights  logical flag governing use of method of lights
c     use_list    logical flag governing use of any neighbor lists
c     use_vlist   logical flag governing use of vdw neighbor list
c     use_clist   logical flag governing use of charge neighbor list
c     use_mlist   logical flag governing use of multipole neighbor list
c     use_ulist   logical flag governing use of preconditioner list
c
c
      module limits2
      implicit none
      real*8 ewaldcutshort
      save
      end
