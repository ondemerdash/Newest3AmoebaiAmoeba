c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module mpole  --  atomic multipoles in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     polsiz    number of multipole components at each atom
c     pollist   multipole site for each atom (0=no multipole)
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c     polaxe    local axis type for each multipole site
c
c
      module totfield
      implicit none
      real*8, allocatable :: fieldnpole(:,:)
      real*8, allocatable :: fieldpnpole(:,:)      
      character*3 embedtyp
      real*8 vxx_perm,vyx_perm,vzx_perm
      real*8 vyy_perm,vzy_perm,vzz_perm
      real*8, allocatable :: fphi_totfield(:,:)
      real*8, allocatable :: fieldnpole_rcp(:,:)
      save
      end
