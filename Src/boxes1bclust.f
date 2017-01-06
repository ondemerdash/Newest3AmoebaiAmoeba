c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module boxes  --  periodic boundary condition parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     xbox        length of a-axis of periodic box in Angstroms
c     ybox        length of b-axis of periodic box in Angstroms
c     zbox        length of c-axis of periodic box in Angstroms
c     alpha       angle between b- and c-axes of box in degrees
c     beta        angle between a- and c-axes of box in degrees
c     gamma       angle between a- and b-axes of box in degrees
c     xbox2       half of the a-axis length of periodic box
c     ybox2       half of the b-axis length of periodic box
c     zbox2       half of the c-axis length of periodic box
c     box34       three-fourths axis length of truncated octahedron
c     volbox      volume in Ang**3 of the periodic box
c     beta_sin    sine of the beta periodic box angle
c     beta_cos    cosine of the beta periodic box angle
c     gamma_sin   sine of the gamma periodic box angle
c     gamma_cos   cosine of the gamma periodic box angle
c     beta_term   term used in generating triclinic box
c     gamma_term  term used in generating triclinic box
c     lvec        real space lattice vectors as matrix rows
c     recip       reciprocal lattice vectors as matrix columns
c     orthogonal  flag to mark periodic box as orthogonal
c     monoclinic  flag to mark periodic box as monoclinic
c     triclinic   flag to mark periodic box as triclinic
c     octahedron  flag to mark box as truncated octahedron
c     spacegrp    space group symbol for the unitcell type
c
c
      module boxes1bclust
      implicit none
      real*8, allocatable :: xbox1b(:)
      real*8, allocatable :: ybox1b(:)
      real*8, allocatable :: zbox1b(:)
      real*8, allocatable :: alpha1b(:)
      real*8, allocatable :: beta1b(:)
      real*8, allocatable :: gamma1b(:)
      real*8, allocatable :: xbox2_1b(:)
      real*8, allocatable :: ybox2_1b(:)
      real*8, allocatable :: zbox2_1b(:)
      real*8, allocatable :: box34_1b(:)
      real*8, allocatable :: volbox1b(:)
      real*8, allocatable :: beta_sin1b(:)
      real*8, allocatable :: beta_cos1b(:)
      real*8, allocatable :: gamma_sin1b(:)
      real*8, allocatable :: gamma_cos1b(:)
      real*8, allocatable :: beta_term1b(:)
      real*8, allocatable :: gamma_term1b(:)
c      real*8 lvec(3,3)
      real*8, allocatable :: lvec1b(:,:,:)
c      real*8 recip(3,3)
      real*8, allocatable :: recip1b(:,:,:)
      logical, allocatable :: orthogonal1b(:)
      logical, allocatable :: monoclinic1b(:)
      logical, allocatable :: triclinic1b(:)
      logical, allocatable :: octahedron1b(:)
      character*10, allocatable :: spacegrp1b(:)
      save
      end
