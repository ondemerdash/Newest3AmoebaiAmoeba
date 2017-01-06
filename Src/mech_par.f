c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic_parallel
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use potent
      use vdwpot
      use molcul
      implicit none
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
ccc      call bonds
ccc      call angles
ccc      call torsions
ccc      call bitors
ccc      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
      call cutoffs_parallel
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
ccc      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
ccc      call orbital
      call cutoffs2_parallel

c
c     assign bond, angle and cross term potential parameters
c
ccc      if (use_bond .or. use_strbnd .or. use_strtor .or.
ccc     &    (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
ccc      if (use_angle .or. use_strbnd .or.
ccc     &    use_angang .or. use_angtor)  call kangle
ccc      if (use_strbnd)  call kstrbnd
ccc      if (use_urey)  call kurey
ccc      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
ccc      if (use_angle .or. use_opbend)  call kopbend
ccc      if (use_angle .or. use_opdist)  call kopdist
ccc      if (use_improp)  call kimprop
ccc      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
ccc      if (use_tors .or. use_strtor .or.
ccc     &    use_angtor .or. use_tortor)  call ktors
ccc      if (use_pitors)  call kpitors
ccc      if (use_strtor)  call kstrtor
ccc      if (use_angtor)  call kangtor
ccc      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
ccc      if (use_charge .or. use_chgdpl .or.
ccc     &    use_solv)  call kcharge
ccc      if (use_dipole .or. use_chgdpl)  call kdipole
      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld)  call kmpole
      if (use_polar .or. use_solv)  call kpolar
ccc      if (use_ewald)  call kewald
c
c     assign solvation, metal, pisystem and restraint parameters
c
ccc      if (use_solv)  call ksolv
ccc      if (use_metal)  call kmetal
ccc      if (use_orbit)  call korbit
ccc      if (use_geom)  call kgeom
ccc      if (use_extra)  call kextra
c
c     set hybrid parameter values for free energy perturbation
c
ccc      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic_parallel2
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use potent
      use vdwpot
      use molcul
      implicit none
c
c
c     set the bonded connectivity lists and active atoms
c
c      call attach
c      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
ccc      call bonds
ccc      call angles
ccc      call torsions
ccc      call bitors
ccc      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     find unit cell type, lattice parameters and cutoff values
c
c      call unitcell
      call lattice
c      call polymer
c      call cutoffs_parallel
c
c     setup needed for potential energy smoothing methods
c
c      call flatten
c
c     assign atom types, classes and other atomic information
c
       ! TESTING REMOVAL
ccc      call katom
c
c     assign atoms to molecules and set the atom groups
c
ccc      call molecule
c      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
ccc      call orbital
      call cutoffs2_parallel

c
c     assign bond, angle and cross term potential parameters
c
ccc      if (use_bond .or. use_strbnd .or. use_strtor .or.
ccc     &    (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
ccc      if (use_angle .or. use_strbnd .or.
ccc     &    use_angang .or. use_angtor)  call kangle
ccc      if (use_strbnd)  call kstrbnd
ccc      if (use_urey)  call kurey
ccc      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
ccc      if (use_angle .or. use_opbend)  call kopbend
ccc      if (use_angle .or. use_opdist)  call kopdist
ccc      if (use_improp)  call kimprop
ccc      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
ccc      if (use_tors .or. use_strtor .or.
ccc     &    use_angtor .or. use_tortor)  call ktors
ccc      if (use_pitors)  call kpitors
ccc      if (use_strtor)  call kstrtor
ccc      if (use_angtor)  call kangtor
ccc      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
c      if (use_vdw .or. use_solv)  call kvdw
ccc      if (use_charge .or. use_chgdpl .or.
ccc     &    use_solv)  call kcharge
ccc      if (use_dipole .or. use_chgdpl)  call kdipole
c      if (use_mpole .or. use_polar .or.
c     &    use_solv .or. use_rxnfld)  call kmpole
c      if (use_polar .or. use_solv)  call kpolar
ccc      if (use_ewald)  call kewald
c
c     assign solvation, metal, pisystem and restraint parameters
c
ccc      if (use_solv)  call ksolv
ccc      if (use_metal)  call kmetal
ccc      if (use_orbit)  call korbit
ccc      if (use_geom)  call kgeom
ccc      if (use_extra)  call kextra
c
c     set hybrid parameter values for free energy perturbation
c
ccc      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic_parallel2_nolattice
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use potent
      use vdwpot
      use molcul
      implicit none
c
c
c     set the bonded connectivity lists and active atoms
c
c      call attach
c      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
ccc      call bonds
ccc      call angles
ccc      call torsions
ccc      call bitors
ccc      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     find unit cell type, lattice parameters and cutoff values
c
c      call unitcell
c      call lattice
c      call polymer
      call cutoffs_parallel
c
c     setup needed for potential energy smoothing methods
c
c      call flatten
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
ccc      call molecule
c      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
ccc      call orbital
      call cutoffs2_parallel

c
c     assign bond, angle and cross term potential parameters
c
c     assign out-of-plane deformation potential parameters
c
c     assign torsion and torsion cross term potential parameters
c
c     assign van der Waals and electrostatic potential parameters
c
c     assign solvation, metal, pisystem and restraint parameters
c
c
c     set hybrid parameter values for free energy perturbation
c
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
