c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine gradient (energy,derivs)
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      use molcul
      use mpole
      use neigh3b
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /
      real*8 ep3bt,virep3bt(3,3)
      real*8, allocatable :: dep3bt(:,:)
      real*8, allocatable :: uindt(:,:)
      real*8 efindiff(3)
      integer ntript
      character*7 switchmode
      allocate (dep3bt(3,npole))
      allocate (uindt(3,npole))
      print*,"npole=",npole
       
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      eat = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
c      ep = 0.0d0
      ep3b=0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(desum))  allocate (desum(3,n))
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deaa))  allocate (deaa(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(deopd))  allocate (deopd(3,n))
         if (.not. allocated(deid))  allocate (deid(3,n))
         if (.not. allocated(deit))  allocate (deit(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(debt))  allocate (debt(3,n))
         if (.not. allocated(deat))  allocate (deat(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
         if (.not. allocated(dec))  allocate (dec(3,n))
         if (.not. allocated(decd))  allocate (decd(3,n))
         if (.not. allocated(ded))  allocate (ded(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
c         if (.not. allocated(dep))  allocate (dep(3,n))
         if (.not. allocated(dep3b))  allocate (dep3b(3,n))
         if (.not. allocated(der))  allocate (der(3,n))
         if (.not. allocated(des))  allocate (des(3,n))
         if (.not. allocated(delf))  allocate (delf(3,n))
         if (.not. allocated(deg))  allocate (deg(3,n))
         if (.not. allocated(dex))  allocate (dex(3,n))
      end if
c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            deat(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
c            dep(j,i) = 0.0d0
            dep3b(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            delf(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
         end do
      end do
c
c     zero out the virial and the intermolecular energy
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 0.0d0
            virep3b(j,i)=0.0d0
         end do
      end do
      einter = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
      print*,"Before allocs"
      if (.not.allocated(nmollst))  allocate (nmollst(nmol))
      if (.not.allocated(mollst))  allocate (mollst(max2blst,nmol))
      if (.not.allocated(nmollst2))  allocate (nmollst2(nmol))
      if (.not.allocated(mollst2))  allocate (mollst2(max2blst2,nmol))
      if (.not.allocated(xmolold))  allocate (xmolold(nmol))
      if (.not.allocated(ymolold))  allocate (ymolold(nmol))
      if (.not.allocated(zmolold))  allocate (zmolold(nmol))
      print*,"After allocs"
        call mollist2bodyOO6_8_par
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy and gradient routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_angtor)  call eangtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
c      if (use_mpole .or. use_polar)  call empole1
       call prep_pole
      print*,"nmol=",nmol
c      call UindAll2bAll3bNew2SmoothInnerloop1_2cut_wskin(1,nmol,
c     &            ep3bt,virep3bt,dep3bt,ntript,0,uindt)
c      call All2bAll3bNew2SmoothInnerloop1_2cut_wskin(
c     &  1,nmol,ep3bt,virep3bt,dep3bt,0)
       
      if(use_ewald) then
        call bspline_fill
        call table_fill
c        call  LoadBalSmoothInnerloop1_ew(ep3bt,virep3bt,
c     &  dep3bt)
      switchmode = 'EWALD3B'
        call switch (switchmode)
c          call NoCutSmoothInnerloop2body_ew(ep3bt,virep3bt,dep3bt,
c     &          1,nmol,0)
          
c          call NoCutSmoothInnerloop3body_ew (ep3bt,virep3bt,dep3bt,
c     &          ntript,1,nmol,0)

      else
            switchmode = 'MPOLE'
            call switch (switchmode)

         call All2bAll3bNew2SmoothInnerloop1_2cut_wskin(
     &     1,nmol,ep3bt,virep3bt,dep3bt,0)
      end if

      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_solv)  call esolv1
      if (use_metal)  call emetal1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
c
c     sum up to get the total energy and first derivatives
c
c      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
c     &          + et + ept + ebt + eat + ett + ev + ec + ecd + ed
c     &          + em + ep + er + es + elf + eg + ex
        do i=1,3
           do j=1,3
              virep3b(j,i)=virep3bt(j,i)
              vir(j,i)=vir(j,i)+virep3b(j,i)
           end do
        end do

       ntriples=ntript
c       print*,"ep3bt=",ep3bt
       ep3b=ep3bt
  
       esum = eb + ea + eba + eub + eopb + et + ept  + ett
     &                   + ev + em+ep3b
      energy = esum
               print*,"em ev in verlet/nose/beeman setup",em,ev
               print*,"ep3b in verlet/nose/beeman setup",ep3b          
               print*,"eb in verlet/nose/beeman setup",eb
               print*,"ea in verlet/nose/beeman setup",ea
               print*,"eba in verlet/nose/beeman setup",eba
               print*,"eub in verlet/nose/beeman setup",eub
               print*,"eopb in verlet/nose/beeman setup",eopb
               print*,"et in verlet/nose/beeman setup",et
               print*,"ept in verlet/nose/beeman setup",ept
               print*,"ett in verlet/nose/beeman setup",ett
               print*,"Num Triples",ntriples

      do i = 1, n
         do j = 1, 3
c            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
c     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
c     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
c     &                      + det(j,i) + dept(j,i) + debt(j,i)
c     &                      + deat(j,i) + dett(j,i) + dev(j,i)
c     &                      + dec(j,i) + decd(j,i) + ded(j,i)
c     &                      + dem(j,i) + dep(j,i) + der(j,i)
c     &                      + des(j,i) + delf(j,i) + deg(j,i)
c     &                      + dex(j,i)
           dep3b(j,i)=dep3bt(j,i)
           desum(j,i)=deb(j,i) + dea(j,i) + deba(j,i) +deub(j,i)
     &                   + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)+dep3b(j,i)

            derivs(j,i) = desum(j,i)
         end do
            print*,"i dep3b_x",i,dep3b(1,i)
            print*,"i dep3b_y",i,dep3b(2,i)
            print*,"i dep3b_z",i,dep3b(3,i)
      end do
                print*,"virep3b(1,1)",virep3b(1,1)
                print*,"virep3b(1,2)",virep3b(1,2)
                print*,"virep3b(1,3)",virep3b(1,3)
                print*,"virep3b(2,1)",virep3b(2,1)
                print*,"virep3b(2,2)",virep3b(2,2)
                print*,"virep3b(2,3)",virep3b(2,3)
                print*,"virep3b(3,1)",virep3b(3,1)
                print*,"virep3b(3,2)",virep3b(3,2)
                print*,"virep3b(3,3)",virep3b(3,3)

c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      deallocate(uindt)
      deallocate(dep3bt)
      return
      end
