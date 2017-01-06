      subroutine bcast_mechanic_wlattice
      use mpidat
      use sizes
      use atoms
      use couple
      use neigh3b
      use molcul
      use usage
      use boxes
      use group
      use mpole
      use polar
      use polgrp
      use bound
      use cell
      use neigh
      implicit none
      include 'mpif.h'
      integer maxn13,maxn14,maxn15
      integer maxp11,maxp12
      integer maxp13,maxp14,ierr

         call mpi_bcast(polycut2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(use_polymer,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(listbcast,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(nmol,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(imol,2*n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(kmol,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(molmass,n,mpi_real8,master,
     &   mpi_comm_world,ierr)

          call mpi_bcast(n12,maxatm,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(n13,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(n14,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(n15,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
         maxn13 = 3 * maxval
         maxn14 = 9 * maxval
         maxn15 = 27 * maxval
          call mpi_bcast(i12,maxatm*maxval,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(i13,n*maxn13,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(i14,n*maxn14,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(i15,n*maxn15,mpi_integer,master,
     &   mpi_comm_world,ierr)

         call mpi_bcast(nuse,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(iuse,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use,n+1,mpi_logical,master,
     &   mpi_comm_world,ierr) 
          
          call mpi_bcast(use_bounds,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use_replica,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use_polymer,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(xbox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ybox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(zbox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(alpha,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(xbox2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ybox2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(zbox2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(box34,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(volbox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta_sin,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta_cos,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma_sin,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma_cos,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta_term,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma_term,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(orthogonal,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(monoclinic,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(triclinic,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(octahedron,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
c Added 1 more var to bcast
          call mpi_bcast(spacegrp,10,mpi_character,master,
     &   mpi_comm_world,ierr)     
c NOW BROADCAST THE LATTICE VARIABLES IN CELL.F
          call mpi_bcast(ncell,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(icell,3*maxcell,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(xcell,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ycell,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(zcell,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(xcell2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ycell2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(zcell2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)

          call mpi_bcast(ngrp,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use_group,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(kgrp,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(grplist,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(igrp,2*(maxgrp+1),mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(grpmass,(maxgrp+1),mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(wgrp,(maxgrp+1)*(maxgrp+1),mpi_real8,
     &   master,mpi_comm_world,ierr)

          call mpi_bcast(npole,1,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ipole,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(zaxis,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(xaxis,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(yaxis,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(rpole,n*maxpole,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(polaxe,8*n,mpi_character,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np11,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np12,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np13,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np14,n,mpi_integer,
     &   master,mpi_comm_world,ierr)

         call mpi_bcast(polarity,n,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(thole,n,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(pdamp,n,mpi_real8,
     &   master,mpi_comm_world,ierr)
      maxp11 = 150
      maxp12 = 50
      maxp13 = 50
      maxp14 = 50
         call mpi_bcast(ip11,n*maxp11,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ip12,n*maxp12,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ip13,n*maxp13,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ip14,n*maxp14,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(listsend_mpole,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
      return
      end

      subroutine bcast_mechanic
      use mpidat
      use sizes
      use atoms
      use atomid
      use couple
      use neigh3b
      use molcul
      use usage
      use boxes
      use group
      use mpole
      use polar
      use polgrp
      use bound
      use neigh
      use cobar
      use limits
      use totfield
      use aprx
      use polpot
      use neigh2clust
      use cho
      use pcg
      use inform
      implicit none
      include 'mpif.h'
      integer maxn13,maxn14,maxn15
      integer maxp11,maxp12
      integer maxp13,maxp14,ierr

      ! NEW TESTING
         call mpi_bcast(mass,n,mpi_real8,master,
     &   mpi_comm_world,ierr)
      ! END TESTING
         call mpi_bcast(longrangepoldir,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(uzepmedirpolz,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(listbcast,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(nmol,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(imol,2*n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(kmol,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(molmass,n,mpi_real8,master,
     &   mpi_comm_world,ierr)

          call mpi_bcast(n12,maxatm,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(n13,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(n14,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(n15,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
         maxn13 = 3 * maxval
         maxn14 = 9 * maxval
         maxn15 = 27 * maxval
          call mpi_bcast(i12,maxatm*maxval,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(i13,n*maxn13,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(i14,n*maxn14,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(i15,n*maxn15,mpi_integer,master,
     &   mpi_comm_world,ierr)

         call mpi_bcast(nuse,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(iuse,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use,n+1,mpi_logical,master,
     &   mpi_comm_world,ierr) 
          
          call mpi_bcast(use_bounds,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use_replica,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(use_polymer,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(xbox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ybox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(zbox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(alpha,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(xbox2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(ybox2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(zbox2,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(box34,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(volbox,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta_sin,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta_cos,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma_sin,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma_cos,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(beta_term,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(gamma_term,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(orthogonal,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(monoclinic,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(triclinic,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(octahedron,1,mpi_logical,master,
     &   mpi_comm_world,ierr)

          call mpi_bcast(ngrp,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(use_group,1,mpi_logical,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(kgrp,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(grplist,n,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(igrp,2*(maxgrp+1),mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(grpmass,(maxgrp+1),mpi_real8,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(wgrp,(maxgrp+1)*(maxgrp+1),mpi_real8,
     &   master,mpi_comm_world,ierr)

          call mpi_bcast(npole,1,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ipole,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(zaxis,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(xaxis,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(yaxis,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(rpole,n*maxpole,mpi_real8,
     &   master,mpi_comm_world,ierr)

         call mpi_bcast(pole,n*maxpole,mpi_real8,
     &   master,mpi_comm_world,ierr)

         call mpi_bcast(polaxe,8*n,mpi_character,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np11,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np12,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np13,n,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(np14,n,mpi_integer,
     &   master,mpi_comm_world,ierr)

         call mpi_bcast(polarity,n,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(thole,n,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(pdamp,n,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(usolvcut,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(poleps,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(udiag,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(poltyp,6,mpi_character,
     &   master,mpi_comm_world,ierr) 
      maxp11 = 150
      maxp12 = 50
      maxp13 = 50
      maxp14 = 50
         call mpi_bcast(ip11,n*maxp11,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ip12,n*maxp12,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ip13,n*maxp13,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ip14,n*maxp14,mpi_integer,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(listsend_mpole,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(rtapr2b_input,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(Rcut2b_input,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(cut2b_input,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(R123cut3b,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(rtapr3b,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(Cobarcut3b,1,mpi_real8,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(embedtyp,3,mpi_character,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(approxmode,9,mpi_character,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(itermode,1,mpi_integer,
     &   master,mpi_comm_world,ierr)

         call mpi_bcast(doclust,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(doclust3,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(doclust4,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(do2waterclustlist,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(usecholesky,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(inputkmeansclust,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(usepcgscf,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(ompouterloop3b,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(mlistclust,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(use_ewaldclust,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(debug,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(all2bclust,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
         call mpi_bcast(useboxclust,1,mpi_logical,
     &   master,mpi_comm_world,ierr)
      return
      end

      subroutine gradient_setup
      use bound 
      use group
      implicit none
      real*8 cutoff
      if (use_bounds .and. .not.use_group)  call bounds
      cutoff = 0.0d0
      call replica (cutoff)
      return
      end

      subroutine prep_pole
      implicit none
       call chkpole
       call rotpole
      return
      end     

      subroutine prep_pole_omp
      implicit none
       call chkpole_omp
       call rotpole_omp
      return
      end

c      subroutine alloc3b_nblist
c      use sizes
c      use atoms
c      use deriv3b
c      implicit none
c      logical first2
c      save first2
c      data first2  / .true. /
c      integer i,j
c
c      ep3b=0.0d0
c
c      if (first2) then
c         first2 = .false.
c         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
c      end if
c        do i=1,3
c          do j=1,3
c            virep3b(i,j)=0.0d0
c          end do
c        end do
c      do i = 1, n
c         do j=1,3
c            dep3b(j,i) = 0.0d0
c         end do
c      end do
c
c       call nblist_perm_2bod
c
c      return
c      end


      subroutine alloc3b_3blist
      use sizes
      use atoms
      use deriv3b
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep3b=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if
        do i=1,3
          do j=1,3
            virep3b(i,j)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
         end do
      end do
       call mollist2bodyOO6_8_par 
      return
      end


      subroutine alloc3b_vdwlist
      use sizes
      use atoms
      use deriv3b
      use limits
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep3b=0.0d0
      ntriples=0
      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if
        do i=1,3
          do j=1,3
            virep3b(i,j)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
         end do
      end do

      if(use_vlist) call vlist
      call mollist2bodyOO6_8_par
      return
      end

      subroutine alloc3b_vdwlist_nolist
      use sizes
      use atoms
      use deriv3b
      use limits
      use aprx
      implicit none
      logical first2,first3
      save first2
      save first3
      data first2  / .true. /
      data first3  / .true. /
      integer i,j

      ep3b=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
         if (.not. allocated(dep3bmut)) allocate (dep3bmut(3,n))
         if (.not. allocated(dep3b1)) allocate (dep3b1(3,n))
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
            virep3b(j,i)=0.0d0
            virep3bmut(j,i)=0.0d0
            virep3b_recip(j,i)=0.0d0
            virep3b1(j,i)=0.0d0
            virep3b_recip_pt1(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
            dep3bmut(j,i) = 0.0d0
            dep3b_recip(j,i) = 0.0d0
            dep3b1(j,i) = 0.0d0
         end do
      end do

c      if (approxmode.ne.'1BODYMODE') then
         ep3bmut=0.0d0
         ep3b_recip=0.0d0
c      end if
c       call nblist_perm_2bod
c      if(use_vlist) call vlist
c      call mollist2bodyOO6_8_par
      return
      end

      subroutine alloc3b_vdwlist_nolist_vac
      use sizes
      use atoms
      use deriv3b
      use limits
      use aprx
      implicit none
      logical first2,first3
      save first2
      save first3
      data first2  / .true. /
      data first3  / .true. /
      integer i,j

      ep3b=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
            virep3b(j,i)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
         end do
      end do

c      if (approxmode.ne.'1BODYMODE') then
c      end if
c       call nblist_perm_2bod
c      if(use_vlist) call vlist
c      call mollist2bodyOO6_8_par
      return
      end

      subroutine alloc3b_vdwlist_erecip
      use sizes
      use atoms
      use deriv3b
      use deriv
      use virial
      use energi
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep3b=0.0d0
      em = 0.0d0
      ntriples=0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
      end if
        do i=1,3
          do j=1,3
            virep3b(j,i)=0.0d0
            viremrecip(j,i) =0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
            dem(j,i) =0.0d0
         end do
      end do

c       call nblist_perm_2bod
      call vlist
      call mollist2bodyOO6_8_par
      return
      end

      subroutine alloc_erecip
      use sizes
      use atoms
      use deriv
      use virial
      use energi
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      em = 0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dem))  allocate (dem(3,n))
      end if
        do i=1,3
          do j=1,3
            viremrecip(j,i) =0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dem(j,i) =0.0d0
         end do
      end do

      return
      end

      subroutine alloc3b_vdwlist_erecip2b
      use sizes
      use atoms
      use deriv3b
      use deriv
      use virial
      use energi
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep3b=0.0d0
      em = 0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
      end if
        do i=1,3
          do j=1,3
            virep3b(j,i)=0.0d0
            viremrecip(j,i) =0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
            dem(j,i) =0.0d0
         end do
      end do

c       call nblist_perm_2bod
      print*,"In new 2b alloc"
      call vlist
c      call mollist2bodyOO6_par2bonly

      return
      end


c mollist2bodyOO6_par2bonly

      subroutine alloc2b_vdwlist
      use sizes
      use atoms
      use deriv3b
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep3b=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if
        do i=1,3
          do j=1,3
            virep3b(i,j)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep3b(j,i) = 0.0d0
         end do
      end do

c       call nblist_perm_2bod
      call vlist
      call mollist2bodyOO6_par2bonly
      return
      end

      subroutine alloc2b
      use sizes
      use atoms
      use deriv3b
      implicit none
      logical first2
      save first2
      data first2  / .true. /
      integer i,j

      ep2b=0.0d0

      if (first2) then
         first2 = .false.
         if (.not. allocated(dep2b)) allocate (dep2b(3,n))
      end if
        do i=1,3
          do j=1,3
            virep2b(i,j)=0.0d0
          end do
        end do
      do i = 1, n
         do j=1,3
            dep2b(j,i) = 0.0d0
         end do
      end do

c       call nblist_perm_2bod
c      call vlist
c      call mollist2bodyOO6_par2bonly
      return
      end


      subroutine gradient_serial (energy,derivs)
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
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

      print*,"Very beginning of gradient_serial before zeroing"

c      allocate (dep3b(3,npole))
c      ep3b2=0.0d0
c      ep3b3=0.0d0
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      em = 0.0d0

      ep3b=0.0d0
      ntriples=0
      if (first) then
         first = .false.
c         if (.not. allocated(dep))  allocate (dep(3,n))
         if (.not. allocated(desum))  allocate (desum(3,n))
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
c            virep3b2(i,j)=0.0d0
c            virep3b3(i,j)=0.0d0
            virep3b(i,j)=0.0d0

            vir(i,j)=0.0d0
c            viremreal(i,j)=0.0d0
c            virev(i,j)=0.0d0
          end do
        end do

      do i = 1, n
         do j=1,3
c            dev_total(j,i)=0.0d0
c            dep(j,i) = 0.0d0
c            dep3b2(j,i) = 0.0d0
c            dep3b3(j,i) = 0.0d0
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep3b(j,i) = 0.0d0

c            demreal(j,i) =0.0d0
         end do
      end do

      print*,"Before nblist_perm_2bod in gradient_serial"
c       call nblist_perm_2bod

c
c     call the van der Waals energy component routines
c
      print*,"Before energies in gradient_serial"
      call ebond1
      call eangle1
      call estrbnd1
      call eurey1
      call eopbend1
      call etors1
      call epitors1
      call etortor1
      call ehal1
c         if (vdwtyp .eq. 'BUFFERED-14-7') call ehal1c_new
      print*,"After vdw",ev

c        call empole1c_3b_Perm_PME
        call empole1d_3b_Perm_PME
c
c     call any miscellaneous energy component routines
c

c      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
c     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
c     &          + ep + er + es + elf + eg + ex
c      energy = esum
c      do i = 1, n
c         do j = 1, 3
c            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
c     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
c     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
c     &                      + det(j,i) + dept(j,i) + debt(j,i)
c     &                      + dett(j,i) + dev(j,i) + dec(j,i)
c     &                      + decd(j,i) + ded(j,i) + dem(j,i)
c     &                      + dep(j,i) + der(j,i) + des(j,i)
c     &                      + delf(j,i) + deg(j,i) + dex(j,i)
c         end do
c      end do

      esum =eb + ea + eba + eub + eopb + et + ept  + ett + ev + em 
      energy = esum
c      print*,"em ev",em,ev
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i) + deopb(j,i)
     &                   + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i) 
            derivs(j,i) = desum(j,i)
         end do
      end do
      return
      end

c      subroutine gradient_polar(start,last,moli1,master,ierr)
      subroutine gradient_polar(start,last,moli1)
      use mpole
      use mpidat
      use deriv3b
      implicit none
      include 'mpif.h'
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
      real*8, allocatable :: dep3bt_tot(:,:)
      real*8, allocatable :: dep3bt(:,:)
      integer start,last,moli1,i,j,ierr,ntript
      integer ntript_tot
c      integer master,ierr
      allocate (dep3bt(3,npole))
      allocate (dep3bt_tot(3,npole))
c      print*,"Beginning of gradient_polar"
                ep3bt_tot=0.0d0
c                ep3bt_tot2=0.0d0
c                ep3bt_tot3=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt_tot(j,i) = 0.0d0
c                      dep3bt_tot2(j,i) = 0.0d0
c                      dep3bt_tot3(j,i) = 0.0d0
                   end do
                end do
                do i=1,3
                   do j=1,3
                      virep3bt_tot(j,i)=0.0d0
c                      virep3bt_tot2(i,j)=0.0d0
c                      virep3bt_tot3(i,j)=0.0d0
                   end do
                end do
              ntript_tot=0

c             call Innerloop1_2cut_wskin(start,last,
c    &              ep3bt,virep3bt,dep3bt,ntript)

             ! print*,"After call to Innerloop1_2cut"
                 do i=1,3
                   do j=1,3
                 virep3bt_tot(i,j)=virep3bt_tot(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot=ep3bt_tot+ep3bt
                 ntript_tot=ntript_tot+ntript

       if(moli1.ne.0) then
c                 call Innerloop1_single_2cut_wskin(
c    &               moli1,ep3bt,virep3bt,dep3bt,ntript)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot(i,j)=virep3bt_tot(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot=ep3bt_tot+ep3bt
                 ntript_tot=ntript_tot+ntript
       end if
c       print*,"In gradient_polar",ep3bt_tot
c       print*,"In gradient_polar Before mpi_reduce ",ep3bt_tot
                  call mpi_reduce(ep3bt_tot,ep3b,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot,dep3b,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot,virep3b,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(ntript_tot,ntriples,1,mpi_integer,
     &          mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(dep3bt_tot)
      deallocate(dep3bt)
      return
      end

      subroutine gradient_serial2 (energy,derivs)
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
c      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /

      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      em = 0.0d0

c      ep3b=0.0d0

      if (first) then
         first = .false.
c         if (.not. allocated(dep))  allocate (dep(3,n))
         if (.not. allocated(desum))  allocate (desum(3,n))
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
c         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
c            virep3b2(i,j)=0.0d0
c            virep3b3(i,j)=0.0d0
c            virep3b(i,j)=0.0d0
            vir(i,j)=0.0d0
c            viremreal(i,j)=0.0d0
c            virev(i,j)=0.0d0
          end do
        end do

      do i = 1, n
         do j=1,3
c            dev_total(j,i)=0.0d0
c            dep(j,i) = 0.0d0
c            dep3b2(j,i) = 0.0d0
c            dep3b3(j,i) = 0.0d0
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dem(j,i) = 0.0d0
c            dep3b(j,i) = 0.0d0

c            demreal(j,i) =0.0d0
         end do
      end do

c       call nblist_perm_2bod

c
c     call the van der Waals energy component routines
c
      call ebond1
      call eangle1
      call estrbnd1
      call eurey1
      call eopbend1
      call etors1
      call epitors1
      call etortor1
      call ehal1
        call empole1d_3b_Perm_PME
c
c     call any miscellaneous energy component routines

      esum =eb + ea + eba + eub + eopb + et + ept  + ett + ev + em
      energy = esum
c      print*,"In gradient_serial2 em ev",em,ev
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i) + deopb(j,i)
     &                   + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i)
            derivs(j,i) = desum(j,i)
         end do
      end do

      return
      end

c      subroutine gradient_polar2(start,last,moli1)
c      use mpole
c      use mpidat
c      use deriv3b
c      implicit none
c      include 'mpif.h'
c      real*8 ep3bt,virep3bt(3,3),ep3bt_tot,virep3bt_tot(3,3)
c      real*8, allocatable :: dep3bt_tot(:,:)
c      real*8, allocatable :: dep3bt(:,:)
c      integer start,last,moli1,i,j,ierr,ntriple_tot,ntript
c      integer master1,ierr
c      allocate (dep3bt(3,npole))
c      allocate (dep3bt_tot(3,npole))
c      print*,"Beginning of gradient_polar"
c                ep3bt_tot=0.0d0
c                ep3bt_tot2=0.0d0
c                ep3bt_tot3=0.0d0
c                do i = 1, npole
c                   do j = 1, 3
c                      dep3bt_tot(j,i) = 0.0d0
c                      dep3bt_tot2(j,i) = 0.0d0
c                      dep3bt_tot3(j,i) = 0.0d0
c                   end do
c                end do
c                do i=1,3
c                   do j=1,3
c                      virep3bt_tot(i,j)=0.0d0
c                      virep3bt_tot2(i,j)=0.0d0
c                      virep3bt_tot3(i,j)=0.0d0
c                   end do
c                end do
c                ntriple_tot=0
c              call Innerloop1_NoCut2bonly(start,last,
c     &              ep3bt,virep3bt,dep3bt,ntript)
c
c              print*,"After call to Innerloop1_2cut"
c                 do i=1,3
c                   do j=1,3
c                 virep3bt_tot(i,j)=virep3bt_tot(i,j)+virep3bt(i,j)
c                   end do
c                 end do
c
c                 do i=1,npole
c                    do j=1,3
c                 dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
c                    end do
c                 end do
c
c                 ep3bt_tot=ep3bt_tot+ep3bt
c                 ntriple_tot=ntriple_tot+ntript
c       if(moli1.ne.0) then
c                  call Innerloop1_single_NoCut2bonly(
c     &               moli1,ep3bt,virep3bt,dep3bt,ntript)
c                 do i=1,3
c                   do j=1,3
c                   virep3bt_tot(i,j)=virep3bt_tot(i,j)+virep3bt(i,j)
c                   end do
c                 end do
c                 do i=1,npole
c                    do j=1,3
c                    dep3bt_tot(j,i)=dep3bt_tot(j,i)+dep3bt(j,i)
c                    end do
c                 end do
c                ep3bt_tot=ep3bt_tot+ep3bt
c                ntriple_tot=ntriple_tot+ntript
c       end if
c       print*,"Before mpi_reduce ep3bt_tot",ep3bt_tot
c                  call mpi_reduce(ep3bt_tot,ep3b,1,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(dep3bt_tot,dep3b,npole*3,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
c                  call mpi_reduce(virep3bt_tot,virep3b,3*3,mpi_real8,
c     &          mpi_sum,master,mpi_comm_world,ierr)
c                   call mpi_reduce(ntriple_tot,ntriples,1,mpi_integer,
c     &          mpi_sum,master,mpi_comm_world,ierr) 
c
c      deallocate(dep3bt_tot)
c      deallocate(dep3bt)
c      return
c      end

      subroutine gradient_covalent(energy,derivs)
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
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /


c      allocate (dep3b(3,npole))
c      ep3b2=0.0d0
c      ep3b3=0.0d0
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      em = 0.0d0

      ep3b=0.0d0

      if (first) then
         first = .false.
c         if (.not. allocated(dep))  allocate (dep(3,n))
         if (.not. allocated(desum))  allocate (desum(3,n))
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
c            virep3b2(i,j)=0.0d0
c            virep3b3(i,j)=0.0d0
            virep3b(i,j)=0.0d0

            vir(i,j)=0.0d0
c            viremreal(i,j)=0.0d0
c            virev(i,j)=0.0d0
          end do
        end do

      do i = 1, n
         do j=1,3
c            dev_total(j,i)=0.0d0
c            dep(j,i) = 0.0d0
c            dep3b2(j,i) = 0.0d0
c            dep3b3(j,i) = 0.0d0
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep3b(j,i) = 0.0d0

c            demreal(j,i) =0.0d0
         end do
      end do

c       call nblist_perm_2bod

c
c     call the van der Waals energy component routines
c
      call ebond1
      call eangle1
      call estrbnd1
      call eurey1
      call eopbend1
      call etors1
      call epitors1
      call etortor1
c      call ehal1
c         if (vdwtyp .eq. 'BUFFERED-14-7') call ehal1c_new
c      print*,"After vdw"

c        call empole1c_3b_Perm_PME
c        call empole1d_3b_Perm_PME
c
c     call any miscellaneous energy component routines
c

c      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
c     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
c     &          + ep + er + es + elf + eg + ex
c      energy = esum
c      do i = 1, n
c         do j = 1, 3
c            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
c     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
c     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
c     &                      + det(j,i) + dept(j,i) + debt(j,i)
c     &                      + dett(j,i) + dev(j,i) + dec(j,i)
c     &                      + decd(j,i) + ded(j,i) + dem(j,i)
c     &                      + dep(j,i) + der(j,i) + des(j,i)
c     &                      + delf(j,i) + deg(j,i) + dex(j,i)
c         end do
c      end do

      esum =eb + ea + eba + eub + eopb + et + ept  + ett + ev + em 
      energy = esum
c      print*,"em ev",em,ev
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i) + deopb(j,i)
     &                   + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i) 
            derivs(j,i) = desum(j,i)
         end do
      end do
      return
      end

      subroutine gradient_pairwiseadd(energy,derivs)
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
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /


c      allocate (dep3b(3,npole))
c      ep3b2=0.0d0
c      ep3b3=0.0d0
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      em = 0.0d0

      ep3b=0.0d0

      if (first) then
         first = .false.
c         if (.not. allocated(dep))  allocate (dep(3,n))
         if (.not. allocated(desum))  allocate (desum(3,n))
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
c            virep3b2(i,j)=0.0d0
c            virep3b3(i,j)=0.0d0
            virep3b(i,j)=0.0d0

            vir(i,j)=0.0d0
c            viremreal(i,j)=0.0d0
c            virev(i,j)=0.0d0
          end do
        end do

      do i = 1, n
         do j=1,3
c            dev_total(j,i)=0.0d0
c            dep(j,i) = 0.0d0
c            dep3b2(j,i) = 0.0d0
c            dep3b3(j,i) = 0.0d0
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep3b(j,i) = 0.0d0

c            demreal(j,i) =0.0d0
         end do
      end do

c       call nblist_perm_2bod
c        call nblist
c
c     call the van der Waals energy component routines
c
c      call ebond1
c      call eangle1
c      call estrbnd1
c      call eurey1
c      call eopbend1
c      call etors1
c      call epitors1
c      call etortor1
      print*,"Before vdw",ev

      call ehal1
c         if (vdwtyp .eq. 'BUFFERED-14-7') call ehal1c_new
      print*,"After vdw",ev

c        call empole1c_3b_Perm_PME
      print*,"Before empole1d_3b_Perm_PME",em

        call empole1d_3b_Perm_PME
      print*,"After empole1d_3b_Perm_PME",em
c
c     call any miscellaneous energy component routines
c

c      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
c     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
c     &          + ep + er + es + elf + eg + ex
c      energy = esum
c      do i = 1, n
c         do j = 1, 3
c            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
c     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
c     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
c     &                      + det(j,i) + dept(j,i) + debt(j,i)
c     &                      + dett(j,i) + dev(j,i) + dec(j,i)
c     &                      + decd(j,i) + ded(j,i) + dem(j,i)
c     &                      + dep(j,i) + der(j,i) + des(j,i)
c     &                      + delf(j,i) + deg(j,i) + dex(j,i)
c         end do
c      end do

      esum =eb + ea + eba + eub + eopb + et + ept  + ett + ev + em 
      energy = esum
c      print*,"em ev",em,ev
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i) + deopb(j,i)
     &                   + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i) 
            derivs(j,i) = desum(j,i)
         end do
      end do
      return
      end
