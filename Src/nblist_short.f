c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2006 by David Gohara & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw, electrostatic and polarization interactions
c
c
c
c
c
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlist  --  get atomic multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for atomic multipoles
c
c
      subroutine mlist2
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh2
      use neigh
      use sizes2
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(mbuf2_short)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (domlst2) then
         domlst2 = .false.
         if (octahedron) then
            call mbuild2
         else
            call mlight2
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xmold2(i)
         yr = yi - ymold2(i)
         zr = zi - zmold2(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            listsend_mpole =.true.
            update(i) = .true.
            xmold2(i) = xi
            ymold2(i) = yi
            zmold2(i) = zi
            nelst2(i) = 0
            do k = i+1, npole
               kk = ipole(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2_short) then
                  nelst2(i) = nelst2(i) + 1
                  elst2(nelst2(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule (guided)
      do i = 1, npole
         if (update(i)) then
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xmold2(k)
                  yr = yi - ymold2(k)
                  zr = zi - zmold2(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. mbuf2_short) then
!$OMP CRITICAL
                     do j = 1, nelst2(k)
                        if (elst2(j,k) .eq. i)  goto 20
                     end do
                     nelst2(k) = nelst2(k) + 1
                     elst2(nelst2(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. mbufx_short) then
!$OMP CRITICAL
                     do j = 1, nelst2(k)
                        if (elst2(j,k) .eq. i) then
                           elst2(j,k) = elst2(nelst2(k),k)
                           nelst2(k) = nelst2(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do i = 1, npole
         if (nelst2(i) .ge. maxelst2) then
            write (iout,40)
   40       format (/,' MLIST  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine mbuild  --  build mpole list for all sites  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mbuild" performs a complete rebuild of the atomic multipole
c     electrostatic neighbor list for all sites
c
c
      subroutine mbuild2
      use sizes
      use atoms
      use bound
      use iounit
      use mpole
      use neigh2
      use neigh
      use sizes2
      implicit none
      integer i,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xmold2(i) = xi
         ymold2(i) = yi
         zmold2(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nelst2(i) = 0
         do k = i+1, npole
            kk = ipole(k)
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2_short) then
               nelst2(i) = nelst2(i) + 1
               elst2(nelst2(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nelst2(i) .ge. maxelst2) then
            write (iout,10)
   10       format (/,' MBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlight  --  get multipole pair list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlight" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine mlight2
      use sizes
      use atoms
      use bound
      use cell
      use iounit
      use light
      use mpole
      use neigh2
      use neigh
      use sizes2
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(npole))
      allocate (ysort(npole))
      allocate (zsort(npole))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, npole
         nelst2(i) = 0
         ii = ipole(i)
         xmold2(i) = x(ii)
         ymold2(i) = y(ii)
         zmold2(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(mbuf2_short)
      call lightn (off,npole,xsort,ysort,zsort)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i)
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            if (k .le. i)  goto 20
            kk = ipole(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2_short) then
               nelst2(i) = nelst2(i) + 1
               elst2(nelst2(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i)
            stop = npole
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nelst2(i) .ge. maxelst2) then
            write (iout,30)
   30       format (/,' MLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c
c
