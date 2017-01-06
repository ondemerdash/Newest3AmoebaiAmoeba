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
c     ##  subroutine ulist  --  get atomic multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ulist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for atomic multipoles
c
c
      subroutine ulist_1bclustboxclust(npole3b,pnum,counter)
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh2
      use neigh
      use sizes2
      use neigh2clust
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
      integer npole3b,pnum(*),counter,l1,l3
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(npole3b))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(ubuf2)
      call replicabox1bclust (radius,counter)
      if (use_replica) then
         write (iout,10)
   10    format (/,' ULIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c

      if (doulstclust1b(counter)) then
         doulstclust1b(counter) = .false.
         if (octahedron) then
            call ubuild_1bclustboxclust(npole3b,pnum,counter)
         else
            call ulight_1bclustboxclust(npole3b,pnum,counter)
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2,i
!$OMP& l1,l3)
!$OMP DO schedule(guided)
      do l1=1,npole3b
         i=pnum(l1)
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xuold1b(l1,counter)
         yr = yi - yuold1b(l1,counter)
         zr = zi - zuold1b(l1,counter)
         call imagen1b (xr,yr,zr,counter)
         r2 = xr*xr + yr*yr + zr*zr
         update(l1) = .false.
         if (r2 .ge. pbuf2) then
c            listsend_mpole =.true.
            update(l1) = .true.
            xuold1b(l1,counter) = xi
            yuold1b(l1,counter) = yi
            zuold1b(l1,counter) = zi
            nulst1b(l1,counter) = 0
           ! do k = i+1, npole
            do l3=l1+1,npole3b
               k=pnum(l3)
               kk = ipole(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen1b (xr,yr,zr,counter)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. ubuf2) then
                  nulst1b(l1,counter) = nulst1b(l1,counter) + 1
                  ulst1b(nulst1b(l1,counter),l1,counter) = l3
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule (guided)
c      do i = 1, npole
      do l1 = 1, npole3b
         i=pnum(l1)
         if (update(l1)) then
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
c            do k = 1, i-1
            do l3 = 1, l1-1
               if (.not. update(l3)) then
                  xr = xi - xuold1b(l3,counter)
                  yr = yi - yuold1b(l3,counter)
                  zr = zi - zuold1b(l3,counter)
                  call imagen1b (xr,yr,zr,counter)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. ubuf2) then
!$OMP CRITICAL
                     do j = 1, nulst1b(l3,counter)
                        if (ulst1b(j,l3,counter) .eq. l1)  goto 20
                     end do
                     nulst1b(l3,counter) = nulst1b(l3,counter) + 1
                     ulst1b(nulst1b(l3,counter),l3,counter) = l1
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. ubufx) then
!$OMP CRITICAL
                     do j = 1, nulst1b(l3,counter)
                        if (ulst1b(j,l3,counter) .eq. l1) then
                           ulst1b(j,l3,counter) = 
     &                       ulst1b(nulst1b(l3,counter),l3,counter)
                           nulst1b(l3,counter) = nulst1b(l3,counter) - 1
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
      do i = 1, npole3b
         if (nulst1b(i,counter) .ge. maxulst2) then
            write (iout,40)
   40       format (/,' ULIST  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
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
c     ##  subroutine ubuild  --  build mpole list for all sites  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ubuild" performs a complete rebuild of the atomic multipole
c     electrostatic neighbor list for all sites
c
c
      subroutine ubuild_1bclustboxclust(npole3b,pnum,counter)
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
      integer npole3b,pnum(*),counter
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
     &                 ' Increase MAXULST')
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
c     ##  subroutine ulight  --  get multipole pair list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ulight" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine ulight_1bclustboxclust(npole3b,pnum,counter)
      use sizes
      use atoms
      use bound
      use cell
      use iounit
      use light
      use mpole
      use neigh2clust
      use neigh
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
      integer npole3b,pnum(*),counter,l1
      logical first
      save first
      data first  / .true. /

c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(npole3b))
      allocate (ysort(npole3b))
      allocate (zsort(npole3b))

c
c     transfer interaction site coordinates to sorting arrays
c
      do l1 = 1, npole3b
         i=pnum(l1)
         nulst1b(l1,counter) = 0
         ii = ipole(i)
         xuold1b(l1,counter) = x(ii)
         yuold1b(l1,counter) = y(ii)
         zuold1b(l1,counter) = z(ii)
         xsort(l1) = x(ii)
         ysort(l1) = y(ii)
         zsort(l1) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(ubuf2)
      call lightnbox1bclust (off,npole3b,xsort,ysort,zsort,counter)
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
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat,l1)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
c      do i = 1, npole
      do l1=1,npole3b
         i=pnum(l1)
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(l1) .le. kex(l1)) then
            repeat = .false.
            start = kbx(l1)
            stop = kex(l1)
         else
            repeat = .true.
            start = 1
            stop = kex(l1)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kk = ipole(pnum(k))
            if (k .le. l1)  goto 20
            kgy = rgy(k)
            if (kby(l1) .le. key(l1)) then
               if (kgy.lt.kby(l1) .or. kgy.gt.key(l1))  goto 20
            else
               if (kgy.lt.kby(l1) .and. kgy.gt.key(l1))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(l1) .le. kez(l1)) then
               if (kgz.lt.kbz(l1) .or. kgz.gt.kez(l1))  goto 20
            else
               if (kgz.lt.kbz(l1) .and. kgz.gt.kez(l1))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen1b (xr,yr,zr,counter)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ubuf2) then
               nulst1b(l1,counter) = nulst1b(l1,counter) + 1
               ulst1b(nulst1b(l1,counter),l1,counter) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(l1)
            stop = npole3b
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nulst1b(l1,counter) .ge. maxulst2) then
            write (iout,30) nulst1b(l1,counter)
   30       format (/,' ULIGHT1B  --  Too many Neighbors;',
     &                 ' Increase MAXULST',i6)
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
