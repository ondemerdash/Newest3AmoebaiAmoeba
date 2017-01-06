c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine vlist  --  get van der Waals neighbor lists  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vlist" performs an update or a complete rebuild of the
c     van der Waals neighbor list
c
c
      subroutine vlist_par
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use neigh
      use vdw
      use mpidat
      implicit none
      integer i,j,k
      integer ii,iv
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius
      real*8 rdn,r2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      logical, allocatable :: update(:)
      integer taskid_offset
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (update(n))
c
c     apply reduction factors to find coordinates for each site
c
      do i = 1, nvdw
         ii = ivdw(i)
         iv = ired(ii)
         rdn = kred(ii)
         xred(i) = rdn*(x(ii)-x(iv)) + x(iv)
         yred(i) = rdn*(y(ii)-y(iv)) + y(iv)
         zred(i) = rdn*(z(ii)-z(iv)) + z(iv)
      end do
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(vbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' VLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (dovlst) then
         dovlst = .false.
         if (octahedron) then
            call vbuild (xred,yred,zred)
         else
            call vlight (xred,yred,zred)
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
      if(numtasks_vdw.lt.numtasks) then
        taskid_offset=taskid-numtasks_emreal
      else
        taskid_offset=taskid
      end if
c      print*,"vlist taskid_offset=",taskid_offset
c      print*,"vlist taskid=",taskid_offset,"start_vdw2(taskid_offset)",
c     &              start_vdw2(taskid_offset)
c      print*,"vlist taskid=",taskid_offset,"last_vdw2(taskid_offset)",
c     &              last_vdw2(taskid_offset)


!$OMP PARALLEL default(shared) private(i,j,k,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c      do i = 1, nvdw
      do i =start_vdw2(taskid_offset),last_vdw2(taskid_offset)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xr = xi - xvold(i)
         yr = yi - yvold(i)
         zr = zi - zvold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            listsend_vdw =.true.
            update(i) = .true.
            xvold(i) = xi
            yvold(i) = yi
            zvold(i) = zi
            nvlst(i) = 0
            do k = i+1, nvdw
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  nvlst(i) = nvlst(i) + 1
                  vlst(nvlst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
c      do i = 1, nvdw
      do i =start_vdw2(taskid_offset),last_vdw2(taskid_offset)
         if (update(i)) then
            xi = xred(i)
            yi = yred(i)
            zi = zred(i)
c            do k = 1, i-1
            do k = start_vdw2(taskid_offset),i-1           
               if (.not. update(k)) then
                  xr = xi - xvold(k)
                  yr = yi - yvold(k)
                  zr = zi - zvold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. vbuf2) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i)  goto 20
                     end do
                     nvlst(k) = nvlst(k) + 1
                     vlst(nvlst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. vbufx) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i) then
                           vlst(j,k) = vlst(nvlst(k),k)
                           nvlst(k) = nvlst(k) - 1
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
c      do i = 1, nvdw
      do i =start_vdw2(taskid_offset),last_vdw2(taskid_offset)
         if (nvlst(i) .ge. maxvlst) then
            write (iout,40)
   40       format (/,' VLIST  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
c      print*,"Successful completion of vlist update!"
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (update)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
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
      subroutine mlist_par
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh
      use mpidat
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
      radius = sqrt(mbuf2)
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
      if (domlst) then
         domlst = .false.
         if (octahedron) then
            call mbuild
         else
            call mlight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c      do i = 1, npole
      do i=start_emreal2(taskid),last_emreal2(taskid)      
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xmold(i)
         yr = yi - ymold(i)
         zr = zi - zmold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            listsend_mpole =.true.
            update(i) = .true.
            xmold(i) = xi
            ymold(i) = yi
            zmold(i) = zi
            nelst(i) = 0
            do k = i+1, npole
               kk = ipole(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
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
      do i=start_emreal2(taskid),last_emreal2(taskid)
         if (update(i)) then
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
c            do k = 1, i-1
            do k = start_emreal2(taskid),i-1
               if (.not. update(k)) then
                  xr = xi - xmold(k)
                  yr = yi - ymold(k)
                  zr = zi - zmold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. mbuf2) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i)  goto 20
                     end do
                     nelst(k) = nelst(k) + 1
                     elst(nelst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. mbufx) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i) then
                           elst(j,k) = elst(nelst(k),k)
                           nelst(k) = nelst(k) - 1
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
c      do i = 1, npole
      do i=start_emreal2(taskid),last_emreal2(taskid)
         if (nelst(i) .ge. maxelst) then
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
c
c
