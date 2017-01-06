c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull2body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mfull" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull2body_NewBodies
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'bound.i'
c      include 'cell.i'
c      include 'iounit.i'
c      include 'light.i'
c      include 'mpole.i'
c      include 'neigh.i'
c      include 'neigh2.i'
c      include 'molcul.i'
c      include 'combo.i'
c      include 'atmtyp.i'
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use neigh3b
      use molcul
      use atomid
      use neigh
      use light
      use cell
      use neigh2clust
      implicit none
      integer i,j,k,ii
      integer kgy,kgz,l1,lenmol
      integer start,stop
      real*8 xi,yi,zi,M1
      real*8 xr,yr,zr,r
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
      real*8 lbuffer2b,molbuf2_x
c
c
c     perform dynamic allocation of some local arrays
c
c      print*, "Before allocation Molfull2body"

      nlight = (ncell+1) * clustcount
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c      domollst2bod = .false.
c      print*, "After allocation Molfull2body"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, clustcount
         nmollst(i) = 0         
         xmolold(i) = clust_cm(1,i)
         ymolold(i) = clust_cm(2,i)
         zmolold(i) = clust_cm(3,i)
         xsort(i) = clust_cm(1,i)
         ysort(i) = clust_cm(2,i)
         zsort(i) = clust_cm(3,i)
      end do
c
c     use the method of lights to generate neighbors
c
c      print*,"Before lights Call in Molfull2body"
      lbuffer2b=lbuffer
      molbuf2=8.0d0+lbuffer2b

      lbuffer2b=lbuffer
c      lbuffer2b=lbuffer
      molbuf2l=cut2b_input+lbuffer2b
      molbuf2s=cut2b_input+lbuffer2b

c      off = sqrt(molbuf2)
      off=molbuf2
      off=molbuf2s

      call lights (off,clustcount,xsort,ysort,zsort)

c      print*,"After lights Call in Molfull2body"
c
c     loop over all atoms computing the interactions
c
      do i = 1, clustcount
         xi = xsort(rgx(i))
         yi = ysort(rgy(i))
         zi = zsort(rgz(i))
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
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
            xr = xi - xsort(j)
            yr = yi - ysort(kgy)
            zr = zi - zsort(kgz)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
c            if (r .le. molbuf2) then
            if (r .le. molbuf2s) then
               if (i .lt. k) then
                  nmollst(i) = nmollst(i) + 1
                  mollst(nmollst(i),i) = k
               else
                  nmollst(k) = nmollst(k) + 1
                  mollst(nmollst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      do i=1,clustcount
         print*,"Clust ind, nmollst(i)",i,nmollst(i)
         do j=1,nmollst(i)
            print*,"Clust ind,mollst(j,i)",mollst(j,i)
         end do 
      end do 
c
c     check to see if the neighbor lists are too long
c
      do i = 1,clustcount 
         if (nmollst(i) .ge. max2blst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')
          print*,"nmollst Too many 2b short nb.Tilt!",nmollst(i)
            call fatal
         end if
      end do
      return
      end

      subroutine mollist2body_NewBodies 
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use neigh3b
      use molcul
      use atomid
      use neigh
      use mpidat
      use neigh2clust
      implicit none
      integer i,j,k,ii,l1,lenmol,j1,j2,k1
      real*8 xi,yi,zi,xcm1,ycm1,zcm1
      real*8 xr,yr,zr,M1,lbuffer2b,lbuf2_2b
      real*8 radius,r2,r,xk,yk,zk
      logical, allocatable :: update(:)
      real*8 molbuf2_x
c      lbuffer2b=lbuffer
c      molbuf2l=cut2b_input+lbuffer2b
c      molbuf2s=cut2b_input+lbuffer2b
c      molbuf2l_x=cut2b_input+2.0d0*lbuffer2b
c      molbuf2s_x=cut2b_input+2.0d0*lbuffer2b
c      lbuf2_2b=(0.5d0*lbuffer2b)**2
      lbuffer2b=lbuffer
      molbuf2=8.0d0+lbuffer2b
      lbuf2_2b=(0.5d0*lbuffer2b)**2
      molbuf2_x=8.0d0+2.0d0*lbuffer2b
      lbuffer2b=lbuffer
      molbuf2l=cut2b_input+lbuffer2b
      molbuf2s=cut2b_input+lbuffer2b
      molbuf2l_x=cut2b_input+2.0d0*lbuffer2b
      molbuf2s_x=cut2b_input+2.0d0*lbuffer2b
      lbuf2_2b=(0.5d0*lbuffer2b)**2

      allocate (update(n))
c      print*,"in mollist2bodyOO6_8 cut2b_input",cut2b_input
c
c
c     neighbor list cannot be used with the replicates method
c

c      radius = sqrt(molbuf2)
      radius = sqrt(molbuf2s)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' MOLLST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
c      print*, "Hello from Mollist2body"

      if (domollst2bod) then
         domollst2bod = .false.
c         if (octahedron) then
c            do i = 1, npole
c               call mbuild (i)
c            end do
c         else
c          print*, "Within Mollist2body if domollst2bod"
            call molfull2body_NewBodies 
c            call molfull2bodyOO6_8
c         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
!$OMP PARALLEL default(shared) private(i,lenmol,l1,j,k,xi,yi,zi,
!$OMP& xr,yr,zr,r2,r,j1,j2,k1,xk,yk,zk)
!$OMP DO schedule(guided)
      do i = 1, clustcount 

         xr = clust_cm(1,i) - xmolold(i)
         yr = clust_cm(2,i) - ymolold(i)
         zr = clust_cm(3,i) - zmolold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2_2b) then
            listbcast =.true.
c            print*,"listbcast=",listbcast
            update(i) = .true.
            xmolold(i) = clust_cm(1,i)
            ymolold(i) = clust_cm(2,i)
            zmolold(i) = clust_cm(3,i)
            nmollst(i) = 0
            do k= i+1,clustcount
c               lenmol=imol(2,k)-imol(1,k)+1
c               do l1 = 1, lenmol
c                  k1=l1-1
c                  j=imol(1,k)+k1
c                  if(name(j).eq.'O') then
c                    xk = x(j)
c                    yk = y(j)
c                    zk = z(j)
c                  end if
c               end do

               xr = clust_cm(1,i)-clust_cm(1,k)
               yr = clust_cm(2,i)-clust_cm(2,k)
               zr = clust_cm(3,i)-clust_cm(3,k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               r=sqrt(r2)
c               if (r .le. molbuf2) then
               if (r .le. molbuf2s) then
                  nmollst(i) = nmollst(i) + 1
                  mollst(nmollst(i),i)=k
               end if
            end do
         end if
      end do
!$OMP END DO

!$OMP DO schedule(guided)
      do i = 1, clustcount
         if (update(i)) then
c            lenmol=imol(2,i)-imol(1,i)+1
c
c            do l1 = 1, lenmol
c             k1=l1-1
c             j=imol(1,i)+k1
c             if(name(j).eq.'O') then
c              xi = x(j)
c              yi = y(j)
c              zi = z(j)
c             end if
c            end do

c
c     adjust lists of lower numbered neighbors of updated sites
c

            do k = 1, i-1
               if (.not. update(k)) then
                  xr = clust_cm(1,i) - xmolold(k)
                  yr = clust_cm(2,i) - ymolold(k)
                  zr = clust_cm(3,i) - zmolold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  r=sqrt(r2)
c                  if (r .le. molbuf2) then
                  if (r .le. molbuf2s) then

!$OMP CRITICAL
                     do j = 1, nmollst(k)
                        if (mollst(j,k) .eq. i)  goto 20
                     end do
                     nmollst(k) = nmollst(k) + 1
                     mollst(nmollst(k),k) = i
   20                continue
!$OMP END CRITICAL
c                  else if (r .le. molbuf2_x) then
                  else if (r .le. molbuf2s_x) then
!$OMP CRITICAL
                     do j = 1, nmollst(k)
                        if (mollst(j,k) .eq. i) then
                           mollst(j,k) = mollst(nmollst(k),k)
                           nmollst(k) = nmollst(k) - 1
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
!$OMP DO
          do i=1,clustcount
c          do i =start_polar(taskid),last_polar(taskid)
             if (nmollst(i).gt.max2blst) then
              print*,"nmollst Too many 2b short nbs.Tilt!",nmollst(i)
              call fatal
             end if
          end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c

      deallocate (update)
c      print*,"End of mollist2bodyOO6_8"
      return
      end

