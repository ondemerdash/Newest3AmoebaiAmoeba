c
      subroutine mollist2bodyOO6_8_par_1list 
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
      implicit none
      integer i,j,k,ii,l1,lenmol,j1,j2,k1
      real*8 xi,yi,zi,xcm1,ycm1,zcm1
      real*8 xr,yr,zr,M1,lbuffer2b,lbuf2_2b
      real*8 radius,r2,r,xk,yk,zk
      logical, allocatable :: update(:)
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
            call Newmolfull2bodyOO6_8_par_1list 
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
      do i = 1, nmol
c      do i = start_polar(taskid),last_polar(taskid)
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k1=l1-1
            j=imol(1,i)+k1
           if(name(j).eq.'O') then
            xi = x(j)
            yi = y(j)
            zi = z(j)
           end if
         end do

         xr = xi - xmolold(i)
         yr = yi - ymolold(i)
         zr = zi - zmolold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2_2b) then
            listbcast =.true.
c            print*,"listbcast=",listbcast
            update(i) = .true.
            xmolold(i) = xi
            ymolold(i) = yi
            zmolold(i) = zi
            nmollst(i) = 0
            do k= i+1,nmol
               lenmol=imol(2,k)-imol(1,k)+1
               do l1 = 1, lenmol
                  k1=l1-1
                  j=imol(1,k)+k1
                  if(name(j).eq.'O') then
                    xk = x(j)
                    yk = y(j)
                    zk = z(j)
                  end if
               end do
               xr = xi - xk
               yr = yi - yk
               zr = zi - zk
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               r=sqrt(r2)
               if (r .le. molbuf2s) then
                  nmollst(i) = nmollst(i) + 1
                  mollst(nmollst(i),i)=k
               end if
            end do
         end if
      end do 
!$OMP END DO

!$OMP DO schedule(guided)
      do i = 1, nmol
c      do i = start_polar(taskid),last_polar(taskid)
         if (update(i)) then
            lenmol=imol(2,i)-imol(1,i)+1

            do l1 = 1, lenmol
             k1=l1-1
             j=imol(1,i)+k1
             if(name(j).eq.'O') then
              xi = x(j)
              yi = y(j)
              zi = z(j)
             end if
            end do
c
c     adjust lists of lower numbered neighbors of updated sites
c

            do k = 1, i-1
c            do k =start_polar(taskid),i-1
               if (.not. update(k)) then
                  xr = xi - xmolold(k)
                  yr = yi - ymolold(k)
                  zr = zi - zmolold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  r=sqrt(r2) 
                  if (r .le. molbuf2s) then
!$OMP CRITICAL
                     do j = 1, nmollst(k)
                        if (mollst(j,k) .eq. i)  goto 20
                     end do
                     nmollst(k) = nmollst(k) + 1
                     mollst(nmollst(k),k) = i
   20                continue
!$OMP END CRITICAL
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
          do i=1,nmol
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


c
      subroutine Oldmollist2bodyOO6_8_par_1list
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use neigh3b
      use molcul
      use atomid 
      use neigh
      implicit none
      integer i,j,k,ii,l1,lenmol,j1,j2,k1
      real*8 xi,yi,zi,xcm1,ycm1,zcm1
      real*8 xr,yr,zr,M1,lbuffer2b,lbuf2_2b
      real*8 radius,r2,r,xk,yk,zk
      logical, allocatable :: update(:)
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
            call Newmolfull2bodyOO6_8_par_1list 
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
      do i = 1, nmol
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k1=l1-1
            j=imol(1,i)+k1
           if(name(j).eq.'O') then
            xi = x(j)
            yi = y(j)
            zi = z(j)
           end if
         end do

         xr = xi - xmolold(i)
         yr = yi - ymolold(i)
         zr = zi - zmolold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2_2b) then
            listbcast =.true.
c            print*,"listbcast=",listbcast
            update(i) = .true.
            xmolold(i) = xi
            ymolold(i) = yi
            zmolold(i) = zi
         end if
      end do 
!$OMP END DO

!$OMP DO schedule(guided)
       do i = 1, nmol
         if (update(i)) then
            lenmol=imol(2,i)-imol(1,i)+1

            do l1 = 1, lenmol
             k1=l1-1
             j=imol(1,i)+k1
             if(name(j).eq.'O') then
              xi = x(j)
              yi = y(j)
              zi = z(j)
             end if
            end do
            j1 = 0
            j2 = 0
            do k = i+1, nmol
               lenmol=imol(2,k)-imol(1,k)+1

               do l1 = 1, lenmol
                  k1=l1-1
                  j=imol(1,k)+k1
                  if(name(j).eq.'O') then
                    xk = x(j)
                    yk = y(j)
                    zk = z(j)
                  end if
               end do
               xr = xi - xk
               yr = yi - yk
               zr = zi - zk
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               r=sqrt(r2)
               if (r .le. molbuf2s) then
                  j1 = j1 + 1
                  mollst(j1,i)=k
               end if
            end do 
            nmollst(i)=j1
c
c     adjust lists of lower numbered neighbors of updated sites
c
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xmolold(k)
                  yr = yi - ymolold(k)
                  zr = zi - zmolold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  r=sqrt(r2) 
                  if (r .le. molbuf2s) then
                     do j = 1, nmollst(k)
                        if (mollst(j,k) .eq. i)  goto 20
                     end do
!$OMP CRITICAL
                     nmollst(k) = nmollst(k) + 1
                     mollst(nmollst(k),k) = i
!$OMP END CRITICAL
   20                continue
                  else if (r .le. molbuf2s_x) then
                     do j = 1, nmollst(k)
                        if (mollst(j,k) .eq. i) then
!$OMP CRITICAL
                           mollst(j,k) = mollst(nmollst(k),k)
                           nmollst(k) = nmollst(k) - 1
!$OMP END CRITICAL
                           goto 30
                        end if
                     end do
   30                continue
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
          do i=1,nmol
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


      subroutine Newmolfull2bodyOO6_8_par_1list
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
      implicit none
      integer i,j,k,ii
      integer kgy,kgz,l1,lenmol
      integer start,stop
      real*8 xi,yi,zi,M1
      real*8 xr,yr,zr,r,x1,y1,z1
      real*8 r2,off,xcm1,ycm1,zcm1,lbuffer2b
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      integer, allocatable :: mollst_tmp(:,:)
      integer, allocatable :: nmollst_tmp(:)
      integer, allocatable :: mollst2_tmp(:,:)
      integer, allocatable :: nmollst2_tmp(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
c      print*, "Before allocation Molfull2body"
c      allocate (mollst_tmp(max2blst,nmol))
c      allocate (nmollst_tmp(nmol))
c
c      allocate (mollst2_tmp(max2blst,nmol))
c      allocate (nmollst2_tmp(nmol))

c      nlight = (ncell+1) * nmol
      allocate (xsort(nmol))
      allocate (ysort(nmol))
      allocate (zsort(nmol))
c      domollst2bod = .false.
c      print*, "After allocation Molfull2body"
c
c     transfer interaction site coordinates to sorting arrays
c

      do i = 1, nmol
         nmollst(i) = 0

         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k=l1-1
            j=imol(1,i)+k
           if(name(j).eq.'O') then
            x1 = x(j)
            y1 = y(j)
            z1 = z(j)
           end if
         end do


         xmolold(i) = x1
         ymolold(i) = y1
         zmolold(i) = z1
         xsort(i) = x1
         ysort(i) = y1
         zsort(i) = z1
      end do

c
      lbuffer2b=lbuffer
c      lbuffer2b=lbuffer
      molbuf2l=cut2b_input+lbuffer2b
      molbuf2s=cut2b_input+lbuffer2b

      print*,"in Newmolfull2bodyOO6_8_par cut2b_input",cut2b_input

c      off = sqrt(molbuf2)
c      off=molbuf2
      off=molbuf2s
      call lightn (off,nmol,xsort,ysort,zsort)

      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)

c      print*,"After lights Call in Molfull2body"
c
c     loop over all atoms computing the interactions
c
!$OMP PARALLEL default(private) shared(nmol,mollst,
!$OMP& nmollst,mollst2,xmolold,ymolold,zmolold,kbx,kex,
!$OMP& nmollst2,kby,key,kbz,kez,rgy,rgz,locx,molbuf2s,
!$OMP& molbuf2l)
!$OMP DO schedule(guided)
      do i = 1, nmol
         xi = xmolold(i)
         yi = ymolold(i)
         zi = zmolold(i)

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
            if (k .le. i) goto 20
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

            xr = xi - xmolold(k)
            yr = yi - ymolold(k)
            zr = zi - zmolold(k)

            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            if (r .le. molbuf2s) then
c               if (i .lt. k) then
c                  nmollst_tmp(i) = nmollst_tmp(i) + 1
c                  mollst_tmp(nmollst_tmp(i),i) = k
c               else
c                  nmollst_tmp(k) = nmollst_tmp(k) + 1
c                  mollst_tmp(nmollst_tmp(k),k) = i
c               end if
                  nmollst(i) = nmollst(i) + 1
                  mollst(nmollst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nmol
            goto 10
         end if

         if (nmollst(i).gt.max2blst) then
          print*,"nmollst Too many 2b short nb.Tilt!",nmollst(i)
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL 

c      print*,"End of molfull2bodyOO6_8_par_new"
      return
      end
