c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mollist -- build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mollist" performs an update or a complete rebuild of the
c     body/molecule neighbor lists for 3-body approx of polariz. energy
c
c
c      subroutine mollist2body
      subroutine mollist2bodyOO6_8_par
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
      print*,"in mollist2bodyOO6_8 cut2b_input",cut2b_input
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
            call Newmolfull2bodyOO6_8_par 
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
               if (r .le. molbuf2l) then
                  j2 = j2 + 1
                  mollst2(j2,i)=k
               end if
            end do 
            nmollst(i)=j1
            nmollst2(i)=j2
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


                  if (r .le. molbuf2l) then
                     do j = 1, nmollst2(k)
                        if (mollst2(j,k) .eq. i)  goto 40
                     end do
!$OMP CRITICAL
                     nmollst2(k) = nmollst2(k) + 1
                     mollst2(nmollst2(k),k) = i
!$OMP END CRITICAL
   40                continue
                  else if (r .le. molbuf2l_x) then
                     do j = 1, nmollst2(k)
                        if (mollst2(j,k) .eq. i) then
!$OMP CRITICAL
                           mollst2(j,k) = mollst2(nmollst2(k),k)
                           nmollst2(k) = nmollst2(k) - 1
!$OMP END CRITICAL
                           goto 50
                        end if
                     end do
   50                continue
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
             if (nmollst2(i).gt.max2blst2) then
              print*,"nmollst2 Too many 2b short nbs.Tilt!",nmollst2(i)
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

      subroutine molfull2bodyOO6_8_par
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
      allocate (mollst_tmp(max2blst,nmol))
      allocate (nmollst_tmp(nmol))

      allocate (mollst2_tmp(max2blst,nmol))
      allocate (nmollst2_tmp(nmol))

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
      domollst2bod = .false.
c      print*, "After allocation Molfull2body"
c
c     transfer interaction site coordinates to sorting arrays
c

      do i = 1, nmol
         nmollst(i) = 0
         nmollst2(i) = 0
         nmollst_tmp(i)=0
         nmollst2_tmp(i)=0

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
      molbuf2l=8.0d0+lbuffer2b
      molbuf2s=6.0d0+lbuffer2b

c      off = sqrt(molbuf2)
c      off=molbuf2
      off=molbuf2l
      call lights (off,nmol,xsort,ysort,zsort)

      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)

c      print*,"After lights Call in Molfull2body"
c
c     loop over all atoms computing the interactions
c
!$OMP PARALLEL default(private) shared(nmol,mollst_tmp,
!$OMP& nmollst_tmp,mollst2_tmp,xmolold,ymolold,zmolold,kbx,kex,
!$OMP& nmollst2_tmp,kby,key,kbz,kez,rgy,rgz,locx,molbuf2s,
!$OMP& molbuf2l,nlight)
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
               if (i .lt. k) then
                  nmollst_tmp(i) = nmollst_tmp(i) + 1
                  mollst_tmp(nmollst_tmp(i),i) = k
               else
                  nmollst_tmp(k) = nmollst_tmp(k) + 1
                  mollst_tmp(nmollst_tmp(k),k) = i
               end if
            end if
            if (r .le. molbuf2l) then
               if (i .lt. k) then
                  nmollst2_tmp(i) = nmollst2_tmp(i) + 1
                  mollst2_tmp(nmollst2_tmp(i),i) = k
               else
                  nmollst2_tmp(k) = nmollst2_tmp(k) + 1
                  mollst2_tmp(nmollst2_tmp(k),k) = i
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
!$OMP END DO
!$OMP END PARALLEL 
      do i = 1,nmol
         if (nmollst_tmp(i).gt.max2blst) then
          print*,"nmollst Too many 2b short nb.Tilt!",nmollst_tmp(i)
            call fatal
         end if
         if (nmollst2_tmp(i).gt.max2blst) then
          print*,"nmollst2 Too many 2b short nb.Tilt!",nmollst2_tmp(i)
            call fatal
         end if

         nmollst(i)=nmollst_tmp(i)
         nmollst2(i)=nmollst2_tmp(i)
         print*,"i nmollst",i,nmollst(i)
         print*,"i nmollst2",i,nmollst2(i)

         do j=1,nmollst2_tmp(i)
            mollst2(j,i)=mollst2_tmp(j,i)
         end do
         do j=1,nmollst_tmp(i)
            mollst(j,i)=mollst_tmp(j,i)
         end do
      end do

      deallocate (nmollst_tmp)
      deallocate (mollst_tmp)
      deallocate (nmollst2_tmp)
      deallocate (mollst2_tmp)
      print*,"End of molfull2bodyOO6_8_par"
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mollist -- build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mollist" performs an update or a complete rebuild of the
c     body/molecule neighbor lists for 3-body approx of polariz. energy
c
c
c      subroutine mollist2body
      subroutine mollist2bodyOO6_8_par_633style
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
      molbuf2l=8.0d0+lbuffer2b
      molbuf2s=6.0d0+lbuffer2b
      molbuf2l_x=8.0d0+2.0d0*lbuffer2b
      molbuf2s_x=6.0d0+2.0d0*lbuffer2b      
      lbuf2_2b=(0.5d0*lbuffer2b)**2
c      allocate (update(n))

c
c
c     neighbor list cannot be used with the replicates method
c
c      radius = sqrt(mbuf2)
c      call replica (radius)
c      if (use_replica) then
c         write (iout,10)
c   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
c     &              ' be used with Replicas')
c         call fatal
c      end if
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
            call molfull2bodyOO6_8_par 
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
         if (r2 .ge. lbuf2_2b) then
            listbcast =.true.
            print*,"listbcast=",listbcast
c            update(i) = .true.
            xmolold(i) = xi
            ymolold(i) = yi
            zmolold(i) = zi

            j1 = 0
            j2 = 0
            do k = i+1, nmol
               xr = xi - xmolold(k)
               yr = yi - ymolold(k)
               zr = zi - zmolold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               r=sqrt(r2)
               if (r .le. molbuf2s) then
                  j1 = j1 + 1
                  mollst(j1,i)=k
               end if
               if (r .le. molbuf2l) then
                  j2 = j2 + 1
                  mollst2(j2,i)=k
               end if
            end do 
            nmollst(i)=j1
            nmollst2(i)=j2
            if (nmollst(i) .ge. max2blst) then
             print*,"Too many nmollst",nmollst(i)
               call fatal
            end if
            if (nmollst2(i) .ge. max2blst) then
             print*,"Too many nmollst2",nmollst2(i)
               call fatal
            end if

c
c     adjust lists of lower numbered neighbors of updated sites
c
            do k = 1, i-1
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

                  if (r .le. molbuf2l) then
                     do j = 1, nmollst2(k)
                        if (mollst2(j,k) .eq. i)  goto 40
                     end do
!$OMP CRITICAL
                     nmollst2(k) = nmollst2(k) + 1
                     mollst2(nmollst2(k),k) = i
!$OMP END CRITICAL
   40                continue
                  else if (r .le. molbuf2l_x) then
                     do j = 1, nmollst2(k)
                        if (mollst2(j,k) .eq. i) then
!$OMP CRITICAL
                           mollst2(j,k) = mollst2(nmollst2(k),k)
                           nmollst2(k) = nmollst2(k) - 1
!$OMP END CRITICAL
                           goto 50
                        end if
                     end do
   50                continue
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
             if (nmollst2(i).gt.max2blst) then
              print*,"nmollst2 Too many 2b short nbs.Tilt!",nmollst2(i)
              call fatal
             end if
          end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      
      return
      end

      subroutine molfull2bodyOO6_8
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
      allocate (mollst_tmp(max2blst,nmol))
      allocate (nmollst_tmp(nmol))

      allocate (mollst2_tmp(max2blst,nmol))
      allocate (nmollst2_tmp(nmol))

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
      domollst2bod = .false.
c      print*, "After allocation Molfull2body"
c
c     transfer interaction site coordinates to sorting arrays
c

      do i = 1, nmol
         nmollst(i) = 0
         nmollst2(i) = 0
         nmollst_tmp(i)=0
         nmollst2_tmp(i)=0

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
      molbuf2l=8.0d0+lbuffer2b
      molbuf2s=6.0d0+lbuffer2b

c      off = sqrt(molbuf2)
c      off=molbuf2
      off=molbuf2l
      call lights (off,nmol,xsort,ysort,zsort)

      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)

c      print*,"After lights Call in Molfull2body"
c
c     loop over all atoms computing the interactions
c
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
               if (i .lt. k) then
                  nmollst_tmp(i) = nmollst_tmp(i) + 1
                  mollst_tmp(nmollst_tmp(i),i) = k
               else
                  nmollst_tmp(k) = nmollst_tmp(k) + 1
                  mollst_tmp(nmollst_tmp(k),k) = i
               end if
            end if
            if (r .le. molbuf2l) then
               if (i .lt. k) then
                  nmollst2_tmp(i) = nmollst2_tmp(i) + 1
                  mollst2_tmp(nmollst2_tmp(i),i) = k
               else
                  nmollst2_tmp(k) = nmollst2_tmp(k) + 1
                  mollst2_tmp(nmollst2_tmp(k),k) = i
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
      do i = 1,nmol
         if (nmollst_tmp(i).gt.max2blst) then
          print*,"nmollst Too many 2b short nb.Tilt!",nmollst_tmp(i)
            call fatal
         end if
         if (nmollst2_tmp(i).gt.max2blst) then
          print*,"nmollst2 Too many 2b short nb.Tilt!",nmollst2_tmp(i)
            call fatal
         end if

         nmollst(i)=nmollst_tmp(i)
         nmollst2(i)=nmollst2_tmp(i)
         print*,"i nmollst",i,nmollst(i)
         print*,"i nmollst2",i,nmollst2(i)

         do j=1,nmollst2_tmp(i)
            mollst2(j,i)=mollst2_tmp(j,i)
         end do
         do j=1,nmollst_tmp(i)
            mollst(j,i)=mollst_tmp(j,i)
         end do
      end do

      deallocate (nmollst_tmp)
      deallocate (mollst_tmp)
      deallocate (nmollst2_tmp)
      deallocate (mollst2_tmp)
      print*,"End of molfull2bodyOO6_8"
      return
      end

      subroutine Newmolfull2bodyOO6_8_par
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
         nmollst2(i) = 0

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
      off=molbuf2l
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
            if (r .le. molbuf2l) then
c               if (i .lt. k) then
c                  nmollst2_tmp(i) = nmollst2_tmp(i) + 1
c                  mollst2_tmp(nmollst2_tmp(i),i) = k
c               else
c                  nmollst2_tmp(k) = nmollst2_tmp(k) + 1
c                  mollst2_tmp(nmollst2_tmp(k),k) = i
c               end if
                  nmollst2(i) = nmollst2(i) + 1
                  mollst2(nmollst2(i),i) = k
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
         if (nmollst2(i).gt.max2blst2) then
          print*,"nmollst2 Too many 2b short nb.Tilt!",nmollst2(i)
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL 

      !do i = 1,nmol
c         if (nmollst_tmp(i).gt.max2blst) then
c          print*,"nmollst Too many 2b short nb.Tilt!",nmollst_tmp(i)
c            call fatal
c         end if
c         if (nmollst2_tmp(i).gt.max2blst) then
c          print*,"nmollst2 Too many 2b short nb.Tilt!",nmollst2_tmp(i)
c            call fatal
c         end if
c
c         nmollst(i)=nmollst_tmp(i)
c         nmollst2(i)=nmollst2_tmp(i)
        ! print*,"i nmollst",i,nmollst(i)
        ! print*,"i nmollst2",i,nmollst2(i)
c
c         do j=1,nmollst2_tmp(i)
c            mollst2(j,i)=mollst2_tmp(j,i)
c         end do
c         do j=1,nmollst_tmp(i)
c            mollst(j,i)=mollst_tmp(j,i)
c         end do
      !end do

c      deallocate (nmollst_tmp)
c      deallocate (mollst_tmp)
c      deallocate (nmollst2_tmp)
c      deallocate (mollst2_tmp)
      print*,"End of molfull2bodyOO6_8_par_new"
      return
      end
