c
      subroutine mollist2bodyOO6_8_par_1list_par_newesttink7
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
c      do i = 1, nmol
      do i = start_polar(taskid),last_polar(taskid)
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
c       do i = 1, nmol
      do i = start_polar(taskid),last_polar(taskid)
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

c            do k = 1, i-1
            do k =start_polar(taskid),i-1
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
c          do i=1,nmol
          do i =start_polar(taskid),last_polar(taskid)
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


