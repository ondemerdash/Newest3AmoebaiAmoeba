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
      subroutine molfull2body_NewBodieskmeans
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
      use cobar
      use aprx
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
      real*8 r1_2,xr1,yr1,zr1,r1
      integer clust1,clust2
      integer clust3
      real*8 r2_2,xr2,yr2,zr2
      real*8 r3_2,xr3,yr3,zr3,r3
      integer k1,k2,counter
      logical in_clust_list
c
c
c     perform dynamic allocation of some local arrays
c

c
c
      do i = 1, clustcount
         nmollst(i) = 0         
      end do

c
c
      listbcast=.true.
      lbuffer2b=lbuffer
      molbuf2s=cut2b_input+lbuffer2b
c      molbuf2s=xcell/2.0d0
c      off = sqrt(molbuf2)
      off=molbuf2
      off=molbuf2s

      if(approxmode.eq.'3BODYMODE') then
        do clust1=1,clustcount
           do clust2=1,clustcount
               save2bfor3b(clust2,clust1)=.false.              
           end do
        end do
      end if

      do clust1=1,clustcount-1
        do clust2=clust1+1,clustcount
           xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
           yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
           zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
           call image(xr1,yr1,zr1)
           r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
           r1=sqrt(r1_2)

          if(r1.le.molbuf2s) then
             nmollst(clust1) = nmollst(clust1) + 1
             mollst(nmollst(clust1),clust1) = clust2
          end if
        end do
      end do   
c      print*,"After lights Call in Molfull2body"
c
c     loop over all atoms computing the interactions
c
c
c     perform deallocation of some local arrays
c
c      do i=1,clustcount
c         print*,"Clust ind, nmollst(i)",i,nmollst(i)
c         do j=1,nmollst(i)
c            print*,"Clust ind,mollst(j,i)",mollst(j,i)
c         end do 
c      end do 
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

      if(approxmode.eq.'3BODYMODE') then
        do i = 1, clustcount
           nmollst3(i) = 0
        end do
        cnt3bwork=0
        num2bsave=0
c        do clust1=1,clustcount
        do clust1=1,clustcount-2
           do clust2=clust1+1,clustcount-1
c           do k1=1,nmollst(clust1)
c              if(nmollst(clust1).ne.0) then
c                clust2=mollst(k1,clust1)
                xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
                yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
                zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
                call image(xr1,yr1,zr1)
                r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
                r1=sqrt(r1_2)
                do clust3=clust2+1,clustcount
c                do k2=1,nmollst(clust2)
c                   if(nmollst(clust2).ne.0) then
c                    clust3=mollst(k2,clust2)
                    xr2=clust_cm(1,clust3)-clust_cm(1,clust1)
                    yr2=clust_cm(2,clust3)-clust_cm(2,clust1)
                    zr2=clust_cm(3,clust3)-clust_cm(3,clust1)
                    call image(xr2,yr2,zr2)
                    r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
                    r2=sqrt(r2_2)

                    xr3=clust_cm(1,clust3)-clust_cm(1,clust2)
                    yr3=clust_cm(2,clust3)-clust_cm(2,clust2)
                    zr3=clust_cm(3,clust3)-clust_cm(3,clust2)
                    call image(xr3,yr3,zr3)
                    r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
                    r3=sqrt(r3_2)

                    if((r1.lt.r3).and.(r2.lt.r3) ) then
                       r=r1+r2
              ! print*,"c1=",clust1,"c2=",clust2,"c3=",clust3,"Cobar=",r
                    else if ((r1.lt.r2).and.(r3.lt.r2)) then
                       r=r1+r3
              ! print*,"c1=",clust1,"c2=",clust2,"c3=",clust3,"Cobar=",r
                    else if ((r2.lt.r1).and.(r3.lt.r1)) then
                       r=r2+r3
              ! print*,"c1=",clust1,"c2=",clust2,"c3=",clust3,"Cobar=",r
                    end if

                    if(r.le.(Cobarcut3b+lbuffer2b)) then 
                      nmollst3(clust1) = nmollst3(clust1) + 1
                      mollst3(nmollst3(clust1),clust1)=clust2
                      nmollst3(clust1) = nmollst3(clust1) + 1
                      mollst3(nmollst3(clust1),clust1)=clust3 
                      cnt3bwork=cnt3bwork+1
                    !  num2bsave=num2bsave+2
cccc                      save2bfor3b(k1,clust1)=.true.
cccc                      save2bfor3b(k2,clust2)=.true.
                       in_clust_list=.false.
                       do k1=1,nmollst(clust1)
                          if(mollst(k1,clust1).eq.clust2) then
                             in_clust_list=.true.
                             if(.not.save2bfor3b(k1,clust1)) then
                               num2bsave=num2bsave+1
                             end if 
                             save2bfor3b(k1,clust1)=.true.
                             goto 11
                          end if
                       end do
  11  continue
                       if(.not.in_clust_list) then
                              num2bsave=num2bsave+1
                          nmollst(clust1)=nmollst(clust1)+1
                          mollst(nmollst(clust1),clust1)=clust2
                          save2bfor3b(nmollst(clust1),clust1)=.true.
                       end if 

                       
                       in_clust_list=.false.
                       do k1=1,nmollst(clust1)
                          if(mollst(k1,clust1).eq.clust3) then
                             in_clust_list=.true.
                             if(.not.save2bfor3b(k1,clust1)) then
                               num2bsave=num2bsave+1
                             end if
                             save2bfor3b(k1,clust1)=.true.
                             goto 21
                          end if
                       end do
  21  continue
                       if(.not.in_clust_list) then
                              num2bsave=num2bsave+1
                          nmollst(clust1)=nmollst(clust1)+1
                          mollst(nmollst(clust1),clust1)=clust3
                          save2bfor3b(nmollst(clust1),clust1)=.true.
                       end if


                       in_clust_list=.false.
                       do k1=1,nmollst(clust2)
                          if(mollst(k1,clust2).eq.clust3) then
                             in_clust_list=.true.
                             if(.not.save2bfor3b(k1,clust2)) then
                               num2bsave=num2bsave+1
                             end if
                             save2bfor3b(k1,clust2)=.true.
                             goto 31
                          end if
                       end do
  31  continue
                       if(.not.in_clust_list) then
                              num2bsave=num2bsave+1
                          nmollst(clust2)=nmollst(clust2)+1
                          mollst(nmollst(clust2),clust2)=clust3
                          save2bfor3b(nmollst(clust2),clust2)=.true.
                       end if
                    end if
c                   end if 
                end do
c              end if              
           end do
        end do

c        if(.not.allocated(counter2b_ind)) 
c     &      allocate(counter2b_ind(2,num2bsave))
        !do i = 1, clustcount
        !   do j=1,nmollst3(i),2
        !    if(nmollst3(i).ne.0) then
        ! print*,"clust=",i,"neig1=",mollst3(j,i),"neig2=",mollst3(j+1,i)
        !    end if
        !   end do
        !end do

        do clust1=1,clustcount
           do clust2=1,clustcount
              ind2b_counter(clust2,clust1)=0
           end do
        end do   
  
        counter=0
        do clust1=1,clustcount
           do k1=1,nmollst(clust1)
              clust2=mollst(k1,clust1)
              if(save2bfor3b(k1,clust1)) then
                 counter=counter+1
c                 counter2b_ind(1,counter)=clust1
c                 counter2b_ind(2,counter)=clust2
                 ind2b_counter(clust2,clust1)=counter
                 
              end if
           end do
        end do
      print*,"3-body mode molfull2body_kmeans counter=",counter
      print*,"3-body mode molfull2body_kmeans num2bsave=",num2bsave
      do i = 1,clustcount
         if (nmollst3(i) .ge. 7*clustcount) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')
          print*,"nmollst3 Too many 3b nb.Tilt!",nmollst3(i)
            call fatal
         end if
      end do

      end if

      return
      end

      subroutine mollist2body_NewBodieskmeans 
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
      !radius = sqrt(molbuf2s)
      !call replica (radius)
      !if (use_replica) then
c         write (iout,10)
c   10    format (/,' MOLLST  --  Pairwise Neighbor List cannot',
c     &              ' be used with Replicas')
c         call fatal
!      end if
c
c     perform a complete list build instead of an update
c
c      print*, "Hello from Mollist2body"

      !if (domollst2bod) then
      !   domollst2bod = .false.
      if (do2waterclustlist) then 
c         if (octahedron) then
c            do i = 1, npole
c               call mbuild (i)
c            end do
c         else
c          print*, "Within Mollist2body if domollst2bod"
            call molfull2body_NewBodieskmeans 
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

