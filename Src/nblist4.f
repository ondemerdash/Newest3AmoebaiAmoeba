c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfullNew4Bodies
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
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol,oldcount
      integer start,stop,start2,stop2,start3,stop3
      integer k1,kgy1,kgz1,j1,k2,j2,kgy2,kgz2
      real*8 xi,yi,zi,M1,xl1,yl1,zl1,xcm,ycm,zcm
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 xr4,yr4,zr4,xr5,yr5,zr5,xr6,yr6,zr6
      real*8 r2,off,x1,y1,z1,r1,r3,r4,r5,r6
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8  molbufshell,r1_2,r2_2,r3_2,r4_2,r5_2,r6_2
      integer lowk1_3bod,lowk_3bod,lowk_2bod
      integer lowk_4bod,lowk1_4bod,lowk2_4bod
      real*8 lowr_2bod,lowr1_3bod,lowr2_3bod
      real*8 lowr3_3bod,lowr1_4bod,lowr2_4bod,lowr3_4bod
      real*8 lowr4_4bod,lowr5_4bod,lowr6_4bod
      logical repeat,repeat2,repeat3,done(nmol),empty2,empty3,done2
      logical empty4
      integer lowk,prox_counter2,prox_counter3,prox_counter4
      real*8 lbuffer2b
      integer moli1,moli2

c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation MolfullNew4Bodies! 4!!"
c      doclust = .false.
      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c      print*, "After xyzsort(nlight) alloction"
c      allocate (done(nmol,nmol,nmol))
c      print*, "After done allocation Molfull3bodymod4"
c
c     transfer interaction site coordinates to sorting arrays
c

c      do i = 1, nmol
c       nmollst3mod(i)=0
c      end do

c      do i=1,nmol
c         do j=1,nmol
c            do k=1,nmol
c               done(i,j,k)=0
c            end do
c         end do
c      end do

      do i = 1, nmol
         done(i)=.false.     
         sizeclust(i)=0
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

c      allocate (done(nmol,nmol,nmol))
c
c     use the method of lights to generate neighbors
c      print*,"Before lights Call in Molfull3bodymod4"
c      molbuf2=4.25+3
c      molbuf2=8.5+2
c      molbuf2=4
c      off = sqrt(molbuf2)
      !molbufshell=4.0d0
      !off =molbufshell
      lbuffer2b=lbuffer
      !molbufshell=4.0d0+lbuffer2b
      molbufshell=clustcut
      print*,"clustcut=",clustcut
      off =molbufshell

c      off=molbuf2
      call lights (off,nmol,xsort,ysort,zsort)

c      print*,"After lights Call in Molfull3bodymod4"
c
c     loop over all atoms computing the interactions
c
      clustcount=0

      do i = 1, nmol
c         print*, "In loop i=",i
        if(done(i).eq..true.) goto 21
            empty2=.true.
            empty3=.true.
            empty4=.true.
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

c            print*, "In inner loop i=",i

   10    continue
            prox_counter2=0
            do j = start, stop
c            if(j .ne. i) then
               k = locx(j)
               if(done(k).eq..true.) goto 20
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

               xr1=xr
               yr1=yr
               zr1=zr
               r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
               r1=sqrt(r1_2)
       
               if(r1.le.molbufshell) then
                 if(prox_counter2.eq.0) then
                   lowr_2bod=r1
                   lowk_2bod=k
                 else
                   if(r1.lt.lowr_2bod) then
                     lowr_2bod=r1
                     lowk_2bod=k  
                   end if 
                 end if
                 empty2=.false.
                 prox_counter2=prox_counter2+1
               end if

               if (kbx(k) .le. kex(k)) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = kex(k)
               else
                  repeat2 = .true.
                  start2 = 1
                  stop2 = kex(k)
               end if

   30    continue
               prox_counter3=0
               do j1 = start2, stop2

                  k1=locx(j1)
                  if(done(k1).eq..true.) goto 40
                  kgy1 =rgy(k1)                  
                  if (kby(k) .le. key(k)) then
                   if(kgy1.lt.kby(k)
     &                .or. kgy1.gt.key(k)) goto 40
                  else
                   if (kgy1.lt.kby(k)
     &               .and. kgy1.gt.key(k)) goto 40
                  end if
                  kgz1=rgz(k1)
                  if (kbz(k) .le. kez(k)) then
                   if (kgz1.lt.kbz(k)
     &                 .or. kgz1.gt.kez(k))  goto 40
                  else
                   if (kgz1.lt.kbz(k)
     &                 .and. kgz1.gt.kez(k))  goto 40
                  end if

                  xr = xi - xsort(j1)
                  yr = yi - ysort(kgy1)
                  zr = zi - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xsort(j) - xsort(j1)
                  yr = ysort(kgy) - ysort(kgy1)
                  zr = zsort(kgz) - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr3 = xr
                  yr3 = yr
                  zr3 = zr
                  
c                r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
                r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
                r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
                r2=sqrt(r2_2)
                r3=sqrt(r3_2)

               if( (r1.le.molbufshell).and.((r2.le.molbufshell).or.
     &          (r3.le.molbufshell)) ) then
                  if(prox_counter3.eq.0) then
                     lowr1_3bod=r1
                     lowr2_3bod=r2
                     lowr3_3bod=r3
                     lowk_3bod=k
                     lowk1_3bod=k1
                  else
                     if(r1.lt.lowr1_3bod) then
                      lowr1_3bod=r1
                      lowk_3bod=k
                     end if
                     if((r2.le.molbufshell).and.(r2.lt.lowr2_3bod)) then
                      lowr2_3bod=r2
                      lowk1_3bod=k1
                     end if
                     if((r3.le.molbufshell).and.(r3.lt.lowr3_3bod)) then
                      lowr3_3bod=r3
                      lowk1_3bod=k1
                     end if 
                  end if
                  empty3=.false.
                  prox_counter3=prox_counter3+1
               end if

                   if (kbx(k1) .le. kex(k1)) then
                      repeat3 = .false.
                      start3 = kbx(k1) + 1
                      stop3 = kex(k1)
                   else
                      repeat3 = .true.
                      start3 = 1
                      stop3 = kex(k1)
                   end if

   50    continue
                   prox_counter4=0
                   do j2 = start3, stop3

                      k2=locx(j2)
                      if(done(k2).eq..true.) goto 60
                      kgy2 =rgy(k2)
                      if (kby(k1) .le. key(k1)) then
                       if(kgy2.lt.kby(k1)
     &                  .or. kgy2.gt.key(k1)) goto 60
                      else
                       if (kgy2.lt.kby(k1)
     &                  .and. kgy2.gt.key(k1)) goto 60
                      end if
                      kgz2=rgz(k2)
                      if (kbz(k1) .le. kez(k1)) then
                       if (kgz2.lt.kbz(k1)
     &                  .or. kgz2.gt.kez(k1))  goto 60
                      else
                       if (kgz2.lt.kbz(k1)
     &                  .and. kgz2.gt.kez(k1))  goto 60
                      end if

                      xr = xi - xsort(j2)
                      yr = yi - ysort(kgy2)
                      zr = zi - zsort(kgz2)

                      call imagen (xr,yr,zr)

                      xr4=xr
                      yr4=yr
                      zr4=zr

                      xr = xsort(j) - xsort(j2)
                      yr = ysort(kgy) - ysort(kgy2)
                      zr = zsort(kgz) - zsort(kgz2)

                      call imagen (xr,yr,zr)

                      xr5 = xr
                      yr5 = yr
                      zr5 = zr

                      xr = xsort(j1) - xsort(j2)
                      yr = ysort(kgy1) - ysort(kgy2)
                      zr = zsort(kgz1) - zsort(kgz2)

                      xr6 = xr
                      yr6 = yr
                      zr6 = zr
                    r4_2=xr4*xr4 + yr4*yr4 + zr4*zr4
                    r4=sqrt(r4_2)
                    r5_2=xr5*xr5 + yr5*yr5 + zr5*zr5
                    r5=sqrt(r5_2)
                    r6_2=xr6*xr6 + yr6*yr6 + zr6*zr6
                    r6=sqrt(r6_2)
                   if( (r1.le.molbufshell).and.((r2.le.molbufshell).or.
     &              (r3.le.molbufshell)).and.((r4.le.molbufshell).or.
     &              (r5.le.molbufshell).or.(r6.le.molbufshell)) )then
                      if(prox_counter4.eq.0) then
                         lowr1_4bod=r1
                         lowr2_4bod=r2
                         lowr3_4bod=r3
                         lowr4_4bod=r4
                         lowr5_4bod=r5
                         lowr6_4bod=r6
                         lowk_4bod=k
                         lowk1_4bod=k1
                         lowk2_4bod=k2
                      else
                         if(r1.lt.lowr1_4bod) then
                          lowr1_4bod=r1
                          lowk_4bod=k
                         end if
                         if((r2.le.molbufshell).and.
     &                   (r2.lt.lowr2_4bod)) then
                          lowr2_4bod=r2
                          lowk1_4bod=k1
                         end if
                         if((r3.le.molbufshell).and.
     &                     (r3.lt.lowr3_4bod)) then
                          lowr3_4bod=r3
                          lowk1_4bod=k1
                         end if
                         if((r4.le.molbufshell).and.
     &                     (r4.lt.lowr4_4bod)) then
                          lowr4_4bod=r4
                          lowk2_4bod=k2
                         end if
                         if((r5.le.molbufshell).and.
     &                     (r5.lt.lowr5_4bod)) then
                          lowr5_4bod=r5
                          lowk2_4bod=k2
                         end if
                         if((r6.le.molbufshell).and.
     &                    (r6.lt.lowr6_4bod)) then
                          lowr6_4bod=r6
                          lowk2_4bod=k2
                         end if
                      end if
                      empty4=.false.
                      prox_counter4=prox_counter4+1
                   end if
   60       continue
                   end do
                   if (repeat3) then
                      repeat3 = .false.
                      start3 = kbx(k1) + 1
                      stop3 = nlight
                      goto 50
                   end if

   40       continue          
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nlight
                  goto 30
               end if

   20       continue
c            print*, "In inner loop i=",i," j=",j
            end do
   
            if (empty4.eq..false.) then
               if((done(i).eq..false.).and.(done(lowk_4bod).eq..false.)
     &           .and.(done(lowk1_4bod).eq..false.).and.
     &           (done(lowk2_4bod).eq..false.) ) then
                 clustcount=clustcount+1
                 clust(1,clustcount)=i
                 clust(2,clustcount)=lowk_4bod
                 clust(3,clustcount)=lowk1_4bod
                 clust(4,clustcount)=lowk2_4bod
                 sizeclust(clustcount)=4
                 done(i)=.true.
                 done(lowk_4bod)=.true.
                 done(lowk1_4bod)=.true. 
                 done(lowk2_4bod)=.true.
               end if     
            end if 
            if ((empty3.eq..false.).and.(empty4.eq..true.) ) then
               if((done(i).eq..false.).and.(done(lowk_3bod).eq..false.)
     &           .and.(done(lowk1_3bod).eq..false.)) then
                 clustcount=clustcount+1
                 clust(1,clustcount)=i
                 clust(2,clustcount)=lowk_3bod
                 clust(3,clustcount)=lowk1_3bod
                 sizeclust(clustcount)=3
                 done(i)=.true.
                 done(lowk_3bod)=.true.
                 done(lowk1_3bod)=.true.
               end if    
            end if
            if ( (empty2.eq..false.).and.(empty3.eq..true.).and.
     &       (empty4.eq..true.) ) then
               if( (done(i).eq..false.).and.
     &          (done(lowk_2bod).eq..false.) ) then
                clustcount=clustcount+1
                clust(1,clustcount)=i
                clust(2,clustcount)=lowk_2bod
                sizeclust(clustcount)=2
                done(i)=.true.
                done(lowk_2bod)=.true.
               end if
            end if


            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 10
            end if

   21      continue
      end do

      oldcount=clustcount
      print*,"Oldcount",oldcount

      do i=1,nmol
         done2=.false.
         do j=1,oldcount
            do k=1,sizeclust(j)
               if(i.eq.clust(k,j)) then
                 done2=.true.
                 goto 32
               end if 
            end do
         end do 
         if(done2.eq..false.) then
           clustcount=clustcount+1
           clust(1,clustcount)=i
           sizeclust(clustcount)=1
         end if      
   32      continue
      end do 

      print*,"Clustcount",clustcount


c     perform deallocation of some local arrays
c
      do i=1,clustcount
         xcm=0.0d0
         ycm=0.0d0
         zcm=0.0d0
         print*,"sizeclust(i) i ",i,sizeclust(i)
         do j=1,sizeclust(i)
            moli1=clust(j,i)
            print*,"clust i, members",i,moli1
            xcm=xcm+xsort(rgx(moli1))
            ycm=ycm+ysort(rgy(moli1))
            zcm=zcm+zsort(rgz(moli1))
         end do
         clust_cm(1,i)=xcm/sizeclust(i)
         clust_cm(2,i)=ycm/sizeclust(i)
         clust_cm(3,i)=zcm/sizeclust(i)
      end do

      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c      deallocate (done)
c
c     check to see if the neighbor lists are too long
c
c      do i = 1, nmol
c         if (nmollst3mod(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')            
c            call fatal
c            print*, "Too many in nmollst3mod!!"
c         end if
c      end do
      return
      end

