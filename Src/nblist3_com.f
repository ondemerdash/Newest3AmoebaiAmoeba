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
      subroutine molfullNew3Bodies_com
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
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
c      integer k1,kgy1,kgz1,j1,done(nmol,nmol,nmol)
      integer k1,kgy1,kgz1,j1,oldcount
      real*8 xi,yi,zi,M1,xl1,yl1,zl1,xcm,ycm,zcm
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,x1,y1,z1,r1,r3
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8  molbufshell,lowr_2bod
      real*8 lowr1,lowr2,lowr3,r1_2,r2_2,r3_2
      logical repeat,repeat2,done(nmol),empty2,empty3,done2
      integer lowk,lowk1,lowk_2bod,prox_counter2,prox_counter3
      real*8 lbuffer2b
      integer moli1,moli2

c
c
c     perform dynamic allocation of some local arrays
c
c      print*, "Before allocation Molfull3bodymod4"
c      doclust = .false.
      nlight = (ncell+1) * nmol
      print*, "Before allocation molfullNew3Bodies! 3!"
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
        if(done(i).eq..true.) goto 31
            empty2=.true.
            empty3=.true.
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
            prox_counter2=0
c            prox_counter3=0
            do j = start, stop
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
                 xcm=(xi+xsort(j))/2.0d0
                 ycm=(yi+ysort(kgy))/2.0d0
                 zcm=(zi+zsort(kgz))/2.0d0
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


                  !xr = xi - xsort(j1)
                  !yr = yi - ysort(kgy1)
                  !zr = zi - zsort(kgz1)

                  !call imagen (xr,yr,zr)

                  !xr2=xr
                  !yr2=yr
                  !zr2=zr

                  !xr = xsort(j) - xsort(j1)
                  !yr = ysort(kgy) - ysort(kgy1)
                  !zr = zsort(kgz) - zsort(kgz1)

                  xr = xcm - xsort(j1)
                  yr = ycm - ysort(kgy1)
                  zr = zcm - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr3 = xr
                  yr3 = yr
                  zr3 = zr
                  
c                r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
c                r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
c                r2=sqrt(r2_2)
                r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
                r3=sqrt(r3_2)

               !if( (r1.le.molbufshell).and.((r2.le.molbufshell).or.
     &         ! (r3.le.molbufshell)) ) then

               if(r3.le.molbufshell) then
                  if(prox_counter3.eq.0) then
                     lowr1=r1
                     !lowr2=r2
                     lowr3=r3
                     lowk=k
                     lowk1=k1
                  else
                     if(r1.lt.lowr1) then
                      lowr1=r1
                      lowk=k
                     end if
                   !  if((r2.le.molbufshell).and.(r2.lt.lowr2)) then
                   !   lowr2=r2
                   !   lowk1=k1
                   !  end if
                     if((r3.le.molbufshell).and.(r3.lt.lowr3)) then
                      lowr3=r3
                      lowk1=k1
                     end if 
                  end if
                  empty3=.false.
                  prox_counter3=prox_counter3+1
               end if
   40       continue
          
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nlight
                  goto 30
               end if

c            else
c             goto 20
c            end if

   20       continue
c            print*, "In inner loop i=",i," j=",j
            end do

            if (empty3.eq..false.) then
               if((done(i).eq..false.).and.(done(lowk).eq..false.)
     &           .and.(done(lowk1).eq..false.)) then
                 clustcount=clustcount+1
                 clust(1,clustcount)=i
                 clust(2,clustcount)=lowk
                 clust(3,clustcount)=lowk1
                 sizeclust(clustcount)=3
                 done(i)=.true.
                 done(lowk)=.true.
                 done(lowk1)=.true.
               end if    
            end if
            if ( (empty2.eq..false.).and.(empty3.eq..true.) ) then
               if((done(i).eq..false.).and.(done(lowk_2bod).eq..false.)
     &          ) then
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

   31      continue
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

