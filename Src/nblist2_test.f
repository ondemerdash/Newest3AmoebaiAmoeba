c
      subroutine molfullNew2BodiesCobar_test
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
      integer k1,kgy1,kgz1,j1
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1,xcm,ycm,zcm
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,x1,y1,z1,r1,r3,r,lowr
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8  molbuf_cobar,shellsum,dist,molbufshell
      integer shell1,shell2,shell3,shellsum_mod,oldcount
      logical repeat,repeat2,done(nmol),empty,done2
      integer lowk,prox_counter
      real*8 lbuffer2b
      integer moli1,moli2
c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation molfullNew2BodiesCobar"
c      doclust = .false.
      nlight = (ncell+1) * nmol
c      nlight =  nmol

      print*,"Nlight",nlight
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c      allocate (done(nmol,nmol,nmol))
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

         xmolold2(i) = x1
         ymolold2(i) = y1
         zmolold2(i) = z1
         xsort(i) = x1
         ysort(i) = y1
         zsort(i) = z1
      end do

      lbuffer2b=lbuffer
      !molbufshell=4.0d0+lbuffer2b
      molbufshell=clustcut
      off =molbufshell
      print*,"clustcut=",clustcut

c      off=molbuf2
      !call lights (off,nmol,xsort,ysort,zsort)
      call lightn (off,nmol,xsort,ysort,zsort)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      print*,"After lights Call in Cobar"
c
c     loop over all atoms computing the interactions
c

      clustcount=0

      do i = 1, nmol
        
        if(done(i).eqv..true.) goto 31
         empty=.true.
         !xi = xsort(rgx(i))
         !yi = ysort(rgy(i))
         !zi = zsort(rgz(i))
          xi = xmolold2(i)
          yi = ymolold2(i)
          zi = zmolold2(i)

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

         prox_counter=0         
         do j = start, stop
            k = locx(j)
c           if(done(k).eq..true.) goto 20
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
            !xr = xi - xsort(j)
            !yr = yi - ysort(kgy)
            !zr = zi - zsort(kgz)
            xr = xi - xmolold2(k)
            yr = yi - ymolold2(k)
            zr = zi - zmolold2(k)

            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            
            if (r .le. molbufshell) then
               if(prox_counter.eq.0) then
                 lowr=r
                 lowk=k
               else
                  if(r.lt.lowr) then
                    lowr=r
                    lowk=k
                  end if
               end if
               empty=.false.
               prox_counter=prox_counter+1
            end if                      
c               if (i .lt. k) then
c                  nmollst(i) = nmollst(i) + 1
c                  mollst(nmollst(i),i) = k
c               else
c                  nmollst(k) = nmollst(k) + 1
c                  mollst(nmollst(k),k) = i
c               end if
   20       continue
         end do
 
            if (empty.eqv..false.) then
c             clustcount=clustcount+1
               if((i .lt. lowk).and.(done(i).eqv..false.)
     &            .and.(done(lowk).eqv..false.)) then
                  clustcount=clustcount+1                 
                  clust(1,clustcount)=i
                  clust(2,clustcount)=lowk
c                  nmollst(i) = nmollst(i) + 1
c                  mollst(nmollst(i),i) = k
                  sizeclust(clustcount)=2
                  done(i)=.true.
                  done(lowk)=.true.
               end if
               if((lowk .lt. i).and.(done(i).eqv..false.)
     &            .and.(done(lowk).eqv..false.)) then
                  clustcount=clustcount+1
                  clust(1,clustcount)=lowk
                  clust(2,clustcount)=i
c                  nmollst(k) = nmollst(k) + 1
c                  mollst(nmollst(k),k) = i
                  sizeclust(clustcount)=2
                  done(i)=.true.
                  done(lowk)=.true.
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
c
c     perform deallocation of some local arrays
c
      oldcount=clustcount
      print*,"Oldcount",oldcount

      do i=1,nmol
         done2=.false.
         do j=1,oldcount
            if((i.eq.clust(1,j)).or.(i.eq.clust(2,j))) then
              done2=.true.
              goto 32
            end if           
         end do 
         if(done2.eqv..false.) then
          clustcount=clustcount+1
          clust(1,clustcount)=i
          sizeclust(clustcount)=1
         end if
   32      continue         
      end do

      print*,"Clustcount",clustcount

      do i=1,clustcount
         if(sizeclust(i).eq.1) then
           distmax(i)=10.0d0
         else
           do j=1,sizeclust(i)-1
              do k=j+1,sizeclust(i)
                moli1=clust(j,i)
                moli2=clust(k,i)
            !     dist=sqrt( (xsort(rgx(moli1))-xsort(rgx(moli2)))
     &      ! *(xsort(rgx(moli1))-xsort(rgx(moli2)))
     &      ! + (ysort(rgy(moli1))-ysort(rgy(moli2)))
     &      ! *(ysort(rgy(moli1))-ysort(rgy(moli2)))
     &      ! + (zsort(rgz(moli1))-zsort(rgz(moli2)))
     &      ! *(zsort(rgz(moli1))-zsort(rgz(moli2))) )     

                 dist=sqrt( (xmolold2(moli1)-xmolold2(moli2))
     &       *(xmolold2(moli1)-xmolold2(moli2))
     &       + (ymolold2(moli1)-ymolold2(moli2))
     &       *(ymolold2(moli1)-ymolold2(moli2))
     &       + (zmolold2(moli1)-zmolold2(moli2))
     &       * (zmolold2(moli1)-zmolold2(moli2)) )

                if (j.eq.1.and.k.eq.2) then
                 distmax(i)=dist
                end if
                if(dist.gt.distmax(i)) then
                  distmax(i)=dist
                end if
              end do 
           end do
         end if 
      end do 

      do i=1,clustcount
         xcm=0.0d0
         ycm=0.0d0
         zcm=0.0d0
         print*,"sizeclust(i) i ",i,sizeclust(i) 
         do j=1,sizeclust(i)
            moli1=clust(j,i)
            print*,"clust i, members",i,moli1
           ! xcm=xcm+xsort(rgx(moli1))
           ! ycm=ycm+ysort(rgy(moli1))
           ! zcm=zcm+zsort(rgz(moli1))
            xcm=xcm+xmolold2(moli1)
            ycm=ycm+ymolold2(moli1)
            zcm=zcm+zmolold2(moli1)
         end do 
         clust_cm(1,i)=xcm/sizeclust(i)
         clust_cm(2,i)=ycm/sizeclust(i)
         clust_cm(3,i)=zcm/sizeclust(i)
      end do 



      !deallocate (xsort)
      !deallocate (ysort)
      !deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
c      do i = 1, nmol
c         if (nmollst(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')
c            call fatal
c         end if
c      end do
      return
      end

