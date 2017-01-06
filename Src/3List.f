      subroutine mollist3bodycobarPar
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use neigh3b
      use molcul
      use atomid
      use neigh
      use cobar    
      implicit none
      integer i,j,k,ii,l1,lenmol,j1,j2,k1,k2
      real*8 xi,yi,zi,xcm1,ycm1,zcm1
      real*8 xr,yr,zr,M1,lbuffer2b,lbuf2_2b
      real*8 radius,r,xk,yk,zk
      real*8 xk1,yk1,zk1,r1,r2,r3,r1_2,r2_2,r3_2
      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3,shellsum
      logical, allocatable :: update(:)
      lbuffer2b=lbuffer
      molbufcobar=Cobarcut3b+lbuffer2b
      molbufcobar_x=Cobarcut3b+2.0d0*lbuffer2b
      lbuf2_2b=(0.5d0*lbuffer2b)**2
      allocate (update(nmol))
      print*,"in mollist3bodycobarPar",Cobarcut3b

      radius = sqrt(molbufcobar)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' MOLLST  --  3Bod Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if

      if (domollst3bod) then
         domollst3bod = .false.
         call molfull3bodycobarPar
         return
      end if

!$OMP PARALLEL default(shared) private(i,lenmol,l1,j,k,xi,yi,zi,
!$OMP& xr,yr,zr,r2,r,j1,j2,k1,k2,xk,yk,zk,xk1,yk1,zk1,xr1,yr1,zr1,
!$OMP& xr2,yr2,zr2,
!$OMP& xr3,yr3,zr3,r1,r3,r1_2,r2_2,r3_2,shellsum)
!$OMP DO schedule(guided)
      do i = 1, nmol
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k2=l1-1
            j=imol(1,i)+k2
           if(name(j).eq.'O') then
            xi = x(j)
            yi = y(j)
            zi = z(j)
           end if
         end do

         xr = xi - xmolold3(i)
         yr = yi - ymolold3(i)
         zr = zi - zmolold3(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2_2b) then
            listbcast3b =.true.
c            print*,"listbcast=",listbcast
            update(i) = .true.
            xmolold3(i) = xi
            ymolold3(i) = yi
            zmolold3(i) = zi
         end if
      end do
!$OMP END DO

!$OMP DO schedule(guided)
       do i = 1, nmol
         if (update(i)) then
            lenmol=imol(2,i)-imol(1,i)+1

            do l1 = 1, lenmol
             k2=l1-1
             j=imol(1,i)+k2
             if(name(j).eq.'O') then
              xi = x(j)
              yi = y(j)
              zi = z(j)
             end if
            end do
            j1 = 0
 
            do k = i+1, nmol-1
               lenmol=imol(2,k)-imol(1,k)+1

               do l1 = 1, lenmol
                  k2=l1-1
                  j=imol(1,k)+k2
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
               xr1=xr
               yr1=yr
               zr1=zr

               do k1=k+1,nmol
                  lenmol=imol(2,k1)-imol(1,k1)+1
 
                  do l1 = 1, lenmol
                     k2=l1-1
                     j=imol(1,k1)+k2
                     if(name(j).eq.'O') then
                       xk1 = x(j)
                       yk1 = y(j)
                       zk1 = z(j)
                     end if
                  end do

                  xr = xi - xk1
                  yr = yi - yk1
                  zr = zi - zk1
                  call imagen (xr,yr,zr)
                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xk - xk1
                  yr = yk - yk1
                  zr = zk - zk1
                  call imagen (xr,yr,zr)
                  xr3=xr
                  yr3=yr
                  zr3=zr

                  r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
                  r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
                  r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
                  r1=sqrt(r1_2)
                  r2=sqrt(r2_2)
                  r3=sqrt(r3_2)

                  if( (r1.lt.r3).and.(r2.lt.r3) ) then
                    shellsum=r1+r2
                  else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                    shellsum=r1+r3
                  else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                    shellsum=r2+r3
                  end if


                  if (shellsum .le. molbufcobar) then
                    ! nmollst3(i)=nmollst3(i)+1
                      j1=j1+1
                    ! mollst3(nmollst3(i),i)=k
                      mollst3(j1,i)=k
                    ! nmollst3(i)=nmollst3(i)+1
                      j1=j1+1
                     mollst3(j1,i)=k1
                  end if
               
               end do 
            end do 
          
            nmollst3(i)=j1

            do k=1,i-2
               if (.not. update(k)) then
                  xr = xi - xmolold3(k)
                  yr = yi - ymolold3(k)
                  zr = zi - zmolold3(k)
                  call imagen (xr,yr,zr)
                  xr1=xr
                  yr1=yr
                  zr1=zr
                  r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1

                  r1=sqrt(r1_2)

                  do k1=k+1,i-1
                     if (.not. update(k1)) then
                        xr = xi - xmolold3(k1)
                        yr = yi - ymolold3(k1)
                        zr = zi - zmolold3(k1)
                        call imagen (xr,yr,zr)
                        xr2=xr
                        yr2=yr
                        zr2=zr
                        r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2

                        r2=sqrt(r2_2)

                        xr = xmolold3(k) - xmolold3(k1)
                        yr = ymolold3(k) - ymolold3(k1)
                        zr = zmolold3(k) - zmolold3(k1)
                        call imagen (xr,yr,zr)
                        xr3=xr
                        yr3=yr
                        zr3=zr
                        r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3

                        r3=sqrt(r3_2)

                        if( (r1.lt.r3).and.(r2.lt.r3) ) then
                          shellsum=r1+r2
                        else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                          shellsum=r1+r3
                        else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                          shellsum=r2+r3
                        end if
              
                        if (shellsum .le. molbufcobar) then
                           do j =1,nmollst3(k),2
                              if(mollst3(j-1,k).eq.k1 .and. 
     &                          mollst3(j,k).eq.i ) goto 20
                           end do

!$OMP CRITICAL
                          nmollst3(k) = nmollst3(k) + 1
                          mollst3(nmollst3(k),k) = k1
                          nmollst3(k) = nmollst3(k) + 1
                          mollst3(nmollst3(k),k) = i
!$OMP END CRITICAL

   20             continue

                        else if (shellsum .le. molbufcobar_x) then
                           do j =1,nmollst3(k),2
                              if(mollst3(j-1,k).eq.k1 .and.
     &                          mollst3(j,k).eq.i ) then
!$OMP CRITICAL
                                mollst3(j-1,k)=
     &                                    mollst3(nmollst3(k)-1,k)
             
                                mollst3(j,k)=
     &                                    mollst3(nmollst3(k),k)
                                nmollst3(k)=nmollst3(k)-2
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

         end if
       end do
!$OMP END DO

!$OMP DO
          do i=1,nmol
             if (nmollst3(i).gt.max3blst) then
            print*,"nmollst3 Too many 3b nbs.Tilt!",nmollst3(i)
              call fatal
             end if
          end do 
!$OMP END DO
!$OMP END PARALLEL
      
      print*,"Done With mollist3bodycobarPar (REBUILD)"
      deallocate (update)
      return
      end

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
      subroutine molfull3bodycobarPar
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
      use cobar
      implicit none
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
      integer k1,kgy1,kgz1,j1
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,x1,y1,z1,r1,r3
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8  shellsum,r1_2,r2_2,r3_2
      integer shell1,shell2,shell3,shellsum_mod
      real*8 lbuffer2b
      integer count_i,count_k,count_k1
      logical repeat,repeat2
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(nmol))
      allocate (ysort(nmol))
      allocate (zsort(nmol))
c
c     transfer interaction site coordinates to sorting arrays
c

      do i = 1, nmol
         nmollst3(i)=0
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


         xmolold3(i) = x1
         ymolold3(i) = y1
         zmolold3(i) = z1
         xsort(i) = x1
         ysort(i) = y1
         zsort(i) = z1
      end do

c
      lbuffer2b=lbuffer
      print*,"Cobarcut3b in molfull3bodycobarPar",Cobarcut3b
      molbufcobar=Cobarcut3b+lbuffer2b
      off =molbufcobar
      call lightn (off,nmol,xsort,ysort,zsort)

      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)

c      print*,"After lights Call in Molfull3bodymod4"
c
c     loop over all atoms computing the interactions
c

!$OMP PARALLEL DO default(private) shared(nmol,mollst3,
!$OMP& nmollst3,xmolold3,ymolold3,zmolold3,kbx,kex,
!$OMP& kby,key,kbz,kez,rgy,rgz,locx,molbufcobar)
!$OMP& schedule(guided)
      do i = 1, nmol
            count_i=0
            xi = xmolold3(i)
            yi = ymolold3(i)
            zi = zmolold3(i)


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
               count_k=0
               
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

               xr = xi - xmolold3(k)
               yr = yi - ymolold3(k)
               zr = zi - zmolold3(k)
               call imagen (xr,yr,zr)
               xr1=xr
               yr1=yr
               zr1=zr

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


               do j1 = start2, stop2
                  count_k1=0

                  k1=locx(j1)
                  if(k1 .le. k) goto 40
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

                  xr = xi - xmolold3(k1)
                  yr = yi - ymolold3(k1)
                  zr = zi - zmolold3(k1)
                  call imagen (xr,yr,zr)
                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xmolold3(k) - xmolold3(k1)
                  yr = ymolold3(k) - ymolold3(k1)
                  zr = zmolold3(k) - zmolold3(k1)
                  call imagen (xr,yr,zr)
                  xr3 = xr
                  yr3 = yr
                  zr3 = zr

                  r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
                  r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
                  r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3      
                  r1=sqrt(r1_2)
                  r2=sqrt(r2_2)
                  r3=sqrt(r3_2)

               if( (r1.lt.r3).and.(r2.lt.r3) ) then
                 shellsum=r1+r2
               else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                 shellsum=r1+r3
               else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                 shellsum=r2+r3
               end if


               if (shellsum .le. molbufcobar) then
c                  if (i.lt.k1 .and. k1.lt.k) then
c                     nmollst3(i)=nmollst3(i)+1
c                     mollst3(nmollst3(i),i)=k1
c                     nmollst3(i)=nmollst3(i)+1
c                     mollst3(nmollst3(i),i)=k                     
c                  end if
c                  if (i.lt.k .and. k.lt.k1) then
                     nmollst3(i)=nmollst3(i)+1
                     mollst3(nmollst3(i),i)=k
                     nmollst3(i)=nmollst3(i)+1                    
                     mollst3(nmollst3(i),i)=k1
c                  end if
c                  if (k.lt.i .and. i.lt.k1) then
c                     nmollst3(k)=nmollst3(k)+1
c                     mollst3(nmollst3(k),k)=i
c                     nmollst3(k)=nmollst3(k)+1
c                     mollst3(nmollst3(k),k)=k1
c                  end if
c                  if (k.lt.k1 .and. k1.lt.i) then
c                     nmollst3(k)=nmollst3(k)+1
c                     mollst3(nmollst3(k),k)=k1
c                     nmollst3(k)=nmollst3(k)+1
c                     mollst3(nmollst3(k),k)=i
c                  end if                     
c                  if (k1.lt.i .and. i.lt.k) then
c                     nmollst3(k1)=nmollst3(k1)+1
c                     mollst3(nmollst3(k1),k1)=i
c                     nmollst3(k1)=nmollst3(k1)+1
c                     mollst3(nmollst3(k1),k1)=k
c                  end if
c                  if (k1.lt.k .and. k.lt.i) then
c                     nmollst3(k1)=nmollst3(k1)+1
c                     mollst3(nmollst3(k1),k1)=k
c                     nmollst3(k1)=nmollst3(k1)+1
c                     mollst3(nmollst3(k1),k1)=i
c                  end if
               end if
   40       continue
          
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nmol
                  goto 30
               end if

   20       continue
c            print*, "In inner loop i=",i," j=",j
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nmol
               goto 10
            end if

      end do
!$OMP END PARALLEL DO


c     perform deallocation of some local arrays
c
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nmol
         if (nmollst3(i) .ge. max3blst) then
            write (iout,31)
   31       format (/,' MOLFULL3  --  Too many Neighbors;',
     &                 ' Increase MAX3BLST')            
      print*, "Too many in nmollst3! Tilt! nmollst3=",nmollst3(i)
            call fatal
         end if
      end do
      return
      end



