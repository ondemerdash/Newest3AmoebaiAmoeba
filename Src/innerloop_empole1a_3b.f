c
c
      subroutine Innerloop1_single_2cut_wskin(
     & moli1,ep3bt,virep3bt,dep3bt,ntript) 
c      include 'sizes3b.i'
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'atoms3b.i'
c      include 'atmtyp.i'
c      include 'molcul3b.i'
c      include 'mpole3b_3.i'
c      include 'neigh3b.i'
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      implicit none
      real*8  ep2moli12
      real*8  dep2moli12(3,30)
      integer i,ii,j,l1,i1,i2,i3,k,ntript
      integer k1,k2,k5,start,last
      real*8 eptemp,deptemp(3,30)
      real*8 vir2moli12(3,3)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3,np1,np2,np3
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      logical do2
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,xr,yr,zr
      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r1,r2,r3,r1_2,r2_2,r3_2
      real*8 shellsum

         ntript=0

                ep3bt=0.0d0
c                ep3b=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
c                      dep3b(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
c                      virep3b(i,j)=0.0d0
                      virep3bt(i,j)=0.0d0
                   end do
                end do

      do i = 1, 30
        do j = 1, 3
c          dep1moli1(j,i)=0.0d0
          dep2moli12(j,i)=0.0d0
c          dep1moli2(j,i)=0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
c           vir1moli1(i,j)=0.0d0
           vir2moli12(i,j)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 
c      print*,"Begin inner1 single"
      if(nmollst(moli1).gt.0) then
cc!$OMP PARALLEL default(private) shared(imol,ep3bt,
cc!$OMP& virep3bt,dep3bt,moli1,nmollst2,mollst2,
cc!$OMP& nmollst,mollst,name,x,y,z,ntript)
cc!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
cc!$OMP& schedule(guided)
c      do moli1=start,last
        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
          np1=3
          np2=6
          np3=9
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
            do l1 = 1, np1
              i = pnum(l1)
              if(name(i).eq.'O') then
                x1 = x(i)
                y1 = y(i)
                z1 = z(i)
              end if
            end do
            do l1 = np1+1, np2
              i = pnum(l1)
              if(name(i).eq.'O') then
                x2 = x(i)
                y2 = y(i)
                z2 = z(i)
              end if
            end do
            xr = x1 - x2
            yr = y1 - y2
            zr = z1 - z2
         !   call image(xr,yr,zr)
            xr1=xr
            yr1=yr
            zr1=zr
            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r1=sqrt(r1_2)
          if(r1.le.6.0d0) then          
           call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

           ep3bt = ep3bt + eptemp
           ep2moli12=eptemp
           do l1 = 1, npole3b
             i = pnum(l1)
             do j = 1, 3
               dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
               dep2moli12(j,l1)=deptemp(j,l1)
             end do
           end do
           do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
                vir2moli12(i,j)=virtemp(i,j)
             end do
           end do
          else
            ep2moli12=0.0d0
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
               dep2moli12(j,l1)= 0.0d0
              end do
            end do
            do i=1,3
              do j=1,3
                vir2moli12(i,j)= 0.0d0
              end do
            end do
          end if
c          do k2=1,nmollst3mod2(moli2,moli1) 
c            done(moli2)=.true.
          do k2=1,nmollst2(moli2)
             moli3=mollst2(k2,moli2)
c            moli3=mollst3mod2(k2,moli2,moli1)

            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)

            do l1 = np2+1, np3
              i = pnum(l1)
              if(name(i).eq.'O') then
                x3 = x(i)
                y3 = y(i)
                z3 = z(i)
              end if
            end do

            xr = x1 - x3
            yr = y1 - y3
            zr = z1 - z3
          !  call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr

            xr = x2 - x3
            yr = y2 - y3
            zr = z2 - z3
          !  call image(xr,yr,zr)

            xr3=xr
            yr3=yr
            zr3=zr
c            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
c            r1=sqrt(r1_2)
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)

            if( (r1.lt.r3).and.(r2.lt.r3) ) then
              shellsum=r1+r2
            else if ( (r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
            else if ( (r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
            end if

            if (shellsum .le. 8.0d0) then
c               print*,"moli1 moli2 moli3",moli1,moli2,moli3
               ntript=ntript+1
               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)
               ep3bt = ep3bt + eptemp
               do l1 = 1, npole3b
                 i = pnum(l1)
                 do j = 1, 3
                   dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
                 end do
               end do
               do i=1,3
                 do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j) 
                 end do
               end do
            
               npole3b=6
               ep3bt = ep3bt - ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)-dep2moli12(j,l1)
                  end do
               end do
               do i=1,3
                 do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)-vir2moli12(i,j)
                 end do
               end do

c               do2=.false.
c               do k5=1,nmollst(moli2)
c                 if(mollst(k5,moli2).eq.moli3) then
c                   do2=.true.
c                 goto 32
c                 end if
c               end do

                         
c   32 continue
c              if(do2.eq..true.) then
              if(r3.le.6.0d0) then
               npole3b=6
               pnum(1)=imol(1,moli2)
               pnum(2)=imol(1,moli2)+1
               pnum(3)=imol(2,moli2)
               pnum(4)=imol(1,moli3)
               pnum(5)=imol(1,moli3)+1
               pnum(6)=imol(2,moli3)
               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)
               ep3bt = ep3bt - eptemp
               do l1 = 1, npole3b
                 i = pnum(l1)
                 do j = 1, 3
                   dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
                 end do
               end do
               do i=1,3
                 do j=1,3
                  virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
                 end do
               end do
              end if

c               do2=.false.
c               do k5=1,nmollst(moli1)
c                 if(mollst(k5,moli1).eq.moli3) then
c                   do2=.true.
c                 goto 31
c                 end if
c               end do

c   31 continue

c              if(do2.eq..true.) then
              if(r2.le.6.0d0) then
               npole3b=6
               pnum(1)=imol(1,moli1)
               pnum(2)=imol(1,moli1)+1
               pnum(3)=imol(2,moli1)
               pnum(4)=imol(1,moli3)
               pnum(5)=imol(1,moli3)+1
               pnum(6)=imol(2,moli3)
               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)
               ep3bt = ep3bt - eptemp
               do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
                end do
               end do
               do i=1,3
                do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
                end do
               end do
              end if

            end if
          end do
        end do 
c      end do  
cc!$OMP END DO
cc!$OMP END PARALLEL
      end if
c      print*,"InnerloopSingle= ep3bt=",moli1,ep3bt
c       print*,"End inner1 single",ep3bt

      return
      end

c
c     ###############################################################
c     ##                                                           ##
c     ##               Subroutine empole1c_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ###############################################################
c
c
c     "empole1c_3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using the 3-body approximation
c
c
      subroutine Innerloop1_2cut_wskin(
     &  start,last,ep3bt,virep3bt,dep3bt,ntript) 
c      implicit none
c      include 'sizes3b.i'
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'atoms3b.i'
c      include 'atmtyp.i'
c      include 'molcul3b.i'
c      include 'mpole3b_3.i'
c      include 'neigh3b.i'
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      implicit none
      real*8  ep2moli12
      real*8  dep2moli12(3,30)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp,deptemp(3,30)
      real*8 vir2moli12(3,3)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3,np1,np2,np3
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      logical do2
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,xr,yr,zr
      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r1,r2,r3,r1_2,r2_2,r3_2
      real*8 shellsum

       ntript=0

                ep3bt=0.0d0
c                ep3b=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
c                      dep3b(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
c                      virep3b(i,j)=0.0d0
                      virep3bt(i,j)=0.0d0
                   end do
                end do

      do i = 1, 30
        do j = 1, 3
c          dep1moli1(j,i)=0.0d0
          dep2moli12(j,i)=0.0d0
c          dep1moli2(j,i)=0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
c           vir1moli1(i,j)=0.0d0
           vir2moli12(i,j)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 

c      print*,"Begin inner1"

cc!$OMP PARALLEL default(private) shared(imol,ep3bt,
cc!$OMP& virep3bt,dep3bt,start,last,nmollst2,mollst2,
cc!$OMP& nmollst,mollst,name,x,y,z,ntript)
cc!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
cc!$OMP& schedule(guided)
      do moli1=start,last
        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
          np1=3
          np2=6
          np3=9
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
            do l1 = 1, np1
              i = pnum(l1)
              if(name(i).eq.'O') then
                x1 = x(i)
                y1 = y(i)
                z1 = z(i)
              end if
            end do
            do l1 = np1+1, np2
              i = pnum(l1)
              if(name(i).eq.'O') then
                x2 = x(i)
                y2 = y(i)
                z2 = z(i)
              end if
            end do
            xr = x1 - x2
            yr = y1 - y2
            zr = z1 - z2
         !   call image(xr,yr,zr)
            xr1=xr
            yr1=yr
            zr1=zr

            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r1=sqrt(r1_2)
          if(r1.le.6.0d0) then
           call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

           ep3bt = ep3bt + eptemp
           ep2moli12=eptemp
           do l1 = 1, npole3b
             i = pnum(l1)
             do j = 1, 3
               dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
               dep2moli12(j,l1)=deptemp(j,l1)
             end do
           end do
           do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
                vir2moli12(i,j)=virtemp(i,j)
             end do
           end do
          else
            ep2moli12=0.0d0
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
               dep2moli12(j,l1)= 0.0d0
              end do
            end do
            do i=1,3
              do j=1,3
                vir2moli12(i,j)= 0.0d0
              end do
            end do
          end if

c          do k2=1,nmollst3mod2(moli2,moli1) 
c            done(moli2)=.true.
          do k2=1,nmollst2(moli2)
             moli3=mollst2(k2,moli2)
c            moli3=mollst3mod2(k2,moli2,moli1)

            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)

            do l1 = np2+1, np3
              i = pnum(l1)
              if(name(i).eq.'O') then
                x3 = x(i)
                y3 = y(i)
                z3 = z(i)
              end if
            end do

            xr = x1 - x3
            yr = y1 - y3
            zr = z1 - z3
         !   call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr

            xr = x2 - x3
            yr = y2 - y3
            zr = z2 - z3
         !   call image(xr,yr,zr)

            xr3=xr
            yr3=yr
            zr3=zr
c            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
c            r1=sqrt(r1_2)
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)

            if( (r1.lt.r3).and.(r2.lt.r3) ) then
              shellsum=r1+r2
            else if ( (r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
            else if ( (r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
            end if

            if (shellsum .le. 8.0d0) then
c               print*,"moli1 moli2 moli3",moli1,moli2,moli3
               ntript=ntript+1
               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)
               ep3bt = ep3bt + eptemp
               do l1 = 1, npole3b
                 i = pnum(l1)
                 do j = 1, 3
                   dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
                 end do
               end do
               do i=1,3
                 do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j) 
                 end do
               end do
            
               npole3b=6
               ep3bt = ep3bt - ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)-dep2moli12(j,l1)
                  end do
               end do
               do i=1,3
                 do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)-vir2moli12(i,j)
                 end do
               end do


c               do2=.false.
c               do k5=1,nmollst(moli2)
c                 if(mollst(k5,moli2).eq.moli3) then
c                   do2=.true.
c                 goto 32
c                 end if
c               end do

c   32 continue
c              if(do2.eq..true.) then
              if(r3.le.6.0d0) then
               npole3b=6
               pnum(1)=imol(1,moli2)
               pnum(2)=imol(1,moli2)+1
               pnum(3)=imol(2,moli2)
               pnum(4)=imol(1,moli3)
               pnum(5)=imol(1,moli3)+1
               pnum(6)=imol(2,moli3)
               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)
               ep3bt = ep3bt - eptemp
               do l1 = 1, npole3b
                 i = pnum(l1)
                 do j = 1, 3
                   dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
                 end do
               end do
               do i=1,3
                 do j=1,3
                  virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
                 end do
               end do
              end if

c               do2=.false.
c               do k5=1,nmollst(moli1)
c                 if(mollst(k5,moli1).eq.moli3) then
c                   do2=.true.
c                 goto 31
c                 end if
c               end do

c   31 continue
c              if(do2.eq..true.) then
              if(r2.le.6.0d0) then
               npole3b=6
               pnum(1)=imol(1,moli1)
               pnum(2)=imol(1,moli1)+1
               pnum(3)=imol(2,moli1)
               pnum(4)=imol(1,moli3)
               pnum(5)=imol(1,moli3)+1
               pnum(6)=imol(2,moli3)
               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)
               ep3bt = ep3bt - eptemp
               do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
                end do
               end do
               do i=1,3
                do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
                end do
               end do
              end if

            end if
          end do
        end do 
      end do  
cc!$OMP END DO
cc!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      return
      end
