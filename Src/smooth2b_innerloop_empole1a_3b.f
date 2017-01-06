c
      subroutine Smooth2bInnerloop1_2cut_wskin(
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
      real*8 tapr2b,tapr3b,dtapr2b_ix,dtapr2b_iy,dtapr2b_iz
      real*8 dtapr2b_jx,dtapr2b_jy,dtapr2b_jz
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 rtapr2b,rtapr3b
      real*8 tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      

      rtapr2b=5.0d0
      rtapr3b=6.0d0

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
                      virep3bt(j,i)=0.0d0
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
           vir2moli12(j,i)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 

c      print*,"Begin inner1"

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start,last,nmollst2,mollst2,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
!$OMP& schedule(guided)
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
            call image(xr,yr,zr)
            xr1=xr
            yr1=yr
            zr1=zr

            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r1=sqrt(r1_2)
          if(r1.le.6.0d0) then
           call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
           if(r1.gt.rtapr2b) then
             tapr2b_2=(r1-rtapr2b)*(r1-rtapr2b)
             tapr2b_3=tapr2b_2*(r1-rtapr2b)
             tapr2b_4=tapr2b_3*(r1-rtapr2b)
             tapr2b_5=tapr2b_4*(r1-rtapr2b)

             tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5

             dtapr2b_ix = (-6.0d0*5.0d0*tapr2b_4 + 
     &                    15.0d0*4.0d0*tapr2b_3 
     &                    -10.0d0*3.0d0*tapr2b_2)*xr1/r1
             dtapr2b_iy = (-6.0d0*5.0d0*tapr2b_4 +      
     &                    15.0d0*4.0d0*tapr2b_3       
     &                    -10.0d0*3.0d0*tapr2b_2)*yr1/r1
             dtapr2b_iz = (-6.0d0*5.0d0*tapr2b_4 +      
     &                    15.0d0*4.0d0*tapr2b_3       
     &                    -10.0d0*3.0d0*tapr2b_2)*zr1/r1
c             dtapr2b_jx = (-6.0d0*-5.0d0*tapr2b_4 +      
c     &                    15.0d0*-4.0d0*tapr2b_3       
c     &                    -10.0d0*-3.0d0*tapr2b_2)*xr/r1
c             dtapr2b_jy = (-6.0d0*-5.0d0*tapr2b_4 +
c     &                    15.0d0*-4.0d0*tapr2b_3
c     &                    -10.0d0*-3.0d0*tapr2b_2)*yr/r1
c             dtapr2b_jz = (-6.0d0*-5.0d0*tapr2b_4 +
c     &                    15.0d0*-4.0d0*tapr2b_3
c     &                    -10.0d0*-3.0d0*tapr2b_2)*zr/r1
             dtapr2b_jx = -dtapr2b_ix
             dtapr2b_jy = -dtapr2b_iy
             dtapr2b_jz = -dtapr2b_iz
  
             ep3bt = ep3bt + eptemp*tapr2b
c             print*,"Taper function,r1",tapr2b,r1,moli1,moli2
c             print*,"Tapered 2body eng,r1",eptemp*tapr2b,r1,moli1,moli2
c             print*,"UNtapered 2body eng,r1",eptemp,r1,moli1,moli2
             ep2moli12=eptemp
             do l1 = 1,np1
               i = pnum(l1)
               dep3bt(1,i) = dep3bt(1,i) + deptemp(1,l1)*tapr2b 
     &                       + dtapr2b_ix*eptemp
               dep3bt(2,i) = dep3bt(2,i) + deptemp(2,l1)*tapr2b
     &                       + dtapr2b_iy*eptemp 
               dep3bt(3,i) = dep3bt(3,i) + deptemp(3,l1)*tapr2b
     &                       + dtapr2b_iz*eptemp
               dep2moli12(1,l1)=deptemp(1,l1)
               dep2moli12(2,l1)=deptemp(2,l1)
               dep2moli12(3,l1)=deptemp(3,l1)
             end do 
             do l1 = np1+1, np2
               i = pnum(l1)
               dep3bt(1,i) = dep3bt(1,i) + deptemp(1,l1)*tapr2b  
     &                       + dtapr2b_jx*eptemp
               dep3bt(2,i) = dep3bt(2,i) + deptemp(2,l1)*tapr2b    
     &                       + dtapr2b_jy*eptemp
               dep3bt(3,i) = dep3bt(3,i) + deptemp(3,l1)*tapr2b
     &                       + dtapr2b_jz*eptemp
               dep2moli12(1,l1)=deptemp(1,l1)
               dep2moli12(2,l1)=deptemp(2,l1)
               dep2moli12(3,l1)=deptemp(3,l1)
             end do             
! NEED TO SMOOTHEN VIRIAL !
             do i=1,3
               do j=1,3
                  vir2moli12(j,i)=virtemp(j,i)
                  virtemp(j,i)=virtemp(j,i)*tapr2b
               end do 
             end do  
             do l1=1,np1
                i = pnum(l1)                
               virtemp(1,1)=virtemp(1,1)+dtapr2b_ix*eptemp*x(i)
               virtemp(2,1)=virtemp(2,1)+dtapr2b_ix*eptemp*y(i)
               virtemp(3,1)=virtemp(3,1)+dtapr2b_ix*eptemp*z(i)

               virtemp(1,2)=virtemp(1,2)+dtapr2b_iy*eptemp*x(i)
               virtemp(2,2)=virtemp(2,2)+dtapr2b_iy*eptemp*y(i)
               virtemp(3,2)=virtemp(3,2)+dtapr2b_iy*eptemp*z(i)

               virtemp(1,3)=virtemp(1,3)+dtapr2b_iz*eptemp*x(i)
               virtemp(2,3)=virtemp(2,3)+dtapr2b_iz*eptemp*y(i)
               virtemp(3,3)=virtemp(3,3)+dtapr2b_iz*eptemp*z(i)
             end do
             do l1=np1+1,npole3b
                i = pnum(l1)
               virtemp(1,1)=virtemp(1,1)+dtapr2b_jx*eptemp*x(i)
               virtemp(2,1)=virtemp(2,1)+dtapr2b_jx*eptemp*y(i)
               virtemp(3,1)=virtemp(3,1)+dtapr2b_jx*eptemp*z(i)

               virtemp(1,2)=virtemp(1,2)+dtapr2b_jy*eptemp*x(i)
               virtemp(2,2)=virtemp(2,2)+dtapr2b_jy*eptemp*y(i)
               virtemp(3,2)=virtemp(3,2)+dtapr2b_jy*eptemp*z(i)         

               virtemp(1,3)=virtemp(1,3)+dtapr2b_jz*eptemp*x(i)
               virtemp(2,3)=virtemp(2,3)+dtapr2b_jz*eptemp*y(i)
               virtemp(3,3)=virtemp(3,3)+dtapr2b_jz*eptemp*z(i)
             end do 
             do i=1,3
               do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
               end do
             end do
           else  
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
                 virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
                 vir2moli12(j,i)=virtemp(j,i)
                end do
              end do
           end if
          end if

        end do 
      end do  
!$OMP END DO
!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      return
      end
c
      subroutine Smooth2bInnerloop1_single_2cut_wskin(
     &  moli1,ep3bt,virep3bt,dep3bt,ntript) 
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
      real*8 tapr2b,tapr3b,dtapr2b_ix,dtapr2b_iy,dtapr2b_iz
      real*8 dtapr2b_jx,dtapr2b_jy,dtapr2b_jz
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 rtapr2b,rtapr3b
      real*8 tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5
      
      rtapr2b=5.0d0
      rtapr3b=6.0d0

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
                      virep3bt(j,i)=0.0d0
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
           vir2moli12(j,i)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 

c      print*,"Begin inner1"
      if(nmollst(moli1).gt.0) then
!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,moli1,nmollst2,mollst2,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
!$OMP& schedule(guided)
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
            call image(xr,yr,zr)
            xr1=xr
            yr1=yr
            zr1=zr

            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r1=sqrt(r1_2)
          if(r1.le.6.0d0) then
           call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
           if(r1.gt.rtapr2b) then
             tapr2b_2=(r1-rtapr2b)*(r1-rtapr2b)
             tapr2b_3=tapr2b_2*(r1-rtapr2b)
             tapr2b_4=tapr2b_3*(r1-rtapr2b)
             tapr2b_5=tapr2b_4*(r1-rtapr2b)

             tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5

             dtapr2b_ix = (-6.0d0*5.0d0*tapr2b_4 + 
     &                    15.0d0*4.0d0*tapr2b_3 
     &                    -10.0d0*3.0d0*tapr2b_2)*xr/r1
             dtapr2b_iy = (-6.0d0*5.0d0*tapr2b_4 +      
     &                    15.0d0*4.0d0*tapr2b_3       
     &                    -10.0d0*3.0d0*tapr2b_2)*yr/r1
             dtapr2b_iz = (-6.0d0*5.0d0*tapr2b_4 +      
     &                    15.0d0*4.0d0*tapr2b_3       
     &                    -10.0d0*3.0d0*tapr2b_2)*zr/r1
c             dtapr2b_jx = (-6.0d0*-5.0d0*tapr2b_4 +      
c     &                    15.0d0*-4.0d0*tapr2b_3       
c     &                    -10.0d0*-3.0d0*tapr2b_2)*xr/r1
c             dtapr2b_jy = (-6.0d0*-5.0d0*tapr2b_4 +
c     &                    15.0d0*-4.0d0*tapr2b_3
c     &                    -10.0d0*-3.0d0*tapr2b_2)*yr/r1
c             dtapr2b_jz = (-6.0d0*-5.0d0*tapr2b_4 +
c     &                    15.0d0*-4.0d0*tapr2b_3
c     &                    -10.0d0*-3.0d0*tapr2b_2)*zr/r1
             dtapr2b_jx = -dtapr2b_ix
             dtapr2b_jy = -dtapr2b_iy
             dtapr2b_jz = -dtapr2b_iz
  
             ep3bt = ep3bt + eptemp*tapr2b
             ep2moli12=eptemp
             do l1 = 1,np1
               i = pnum(l1)
               dep3bt(1,i) = dep3bt(1,i) + deptemp(1,l1)*tapr2b
     &                       + dtapr2b_ix*eptemp
               dep3bt(2,i) = dep3bt(2,i) + deptemp(2,l1)*tapr2b
     &                       + dtapr2b_iy*eptemp
               dep3bt(3,i) = dep3bt(3,i) + deptemp(3,l1)*tapr2b
     &                       + dtapr2b_iz*eptemp
               dep2moli12(1,l1)=deptemp(1,l1)
               dep2moli12(2,l1)=deptemp(2,l1)
               dep2moli12(3,l1)=deptemp(3,l1)
             end do
             do l1 = np1+1, np2
               i = pnum(l1)
               dep3bt(1,i) = dep3bt(1,i) + deptemp(1,l1)*tapr2b
     &                       + dtapr2b_jx*eptemp
               dep3bt(2,i) = dep3bt(2,i) + deptemp(2,l1)*tapr2b
     &                       + dtapr2b_jy*eptemp
               dep3bt(3,i) = dep3bt(3,i) + deptemp(3,l1)*tapr2b
     &                       + dtapr2b_jz*eptemp
               dep2moli12(1,l1)=deptemp(1,l1)
               dep2moli12(2,l1)=deptemp(2,l1)
               dep2moli12(3,l1)=deptemp(3,l1)
             end do
! NEED TO SMOOTHEN VIRIAL !
             do i=1,3
               do j=1,3
                  vir2moli12(j,i)=virtemp(j,i)
                  virtemp(j,i)=virtemp(j,i)*tapr2b
               end do
             end do

             do l1=1,np1
                i = pnum(l1)                
               virtemp(1,1)=virtemp(1,1)+dtapr2b_ix*eptemp*x(i)
               virtemp(2,1)=virtemp(2,1)+dtapr2b_ix*eptemp*y(i)
               virtemp(3,1)=virtemp(3,1)+dtapr2b_ix*eptemp*z(i)

               virtemp(1,2)=virtemp(1,2)+dtapr2b_iy*eptemp*x(i)
               virtemp(2,2)=virtemp(2,2)+dtapr2b_iy*eptemp*y(i)
               virtemp(3,2)=virtemp(3,2)+dtapr2b_iy*eptemp*z(i)

               virtemp(1,3)=virtemp(1,3)+dtapr2b_iz*eptemp*x(i)
               virtemp(2,3)=virtemp(2,3)+dtapr2b_iz*eptemp*y(i)
               virtemp(3,3)=virtemp(3,3)+dtapr2b_iz*eptemp*z(i)
             end do
             do l1=np1+1,npole3b
                i = pnum(l1)
               virtemp(1,1)=virtemp(1,1)+dtapr2b_jx*eptemp*x(i)
               virtemp(2,1)=virtemp(2,1)+dtapr2b_jx*eptemp*y(i)
               virtemp(3,1)=virtemp(3,1)+dtapr2b_jx*eptemp*z(i)

               virtemp(1,2)=virtemp(1,2)+dtapr2b_jy*eptemp*x(i)
               virtemp(2,2)=virtemp(2,2)+dtapr2b_jy*eptemp*y(i)
               virtemp(3,2)=virtemp(3,2)+dtapr2b_jy*eptemp*z(i)         

               virtemp(1,3)=virtemp(1,3)+dtapr2b_jz*eptemp*x(i)
               virtemp(2,3)=virtemp(2,3)+dtapr2b_jz*eptemp*y(i)
               virtemp(3,3)=virtemp(3,3)+dtapr2b_jz*eptemp*z(i)
             end do 
             do i=1,3
               do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
               end do
             end do
           else  
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
                 virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
                 vir2moli12(j,i)=virtemp(j,i)
                end do
              end do
           end if
          end if

        end do 
c      end do  
!$OMP END DO
!$OMP END PARALLEL
      end if

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      return
      end
