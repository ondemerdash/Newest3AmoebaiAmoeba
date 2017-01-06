c
c      subroutine R123Smooth3b2bSmoothOnly2(
c     &  start,last,ep3bt,virep3bt,dep3bt,ntript) 
      subroutine COMSmooth2bInnerloop1_2cut_wskin( 
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
      real*8 tapr2b,tapr3b
      real*8 dtapr2b_x(30),dtapr2b_y(30),dtapr2b_z(30)
      real*8 dtapr3b_x(30),dtapr3b_y(30),dtapr3b_z(30)
      real*8 tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 xcm,ycm,zcm,totmas,totmas1,totmas2,totmas3
      real*8 r2b_12,r2b_13,r2b_23,r2b_12_2,r2b_13_2,r2b_23_2
      real*8 SwitchR3b,SwitchR2b
      real*8 ep3moli123,ep2moli13,ep2moli23
      real*8 untapep3moli123,untapep2moli13,untapep2moli23
      real*8 untapep2moli12
      real*8 dep3moli123(3,30),vir3moli123(3,3)
      integer l2
c      rtapr3b=10.0d0
c      rtapr3b=13.0d0
c      R123cut3b=15.0d0
c      SwitchR3b=R123cut3b-rtapr3b

c      rtapr2b_input=5.5d0
c      rtapr2b_input=4.0d0
c      Rcut2b_input=7.0d0
      Rcut2b_input=1.0d12
      rtapr2b_input=1.0d12
      SwitchR2b=Rcut2b_input-rtapr2b_input
      !print*,"rtapr2b_input",rtapr2b_input
      
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
          dep3moli123(j,i)=0.0d0
c          dep1moli2(j,i)=0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
c           vir1moli1(i,j)=0.0d0
           vir2moli12(j,i)=0.0d0
           vir3moli123(j,i)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 
      ep3moli123=0.0d0

c      print*,"Begin L3ong New Taper"

!$OMP PARALLEL default(private) shared(imol,ep3bt,Rcut2b_input,
!$OMP& SwitchR2b,virep3bt,dep3bt,start,last,
!$OMP& nmollst,mollst,mass,x,y,z,rtapr2b_input)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
!$OMP& schedule(guided)
      do moli1=start,last
        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
          np1=3
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
            x1=0.0d0
            y1=0.0d0
            z1=0.0d0
            totmas1=0.0d0            
            xcm=0.0d0
            ycm=0.0d0
            zcm=0.0d0
            totmas=0.0d0
            do l1 = 1, np1
              i = pnum(l1)
                x1 = x1+mass(i)*x(i)
                y1 = y1+mass(i)*y(i)
                z1 = z1+mass(i)*z(i)
                xcm=xcm+mass(i)*x(i)
                ycm=ycm+mass(i)*y(i)
                zcm=zcm+mass(i)*z(i)
                totmas1=totmas1+mass(i)
                totmas=totmas+mass(i)
            end do
            x1=x1/totmas1
            y1=y1/totmas1
            z1=z1/totmas1
            x2=0.0d0
            y2=0.0d0
            z2=0.0d0
            totmas2=0.0d0
            do l1 = np1+1, npole3b
              i = pnum(l1)
                x2 = x2+mass(i)*x(i)
                y2 = y2+mass(i)*y(i)
                z2 = z2+mass(i)*z(i)
                xcm=xcm+mass(i)*x(i)
                ycm=ycm+mass(i)*y(i)
                zcm=zcm+mass(i)*z(i)
                totmas2=totmas2+mass(i)
                totmas=totmas+mass(i)
            end do
            x2=x2/totmas2
            y2=y2/totmas2
            z2=z2/totmas2
            xr1 = x1 - x2
            yr1 = y1 - y2
            zr1 = z1 - z2
            r2b_12_2=xr1*xr1+yr1*yr1+zr1*zr1
            r2b_12=sqrt(r2b_12_2)

          if(r2b_12.le.Rcut2b_input) then

               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)

              ep2moli12=eptemp
              
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,l1)=deptemp(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                end do
              end do


            if (r2b_12.gt.rtapr2b_input) then
               tapr2b_2=(r2b_12-rtapr2b_input)/SwitchR2b
     &          *(r2b_12-rtapr2b_input)/SwitchR2b
               tapr2b_3=tapr2b_2*(r2b_12-rtapr2b_input)/SwitchR2b
               tapr2b_4=tapr2b_3*(r2b_12-rtapr2b_input)/SwitchR2b
               tapr2b_5=tapr2b_4*(r2b_12-rtapr2b_input)/SwitchR2b

               tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5

                 do l1=1,np1
                    i=pnum(l1)
                 dtapr2b_x(l1) = (-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                  15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &               -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(x1-x2)/r2b_12*
     &                    mass(i)/totmas1
                 dtapr2b_y(l1) = (-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                    15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &            -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(y1-y2)/r2b_12*
     &                    mass(i)/totmas1
                 dtapr2b_z(l1) = (-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(z1-z2)/r2b_12*
     &                    mass(i)/totmas1
                 end do
                 do l1=np1+1,npole3b
                    i=pnum(l1)
                 dtapr2b_x(l1) = -(-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                    15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(x1-x2)/r2b_12*
     &                    mass(i)/totmas2
                 dtapr2b_y(l1) = -(-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(y1-y2)/r2b_12*
     &                    mass(i)/totmas2
                 dtapr2b_z(l1) = -(-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                 15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(z1-z2)/r2b_12*
     &                    mass(i)/totmas2
                 end do

                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep2moli12(1,l1)=dep2moli12(1,l1)*tapr2b
     &                            + ep2moli12*dtapr2b_x(l1)
                    dep2moli12(2,l1)=dep2moli12(2,l1)*tapr2b
     &                            + ep2moli12*dtapr2b_y(l1)
                    dep2moli12(3,l1)=dep2moli12(3,l1)*tapr2b
     &                            + ep2moli12*dtapr2b_z(l1)
                  end do
                  do i=1,3
                    do j=1,3
                      vir2moli12(j,i)=vir2moli12(j,i)*tapr2b
                    end do
                  end do
                  do l1=1,npole3b
                     i = pnum(l1)
                    vir2moli12(1,1)=vir2moli12(1,1)
     &                              +dtapr2b_x(l1)*ep2moli12*x(i)
                    vir2moli12(2,1)=vir2moli12(2,1)
     &                              +dtapr2b_x(l1)*ep2moli12*y(i)
                    vir2moli12(3,1)=vir2moli12(3,1)
     &                              +dtapr2b_x(l1)*ep2moli12*z(i)
                    vir2moli12(1,2)=vir2moli12(1,2)
     &                              +dtapr2b_y(l1)*ep2moli12*x(i)
                    vir2moli12(2,2)=vir2moli12(2,2)
     &                              +dtapr2b_y(l1)*ep2moli12*y(i)
                    vir2moli12(3,2)=vir2moli12(3,2)
     &                              +dtapr2b_y(l1)*ep2moli12*z(i)
                    vir2moli12(1,3)=vir2moli12(1,3)
     &                              +dtapr2b_z(l1)*ep2moli12*x(i)
                    vir2moli12(2,3)=vir2moli12(2,3)
     &                              +dtapr2b_z(l1)*ep2moli12*y(i)
                    vir2moli12(3,3)=vir2moli12(3,3)
     &                              +dtapr2b_z(l1)*ep2moli12*z(i)
                  end do
                  ep2moli12=tapr2b*ep2moli12
            end if            
              ep3bt=ep3bt+ep2moli12
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do
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
c
c
c
c
c      subroutine R123Smooth3b2bSmoothOnly2(
c     &  start,last,ep3bt,virep3bt,dep3bt,ntript) 
      subroutine COMSmooth2bInnerloop1_2cut_wskin_single(
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
      real*8 tapr2b,tapr3b
      real*8 dtapr2b_x(30),dtapr2b_y(30),dtapr2b_z(30)
      real*8 dtapr3b_x(30),dtapr3b_y(30),dtapr3b_z(30)
      real*8 tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 xcm,ycm,zcm,totmas,totmas1,totmas2,totmas3
      real*8 r2b_12,r2b_13,r2b_23,r2b_12_2,r2b_13_2,r2b_23_2
      real*8 SwitchR3b,SwitchR2b
      real*8 ep3moli123,ep2moli13,ep2moli23
      real*8 untapep3moli123,untapep2moli13,untapep2moli23
      real*8 untapep2moli12
      real*8 dep3moli123(3,30),vir3moli123(3,3)
      integer l2
c      rtapr3b=10.0d0

c      rtapr2b_input=4.0d0
c      Rcut2b_input=7.0d0
      Rcut2b_input=1.0d12
      rtapr2b_input=1.0d12

      SwitchR2b=Rcut2b_input-rtapr2b_input
      !print*,"rtapr2b_input",rtapr2b_input

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
          dep3moli123(j,i)=0.0d0
c          dep1moli2(j,i)=0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
c           vir1moli1(i,j)=0.0d0
           vir2moli12(j,i)=0.0d0
           vir3moli123(j,i)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 
      ep3moli123=0.0d0

      print*,"Begin L3ong New Taper"
      if(nmollst(moli1).gt.0) then
!$OMP PARALLEL default(private) shared(imol,ep3bt,Rcut2b_input,
!$OMP& SwitchR2b,virep3bt,dep3bt,moli1,
!$OMP& nmollst,mollst,mass,x,y,z,rtapr2b_input)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
!$OMP& schedule(guided)
c      do moli1=start,last
        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
          np1=3
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
            x1=0.0d0
            y1=0.0d0
            z1=0.0d0
            totmas1=0.0d0            
            xcm=0.0d0
            ycm=0.0d0
            zcm=0.0d0
            totmas=0.0d0
            do l1 = 1, np1
              i = pnum(l1)
                x1 = x1+mass(i)*x(i)
                y1 = y1+mass(i)*y(i)
                z1 = z1+mass(i)*z(i)
                xcm=xcm+mass(i)*x(i)
                ycm=ycm+mass(i)*y(i)
                zcm=zcm+mass(i)*z(i)
                totmas1=totmas1+mass(i)
                totmas=totmas+mass(i)
            end do
            x1=x1/totmas1
            y1=y1/totmas1
            z1=z1/totmas1
            x2=0.0d0
            y2=0.0d0
            z2=0.0d0
            totmas2=0.0d0
            do l1 = np1+1, npole3b
              i = pnum(l1)
                x2 = x2+mass(i)*x(i)
                y2 = y2+mass(i)*y(i)
                z2 = z2+mass(i)*z(i)
                xcm=xcm+mass(i)*x(i)
                ycm=ycm+mass(i)*y(i)
                zcm=zcm+mass(i)*z(i)
                totmas2=totmas2+mass(i)
                totmas=totmas+mass(i)
            end do
            x2=x2/totmas2
            y2=y2/totmas2
            z2=z2/totmas2
            xr1 = x1 - x2
            yr1 = y1 - y2
            zr1 = z1 - z2
            r2b_12_2=xr1*xr1+yr1*yr1+zr1*zr1
            r2b_12=sqrt(r2b_12_2)

          if(r2b_12.le.Rcut2b_input) then

               call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &         virtemp)

              ep2moli12=eptemp
              untapep2moli12=eptemp
              
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,l1)=deptemp(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                end do
              end do


            if (r2b_12.gt.rtapr2b_input) then
               tapr2b_2=(r2b_12-rtapr2b_input)/SwitchR2b
     &          *(r2b_12-rtapr2b_input)/SwitchR2b
               tapr2b_3=tapr2b_2*(r2b_12-rtapr2b_input)/SwitchR2b
               tapr2b_4=tapr2b_3*(r2b_12-rtapr2b_input)/SwitchR2b
               tapr2b_5=tapr2b_4*(r2b_12-rtapr2b_input)/SwitchR2b

               tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5

                 do l1=1,np1
                    i=pnum(l1)
                 dtapr2b_x(l1) = (-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                  15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &               -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(x1-x2)/r2b_12*
     &                    mass(i)/totmas1
                 dtapr2b_y(l1) = (-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                    15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &            -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(y1-y2)/r2b_12*
     &                    mass(i)/totmas1
                 dtapr2b_z(l1) = (-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(z1-z2)/r2b_12*
     &                    mass(i)/totmas1
                 end do
                 do l1=np1+1,npole3b
                    i=pnum(l1)
                 dtapr2b_x(l1) = -(-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                    15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(x1-x2)/r2b_12*
     &                    mass(i)/totmas2
                 dtapr2b_y(l1) = -(-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(y1-y2)/r2b_12*
     &                    mass(i)/totmas2
                 dtapr2b_z(l1) = -(-6.0d0*5.0d0*tapr2b_4/SwitchR2b +
     &                 15.0d0*4.0d0*tapr2b_3/SwitchR2b
     &                -10.0d0*3.0d0*tapr2b_2/SwitchR2b)*(z1-z2)/r2b_12*
     &                    mass(i)/totmas2
                 end do

                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep2moli12(1,l1)=dep2moli12(1,l1)*tapr2b
     &                            + ep2moli12*dtapr2b_x(l1)
                    dep2moli12(2,l1)=dep2moli12(2,l1)*tapr2b
     &                            + ep2moli12*dtapr2b_y(l1)
                    dep2moli12(3,l1)=dep2moli12(3,l1)*tapr2b
     &                            + ep2moli12*dtapr2b_z(l1)
                  end do
                  do i=1,3
                    do j=1,3
                      vir2moli12(j,i)=vir2moli12(j,i)*tapr2b
                    end do
                  end do
                  do l1=1,npole3b
                     i = pnum(l1)
                    vir2moli12(1,1)=vir2moli12(1,1)
     &                              +dtapr2b_x(l1)*ep2moli12*x(i)
                    vir2moli12(2,1)=vir2moli12(2,1)
     &                              +dtapr2b_x(l1)*ep2moli12*y(i)
                    vir2moli12(3,1)=vir2moli12(3,1)
     &                              +dtapr2b_x(l1)*ep2moli12*z(i)
                    vir2moli12(1,2)=vir2moli12(1,2)
     &                              +dtapr2b_y(l1)*ep2moli12*x(i)
                    vir2moli12(2,2)=vir2moli12(2,2)
     &                              +dtapr2b_y(l1)*ep2moli12*y(i)
                    vir2moli12(3,2)=vir2moli12(3,2)
     &                              +dtapr2b_y(l1)*ep2moli12*z(i)
                    vir2moli12(1,3)=vir2moli12(1,3)
     &                              +dtapr2b_z(l1)*ep2moli12*x(i)
                    vir2moli12(2,3)=vir2moli12(2,3)
     &                              +dtapr2b_z(l1)*ep2moli12*y(i)
                    vir2moli12(3,3)=vir2moli12(3,3)
     &                              +dtapr2b_z(l1)*ep2moli12*z(i)
                  end do
                  ep2moli12=tapr2b*ep2moli12
            end if            
              ep3bt=ep3bt+ep2moli12
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do
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
c
c
c
