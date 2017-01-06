c
c      subroutine totfieldNpolesmoothInner_1c_3bPolar(
c     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
      subroutine totfieldNpolesmoothInner_1c_3bPolar(
     &  ep3bt,virep3bt,dep3bt)
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use mpidat
      implicit none
      real*8 ep2moli12
      real*8 dep2moli12(3,npole)
      real*8 ep1moli1,ep1moli2
      real*8 dep1moli1(3,npole),dep1moli2(3,npole)      
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp,deptemp(3,npole)
      real*8 vir2moli12(3,3),vir1moli1(3,3),vir1moli2(3,3)
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
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 rtapr2b
      real*8 tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,npole) 
      real*8 Cobarcut3b,Rcut2b,SwitchR2b,SwitchR3b,Rcut2b2
      integer moli1rmndr,l2
      real*8 term,term2,term3,term4
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      character*6 mode

c      rtapr2b=1.0d12
c      Rcut2b=1.0d12
c      rtapr2b=6.0d0
c      Rcut2b=8.0d0
      rtapr2b=rtapr2b_input
      Rcut2b=cut2b_input

      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b

c      rtapr3b=6.0d0
c      Cobarcut3b=8.0d0
      SwitchR3b=Cobarcut3b-rtapr3b


                ep3bt=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(j,i)=0.0d0
                   end do
                end do

      mode = 'EWALD'
      call switch (mode)

        print*,"TOTFIELDEWALD NOCut3b NoCut2b"
       ! print*,"task start last",taskid,start_polar(taskid),
     & !   last_polar(taskid)
        print*,"rtapr2b",rtapr2b
        print*,"rtapr3b",rtapr3b

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,npole,
c!$OMP& nmol)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
c!$OMP& schedule(guided)
c      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
!$OMP& Rcut2b2,npole)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
!$OMP& schedule(guided)
c      do kouter=start_polar(taskid),last_polar(taskid)
c        moli1=molnew(kouter)
      do moli1=start_polar(taskid),last_polar(taskid)

          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=3
           call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &          deptemp,virtemp)
          ep1moli1=eptemp
c          print*,"moli1 ep1moli1=",moli1,ep1moli1
          ep3bt=ep3bt+ep1moli1
          !do i = 1, npole
          !  !i = pnum(l1)
          do i=1,npole
            do j = 1, 3
              dep1moli1(j,i)=deptemp(j,i)
c              print*,"i j dep1moli1(j,i)",i,j,dep1moli1(j,i)
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,i)
            end do
          end do
              do i=1,3
                do j=1,3
                 vir1moli1(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)+vir1moli1(j,i)
                end do
              end do


        !do moli2=moli1+1,nmol
        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
           npole3b=3
           pnum(1)=imol(1,moli2)
           pnum(2)=imol(1,moli2)+1
           pnum(3)=imol(2,moli2)
           call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &             deptemp,virtemp)
           ep1moli2=eptemp
           do i = 1, npole
              do j = 1, 3
               dep1moli2(j,i)=deptemp(j,i)
              end do
           end do
           do i=1,3
              do j=1,3
                vir1moli2(j,i)=virtemp(j,i)
              end do
           end do

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
            if (name(i).eq.'O') then
              x1 = x(i)
              y1 = y(i)
              z1 = z(i)
            end if
          end do
          do l1 = np1+1, np2
            i = pnum(l1)
            if (name(i).eq.'O') then
              x2 = x(i)
              y2 = y(i)
              z2 = z(i)
            end if
          end do
          xr1 = x1 - x2
          yr1 = y1 - y2
          zr1 = z1 - z2
          call image(xr1,yr1,zr1)
          r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1


          if(r1_2.le.Rcut2b2) then
             r1=sqrt(r1_2)
 
             call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &             deptemp,virtemp)
              ep2moli12=eptemp
              do i=1,npole
                do j = 1, 3
                  dep2moli12(j,i)=deptemp(j,i)
                end do
              end do
              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                end do
              end do

              if (r1.gt.rtapr2b) then
                 term=(r1-rtapr2b)/SwitchR2b
                 tapr2b_2=term*term
                 tapr2b_3=tapr2b_2*term
                 tapr2b_4=tapr2b_3*term
                 tapr2b_5=tapr2b_4*term
                 tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5
                 npole3b=6
                 np1=3
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli2)
                 pnum(5)=imol(1,moli2)+1
                 pnum(6)=imol(2,moli2)
                 do l1=1,np1
                   i = pnum(l1)
                   if (name(i).eq.'O') then
                     term2=(-30.0d0*tapr2b_4+60.0d0*tapr2b_3
     &                 -30.0d0*tapr2b_2)/r1/SwitchR2b
                     dtapr2b_x(l1) = xr1*term2
                     dtapr2b_y(l1) = yr1*term2
                     dtapr2b_z(l1) = zr1*term2
                   else
                     dtapr2b_x(l1) =0.0d0
                     dtapr2b_y(l1) =0.0d0
                     dtapr2b_z(l1) =0.0d0
                   end if
                 end do
                 do l1=np1+1,npole3b
                   i = pnum(l1)
                   if (name(i).eq.'O') then
                     term2=(-30.0d0*tapr2b_4+60.0d0*tapr2b_3
     &                -30.0d0*tapr2b_2)/r1/SwitchR2b
                     dtapr2b_x(l1) = -xr1*term2
                     dtapr2b_y(l1) = -yr1*term2
                     dtapr2b_z(l1) = -zr1*term2
                   else
                     dtapr2b_x(l1) =0.0d0
                     dtapr2b_y(l1) =0.0d0
                     dtapr2b_z(l1) =0.0d0
                   end if
                 end do

                 ep3bt= ep3bt + tapr2b*(ep2moli12-ep1moli1-ep1moli2)
                 do i =1,npole
                  dep3bt(1,i)=dep3bt(1,i) + tapr2b*(dep2moli12(1,i)
     &                       -dep1moli1(1,i)-dep1moli2(1,i))
                  dep3bt(2,i)=dep3bt(2,i) + tapr2b*(dep2moli12(2,i)
     &                       -dep1moli1(2,i)-dep1moli2(2,i))
                  dep3bt(3,i)=dep3bt(3,i) + tapr2b*(dep2moli12(3,i)
     &                       -dep1moli1(3,i)-dep1moli2(3,i))
                 end do
                 do l1 = 1,npole3b
                  i = pnum(l1)
                  dep3bt(1,i)=dep3bt(1,i) 
     &                   + (ep2moli12-ep1moli1-ep1moli2)*dtapr2b_x(l1)
                  dep3bt(2,i)=dep3bt(2,i) 
     &                   + (ep2moli12-ep1moli1-ep1moli2)*dtapr2b_y(l1)
                  dep3bt(3,i)=dep3bt(3,i) 
     &                   + (ep2moli12-ep1moli1-ep1moli2)*dtapr2b_z(l1)
                 end do
                 
                 do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i) + tapr2b*(vir2moli12(j,i)
     &               -vir1moli1(j,i)-vir1moli2(j,i))
                  end do
                 end do

                do l1=1,np1
                  i=pnum(l1)
                  if(name(i).eq.'O') then

                   virep3bt(1,1)=virep3bt(1,1)
     &               +dtapr2b_x(l1)*(ep2moli12-ep1moli1-ep1moli2)*xr1
                   virep3bt(2,1)=virep3bt(2,1)
     &                  +dtapr2b_x(l1)*(ep2moli12-ep1moli1-ep1moli2)*yr1
                   virep3bt(3,1)=virep3bt(3,1)
     &                  +dtapr2b_x(l1)*(ep2moli12-ep1moli1-ep1moli2)*zr1
                   virep3bt(1,2)=virep3bt(2,1)
                   virep3bt(2,2)=virep3bt(2,2)
     &                  +dtapr2b_y(l1)*(ep2moli12-ep1moli1-ep1moli2)*yr1
                   virep3bt(3,2)=virep3bt(3,2)
     &                  +dtapr2b_y(l1)*(ep2moli12-ep1moli1-ep1moli2)*zr1
                   virep3bt(1,3)=virep3bt(3,1)
                   virep3bt(2,3)=virep3bt(3,2)
                   virep3bt(3,3)=virep3bt(3,3)
     &                  +dtapr2b_z(l1)*(ep2moli12-ep1moli1-ep1moli2)*zr1

                  end if
                end do
              else
                ep3bt=ep3bt+ep2moli12-ep1moli1-ep1moli2
                 do i =1,npole
                  dep3bt(1,i)=dep3bt(1,i) + (dep2moli12(1,i)
     &                       -dep1moli1(1,i)-dep1moli2(1,i))
                  dep3bt(2,i)=dep3bt(2,i) + (dep2moli12(2,i)
     &                       -dep1moli1(2,i)-dep1moli2(2,i))
                  dep3bt(3,i)=dep3bt(3,i) + (dep2moli12(3,i)
     &                       -dep1moli1(3,i)-dep1moli2(3,i))
                 end do
                 do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i) + (vir2moli12(j,i)
     &               -vir1moli1(j,i)-vir1moli2(j,i))
                  end do
                 end do

              end if
          else
              ep2moli12=0.0d0
              do i=1,npole
                do j = 1, 3
                  dep2moli12(j,i)=0.0d0
                end do
              end do
              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=0.0d0
                end do
              end do
          end if 
       

          !do moli3=moli2+1,nmol
          do k2=1,nmollst(moli2)
            moli3=mollst(k2,moli2)

            npole3b=9
            np2 = 6
            np3 =9
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
              if (name(i).eq.'O') then
                x3 = x(i)
                y3 = y(i)
                z3 = z(i)
              end if
            end do

            xr2 = x1 - x3
            yr2 = y1 - y3
            zr2 = z1 - z3
            call image(xr2,yr2,zr2)
            xr3=x2-x3
            yr3=y2-y3
            zr3=z2-z3
            call image(xr3,yr3,zr3)

            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
            r1=sqrt(r1_2)
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)

            if ((r1.lt.r3).and.(r2.lt.r3) ) then
              shellsum=r1+r2
              if (shellsum.le.Cobarcut3b) then            
               call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                 ep3moli123=eptemp
                 do i =1,npole
                   dep3moli123(1,i) = deptemp(1,i)
                   dep3moli123(2,i) = deptemp(2,i)
                   dep3moli123(3,i) = deptemp(3,i)
                 end do
                 do i=1,3
                   do j=1,3
                    vir3moli123(j,i)=virtemp(j,i)
                   end do
                 end do
            
                 if (shellsum.gt.rtapr3b) then
                  term3=(shellsum-rtapr3b)/SwitchR3b
                  tapr3b_2=term3*term3
                  tapr3b_3=tapr3b_2*term3
                  tapr3b_4=tapr3b_3*term3
                  tapr3b_5=tapr3b_4*term3
                  tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
     &                       - 6.0d0*tapr3b_5
                  term4 = (-30.0d0*tapr3b_4 + 60.0d0*tapr3b_3 -
     &                       30.0d0*tapr3b_2)/SwitchR3b
                  dtapr3b_ix = term4*(xr1/r1+xr2/r2)
                  dtapr3b_iy = term4*(yr1/r1+yr2/r2)
                  dtapr3b_iz = term4*(zr1/r1+zr2/r2)
                  dtapr3b_jx = term4*(-xr1/r1)
                  dtapr3b_jy = term4*(-yr1/r1)
                  dtapr3b_jz = term4*(-zr1/r1)
                  dtapr3b_kx = term4*(-xr2/r2)
                  dtapr3b_ky = term4*(-yr2/r2)
                  dtapr3b_kz = term4*(-zr2/r2)
                  virtapr3b_xx=term4*( xr1/r1*xr1 + xr2/r2*xr2)
                  virtapr3b_xy=term4*( xr1/r1*yr1 + xr2/r2*yr2)
                  virtapr3b_xz=term4*( xr1/r1*zr1 + xr2/r2*zr2)
                  virtapr3b_yy=term4*( yr1/r1*yr1 + yr2/r2*yr2)
                  virtapr3b_yz=term4*( yr1/r1*zr1 + yr2/r2*zr2)
                  virtapr3b_zz=term4*( zr1/r1*zr1 + zr2/r2*zr2)
                 end if
              end if
            else if ((r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
              if (shellsum.le.Cobarcut3b) then
               call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                 ep3moli123=eptemp
                 do i =1,npole
                   dep3moli123(1,i) = deptemp(1,i)
                   dep3moli123(2,i) = deptemp(2,i)
                   dep3moli123(3,i) = deptemp(3,i)
                 end do
                 do i=1,3
                   do j=1,3
                    vir3moli123(j,i)=virtemp(j,i)
                   end do
                 end do

                 if (shellsum.gt.rtapr3b) then
                   term3=(shellsum-rtapr3b)/SwitchR3b
                   tapr3b_2=term3*term3
                   tapr3b_3=tapr3b_2*term3
                   tapr3b_4=tapr3b_3*term3
                   tapr3b_5=tapr3b_4*term3
                   tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
     &                - 6.0d0*tapr3b_5
                   term4 = (-30.0d0*tapr3b_4 + 60.0d0*tapr3b_3 -
     &                     30.0d0*tapr3b_2)/SwitchR3b
                   dtapr3b_ix = term4*(xr1/r1)
                   dtapr3b_iy = term4*(yr1/r1)
                   dtapr3b_iz = term4*(zr1/r1)
                   dtapr3b_jx = term4*(-xr1/r1+xr3/r3)
                   dtapr3b_jy = term4*(-yr1/r1+yr3/r3)
                   dtapr3b_jz = term4*(-zr1/r1+zr3/r3)
                   dtapr3b_kx = term4*(-xr3/r3)
                   dtapr3b_ky = term4*(-yr3/r3)
                   dtapr3b_kz = term4*(-zr3/r3)
                   virtapr3b_xx=term4*((xr1/r1)*xr1+(xr3/r3)*xr3)
                   virtapr3b_xy=term4*((xr1/r1)*yr1+(xr3/r3)*yr3)
                   virtapr3b_xz=term4*((xr1/r1)*zr1+(xr3/r3)*zr3)
                   virtapr3b_yy=term4*((yr1/r1)*yr1+(yr3/r3)*yr3)
                   virtapr3b_yz=term4*((yr1/r1)*zr1+(yr3/r3)*zr3)
                   virtapr3b_zz=term4*((zr1/r1)*zr1+(zr3/r3)*zr3)
                 end if
              end if
            else if ((r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
              if (shellsum.le.Cobarcut3b) then
               call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                 ep3moli123=eptemp
                 do i =1,npole
                   dep3moli123(1,i) = deptemp(1,i)
                   dep3moli123(2,i) = deptemp(2,i)
                   dep3moli123(3,i) = deptemp(3,i)
                 end do
                 do i=1,3
                   do j=1,3
                    vir3moli123(j,i)=virtemp(j,i)
                   end do
                 end do

                 if (shellsum.gt.rtapr3b) then
                   term3=(shellsum-rtapr3b)/SwitchR3b
                   tapr3b_2=term3*term3
                   tapr3b_3=tapr3b_2*term3
                   tapr3b_4=tapr3b_3*term3
                   tapr3b_5=tapr3b_4*term3
                   tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
     &                - 6.0d0*tapr3b_5
                   term4 = (-30.0d0*tapr3b_4 + 60.0d0*tapr3b_3 -
     &                     30.0d0*tapr3b_2)/SwitchR3b
                   dtapr3b_ix = term4*(xr2/r2)
                   dtapr3b_iy = term4*(yr2/r2)
                   dtapr3b_iz = term4*(zr2/r2)
                   dtapr3b_jx = term4*(xr3/r3)
                   dtapr3b_jy = term4*(yr3/r3)
                   dtapr3b_jz = term4*(zr3/r3)
                   dtapr3b_kx = term4*(-xr2/r2-xr3/r3)
                   dtapr3b_ky = term4*(-yr2/r2-yr3/r3)
                   dtapr3b_kz = term4*(-zr2/r2-zr3/r3)
                   virtapr3b_xx=term4*((xr2/r2)*xr2+(xr3/r3)*xr3)
                   virtapr3b_xy=term4*((xr2/r2)*yr2+(xr3/r3)*yr3)
                   virtapr3b_xz=term4*((xr2/r2)*zr2+(xr3/r3)*zr3)
                   virtapr3b_yy=term4*((yr2/r2)*yr2+(yr3/r3)*yr3)
                   virtapr3b_yz=term4*((yr2/r2)*zr2+(yr3/r3)*zr3)
                   virtapr3b_zz=term4*((zr2/r2)*zr2+(zr3/r3)*zr3)
                 end if
              end if
            end if


            if (shellsum .le. Cobarcut3b) then
                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli2)
                 pnum(5)=imol(1,moli2)+1
                 pnum(6)=imol(2,moli2)

               if(r1_2.le.Rcut2b2) then
                 ep3moli123=ep3moli123-ep2moli12
                 do i = 1, npole
                   do j = 1, 3
                    dep3moli123(j,i)=dep3moli123(j,i)-dep2moli12(j,i)
                   end do
                 end do
                 do i=1,3
                   do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)-vir2moli12(j,i)
                   end do
                 end do
               else
               call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &            deptemp,virtemp)
                  ep3moli123=ep3moli123-eptemp
                 do i = 1, npole
                   do j = 1, 3
                    dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                   end do
                 end do
                 do i=1,3
                   do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                   end do
                 end do
               end if

                  npole3b=6
                  np1 =3
                  pnum(1)=imol(1,moli2)
                  pnum(2)=imol(1,moli2)+1
                  pnum(3)=imol(2,moli2)
                  pnum(4)=imol(1,moli3)
                  pnum(5)=imol(1,moli3)+1
                  pnum(6)=imol(2,moli3)              
               call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &            deptemp,virtemp)
                  ep3moli123=ep3moli123-eptemp
                  do i =1,npole
                    dep3moli123(1,i)= dep3moli123(1,i)-deptemp(1,i)
                    dep3moli123(2,i)= dep3moli123(2,i)-deptemp(2,i)
                    dep3moli123(3,i)= dep3moli123(3,i)-deptemp(3,i)
                  end do
                  do i=1,3
                    do j=1,3
                     vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                    end do
                  end do


                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
               call empole1a_3b_Polar_totfieldnpole(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                  ep3moli123=ep3moli123-eptemp
                 do i =1,npole
                    dep3moli123(1,i)= dep3moli123(1,i)-deptemp(1,i)
                    dep3moli123(2,i)= dep3moli123(2,i)-deptemp(2,i)
                    dep3moli123(3,i)= dep3moli123(3,i)-deptemp(3,i)
                 end do
                 do i=1,3
                   do j=1,3
                    vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                   end do
                 end do


                 npole3b=3
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 call empole1a_3b_Polar_totfieldnpole(
     &           npole3b,pnum,eptemp,deptemp,virtemp)
                 ep3moli123=ep3moli123+eptemp+ep1moli1+ep1moli2

                do i =1,npole
                  dep3moli123(1,i) = dep3moli123(1,i)+deptemp(1,i)
     &                     +dep1moli1(1,i)+dep1moli2(1,i)
                  dep3moli123(2,i) = dep3moli123(2,i)+deptemp(2,i)
     &                     +dep1moli1(2,i)+dep1moli2(2,i)
                  dep3moli123(3,i) = dep3moli123(3,i)+deptemp(3,i)
     &                     +dep1moli1(3,i)+dep1moli2(3,i)
                end do
                do i=1,3
                  do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)+virtemp(j,i)
     &               +vir1moli1(j,i)+vir1moli2(j,i)
                  end do
                end do

c               print*,"3body Cobar ep3body",shellsum,ep3moli123
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
               if(shellsum.gt.rtapr3b) then
                  ep3bt=ep3bt+ep3moli123*tapr3b
                  do i =1,npole
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,i)*tapr3b
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,i)*tapr3b
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,i)*tapr3b
                  end do

                  do l1 = 1,np1
                    i = pnum(l1)
                    if(name(i).eq.'O') then
                    dep3bt(1,i) = dep3bt(1,i)
     &                            + ep3moli123*dtapr3b_ix
                    dep3bt(2,i) = dep3bt(2,i)
     &                            + ep3moli123*dtapr3b_iy
                    dep3bt(3,i) = dep3bt(3,i)
     &                            + ep3moli123*dtapr3b_iz
                    end if
                  end do
                  do l1 = np1+1,np2
                    i = pnum(l1)
                    if(name(i).eq.'O') then
                    dep3bt(1,i) = dep3bt(1,i)
     &                            + ep3moli123*dtapr3b_jx
                    dep3bt(2,i) = dep3bt(2,i)
     &                            + ep3moli123*dtapr3b_jy
                    dep3bt(3,i) = dep3bt(3,i)
     &                            + ep3moli123*dtapr3b_jz
                    end if
                  end do
                  do l1 = np2+1,npole3b
                    i = pnum(l1)
                    if(name(i).eq.'O') then
                    dep3bt(1,i) = dep3bt(1,i)
     &                            + ep3moli123*dtapr3b_kx
                    dep3bt(2,i) = dep3bt(2,i)
     &                            + ep3moli123*dtapr3b_ky
                    dep3bt(3,i) = dep3bt(3,i)
     &                            + ep3moli123*dtapr3b_kz
                    end if
                  end do
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)*tapr3b
                    end do
                  end do
                  vir3moli123(1,1)=vir3moli123(1,1)
     &                               +virtapr3b_xx*ep3moli123
                  vir3moli123(2,1)=vir3moli123(2,1)
     &                               +virtapr3b_xy*ep3moli123
                  vir3moli123(3,1)=vir3moli123(3,1)
     &                              +virtapr3b_xz*ep3moli123
                  vir3moli123(1,2)=vir3moli123(1,2)
     &                               +virtapr3b_xy*ep3moli123
                  vir3moli123(2,2)=vir3moli123(2,2)
     &                                +virtapr3b_yy*ep3moli123
                  vir3moli123(3,2)=vir3moli123(3,2)
     &                                +virtapr3b_yz*ep3moli123
                  vir3moli123(1,3)=vir3moli123(1,3)
     &                                +virtapr3b_xz*ep3moli123
                  vir3moli123(2,3)=vir3moli123(2,3)
     &                                  +virtapr3b_yz*ep3moli123
                  vir3moli123(3,3)=vir3moli123(3,3)
     &                                  +virtapr3b_zz*ep3moli123
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
               else
                  ep3bt = ep3bt + ep3moli123
                  do i = 1, npole
                    do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,i)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
               end if
            end if
          end do
        end do 
      end do  
!$OMP END DO
!$OMP END PARALLEL
      
      return
      end

