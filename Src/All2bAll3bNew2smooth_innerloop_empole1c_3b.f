c
      subroutine All2bAll3bNew2SmoothInnerloop1_empole1c(
     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      implicit none
      real*8 ep2moli12
      real*8 dep2moli12(3,30)
      real*8 ep1moli1,ep1moli2
      real*8 dep1moli1(3,30),dep1moli2(3,30)      
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
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 rtapr2b
      real*8 tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,30) 
      real*8 Cobarcut3b,Rcut2b,SwitchR2b,SwitchR3b
      integer moli1rmndr,l2
      rtapr2b=1.0d12
      Rcut2b=1.0d12
      SwitchR2b=Rcut2b-rtapr2b

      rtapr3b=1.0d12
      Cobarcut3b=1.0d12
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

        print*,"EWALD NOCut3b NoCut2b"
        print*,"start last",start,last,moli1rmndr

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,
c!$OMP& nmol)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
c!$OMP& schedule(guided)
      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)

          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=3
           call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
          ep1moli1=eptemp
c          print*,"moli1 ep1moli1=",moli1,ep1moli1
          ep3bt=ep3bt+ep1moli1
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep1moli1(j,l1)=deptemp(j,l1)
c              print*,"i j dep1moli1(j,i)",i,j,dep1moli1(j,i)
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,l1)
            end do
          end do


        do moli2=moli1+1,nmol
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


           call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,l1)=deptemp(j,l1)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do

              npole3b=3
              pnum(1)=imol(1,moli2)
              pnum(2)=imol(1,moli2)+1
              pnum(3)=imol(2,moli2)
             call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
             ep1moli2=eptemp
             ep3bt=ep3bt-ep1moli2
c             do i=1,npole
              do l1 = 1, npole3b
                i = pnum(l1)
               do j = 1, 3
                 dep1moli2(j,l1)=deptemp(j,l1)
                dep3bt(j,i)=dep3bt(j,i)-dep1moli2(j,l1)
               end do
             end do

             pnum(1)=imol(1,moli1)
             pnum(2)=imol(1,moli1)+1
             pnum(3)=imol(2,moli1)
             ep3bt=ep3bt-ep1moli1
              do l1 = 1, npole3b
                i = pnum(l1)
c              do i=1,npole
                do j = 1, 3
                   dep3bt(j,i)=dep3bt(j,i)-dep1moli1(j,l1)
                end do
              end do

           ! print*,"2body r1 ep2body",r1,ep2moli12-ep1moli1-ep1moli2

          do moli3=moli2+1,nmol
            npole3b=9
            np3=9
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
            call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr

            xr = x2 - x3
            yr = y2 - y3
            zr = z2 - z3
            call image(xr,yr,zr)

            xr3=xr
            yr3=yr
            zr3=zr
c            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
c            r1=sqrt(r1_2)
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)

            if((r1.lt.r3).and.(r2.lt.r3) ) then
               shellsum=r1+r2
            else if ( (r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
            else if ( (r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
            end if

            
                 call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ep3bt=ep3bt+eptemp
                ep3moli123=eptemp

                do l1 = 1,npole3b
                  i = pnum(l1)
c                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,l1)
                end do
            


               npole3b=6
               ep3bt=ep3bt-ep2moli12
               ep3moli123=ep3moli123-ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3bt(j,i)=dep3bt(j,i)-dep2moli12(j,l1)
                  end do
               end do

                  npole3b=6
                  np1 =3
                  pnum(1)=imol(1,moli2)
                  pnum(2)=imol(1,moli2)+1
                  pnum(3)=imol(2,moli2)
                  pnum(4)=imol(1,moli3)
                  pnum(5)=imol(1,moli3)+1
                  pnum(6)=imol(2,moli3)              
                  call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
                  ep3bt=ep3bt-eptemp
                  ep3moli123=ep3moli123-eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
               ! do i=1,npole
                    dep3bt(1,i) = dep3bt(1,i)-deptemp(1,l1)
                    dep3bt(2,i) = dep3bt(2,i)-deptemp(2,l1)
                    dep3bt(3,i) = dep3bt(3,i)-deptemp(3,l1)
                  end do


                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ep3bt = ep3bt - eptemp
                  ep3moli123=ep3moli123-eptemp
                do l1 = 1,npole3b
                  i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)-deptemp(3,l1)
                end do

                 np1=3
                 npole3b=3
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 ep3bt=ep3bt+ep1moli1
                 ep3moli123=ep3moli123+ep1moli1
                do l1 = 1,npole3b
                  i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli1(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli1(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+dep1moli1(3,l1)
                end do

                 npole3b=3
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                ep3bt=ep3bt+ep1moli2
                   ep3moli123=ep3moli123+ep1moli2
                do l1 = 1,npole3b
                  i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli2(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli2(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+dep1moli2(3,l1)
                end do


                 npole3b=3
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 call empole1c_3b_Polar(
     &           npole3b,pnum,eptemp,deptemp,virtemp)
                ep3bt=ep3bt+eptemp
                   ep3moli123=ep3moli123+eptemp

                do l1 = 1,npole3b
                  i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,l1)
                end do

               !print*,"3body Cobar ep3body",shellsum,ep3moli123


          end do
        end do 
      end do  
c!$OMP END DO
c!$OMP END PARALLEL

      print*,"Ewald Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      if(moli1rmndr.eq.0) then
        return
      else
        goto 31
      end if

   31 continue

          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
          npole3b=3
          call empole1c_3b_Polar(npole3b,
     &           pnum,eptemp,deptemp,virtemp)
          ep1moli1=eptemp
        !  print*,"moli1 ep1moli1=",moli1rmndr,ep1moli1
          ep3bt=ep3bt+ep1moli1
          do l1 = 1, npole3b
            i = pnum(l1)
          !do i=1,npole
            do j = 1, 3
              dep1moli1(j,l1)=deptemp(j,l1)
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,l1)
            end do
          end do

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,moli1rmndr,
c!$OMP& nmol,ep1moli1,dep1moli1)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
c!$OMP& schedule(guided)
c      do moli1=start,last
c        do k1=1,nmollst(moli1rmndr)
c          moli2=mollst(k1,moli1rmndr)
        do moli2=moli1rmndr+1,nmol
          np1=3
          np2=6
          np3=9
          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
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

           call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,l1)=deptemp(j,l1)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                end do
              end do

              npole3b=3
              pnum(1)=imol(1,moli2)
              pnum(2)=imol(1,moli2)+1
              pnum(3)=imol(2,moli2)
            call empole1c_3b_Polar(npole3b,
     &           pnum,eptemp,deptemp,virtemp)
          ep1moli2=eptemp
          ep3bt=ep3bt-ep1moli2
          do l1 = 1, npole3b
            i = pnum(l1)
         ! do i=1,npole
            do j = 1, 3
              dep1moli2(j,l1)=deptemp(j,l1)
              dep3bt(j,i)=dep3bt(j,i)-dep1moli2(j,l1)
            end do
          end do

          ep3bt=ep3bt-ep1moli1
          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
              do l1 = 1, npole3b
                i = pnum(l1)
             ! do i=1,npole
                do j = 1, 3
                   dep3bt(j,i)=dep3bt(j,i)-dep1moli1(j,l1)
                end do
              end do

            !print*,"2body r1 ep2body",r1,ep2moli12-ep1moli1-ep1moli2


          do moli3=moli2+1,nmol
            npole3b=9
            pnum(1)=imol(1,moli1rmndr)
            pnum(2)=imol(1,moli1rmndr)+1
            pnum(3)=imol(2,moli1rmndr)
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
            call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr

            xr = x2 - x3
            yr = y2 - y3
            zr = z2 - z3
            call image(xr,yr,zr)

            xr3=xr
            yr3=yr
            zr3=zr
c            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
c            r1=sqrt(r1_2)
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)

            if((r1.lt.r3).and.(r2.lt.r3) ) then
               shellsum=r1+r2
            else if ( (r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
            else if ( (r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
            end if

                 call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ep3bt=ep3bt+eptemp
                ep3moli123=eptemp
                do l1 = 1,npole3b
                  i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,l1)
                end do
            

                  npole3b=6
                   ep3bt=ep3bt-ep2moli12
                  ep3moli123=ep3moli123-ep2moli12

                   do l1 = 1, npole3b
                    i = pnum(l1)
                  ! do i=1,npole
                    do j = 1, 3
                      dep3bt(j,i)=dep3bt(j,i)-dep2moli12(j,l1)
                    end do
                   end do

                  npole3b=6
                  np1 =3
                  pnum(1)=imol(1,moli2)
                  pnum(2)=imol(1,moli2)+1
                  pnum(3)=imol(2,moli2)
                  pnum(4)=imol(1,moli3)
                  pnum(5)=imol(1,moli3)+1
                  pnum(6)=imol(2,moli3)              
                  call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
                  ep3bt=ep3bt-eptemp
                  ep3moli123=ep3moli123-eptemp
                 do l1 = 1,npole3b
                  i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)-deptemp(3,l1)
                 end do


                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                   ep3moli123 = ep3moli123 - eptemp
                ep3bt=ep3bt-eptemp
                do l1 = 1,npole3b
                  i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)-deptemp(3,l1)
                end do

                 np1=3
                 npole3b=3
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                ep3bt=ep3bt+ep1moli1
                  ep3moli123=ep3moli123+ep1moli1
                do l1 = 1,npole3b
                  i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli1(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli1(2,l1)
                 dep3bt(3,i) = dep3bt(3,i)+dep1moli1(3,l1)
                end do

                 npole3b=3
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 ep3bt=ep3bt+ep1moli2
                  ep3moli123=ep3moli123+ep1moli2
                do l1 = 1,npole3b
                  i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli2(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli2(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+dep1moli2(3,l1)
                end do


                 npole3b=3
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 call empole1c_3b_Polar(
     &           npole3b,pnum,eptemp,deptemp,virtemp)
                ep3bt=ep3bt+eptemp
                  ep3moli123=ep3moli123+eptemp
                do l1 = 1,npole3b
                  i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,l1)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,l1)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,l1)
                end do
               !print*,"3body Cobar ep3body",shellsum,ep3moli123


          end do
        end do 
c      end do  
c!$OMP END DO
c!$OMP END PARALLEL
c      end if
      !print*,"Ewald Innerloop= ep3bt=",ep3bt
c       print*,"End inner1 single",ep3bt
      return
      end

