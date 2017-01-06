
      subroutine s1btotfieldLoadBalNew2SmoothInnerloop1_rtrndeptemp2
     &  (start,last,ep3bt,virep3bt,dep3bt,moli1rmndr) 
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use mpidat
      use cobar
c      use cell
      implicit none
      real*8  ep2moli12
      real*8  dep2moli12(3,npole)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript,len1,len2,len3
c      real*8 eptemp,deptemp(3,30)
      real*8 eptemp,deptemp2(3,npole)
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
      real*8 dtapr2b_x(npole),dtapr2b_y(npole),dtapr2b_z(npole)
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 term,term2,term3,term4
      real*8 rtapr2b,tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,npole)
      real*8 ep1moli1,dep1moli1(3,npole)
      real*8 ep1moli2,dep1moli2(3,npole)
      real*8 Rcut2b,Rcut2b2,SwitchR2b,SwitchR3b
      integer moli1rmndr,l2,kouter
      real*8 x2min,y2min,z2min,x3min,y3min,z3min
      real*8 xr3test,yr3test,zr3test,xr1t,yr1t,zr1t
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      rtapr2b=rtapr2b_input
      Rcut2b=cut2b_input
c      rtapr2b=6.0d0
c      Rcut2b=7.0d0
      rtapr2b=1.0d12
      Rcut2b=1.0d12
      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b
c      rtapr3b=7.0d0
c      Cobarcut3b=8.5d0
      rtapr3b=1.0d12
      Cobarcut3b=1.0d12
      SwitchR3b=Cobarcut3b-rtapr3b
c      print*,"xcell2=",xcell2

      ntript=0
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
c      print*,"12/17 Original Version!"
      print*,"Totfieldversion taskid",taskid
c      print*,"Rcut2b=",Rcut2b
c      print*,"rtapr2b=",rtapr2b
c      print*,"Cobarcut3b=",Cobarcut3b
c      print*,"rtapr3b=",rtapr3b
        print*,"start last",start,last,moli1rmndr

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,nmollst2,mollst2,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
!$OMP& Rcut2b2,nmol,npole,start,last)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
!$OMP& schedule(guided)
c      do kouter=start_polar(taskid),last_polar(taskid)
c        moli1=molnew(kouter)
c      do moli1=1,nmol
      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
            len1=imol(2,moli1)-imol(1,moli1)+1

            if(len1 .eq. 3) then
              npole3b=3
              pnum(1)=imol(1,moli1)
              pnum(2)=imol(1,moli1)+1
              pnum(3)=imol(2,moli1)
              np1=3
            end if

            if(len1 .eq. 5) then
              npole3b=5
              pnum(1)=imol(1,moli1)
              pnum(2)=imol(1,moli1)+1
              pnum(3)=imol(1,moli1)+2
              pnum(4)=imol(1,moli1)+3
              pnum(5)=imol(2,moli1)
              np1=5
            end if

          call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(npole3b,
     &           pnum,eptemp,deptemp2,virtemp)
          ep1moli1=eptemp
c          print*,"moli1 ep1moli1=",moli1,ep1moli1
          ep3bt=ep3bt+ep1moli1
c          do l1 = 1, npole3b
c            i = pnum(l1)
          do i=1,npole
            do j = 1, 3
              dep1moli1(j,i)=deptemp2(j,i)
c              print*,"i j dep1moli1(j,i)",i,j,dep1moli1(j,i)
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,i)
            end do
          end do


      end do  
!$OMP END DO
!$OMP END PARALLEL
c      print*,"End testvirpolar",ep3bt 

      if(moli1rmndr.eq.0) then
        return
      else
        goto 31
      end if

   31 continue

c      do kouter=start_polar(taskid),last_polar(taskid)
c        moli1=molnew(kouter)
c      do moli1=1,nmol
c      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
          npole3b=3
          call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(npole3b,
     &           pnum,eptemp,deptemp2,virtemp)
          ep1moli1=eptemp
        !  print*,"moli1 ep1moli1=",moli1rmndr,ep1moli1
          ep3bt=ep3bt+ep1moli1
c          do l1 = 1, npole3b
c            i = pnum(l1)
          do i=1,npole
            do j = 1, 3
              dep1moli1(j,i)=deptemp2(j,i)
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,i)
            end do
          end do

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,nmollst2,mollst2,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
!$OMP& Rcut2b2,nmol,npole,start,last,moli1rmndr,ep1moli1,dep1moli1)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
!$OMP& schedule(guided)


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

          call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(npole3b,
     &           pnum,eptemp,deptemp2,virtemp)
          ep2moli12=eptemp
          ep3bt=ep3bt+ep2moli12
c          do l1 = 1, npole3b
c            i = pnum(l1)
          do i=1,npole
            do j = 1, 3
              dep2moli12(j,i)=deptemp2(j,i)
              dep3bt(j,i)=dep3bt(j,i)+dep2moli12(j,i)  
            end do
          end do
          do i=1,3
            do j=1,3
              vir2moli12(j,i)=virtemp(j,i)
            end do
          end do

          npole3b=3
          pnum(1)=imol(1,moli2)
          pnum(2)=imol(1,moli2)+1
          pnum(3)=imol(2,moli2)
          call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(npole3b,
     &           pnum,eptemp,deptemp2,virtemp)
          ep1moli2=eptemp
          ep3bt=ep3bt-ep1moli2
c          do l1 = 1, npole3b
c            i = pnum(l1)
          do i=1,npole
            do j = 1, 3
              dep1moli2(j,i)=deptemp2(j,i)              
              dep3bt(j,i)=dep3bt(j,i)-dep1moli2(j,i)
            end do
          end do
        
          ep3bt=ep3bt-ep1moli1
          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
c              do l1 = 1, npole3b
c                i = pnum(l1)
              do i=1,npole
                do j = 1, 3
                   dep3bt(j,i)=dep3bt(j,i)-dep1moli1(j,i)
                end do
              end do

            print*,"2body r1 ep2body",r1,ep2moli12-ep1moli1-ep1moli2

         !   print*,"moli1 moli2 ep2moli12",moli1,moli2,ep2moli12
          do moli3=moli2+1,nmol
            npole3b=9
            np3=9
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

                call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(
     &           npole3b,pnum,eptemp,deptemp2,virtemp)
                ntript=ntript+1
          !   print*,"moli1moli2moli3 eptemp",moli1,moli2,moli3,eptemp
                ep3bt=ep3bt+eptemp
                ep3moli123=eptemp

c                do l1 = 1,npole3b
c                  i = pnum(l1)
                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp2(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp2(3,i) 
                end do
            
                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli2)
                 pnum(5)=imol(1,moli2)+1
                 pnum(6)=imol(2,moli2)
                   ep3bt=ep3bt-ep2moli12
                  ep3moli123=ep3moli123-ep2moli12

                  ! do l1 = 1, npole3b
                  !  i = pnum(l1)
                   do i=1,npole
                    do j = 1, 3
                      dep3bt(j,i)=dep3bt(j,i)-dep2moli12(j,i)
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
                  call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(
     &            npole3b,pnum,eptemp,deptemp2,virtemp)
                ep3bt=ep3bt-eptemp
                  ep3moli123=ep3moli123-eptemp

                ! do l1 = 1,npole3b
                !  i = pnum(l1)
                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp2(2,i)
                 dep3bt(3,i) = dep3bt(3,i)-deptemp2(3,i)
                end do


                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(
     &            npole3b,pnum,eptemp,deptemp2,virtemp)
                ep3bt=ep3bt-eptemp
                  ep3moli123=ep3moli123-eptemp

                !do l1 = 1,npole3b
                !  i = pnum(l1)
                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp2(2,i)
                 dep3bt(3,i) = dep3bt(3,i)-deptemp2(3,i)
                end do

                 np1=3
                 npole3b=3
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                ep3bt=ep3bt+ep1moli1
                  ep3moli123=ep3moli123+ep1moli1

                !do l1 = 1,npole3b
                !  i = pnum(l1)
                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli1(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli1(2,i)
                 dep3bt(3,i) = dep3bt(3,i)+dep1moli1(3,i)
                end do

                 npole3b=3
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                ep3bt=ep3bt+ep1moli2
                  ep3moli123=ep3moli123+ep1moli2

                !do l1 = 1,npole3b
                !  i = pnum(l1)
                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli2(2,i)
                 dep3bt(3,i) = dep3bt(3,i)+dep1moli2(3,i)
                end do


                 npole3b=3
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 call empole1a_3b_Polar_totfieldnpole3b_rtrndeptemp2(
     &           npole3b,pnum,eptemp,deptemp2,virtemp)
                ep3bt=ep3bt+eptemp
                  ep3moli123=ep3moli123+eptemp

                !do l1 = 1,npole3b
                !  i = pnum(l1)
                do i=1,npole 
                  dep3bt(1,i) = dep3bt(1,i)+deptemp2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp2(2,i)
                 dep3bt(3,i) = dep3bt(3,i)+deptemp2(3,i)
                end do
               print*,"3body Cobar ep3body",shellsum,ep3moli123


          end do
        end do 
c      end do  
!$OMP END DO
!$OMP END PARALLEL
c      print*,"End testvirpolar",ep3bt 
      return
      end

