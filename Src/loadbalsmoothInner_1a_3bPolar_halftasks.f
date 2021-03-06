
      subroutine LoadBalNew2SmoothInnerloop1_2cut_wskin_half(
     &  ep3bt,virep3bt,dep3bt,ntript) 
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
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 term,term2,term3,term4
      real*8 rtapr2b,tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,30) 
      real*8 Rcut2b,Rcut2b2,SwitchR2b,SwitchR3b
      integer moli1rmndr,l2,kouter
      real*8 x2min,y2min,z2min,x3min,y3min,z3min
      real*8 xr3test,yr3test,zr3test,xr1t,yr1t,zr1t
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      integer taskid_offset

      rtapr2b=rtapr2b_input
      Rcut2b=cut2b_input
c      rtapr2b=6.0d0
c      Rcut2b=7.0d0
      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b
c      rtapr3b=7.0d0
c      Cobarcut3b=8.5d0
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
c      print*,"CorrectNewSmoothtaskid",taskid,start_polar(taskid),
c    & last_polar(taskid)
c     print*,"Rcut2b=",Rcut2b
c     print*,"rtapr2b=",rtapr2b
c     print*,"Cobarcut3b=",Cobarcut3b
c     print*,"rtapr3b=",rtapr3b

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start_polar,last_polar,nmollst2,mollst2,
c!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
c!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
c!$OMP& Rcut2b2)
c      if(numtasks_vdw.lt.numtasks) then
c        if(numtasks_polz.eq.int((numtasks-2)/2)) then
         taskid_offset=taskid-(numtasks_emreal+2)
c        else
c        taskid_offset=taskid-numtasks_emreal
c        end if
c      else
c        taskid_offset=taskid
c      end if

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid_offset,
!$OMP& Rcut2b2)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
!$OMP& schedule(guided)
      do kouter=start_polar(taskid_offset),last_polar(taskid_offset)
        moli1=molnew(kouter)
c      do moli1=start_polar(taskid_offset),last_polar(taskid_offset)

c       print*,"moli1,kouter,taskid",moli1,kouter,taskid
c       print*,"moli1,taskid,nmollst(moli1)",moli1,taskid,nmollst(moli1)
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
c          print*,"moli1 moli2 xr1 yr1 zr1",moli1,moli2,xr1,yr1,zr1

          xr1 = x1 - x2
          yr1 = y1 - y2
          zr1 = z1 - z2
          call image(xr1,yr1,zr1)
c          print*,"img moli1 moli2 xr1 yr1 zr1",moli1,moli2,xr1,yr1,zr1

c          x1=x1-xcell2*anint(x1/xcell2)
c          y1=y1-ycell2*anint(y1/ycell2)
c          z1=z1-zcell2*anint(z1/zcell2)

c          x2=x2-xcell2*anint(x2/xcell2)
c          y2=y2-ycell2*anint(y2/ycell2)
c          z2=z2-zcell2*anint(z2/zcell2)
c          xr1t=x1-x2
c          yr1t=y1-y2
c          zr1t=z1-z2
c        print*,"Anint moli12 xr1t yr1t zr1t",moli1,moli2,xr1t,yr1t,zr1t
         
c           x2min=x1-xr1
c          y2min=y1-yr1
c          z2min=z1-zr1
          r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
          

C
C   Pair cutoffs
C
          if(r1_2.le.Rcut2b2) then
          r1=sqrt(r1_2)
          call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
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
           if (r1.gt.rtapr2b) then
            term=(r1-rtapr2b)/SwitchR2b
            tapr2b_2=term*term
            tapr2b_3=tapr2b_2*term
            tapr2b_4=tapr2b_3*term
            tapr2b_5=tapr2b_4*term
            tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5
            do l1=1,np1
              i = pnum(l1)
              if (name(i).eq.'O') then
                term2=(-30.0d0*tapr2b_4+60.0d0*tapr2b_3-30.0d0*tapr2b_2)
     &                 /r1/SwitchR2b
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
                term2=(-30.0d0*tapr2b_4+60.0d0*tapr2b_3-30.0d0*tapr2b_2)
     &                 /r1/SwitchR2b
                dtapr2b_x(l1) = -xr1*term2
                dtapr2b_y(l1) = -yr1*term2
                dtapr2b_z(l1) = -zr1*term2
              else
                dtapr2b_x(l1) =0.0d0
                dtapr2b_y(l1) =0.0d0
                dtapr2b_z(l1) =0.0d0
              end if
            end do 
            

            do l1 = 1,npole3b
              i = pnum(l1)
              dep3bt(1,i)=dep3bt(1,i) + dep2moli12(1,l1)*tapr2b
     &                       + ep2moli12*dtapr2b_x(l1)
              dep3bt(2,i)=dep3bt(2,i) + dep2moli12(2,l1)*tapr2b
     &                       + ep2moli12*dtapr2b_y(l1)
              dep3bt(3,i)=dep3bt(3,i) + dep2moli12(3,l1)*tapr2b
     &                       + ep2moli12*dtapr2b_z(l1)
            end do
           
! NEED TO SMOOTHEN VIRIAL !

            do i=1,3
              do j=1,3
                virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)*tapr2b
              end do
            end do
            do l1=1,np1
               i=pnum(l1)
               if(name(i).eq.'O') then
c              virep3bt(1,1)=virep3bt(1,1)
c     &                          +dtapr2b_x(l1)*ep2moli12*x1
c              virep3bt(2,1)=virep3bt(2,1)
c     &                          +dtapr2b_x(l1)*ep2moli12*y1
c              virep3bt(3,1)=virep3bt(3,1)
c     &                           +dtapr2b_x(l1)*ep2moli12*z1
c              virep3bt(1,2)=virep3bt(1,2)
c     &                           +dtapr2b_y(l1)*ep2moli12*x1
c              virep3bt(2,2)=virep3bt(2,2)
c     &                           +dtapr2b_y(l1)*ep2moli12*y1
c              virep3bt(3,2)=virep3bt(3,2)
c     &                           +dtapr2b_y(l1)*ep2moli12*z1
c              virep3bt(1,3)=virep3bt(1,3)
c     &                           +dtapr2b_z(l1)*ep2moli12*x1
c              virep3bt(2,3)=virep3bt(2,3)
c     &                           +dtapr2b_z(l1)*ep2moli12*y1
c              virep3bt(3,3)=virep3bt(3,3)
c     &                           +dtapr2b_z(l1)*ep2moli12*z1

              virep3bt(1,1)=virep3bt(1,1)
     &                          +dtapr2b_x(l1)*ep2moli12*xr1
              virep3bt(2,1)=virep3bt(2,1)
     &                          +dtapr2b_x(l1)*ep2moli12*yr1
              virep3bt(3,1)=virep3bt(3,1)
     &                           +dtapr2b_x(l1)*ep2moli12*zr1
              virep3bt(1,2)=virep3bt(2,1)
              virep3bt(2,2)=virep3bt(2,2)
     &                           +dtapr2b_y(l1)*ep2moli12*yr1
              virep3bt(3,2)=virep3bt(3,2)
     &                           +dtapr2b_y(l1)*ep2moli12*zr1
              virep3bt(1,3)=virep3bt(3,1)
              virep3bt(2,3)=virep3bt(3,2)
              virep3bt(3,3)=virep3bt(3,3)
     &                           +dtapr2b_z(l1)*ep2moli12*zr1

               end if
            end do
c            do l1=np1+1,npole3b
c               i=pnum(l1)
c               if(name(i).eq.'O') then
c              virep3bt(1,1)=virep3bt(1,1)
c     &                          +dtapr2b_x(l1)*ep2moli12*x2min
c              virep3bt(2,1)=virep3bt(2,1)
c     &                          +dtapr2b_x(l1)*ep2moli12*y2min
c              virep3bt(3,1)=virep3bt(3,1)
c     &                           +dtapr2b_x(l1)*ep2moli12*z2min
c              virep3bt(1,2)=virep3bt(1,2)
c     &                           +dtapr2b_y(l1)*ep2moli12*x2min
c              virep3bt(2,2)=virep3bt(2,2)
c     &                           +dtapr2b_y(l1)*ep2moli12*y2min
c              virep3bt(3,2)=virep3bt(3,2)
c     &                           +dtapr2b_y(l1)*ep2moli12*z2min
c              virep3bt(1,3)=virep3bt(1,3)
c     &                           +dtapr2b_z(l1)*ep2moli12*x2min
c              virep3bt(2,3)=virep3bt(2,3)
c     &                           +dtapr2b_z(l1)*ep2moli12*y2min
c              virep3bt(3,3)=virep3bt(3,3)
c     &                           +dtapr2b_z(l1)*ep2moli12*z2min
c               end if
c            end do
             ep3bt=ep3bt+tapr2b*ep2moli12
           else
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
                vir2moli12(j,i)= 0.0d0
              end do
            end do
          end if

c          do k2=1,nmollst2(moli2)
c            moli3=mollst2(k2,moli2)
          do k2=1,nmollst(moli2)
            moli3=mollst(k2,moli2)

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
            x3min=x1-xr1
            y3min=y1-yr1
            z3min=z1-zr1

            xr3=x2-x3
            yr3=y2-y3
            zr3=z2-z3
            call image(xr3,yr3,zr3)
c        print*,"moli123 xr3 yr3 zr3",moli1,moli2,moli3,xr3,yr3,z3
c            xr3test=x2min-y3min
c            yr3test=y2min-y3min
c            zr3test=z2min-z3min
c            call image(xr3test,yr3test,zr3test)
c        print*,"Testxr3 etc",moli1,moli2,moli3,xr3test,yr3test,zr3test

            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
            r1=sqrt(r1_2)         
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)
            
c          x1=x1-xcell2*anint(x1/xcell2)
c          y1=y1-ycell2*anint(y1/ycell2)
c          z1=z1-zcell2*anint(z1/zcell2)



c          x2=x2-xcell2*anint(x2/xcell2)
c          y2=y2-ycell2*anint(y2/ycell2)
c          z2=z2-zcell2*anint(z2/zcell2)
    

c          x3=x3-xcell2*anint(x3/xcell2)
c          y3=y3-ycell2*anint(y3/ycell2)
c          z3=z3-zcell2*anint(z3/zcell2)
c          print*,"xr2 anint",x1-x3

           
            if ((r1.lt.r3).and.(r2.lt.r3) ) then 
              shellsum=r1+r2
              if (shellsum.le.Cobarcut3b) then
                call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ntript=ntript+1
                ep3moli123=eptemp
                do l1 = 1,npole3b
                  i = pnum(l1)
                  dep3moli123(1,l1) = deptemp(1,l1)
                  dep3moli123(2,l1) = deptemp(2,l1)
                  dep3moli123(3,l1) = deptemp(3,l1)
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
c                  dtapr3b_xr_xr = term4*((xr1/r1)*xr1+(xr2/r2)*xr2)
c                  dtapr3b_xr_yr = term4*((xr1/r1)*yr1+(xr2/r2)*yr2)
c                  dtapr3b_xr_zr = term4*((xr1/r1)*zr1+(xr2/r2)*zr2)
c                  dtapr3b_yr_yr = term4*((yr1/r1)*yr1+(yr2/r2)*yr2)
c                  dtapr3b_yr_zr = term4*((yr1/r1)*zr1+(yr2/r2)*zr2)
c                  dtapr3b_zr_zr = term4*((zr1/r1)*zr1+(zr2/r2)*zr2)
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
                call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ntript=ntript+1
                ep3moli123=eptemp
                do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,l1) = deptemp(1,l1)
                    dep3moli123(2,l1) = deptemp(2,l1)
                    dep3moli123(3,l1) = deptemp(3,l1)
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
                call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,l1) = deptemp(1,l1)
                    dep3moli123(2,l1) = deptemp(2,l1)
                    dep3moli123(3,l1) = deptemp(3,l1)
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
                   do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                      dep3moli123(j,l1)=dep3moli123(j,l1)
     &                                 -dep2moli12(j,l1)
                    end do
                   end do
                   do i=1,3
                    do j=1,3
                     vir3moli123(j,i)=vir3moli123(j,i)-vir2moli12(j,i)
                    end do
                   end do
                 else
                  call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
                  ep3moli123=ep3moli123-eptemp
                  do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                     dep3moli123(j,l1)=dep3moli123(j,l1)-deptemp(j,l1)
                    end do
                  end do
                  do i=1,3
                   do j=1,3
                    vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                   end do
                  end do
                 end if

c               if (r3.le.Rcut2b) then
c               if (r3.le.1.0d12) then
                  npole3b=6
                  np1 =3
                  pnum(1)=imol(1,moli2)
                  pnum(2)=imol(1,moli2)+1
                  pnum(3)=imol(2,moli2)
                  pnum(4)=imol(1,moli3)
                  pnum(5)=imol(1,moli3)+1
                  pnum(6)=imol(2,moli3)              
                  call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
c                  if (r3.gt.rtapr2b) then               
c                  if (r3.gt.1.0d12) then
c                    term3=(r3-rtapr2b)/SwitchR2b
c                    tapr2b_2=term3*term3
c                    tapr2b_3=tapr2b_2*term3
c                    tapr2b_4=tapr2b_3*term3
c                    tapr2b_5=tapr2b_4*term3
c                    tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
c     &                 - 6.0d0*tapr2b_5
c
c                     do l1=1,np1
c                      i=pnum(l1)
c                      if (name(i).eq.'O') then
c                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
c     &                            30.0d0*tapr2b_2)/SwitchR2b/r3
c                        dtapr2b_x(l1) = term4*xr3
c                        dtapr2b_y(l1) = term4*yr3
c                        dtapr2b_z(l1) = term4*zr3
c                      else
c                        dtapr2b_x(l1) =0.0d0
c                        dtapr2b_y(l1) =0.0d0
c                        dtapr2b_z(l1) =0.0d0
c                      end if
c                     end do
c                     do l1=np1+1,npole3b
c                       i=pnum(l1)
c                      if (name(i).eq.'O') then
c                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
c     &                            30.0d0*tapr2b_2)/SwitchR2b/r3
c                        dtapr2b_x(l1)= -term4*xr3
c                        dtapr2b_y(l1)= -term4*yr3
c                        dtapr2b_z(l1)= -term4*zr3
c                      else
c                        dtapr2b_x(l1) =0.0d0
c                        dtapr2b_y(l1) =0.0d0
c                        dtapr2b_z(l1) =0.0d0
c                      end if
c                     end do
c                     ep3moli123=ep3moli123-tapr2b*eptemp
c                     do l1 = 1, npole3b
c                       i = pnum(l1)
c                       l2=l1+np1
c                       dep3moli123(1,l2)=dep3moli123(1,l2)
c     &                  -deptemp(1,l1)*tapr2b- eptemp*dtapr2b_x(l1)
c                       dep3moli123(2,l2)=dep3moli123(2,l2)
c     &                  -deptemp(2,l1)*tapr2b- eptemp*dtapr2b_y(l1)
c                       dep3moli123(3,l2)=dep3moli123(3,l2)
c     &                  -deptemp(3,l1)*tapr2b- eptemp*dtapr2b_z(l1)
c                     end do
c                     do i=1,3
c                       do j=1,3
c                        virtemp(j,i)=virtemp(j,i)*tapr2b
c                       end do
c                     end do
c                     do l1=1,npole3b
c                       i = pnum(l1)
c                     virtemp(1,1)=virtemp(1,1)+dtapr2b_x(l1)*eptemp*x(i)
c                     virtemp(2,1)=virtemp(2,1)+dtapr2b_x(l1)*eptemp*y(i)
c                     virtemp(3,1)=virtemp(3,1)+dtapr2b_x(l1)*eptemp*z(i)
c                     virtemp(1,2)=virtemp(1,2)+dtapr2b_y(l1)*eptemp*x(i)
c                     virtemp(2,2)=virtemp(2,2)+dtapr2b_y(l1)*eptemp*y(i)
c                     virtemp(3,2)=virtemp(3,2)+dtapr2b_y(l1)*eptemp*z(i)
c                     virtemp(1,3)=virtemp(1,3)+dtapr2b_z(l1)*eptemp*x(i)
c                     virtemp(2,3)=virtemp(2,3)+dtapr2b_z(l1)*eptemp*y(i)
c                     virtemp(3,3)=virtemp(3,3)+dtapr2b_z(l1)*eptemp*z(i)
c                     end do
c                     do i=1,3
c                       do j=1,3
c                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
c                       end do
c                     end do
c                  else
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       end do
                     end do
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       end do
                     end do
c                  end if
c               end if

c               if(r2.le.Rcut2b) then
c               if(r2.le.1.0d12) then
                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
c                 if (r2.gt.rtapr2b) then
c                 if (r2.gt.1.0d12) then
c                   term3=(r2-rtapr2b)/SwitchR2b
c                   tapr2b_2=term3*term3
c                   tapr2b_3=tapr2b_2*term3
c                   tapr2b_4=tapr2b_3*term3
c                   tapr2b_5=tapr2b_4*term3
c                   tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
c     &                - 6.0d0*tapr2b_5
c
c                     do l1=1,np1
c                      i=pnum(l1)
c                      if (name(i).eq.'O') then
c                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
c     &                            30.0d0*tapr2b_2)/SwitchR2b/r2
c                        dtapr2b_x(l1) = term4*xr2
c                        dtapr2b_y(l1) = term4*yr2
c                        dtapr2b_z(l1) = term4*zr2
c                      else
c                        dtapr2b_x(l1) = 0.0d0
c                        dtapr2b_y(l1) = 0.0d0
c                        dtapr2b_z(l1) = 0.0d0
c                      end if
c                     end do
c                     do l1=np1+1,npole3b
c                       i=pnum(l1)
c                      if (name(i).eq.'O') then
c                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
c     &                            30.0d0*tapr2b_2)/SwitchR2b/r2
c                        dtapr2b_x(l1) = -term4*xr2
c                        dtapr2b_y(l1) = -term4*yr2
c                        dtapr2b_z(l1) = -term4*zr2
c                      else
c                       dtapr2b_x(l1)= 0.0d0
c                       dtapr2b_y(l1)= 0.0d0
c                       dtapr2b_z(l1)= 0.0d0
c                      end if
c                     end do
c                     ep3moli123 = ep3moli123 -tapr2b*eptemp
c                     do l1 = 1, np1
c                       i = pnum(l1)
c                       dep3moli123(1,l1)=dep3moli123(1,l1)
c     &                   -deptemp(1,l1)*tapr2b- eptemp*dtapr2b_x(l1)
c                       dep3moli123(2,l1)=dep3moli123(2,l1)
c     &                   -deptemp(2,l1)*tapr2b- eptemp*dtapr2b_y(l1)
c                       dep3moli123(3,l1)=dep3moli123(3,l1)
c     &                   -deptemp(3,l1)*tapr2b- eptemp*dtapr2b_z(l1)
c                     end do
c                     do l1=np1+1,npole3b
c                       i = pnum(l1)
c                       l2=l1+np1
c                       dep3moli123(1,l2)=dep3moli123(1,l2)
c     &                  -deptemp(1,l1)*tapr2b-eptemp*dtapr2b_x(l1)
c                       dep3moli123(2,l2)=dep3moli123(2,l2)
c     &                  -deptemp(2,l1)*tapr2b-eptemp*dtapr2b_y(l1)
c                       dep3moli123(3,l2)=dep3moli123(3,l2)
c     &                  -deptemp(3,l1)*tapr2b-eptemp*dtapr2b_z(l1)
c                     end do
c                     do i=1,3
c                       do j=1,3
c                        virtemp(j,i)=virtemp(j,i)*tapr2b
c                       end do
c                     end do
c                     do l1 = 1, npole3b
c                      i = pnum(l1)
c                     virtemp(1,1)=virtemp(1,1)+dtapr2b_x(l1)*eptemp*x(i)
c                     virtemp(2,1)=virtemp(2,1)+dtapr2b_x(l1)*eptemp*y(i)
c                     virtemp(3,1)=virtemp(3,1)+dtapr2b_x(l1)*eptemp*z(i)
c                     virtemp(1,2)=virtemp(1,2)+dtapr2b_y(l1)*eptemp*x(i)
c                     virtemp(2,2)=virtemp(2,2)+dtapr2b_y(l1)*eptemp*y(i)
c                     virtemp(3,2)=virtemp(3,2)+dtapr2b_y(l1)*eptemp*z(i)
c                     virtemp(1,3)=virtemp(1,3)+dtapr2b_z(l1)*eptemp*x(i)
c                     virtemp(2,3)=virtemp(2,3)+dtapr2b_z(l1)*eptemp*y(i)
c                     virtemp(3,3)=virtemp(3,3)+dtapr2b_z(l1)*eptemp*z(i)
c                     end do
c                     do i=1,3
c                       do j=1,3
c                        vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
c                       end do
c                     end do                  
c                  else
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, np1
                      i = pnum(l1)
                      do j = 1, 3
                       dep3moli123(j,l1)=dep3moli123(j,l1)-deptemp(j,l1)
                      end do
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1,3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                      end do
                     end do
c                  end if
c               end if

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

                  do l1 = 1,np1                  
                    i = pnum(l1)
                    if(name(i).eq.'O') then
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_ix
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_iy
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_iz
                    else
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,l1)*tapr3b
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,l1)*tapr3b
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,l1)*tapr3b
                    end if
                  end do
                  do l1 = np1+1,np2
                    i = pnum(l1)
                    if(name(i).eq.'O') then
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_jx
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_jy
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_jz
                    else
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,l1)*tapr3b
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,l1)*tapr3b
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,l1)*tapr3b
                    end if
                  end do
                  do l1 = np2+1,npole3b
                    i = pnum(l1)
                    if(name(i).eq.'O') then
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_kx
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_ky
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,l1)*tapr3b
     &                            + ep3moli123*dtapr3b_kz
                    else
                    dep3bt(1,i) = dep3bt(1,i)+dep3moli123(1,l1)*tapr3b
                    dep3bt(2,i) = dep3bt(2,i)+dep3moli123(2,l1)*tapr3b
                    dep3bt(3,i) = dep3bt(3,i)+dep3moli123(3,l1)*tapr3b
                    end if
                  end do 
!  NEED TO CORRECT VIRIAL W/ TAPERING FUNC
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
c                  do l1=1,np1
c                     i = pnum(l1)
c                    if(name(i).eq.'O') then
c
c                    vir3moli123(1,1)=vir3moli123(1,1)
c     &                 +dtapr3b_ix*ep3moli123*x1
c                    vir3moli123(2,1)=vir3moli123(2,1)
c     &                 +dtapr3b_ix*ep3moli123*y1
c                    vir3moli123(3,1)=vir3moli123(3,1)
c     &                 +dtapr3b_ix*ep3moli123*z1
c                    vir3moli123(1,2)=vir3moli123(1,2)
c     &                 +dtapr3b_iy*ep3moli123*x1
c                    vir3moli123(2,2)=vir3moli123(2,2)
c     &                 +dtapr3b_iy*ep3moli123*y1
c                    vir3moli123(3,2)=vir3moli123(3,2)
c     &                 +dtapr3b_iy*ep3moli123*z1
c                    vir3moli123(1,3)=vir3moli123(1,3)
c     &                 +dtapr3b_iz*ep3moli123*x1
c                    vir3moli123(2,3)=vir3moli123(2,3)
c     &                 +dtapr3b_iz*ep3moli123*y1
c                    vir3moli123(3,3)=vir3moli123(3,3)
c     &                 +dtapr3b_iz*ep3moli123*z1
c                    end if
c                  end do
c                  do l1=np1+1,np2
c                     i = pnum(l1)
c                    if(name(i).eq.'O') then
c
c                    vir3moli123(1,1)=vir3moli123(1,1)
c     &                   +dtapr3b_jx*ep3moli123*x2
c                    vir3moli123(2,1)=vir3moli123(2,1)
c     &                   +dtapr3b_jx*ep3moli123*y2
c                    vir3moli123(3,1)=vir3moli123(3,1)
c     &                   +dtapr3b_jx*ep3moli123*z2
c                    vir3moli123(1,2)=vir3moli123(1,2)
c     &                   +dtapr3b_jy*ep3moli123*x2
c                    vir3moli123(2,2)=vir3moli123(2,2)
c     &                   +dtapr3b_jy*ep3moli123*y2
c                    vir3moli123(3,2)=vir3moli123(3,2)
c     &                   +dtapr3b_jy*ep3moli123*z2
c                    vir3moli123(1,3)=vir3moli123(1,3)
c     &                   +dtapr3b_jz*ep3moli123*x2
c                    vir3moli123(2,3)=vir3moli123(2,3)
c     &                   +dtapr3b_jz*ep3moli123*y2
c                    vir3moli123(3,3)=vir3moli123(3,3)
c     &                   +dtapr3b_jz*ep3moli123*z2
c                    end if
c                  end do
c                  do l1=np2+1,npole3b
c                     i = pnum(l1)
c                    if(name(i).eq.'O') then
c                    vir3moli123(1,1)=vir3moli123(1,1)
c     &                  +dtapr3b_kx*ep3moli123*x3
c                    vir3moli123(2,1)=vir3moli123(2,1)
c     &                  +dtapr3b_kx*ep3moli123*y3
c                    vir3moli123(3,1)=vir3moli123(3,1)
c     &                  +dtapr3b_kx*ep3moli123*z3
c                    vir3moli123(1,2)=vir3moli123(1,2)
c     &                  +dtapr3b_ky*ep3moli123*x3
c                    vir3moli123(2,2)=vir3moli123(2,2)
c     &                  +dtapr3b_ky*ep3moli123*y3
c                    vir3moli123(3,2)=vir3moli123(3,2)
c     &                  +dtapr3b_ky*ep3moli123*z3
c                    vir3moli123(1,3)=vir3moli123(1,3)
c     &                  +dtapr3b_kz*ep3moli123*x3
c                    vir3moli123(2,3)=vir3moli123(2,3)
c     &                  +dtapr3b_kz*ep3moli123*y3
c                    vir3moli123(3,3)=vir3moli123(3,3)
c     &                  +dtapr3b_kz*ep3moli123*z3
c                    end if
c                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
               else
                  ep3bt = ep3bt + ep3moli123
                  do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,l1)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
               end if 
c             print*,"3body shellsum",shellsum,ep3moli123
            end if
          end do
        end do 
      end do  
!$OMP END DO
!$OMP END PARALLEL
      !print*,"Energy,taskid,loadbalsmoothInner_1a_3bPolar",ep3bt,taskid
      return
      end

