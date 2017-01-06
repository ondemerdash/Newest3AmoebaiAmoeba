
      subroutine ClustLoadBalNew2SmoothInnerloop1_2cut_wskinOrig 
     &  (ep3bt,virep3bt,dep3bt,ntript) 
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use mpidat
      use cobar
      use neigh2clust
c      use cell
      implicit none
      real*8  ep2moli12
      real*8  dep2moli12(3,36)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp,deptemp(3,36)
      real*8 vir2moli12(3,3)
      real*8 virtemp(3,3)
      integer pnum(36),npole3b,moli1,moli2,moli3,np1,np2,np3
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      logical do2
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,xr,yr,zr
      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r1,r2,r3,r1_2,r2_2,r3_2
      real*8 shellsum
      real*8 tapr2b,tapr3b
      real*8 dtapr2b_x(36),dtapr2b_y(36),dtapr2b_z(36)
      real*8 dtapr3b_ix,dtapr3b_iy,dtapr3b_iz
      real*8 dtapr3b_jx,dtapr3b_jy,dtapr3b_jz
      real*8 dtapr3b_kx,dtapr3b_ky,dtapr3b_kz
      real*8 term,term2,term3,term4
      real*8 rtapr2b,tapr2b_2,tapr2b_3,tapr2b_4,tapr2b_5
      real*8 tapr3b_2,tapr3b_3,tapr3b_4,tapr3b_5      
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,36) 
      real*8 Rcut2b,Rcut2b2,SwitchR2b,SwitchR3b
      integer moli1rmndr,l2,kouter
      real*8 x2min,y2min,z2min,x3min,y3min,z3min
      real*8 xr3test,yr3test,zr3test,xr1t,yr1t,zr1t
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      integer clust1,clust2,clust3,npole3b12
      real*8 ep1moli1,ep1moli2,dep1moli1(3,36),dep1moli2(3,36)
      real*8 vir1moli1(3,3),vir1moli2(3,3)
      real*8 ep2moli12nosubtr,dep2moli12nosubtr(3,36)
      real*8 vir2moli12nosubtr(3,3)

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
      print*,"ClUSTUSINGLIST!taskid,start,last",taskid,
     &  start_polar(taskid),last_polar(taskid)

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
c      print*,"ClustClust Version!"
c      print*,"CorrectNewSmoothtaskid",taskid,start_polar(taskid),
c    & last_polar(taskid)
c     print*,"Rcut2b=",Rcut2b
c     print*,"rtapr2b=",rtapr2b
c     print*,"Cobarcut3b=",Cobarcut3b
c     print*,"rtapr3b=",rtapr3b


!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
!$OMP& Rcut2b2,sizeclust,clust,clust_cm)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
!$OMP& schedule(guided)
c      do kouter=start_polar(taskid),last_polar(taskid)
c        moli1=molnew(kouter)
      do clust1=start_polar(taskid),last_polar(taskid)

c       print*,"moli1,kouter,taskid",moli1,kouter,taskid
c       print*,"moli1,taskid,nmollst(moli1)",moli1,taskid,nmollst(moli1)
           npole3b=0
           np1=0
           do i=1,sizeclust(clust1)
             moli1=clust(i,clust1)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           npole3b=np1
          call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,ep1moli1,
     &         dep1moli1,vir1moli1)
          ep3bt = ep3bt + ep1moli1
          do l1 = 1, npole3b
             i = pnum(l1)
             do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)+dep1moli1(j,l1)
             end do
          end do
          do i=1,3
            do j=1,3
              virep3bt(j,i) = virep3bt(j,i)+vir1moli1(j,i)
            end do
          end do


        do k1=1,nmollst(clust1)
           clust2=mollst(k1,clust1)
           xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
           yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
           zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
           call image(xr1,yr1,zr1)
           r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
           r1=sqrt(r1_2)

          if(r1_2.le.Rcut2b2) then
          np2=0
          do i=1,sizeclust(clust2)
             moli1=clust(i,clust2)
             pnum(np1+3*(i-1)+1)=imol(1,moli1)
             pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
             pnum(np1+3*(i-1)+3)=imol(2,moli1)
             np2=np2+3
          end do

          npole3b=np1+np2
          call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &         ep2moli12nosubtr,dep2moli12nosubtr,vir2moli12nosubtr)

          ep2moli12 = ep2moli12nosubtr - ep1moli1
          do l1 = 1, np1
             i = pnum(l1)
             do j = 1, 3
              dep2moli12(j,l1) = dep2moli12nosubtr(j,l1)-dep1moli1(j,l1)
             end do
          end do
          do i=1,3
            do j=1,3
              vir2moli12(j,i) = vir2moli12nosubtr(j,i)-vir1moli1(j,i)
            end do
          end do


          npole3b=np2
          do i=1,sizeclust(clust2)
             moli1=clust(i,clust2)
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
          end do

          call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,ep1moli2,
     &         dep1moli2,vir1moli2)
          ep2moli12 = ep2moli12 - ep1moli2
          do l1 =1, npole3b
             l2=l1+np1
             do j=1,3
              dep2moli12(j,l2) = dep2moli12(j,l2)-dep1moli2(j,l1)
             end do
          end do
          do i=1,3
            do j=1,3
              vir2moli12(j,i) = vir2moli12(j,i)-vir1moli2(j,i)
            end do
          end do

C
C   Pair cutoffs
C
           np1=0
           do i=1,sizeclust(clust1)
             moli1=clust(i,clust1)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
          np2=0
          do i=1,sizeclust(clust2)
             moli1=clust(i,clust2)
             pnum(np1+3*(i-1)+1)=imol(1,moli1)
             pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
             pnum(np1+3*(i-1)+3)=imol(2,moli1)
             np2=np2+3
          end do

          npole3b=np1+np2

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
          np2=0
          do i=1,sizeclust(mollst(k1,clust1))
             moli1=clust(i,mollst(k1,clust1))
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
             np2=np2+3
          end do
          npole3b=np2

          call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,ep1moli2,
     &         dep1moli2,vir1moli2)
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
          do k2=1,nmollst(clust2)
            clust3=mollst(k2,clust2)

            npole3b=0
            np1=0
            np2=0
            np3=0
            do i=1,sizeclust(clust1)
              moli1=clust(i,clust1)
              npole3b=npole3b+3
              pnum(3*(i-1)+1)=imol(1,moli1)
              pnum(3*(i-1)+2)=imol(1,moli1)+1
              pnum(3*(i-1)+3)=imol(2,moli1)
              np1=np1+3
            end do
            do i=1,sizeclust(clust2)
              moli1=clust(i,clust2)
              npole3b=npole3b+3
              pnum(np1+3*(i-1)+1)=imol(1,moli1)
              pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np1+3*(i-1)+3)=imol(2,moli1)
              np2=np2+3
            end do
            do i=1,sizeclust(clust3)
              moli1=clust(i,clust3)
              npole3b=npole3b+3
              pnum(np1+np2+3*(i-1)+1)=imol(1,moli1)
              pnum(np1+np2+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np1+np2+3*(i-1)+3)=imol(2,moli1)
              np3=np3+3
            end do

           xr2=clust_cm(1,clust3)-clust_cm(1,clust1)
           yr2=clust_cm(2,clust3)-clust_cm(2,clust1)
           zr2=clust_cm(3,clust3)-clust_cm(3,clust1)
           call image(xr2,yr2,zr2)
           r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
           r2=sqrt(r2_2)

           xr3=clust_cm(1,clust3)-clust_cm(1,clust2)
           yr3=clust_cm(2,clust3)-clust_cm(2,clust2)
           zr3=clust_cm(3,clust3)-clust_cm(3,clust2)
           call image(xr3,yr3,zr3)
           r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
           r3=sqrt(r3_2)
           shellsum=r1+r2+r3           
c            if ((r1.lt.r3).and.(r2.lt.r3) ) then 
c              shellsum=r1+r2
            if (shellsum.le.Cobarcut3b) then

c              if (shellsum.le.Cobarcut3b) then
                call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp,deptemp,virtemp)
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
                virtapr3b_xx=term4*( xr1/r1*xr1 + xr2/r2*xr2)
                virtapr3b_xy=term4*( xr1/r1*yr1 + xr2/r2*yr2)
                virtapr3b_xz=term4*( xr1/r1*zr1 + xr2/r2*zr2)
                virtapr3b_yy=term4*( yr1/r1*yr1 + yr2/r2*yr2)
                virtapr3b_yz=term4*( yr1/r1*zr1 + yr2/r2*zr2)
                virtapr3b_zz=term4*( zr1/r1*zr1 + zr2/r2*zr2)
                end if
c              end if             
            end if 
c            else if ((r1.lt.r2).and.(r3.lt.r2)) then
c              shellsum=r1+r3
c              if (shellsum.le.Cobarcut3b) then
c                call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
c     &           eptemp,deptemp,virtemp)
c                ntript=ntript+1
c                ep3moli123=eptemp
c                do l1 = 1,npole3b
c                    i = pnum(l1)
c                    dep3moli123(1,l1) = deptemp(1,l1)
c                    dep3moli123(2,l1) = deptemp(2,l1)
c                    dep3moli123(3,l1) = deptemp(3,l1)
c                end do
c                  do i=1,3
c                    do j=1,3
c                       vir3moli123(j,i)=virtemp(j,i)
c                    end do
c                  end do
c
c                 if (shellsum.gt.rtapr3b) then
c                   term3=(shellsum-rtapr3b)/SwitchR3b
c                   tapr3b_2=term3*term3
c                   tapr3b_3=tapr3b_2*term3
c                   tapr3b_4=tapr3b_3*term3
c                   tapr3b_5=tapr3b_4*term3
c                   tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
c     &                - 6.0d0*tapr3b_5
c                   term4 = (-30.0d0*tapr3b_4 + 60.0d0*tapr3b_3 - 
c     &                     30.0d0*tapr3b_2)/SwitchR3b
c                   dtapr3b_ix = term4*(xr1/r1)
c                   dtapr3b_iy = term4*(yr1/r1)
c                   dtapr3b_iz = term4*(zr1/r1)
c                   dtapr3b_jx = term4*(-xr1/r1+xr3/r3)
c                   dtapr3b_jy = term4*(-yr1/r1+yr3/r3)
c                   dtapr3b_jz = term4*(-zr1/r1+zr3/r3)
c                   dtapr3b_kx = term4*(-xr3/r3)
c                   dtapr3b_ky = term4*(-yr3/r3)
c                   dtapr3b_kz = term4*(-zr3/r3) 
c                   virtapr3b_xx=term4*((xr1/r1)*xr1+(xr3/r3)*xr3)
c                   virtapr3b_xy=term4*((xr1/r1)*yr1+(xr3/r3)*yr3)
c                   virtapr3b_xz=term4*((xr1/r1)*zr1+(xr3/r3)*zr3)
c                   virtapr3b_yy=term4*((yr1/r1)*yr1+(yr3/r3)*yr3)
c                   virtapr3b_yz=term4*((yr1/r1)*zr1+(yr3/r3)*zr3)
c                   virtapr3b_zz=term4*((zr1/r1)*zr1+(zr3/r3)*zr3)
c                 end if
c              end if
c            else if ((r2.lt.r1).and.(r3.lt.r1)) then
c              shellsum=r2+r3
c               if (shellsum.le.Cobarcut3b) then
c                call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
c     &           eptemp,deptemp,virtemp)
c                  ntript=ntript+1
c                  ep3moli123=eptemp
c                  do l1 = 1,npole3b
c                    i = pnum(l1)
c                    dep3moli123(1,l1) = deptemp(1,l1)
c                    dep3moli123(2,l1) = deptemp(2,l1)
c                    dep3moli123(3,l1) = deptemp(3,l1)
c                  end do
c                  do i=1,3
c                    do j=1,3
c                       vir3moli123(j,i)=virtemp(j,i)
c                    end do
c                  end do
c
c                 if (shellsum.gt.rtapr3b) then
c                   term3=(shellsum-rtapr3b)/SwitchR3b
c                   tapr3b_2=term3*term3
c                   tapr3b_3=tapr3b_2*term3
c                   tapr3b_4=tapr3b_3*term3
c                   tapr3b_5=tapr3b_4*term3
c                   tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
c     &                - 6.0d0*tapr3b_5
c                   term4 = (-30.0d0*tapr3b_4 + 60.0d0*tapr3b_3 - 
c     &                     30.0d0*tapr3b_2)/SwitchR3b
c                   dtapr3b_ix = term4*(xr2/r2)
c                   dtapr3b_iy = term4*(yr2/r2)
c                   dtapr3b_iz = term4*(zr2/r2)
c                   dtapr3b_jx = term4*(xr3/r3)
c                   dtapr3b_jy = term4*(yr3/r3)
c                   dtapr3b_jz = term4*(zr3/r3)
c                   dtapr3b_kx = term4*(-xr2/r2-xr3/r3)
c                   dtapr3b_ky = term4*(-yr2/r2-yr3/r3)
c                   dtapr3b_kz = term4*(-zr2/r2-zr3/r3) 
c                   virtapr3b_xx=term4*((xr2/r2)*xr2+(xr3/r3)*xr3)
c                   virtapr3b_xy=term4*((xr2/r2)*yr2+(xr3/r3)*yr3)
c                   virtapr3b_xz=term4*((xr2/r2)*zr2+(xr3/r3)*zr3)
c                   virtapr3b_yy=term4*((yr2/r2)*yr2+(yr3/r3)*yr3)
c                   virtapr3b_yz=term4*((yr2/r2)*zr2+(yr3/r3)*zr3)
c                   virtapr3b_zz=term4*((zr2/r2)*zr2+(zr3/r3)*zr3)
c                 end if
c               end if
c            end if

c            if (shellsum .le. Cobarcut3b) then
C     LEFT OFF HERE! ACK! ACK!  

            if (shellsum.le.Cobarcut3b) then

            npole3b=0
            np1=0
            np2=0
            do i=1,sizeclust(clust1)
              moli1=clust(i,clust1)
c              npole3b=npole3b+3
              pnum(3*(i-1)+1)=imol(1,moli1)
              pnum(3*(i-1)+2)=imol(1,moli1)+1
              pnum(3*(i-1)+3)=imol(2,moli1)
              np1=np1+3
            end do
            do i=1,sizeclust(clust2)
              moli1=clust(i,clust2)
c              npole3b=npole3b+3
              pnum(np1+3*(i-1)+1)=imol(1,moli1)
              pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np1+3*(i-1)+3)=imol(2,moli1)
              np2=np2+3
            end do
            npole3b=np1+np2
            npole3b12=npole3b
                 if(r1_2.le.Rcut2b2) then
                   ep3moli123=ep3moli123-ep2moli12nosubtr
                   do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                      dep3moli123(j,l1)=dep3moli123(j,l1)
     &                               -dep2moli12nosubtr(j,l1)
                    end do
                   end do
                   do i=1,3
                    do j=1,3
                     vir3moli123(j,i)=vir3moli123(j,i)
     &                         -vir2moli12nosubtr(j,i)
                    end do
                   end do
                 else
                  call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &            eptemp,deptemp,virtemp)
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

                   ep3moli123 = ep3moli123 + ep1moli1
                   do l1 = 1, np1
                    i = pnum(l1)
                     do j = 1, 3
                    dep3moli123(j,l1)=dep3moli123(j,l1)+dep1moli1(j,l1)
                     end do
                   end do
                   do i=1,3
                    do j=1,3
                     vir3moli123(j,i)=vir3moli123(j,i)+vir1moli1(j,i)
                    end do
                   end do

                   ep3moli123 = ep3moli123 + ep1moli2
                   do l1 = 1, np2
                    i = pnum(l1)
                     l2=l1+np1
                     do j = 1, 3
                    dep3moli123(j,l2)=dep3moli123(j,l2)+dep1moli2(j,l1)
                     end do
                   end do
                   do i=1,3
                    do j=1,3
                     vir3moli123(j,i)=vir3moli123(j,i)+vir1moli2(j,i)
                    end do
                   end do

c               if (r3.le.Rcut2b) then
c               if (r3.le.1.0d12) then
            npole3b=0
            np2=0
            np3=0
            do i=1,sizeclust(clust2)
              moli1=clust(i,clust2)
c              npole3b=npole3b+3
              pnum(3*(i-1)+1)=imol(1,moli1)
              pnum(3*(i-1)+2)=imol(1,moli1)+1
              pnum(3*(i-1)+3)=imol(2,moli1)
              np2=np2+3
            end do
            do i=1,sizeclust(clust3)
              moli1=clust(i,clust3)
c              npole3b=npole3b+3
              pnum(np2+3*(i-1)+1)=imol(1,moli1)
              pnum(np2+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np2+3*(i-1)+3)=imol(2,moli1)
              np3=np3+3
            end do
            npole3b=np2+np3

                  call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &            eptemp,deptemp,virtemp)
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


            npole3b=0
            np1=0
            np3=0
            do i=1,sizeclust(clust1)
              moli1=clust(i,clust1)
c              npole3b=npole3b+3
              pnum(3*(i-1)+1)=imol(1,moli1)
              pnum(3*(i-1)+2)=imol(1,moli1)+1
              pnum(3*(i-1)+3)=imol(2,moli1)
              np1=np1+3
            end do
            do i=1,sizeclust(clust3)
              moli1=clust(i,clust3)
c              npole3b=npole3b+3
              pnum(np1+3*(i-1)+1)=imol(1,moli1)
              pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np1+3*(i-1)+3)=imol(2,moli1)
              np3=np3+3
            end do
            npole3b=np1+np3

                 call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp,deptemp,virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, np1
                      i = pnum(l1)
                      do j = 1, 3
                       dep3moli123(j,l1)=dep3moli123(j,l1)-deptemp(j,l1)
                      end do
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np2
                       do j = 1,3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                      end do
                     end do

            npole3b=0
            np3=0
            do i=1,sizeclust(clust3)
              moli1=clust(i,clust3)
c              npole3b=npole3b+3
              pnum(3*(i-1)+1)=imol(1,moli1)
              pnum(3*(i-1)+2)=imol(1,moli1)+1
              pnum(3*(i-1)+3)=imol(2,moli1)
              np3=np3+3
            end do
            npole3b=np3
                 call empole1a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp,deptemp,virtemp)
                     ep3moli123 = ep3moli123 +eptemp

                     do l1 = 1, npole3b
                      i = pnum(l1)
                      l2=l1+npole3b12
                      do j = 1, 3
                       dep3moli123(j,l2)=dep3moli123(j,l2)+deptemp(j,l1)
                      end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)+virtemp(j,i)
                      end do
                     end do



            npole3b=0
            np1=0
            np2=0
            np3=0
            do i=1,sizeclust(clust1)
              moli1=clust(i,clust1)
              npole3b=npole3b+3
              pnum(3*(i-1)+1)=imol(1,moli1)
              pnum(3*(i-1)+2)=imol(1,moli1)+1
              pnum(3*(i-1)+3)=imol(2,moli1)
              np1=np1+3
            end do
            do i=1,sizeclust(clust2)
              moli1=clust(i,clust2)
              npole3b=npole3b+3
              pnum(np1+3*(i-1)+1)=imol(1,moli1)
              pnum(np1+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np1+3*(i-1)+3)=imol(2,moli1)
              np2=np2+3
            end do
            do i=1,sizeclust(clust3)
              moli1=clust(i,clust3)
              npole3b=npole3b+3
              pnum(np1+np2+3*(i-1)+1)=imol(1,moli1)
              pnum(np1+np2+3*(i-1)+2)=imol(1,moli1)+1
              pnum(np1+np2+3*(i-1)+3)=imol(2,moli1)
              np3=np3+3
            end do

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
                  do l1 = np1+1,np1+np2
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
                  do l1 = np1+np2+1,npole3b
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

