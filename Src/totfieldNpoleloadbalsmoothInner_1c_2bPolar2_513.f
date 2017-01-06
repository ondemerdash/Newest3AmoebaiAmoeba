c
c      subroutine totfieldNpoleloadbalsmoothInner_1c_3bPolar(
c     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
      subroutine totfieldNpoleloadbalsmoothInner_1c_2bPolar2_513(
     &  ep3bt,virep3bt,dep3bt)
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use mpidat
      use limits
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
      integer moli1rmndr,l2,kk
      real*8 term,term2,term3,term4
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      integer, allocatable :: embedlst(:,:)
      integer, allocatable :: nembedlst(:)
      real*8 ewaldcut2
      allocate (nembedlst(9))
      allocate (embedlst(maxelst,9))

      ewaldcut2=ewaldcut*ewaldcut

      rtapr2b=1.0d12
      Rcut2b=1.0d12
c      rtapr2b=6.0d0
c      Rcut2b=8.0d0

      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b

c      rtapr3b=6.0d0
c      Cobarcut3b=8.0d0
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

        !print*,"TOTFIELDEWALD List"
        print*,"startlast",taskid,start_polar(taskid),last_polar(taskid)
     
c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start_polar,last_polar,
c!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
c!$OMP& Rcut2b,SwitchR2b,SwitchR3b,taskid,
c!$OMP& Rcut2b2)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
c!$OMP& schedule(guided)

      !do moli1=start,last
      !do kouter=start_polar(taskid),last_polar(taskid)
      !  moli1=molnew(kouter)
      do moli1=start_polar(taskid),last_polar(taskid)

          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=3


          do l1 = 1, npole3b
             i=pnum(l1)
             nembedlst(l1) = 0
             do k = 1,npole
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if ((r2 .le. ewaldcut2).and.(i.ne.k)) then
                nembedlst(l1) = nembedlst(l1)+1
                embedlst(nembedlst(l1),l1) = k
               end if
             end do
          end do          
           call empole1c_3b_Polar_totfield_513(embedlst,nembedlst,
     &          npole3b,pnum,ep1moli1,dep1moli1,vir1moli1)
c          print*,"moli1 ep1moli1=",moli1,ep1moli1
          ep3bt=ep3bt+ep1moli1
        !  do l1=1,npole3b
        !    i=pnum(l1)            
        !    dep3bt(1,i)=dep3bt(1,i)+dep1moli1(1,i)
        !    dep3bt(2,i)=dep3bt(2,i)+dep1moli1(2,i)
        !    dep3bt(3,i)=dep3bt(3,i)+dep1moli1(3,i)
        !    do kk=1,nembedlst(l1)
        !          k=embedlst(kk,l1)
        !       dep3bt(1,k)=dep3bt(1,k)+dep1moli1(1,k)
        !       dep3bt(2,k)=dep3bt(2,k)+dep1moli1(2,k)
        !       dep3bt(3,k)=dep3bt(3,k)+dep1moli1(3,k)
        !    end do
        !  end do

          do i=1,npole
            do j = 1, 3
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,i)
            end do
          end do
  
            do i=1,3
              do j=1,3
               virep3bt(j,i)=virep3bt(j,i)+vir1moli1(j,i)
              end do
            end do


        !do moli2=moli1+1,nmol
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
          xr1 = x1 - x2
          yr1 = y1 - y2
          zr1 = z1 - z2
          call imagen(xr1,yr1,zr1)
          r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1


          if(r1_2.le.Rcut2b2) then
             r1=sqrt(r1_2)

             do l1 = np1+1, npole3b
                i=pnum(l1)
                nembedlst(l1) = 0
                do k = 1,npole
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if ((r2 .le. ewaldcut2).and.(i.ne.k)) then
                   nembedlst(l1) = nembedlst(l1)+1
                   embedlst(nembedlst(l1),l1) = k
                  end if
                end do
             end do

 
             call empole1c_3b_Polar_totfield_513(embedlst,nembedlst,
     &             npole3b,pnum,ep2moli12,dep2moli12,vir2moli12)


             call empole1c_3b_Polar_totfield_513offset(np1,embedlst,
     &           nembedlst,npole3b,pnum,ep1moli2,dep1moli2,vir1moli2)
          !  pnum(1)=imol(1,moli2)
          !  pnum(2)=imol(1,moli2)+1
          !  pnum(3)=imol(2,moli2)
          !  npole3b=3
          !  call empole1c_3b_Polar_totfield(npole3b,pnum,ep1moli2,
     &    !         dep1moli2,vir1moli2)

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
                 !do i =1,npole
                 do l1 = 1,npole3b
                    i =pnum(l1)
                  dep3bt(1,i)=dep3bt(1,i) + tapr2b*(dep2moli12(1,i)
     &                       -dep1moli1(1,i)-dep1moli2(1,i))
                  dep3bt(2,i)=dep3bt(2,i) + tapr2b*(dep2moli12(2,i)
     &                       -dep1moli1(2,i)-dep1moli2(2,i))
                  dep3bt(3,i)=dep3bt(3,i) + tapr2b*(dep2moli12(3,i)
     &                       -dep1moli1(3,i)-dep1moli2(3,i))
                 end do
                 do l1 = 1,npole3b
                    do kk=1,nembedlst(l1)
                        k=embedlst(kk,l1)
                  dep3bt(1,k)=dep3bt(1,k) + tapr2b*(dep2moli12(1,k)
     &                       -dep1moli1(1,k)-dep1moli2(1,k))
                  dep3bt(2,k)=dep3bt(2,k) + tapr2b*(dep2moli12(2,k)
     &                       -dep1moli1(2,k)-dep1moli2(2,k))
                  dep3bt(3,k)=dep3bt(3,k) + tapr2b*(dep2moli12(3,k)
     &                       -dep1moli1(3,k)-dep1moli2(3,k)) 
                    end do
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
           !      do l1=1,npole3b
           !         i=pnum(l1)
           !       dep3bt(1,i)=dep3bt(1,i) + (dep2moli12(1,i)
     &     !                  -dep1moli1(1,i)-dep1moli2(1,i))
           !       dep3bt(2,i)=dep3bt(2,i) + (dep2moli12(2,i)
     &     !                  -dep1moli1(2,i)-dep1moli2(2,i))
           !       dep3bt(3,i)=dep3bt(3,i) + (dep2moli12(3,i)
     &     !                  -dep1moli1(3,i)-dep1moli2(3,i))
           !      end do
           !      do l1 =1,npole3b
           !         do kk=1,nembedlst(l1)
           !             k=embedlst(kk,l1)
           !       dep3bt(1,k)=dep3bt(1,k) + tapr2b*(dep2moli12(1,k)
     &     !                  -dep1moli1(1,k)-dep1moli2(1,k))
           !       dep3bt(2,k)=dep3bt(2,k) + tapr2b*(dep2moli12(2,k)
     &     !                  -dep1moli1(2,k)-dep1moli2(2,k))
           !       dep3bt(3,k)=dep3bt(3,k) + tapr2b*(dep2moli12(3,k)
     &     !                  -dep1moli1(3,k)-dep1moli2(3,k))            
           !         end do
           !      end do

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
          !else
              !ep2moli12=0.0d0
              !do i=1,npole
              !  do j = 1, 3
              !    dep2moli12(j,i)=0.0d0
              !  end do
              !end do
              !do i=1,3
              !  do j=1,3
              !   vir2moli12(j,i)=0.0d0
              !  end do
              !end do
          end if 
       

        end do 
      end do  
c!$OMP END DO
c!$OMP END PARALLEL

      !print*,"Ewald Innerloop= ep3bt=",ep3bt
c       print*,"End inner1 single",ep3bt
      deallocate(nembedlst)
      deallocate(embedlst)
      return
      end

