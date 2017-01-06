c
      subroutine UindLoadBalNew2SmoothInnerloop1_2cut_wskin(
     &  ep3bt,virep3bt,dep3bt,ntript,uindt) 
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use mpidat
      use cobar
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
      real*8 uind_local(3,30),uindt(3,npole)
      real*8 uind2moli12(3,30),uind3moli123(3,30)
      integer moli1rmndr,l2,kouter

      rtapr2b=5.0d0
      Rcut2b=6.0d0
      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b
c      rtapr3b=7.0d0
c      Cobarcut3b=8.5d0
      SwitchR3b=Cobarcut3b-rtapr3b

      ntript=0
      ep3bt=0.0d0
      do i = 1, npole
         do j = 1, 3
            dep3bt(j,i) = 0.0d0
            uindt(j,i)=0.0d0
         end do
      end do
      do i=1,3
         do j=1,3
            virep3bt(j,i)=0.0d0
         end do
      end do

c      print*,"LOAD-BALANCED Correct NewSmooth Cobar taskid",taskid

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start_polar,last_polar,nmollst2,mollst2,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
!$OMP& Rcut2b2,uindt)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript,uindt)
!$OMP& schedule(guided)
      do kouter=start_polar(taskid),last_polar(taskid)
        moli1=molnew(kouter)
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
          xr = x1 - x2
          yr = y1 - y2
          zr = z1 - z2
          !call image(xr,yr,zr)
          xr1=xr
          yr1=yr
          zr1=zr
          r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
C
C   Pair cutoffs
C
          if(r1_2.le.Rcut2b2) then
          r1=sqrt(r1_2)
          call empole1a_3b_Polar_uind(npole3b,pnum,eptemp,deptemp,
     &          virtemp,uind_local)
          ep2moli12=eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep2moli12(j,l1)=deptemp(j,l1)
              uind2moli12(j,l1)=uind_local(j,l1)
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
              dep2moli12(1,l1)=dep2moli12(1,l1)*tapr2b
     &                       + ep2moli12*dtapr2b_x(l1)
              dep2moli12(2,l1)=dep2moli12(2,l1)*tapr2b
     &                       + ep2moli12*dtapr2b_y(l1)
              dep2moli12(3,l1)=dep2moli12(3,l1)*tapr2b
     &                       + ep2moli12*dtapr2b_z(l1)
            end do
           
! NEED TO SMOOTHEN VIRIAL !
            do i=1,3
              do j=1,3
                vir2moli12(j,i)=vir2moli12(j,i)*tapr2b
              end do
            end do
            do l1=1,npole3b
              i = pnum(l1)
              vir2moli12(1,1)=vir2moli12(1,1)
     &                          +dtapr2b_x(l1)*ep2moli12*x(i)
              vir2moli12(2,1)=vir2moli12(2,1)
     &                          +dtapr2b_x(l1)*ep2moli12*y(i)
              vir2moli12(3,1)=vir2moli12(3,1)
     &                           +dtapr2b_x(l1)*ep2moli12*z(i)
              vir2moli12(1,2)=vir2moli12(1,2)
     &                           +dtapr2b_y(l1)*ep2moli12*x(i)
              vir2moli12(2,2)=vir2moli12(2,2)
     &                           +dtapr2b_y(l1)*ep2moli12*y(i)
              vir2moli12(3,2)=vir2moli12(3,2)
     &                           +dtapr2b_y(l1)*ep2moli12*z(i)
              vir2moli12(1,3)=vir2moli12(1,3)
     &                           +dtapr2b_z(l1)*ep2moli12*x(i)
              vir2moli12(2,3)=vir2moli12(2,3)
     &                           +dtapr2b_z(l1)*ep2moli12*y(i)
              vir2moli12(3,3)=vir2moli12(3,3)
     &                           +dtapr2b_z(l1)*ep2moli12*z(i)
            end do
            ep2moli12=tapr2b*ep2moli12
           end if  
              ep3bt=ep3bt+ep2moli12
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                  uindt(j,i) = uindt(j,i)+uind2moli12(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                  virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do
          else
            ep2moli12=0.0d0
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep2moli12(j,l1)= 0.0d0
                uind2moli12(j,l1)= 0.0d0
              end do
            end do
            do i=1,3
              do j=1,3
                vir2moli12(j,i)= 0.0d0
              end do
            end do
          end if

          do k2=1,nmollst2(moli2)
            moli3=mollst2(k2,moli2)
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

            xr = x1 - x3
            yr = y1 - y3
            zr = z1 - z3
            !call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr

            xr = x2 - x3
            yr = y2 - y3
            zr = z2 - z3
            !call image(xr,yr,zr)

            xr3=xr
            yr3=yr
            zr3=zr
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
            r1=sqrt(r1_2)         
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)
            
            if ((r1.lt.r3).and.(r2.lt.r3) ) then 
              shellsum=r1+r2
              if (shellsum.le.Cobarcut3b) then
                call empole1a_3b_Polar_uind(npole3b,pnum,eptemp,deptemp,
     &           virtemp,uind_local)
               


                ntript=ntript+1
                ep3moli123=eptemp
                do l1 = 1,npole3b
                  i = pnum(l1)
                  uind3moli123(1,l1)=uind_local(1,l1)
                  uind3moli123(2,l1)=uind_local(2,l1)
                  uind3moli123(3,l1)=uind_local(3,l1)
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
                end if
              end if             
            else if ((r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
              if (shellsum.le.Cobarcut3b) then
                call empole1a_3b_Polar_uind(npole3b,pnum,eptemp,deptemp,
     &           virtemp,uind_local)
                ntript=ntript+1


                ep3moli123=eptemp
                do l1 = 1,npole3b
                    i = pnum(l1)
                    uind3moli123(1,l1)=uind_local(1,l1)
                    uind3moli123(2,l1)=uind_local(2,l1)
                    uind3moli123(3,l1)=uind_local(3,l1)
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
                 end if
              end if
            else if ((r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
               if (shellsum.le.Cobarcut3b) then
                call empole1a_3b_Polar_uind(npole3b,pnum,eptemp,deptemp,
     &           virtemp,uind_local)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    uind3moli123(1,l1)=uind_local(1,l1)
                    uind3moli123(2,l1)=uind_local(2,l1)
                    uind3moli123(3,l1)=uind_local(3,l1)
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
                 end if
               end if
            end if

            if (shellsum .le. Cobarcut3b) then

               npole3b=6
               ep3moli123=ep3moli123-ep2moli12
                 
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3moli123(j,l1)=dep3moli123(j,l1)-dep2moli12(j,l1)
                    uind3moli123(j,l1)=uind3moli123(j,l1)
     &                                   -uind2moli12(j,l1)
                  end do
               end do
               do i=1,3
                 do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)-vir2moli12(j,i)
                 end do
               end do

               if (r3.le.Rcut2b) then
                  npole3b=6
                  np1 =3
                  pnum(1)=imol(1,moli2)
                  pnum(2)=imol(1,moli2)+1
                  pnum(3)=imol(2,moli2)
                  pnum(4)=imol(1,moli3)
                  pnum(5)=imol(1,moli3)+1
                  pnum(6)=imol(2,moli3)              
                  call empole1a_3b_Polar_uind(npole3b,pnum,eptemp,
     &            deptemp,virtemp,uind_local)
                  if (r3.gt.rtapr2b) then               
                    term3=(r3-rtapr2b)/SwitchR2b
                    tapr2b_2=term3*term3
                    tapr2b_3=tapr2b_2*term3
                    tapr2b_4=tapr2b_3*term3
                    tapr2b_5=tapr2b_4*term3
                    tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                 - 6.0d0*tapr2b_5

                     do l1=1,np1
                      i=pnum(l1)
                      if (name(i).eq.'O') then
                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
     &                            30.0d0*tapr2b_2)/SwitchR2b/r3
                        dtapr2b_x(l1) = term4*xr3
                        dtapr2b_y(l1) = term4*yr3
                        dtapr2b_z(l1) = term4*zr3
                      else
                        dtapr2b_x(l1) =0.0d0
                        dtapr2b_y(l1) =0.0d0
                        dtapr2b_z(l1) =0.0d0
                      end if
                     end do
                     do l1=np1+1,npole3b
                       i=pnum(l1)
                      if (name(i).eq.'O') then
                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
     &                            30.0d0*tapr2b_2)/SwitchR2b/r3
                        dtapr2b_x(l1)= -term4*xr3
                        dtapr2b_y(l1)= -term4*yr3
                        dtapr2b_z(l1)= -term4*zr3
                      else
                        dtapr2b_x(l1) =0.0d0
                        dtapr2b_y(l1) =0.0d0
                        dtapr2b_z(l1) =0.0d0
                      end if
                     end do
                     ep3moli123=ep3moli123-tapr2b*eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       uind3moli123(1,l2)=uind3moli123(1,l2)
     &                  -uind_local(1,l1)
                       uind3moli123(2,l2)=uind3moli123(2,l2)
     &                  -uind_local(2,l1)
                       uind3moli123(3,l2)=uind3moli123(3,l2)
     &                  -uind_local(3,l1)
                       dep3moli123(1,l2)=dep3moli123(1,l2)
     &                  -deptemp(1,l1)*tapr2b- eptemp*dtapr2b_x(l1)
                       dep3moli123(2,l2)=dep3moli123(2,l2)
     &                  -deptemp(2,l1)*tapr2b- eptemp*dtapr2b_y(l1)
                       dep3moli123(3,l2)=dep3moli123(3,l2)
     &                  -deptemp(3,l1)*tapr2b- eptemp*dtapr2b_z(l1)
                     end do
                     do i=1,3
                       do j=1,3
                        virtemp(j,i)=virtemp(j,i)*tapr2b
                       end do
                     end do
                     do l1=1,npole3b
                       i = pnum(l1)
                     virtemp(1,1)=virtemp(1,1)+dtapr2b_x(l1)*eptemp*x(i)
                     virtemp(2,1)=virtemp(2,1)+dtapr2b_x(l1)*eptemp*y(i)
                     virtemp(3,1)=virtemp(3,1)+dtapr2b_x(l1)*eptemp*z(i)
                     virtemp(1,2)=virtemp(1,2)+dtapr2b_y(l1)*eptemp*x(i)
                     virtemp(2,2)=virtemp(2,2)+dtapr2b_y(l1)*eptemp*y(i)
                     virtemp(3,2)=virtemp(3,2)+dtapr2b_y(l1)*eptemp*z(i)
                     virtemp(1,3)=virtemp(1,3)+dtapr2b_z(l1)*eptemp*x(i)
                     virtemp(2,3)=virtemp(2,3)+dtapr2b_z(l1)*eptemp*y(i)
                     virtemp(3,3)=virtemp(3,3)+dtapr2b_z(l1)*eptemp*z(i)
                     end do
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       end do
                     end do
                  else
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       uind3moli123(j,l2)=uind3moli123(j,l2)
     &                  -uind_local(j,l1)
                       end do
                     end do
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       end do
                     end do
                  end if
               end if

               if(r2.le.Rcut2b) then
                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole1a_3b_Polar_uind(npole3b,pnum,eptemp,
     &           deptemp,virtemp,uind_local)
                 if (r2.gt.rtapr2b) then
                   term3=(r2-rtapr2b)/SwitchR2b
                   tapr2b_2=term3*term3
                   tapr2b_3=tapr2b_2*term3
                   tapr2b_4=tapr2b_3*term3
                   tapr2b_5=tapr2b_4*term3
                   tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5

                     do l1=1,np1
                      i=pnum(l1)
                      if (name(i).eq.'O') then
                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
     &                            30.0d0*tapr2b_2)/SwitchR2b/r2
                        dtapr2b_x(l1) = term4*xr2
                        dtapr2b_y(l1) = term4*yr2
                        dtapr2b_z(l1) = term4*zr2
                      else
                        dtapr2b_x(l1) = 0.0d0
                        dtapr2b_y(l1) = 0.0d0
                        dtapr2b_z(l1) = 0.0d0
                      end if
                     end do
                     do l1=np1+1,npole3b
                       i=pnum(l1)
                      if (name(i).eq.'O') then
                        term4 = (-30.0d0*tapr2b_4 + 60.0d0*tapr2b_3 - 
     &                            30.0d0*tapr2b_2)/SwitchR2b/r2
                        dtapr2b_x(l1) = -term4*xr2
                        dtapr2b_y(l1) = -term4*yr2
                        dtapr2b_z(l1) = -term4*zr2
                      else
                       dtapr2b_x(l1)= 0.0d0
                       dtapr2b_y(l1)= 0.0d0
                       dtapr2b_z(l1)= 0.0d0
                      end if
                     end do
                     ep3moli123 = ep3moli123 -tapr2b*eptemp
                     do l1 = 1, np1
                       i = pnum(l1)
                       dep3moli123(1,l1)=dep3moli123(1,l1)
     &                   -deptemp(1,l1)*tapr2b- eptemp*dtapr2b_x(l1)
                       dep3moli123(2,l1)=dep3moli123(2,l1)
     &                   -deptemp(2,l1)*tapr2b- eptemp*dtapr2b_y(l1)
                       dep3moli123(3,l1)=dep3moli123(3,l1)
     &                   -deptemp(3,l1)*tapr2b- eptemp*dtapr2b_z(l1)
                       uind3moli123(1,l1)=uind3moli123(1,l1)
     &                   -uind_local(1,l1)
                       uind3moli123(2,l1)=uind3moli123(2,l1)
     &                   -uind_local(2,l1)
                       uind3moli123(3,l1)=uind3moli123(3,l1)
     &                   -uind_local(3,l1)
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       dep3moli123(1,l2)=dep3moli123(1,l2)
     &                  -deptemp(1,l1)*tapr2b-eptemp*dtapr2b_x(l1)
                       dep3moli123(2,l2)=dep3moli123(2,l2)
     &                  -deptemp(2,l1)*tapr2b-eptemp*dtapr2b_y(l1)
                       dep3moli123(3,l2)=dep3moli123(3,l2)
     &                  -deptemp(3,l1)*tapr2b-eptemp*dtapr2b_z(l1)
                       uind3moli123(1,l2)=uind3moli123(1,l2)
     &                   -uind_local(1,l1)
                       uind3moli123(2,l2)=uind3moli123(2,l2)
     &                   -uind_local(2,l1)
                       uind3moli123(3,l2)=uind3moli123(3,l2)
     &                   -uind_local(3,l1)
                     end do
                     do i=1,3
                       do j=1,3
                        virtemp(j,i)=virtemp(j,i)*tapr2b
                       end do
                     end do
                     do l1 = 1, npole3b
                      i = pnum(l1)
                     virtemp(1,1)=virtemp(1,1)+dtapr2b_x(l1)*eptemp*x(i)
                     virtemp(2,1)=virtemp(2,1)+dtapr2b_x(l1)*eptemp*y(i)
                     virtemp(3,1)=virtemp(3,1)+dtapr2b_x(l1)*eptemp*z(i)
                     virtemp(1,2)=virtemp(1,2)+dtapr2b_y(l1)*eptemp*x(i)
                     virtemp(2,2)=virtemp(2,2)+dtapr2b_y(l1)*eptemp*y(i)
                     virtemp(3,2)=virtemp(3,2)+dtapr2b_y(l1)*eptemp*z(i)
                     virtemp(1,3)=virtemp(1,3)+dtapr2b_z(l1)*eptemp*x(i)
                     virtemp(2,3)=virtemp(2,3)+dtapr2b_z(l1)*eptemp*y(i)
                     virtemp(3,3)=virtemp(3,3)+dtapr2b_z(l1)*eptemp*z(i)
                     end do
                     do i=1,3
                       do j=1,3
                        vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       end do
                     end do                  
                  else
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, np1
                      i = pnum(l1)
                      do j = 1, 3
                       dep3moli123(j,l1)=dep3moli123(j,l1)-deptemp(j,l1)
                       uind3moli123(j,l1)=uind3moli123(j,l1)
     &                   -uind_local(j,l1)
                      end do
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1,3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       uind3moli123(j,l2)=uind3moli123(j,l2)
     &                   -uind_local(j,l1)
                       end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                      end do
                     end do
                  end if
               end if

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
                    uindt(1,i)=uindt(1,i)+uind3moli123(1,l1)
                    uindt(2,i)=uindt(2,i)+uind3moli123(2,l1)
                    uindt(3,i)=uindt(3,i)+uind3moli123(3,l1)
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
                    uindt(1,i)=uindt(1,i)+uind3moli123(1,l1)
                    uindt(2,i)=uindt(2,i)+uind3moli123(2,l1)
                    uindt(3,i)=uindt(3,i)+uind3moli123(3,l1)
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
                    uindt(1,i)=uindt(1,i)+uind3moli123(1,l1)
                    uindt(2,i)=uindt(2,i)+uind3moli123(2,l1)
                    uindt(3,i)=uindt(3,i)+uind3moli123(3,l1)

                  end do 
!  NEED TO CORRECT VIRIAL W/ TAPERING FUNC
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)*tapr3b
                    end do
                  end do
                  do l1=1,np1
                     i = pnum(l1)
                    if(name(i).eq.'O') then
                    vir3moli123(1,1)=vir3moli123(1,1)
     &                 +dtapr3b_ix*ep3moli123*x(i)
                    vir3moli123(2,1)=vir3moli123(2,1)
     &                 +dtapr3b_ix*ep3moli123*y(i)
                    vir3moli123(3,1)=vir3moli123(3,1)
     &                 +dtapr3b_ix*ep3moli123*z(i)
                    vir3moli123(1,2)=vir3moli123(1,2)
     &                 +dtapr3b_iy*ep3moli123*x(i)
                    vir3moli123(2,2)=vir3moli123(2,2)
     &                 +dtapr3b_iy*ep3moli123*y(i)
                    vir3moli123(3,2)=vir3moli123(3,2)
     &                 +dtapr3b_iy*ep3moli123*z(i)
                    vir3moli123(1,3)=vir3moli123(1,3)
     &                 +dtapr3b_iz*ep3moli123*x(i)
                    vir3moli123(2,3)=vir3moli123(2,3)
     &                 +dtapr3b_iz*ep3moli123*y(i)
                    vir3moli123(3,3)=vir3moli123(3,3)
     &                 +dtapr3b_iz*ep3moli123*z(i)
                    end if
                  end do
                  do l1=np1+1,np2
                     i = pnum(l1)
                    if(name(i).eq.'O') then
                    vir3moli123(1,1)=vir3moli123(1,1)
     &                   +dtapr3b_jx*ep3moli123*x(i)
                    vir3moli123(2,1)=vir3moli123(2,1)
     &                   +dtapr3b_jx*ep3moli123*y(i)
                    vir3moli123(3,1)=vir3moli123(3,1)
     &                   +dtapr3b_jx*ep3moli123*z(i)
                    vir3moli123(1,2)=vir3moli123(1,2)
     &                   +dtapr3b_jy*ep3moli123*x(i)
                    vir3moli123(2,2)=vir3moli123(2,2)
     &                   +dtapr3b_jy*ep3moli123*y(i)
                    vir3moli123(3,2)=vir3moli123(3,2)
     &                   +dtapr3b_jy*ep3moli123*z(i)
                    vir3moli123(1,3)=vir3moli123(1,3)
     &                   +dtapr3b_jz*ep3moli123*x(i)
                    vir3moli123(2,3)=vir3moli123(2,3)
     &                   +dtapr3b_jz*ep3moli123*y(i)
                    vir3moli123(3,3)=vir3moli123(3,3)
     &                   +dtapr3b_jz*ep3moli123*z(i)
                    end if
                  end do
                  do l1=np2+1,npole3b
                     i = pnum(l1)
                    if(name(i).eq.'O') then
                    vir3moli123(1,1)=vir3moli123(1,1)
     &                  +dtapr3b_kx*ep3moli123*x(i)
                    vir3moli123(2,1)=vir3moli123(2,1)
     &                  +dtapr3b_kx*ep3moli123*y(i)
                    vir3moli123(3,1)=vir3moli123(3,1)
     &                  +dtapr3b_kx*ep3moli123*z(i)
                    vir3moli123(1,2)=vir3moli123(1,2)
     &                  +dtapr3b_ky*ep3moli123*x(i)
                    vir3moli123(2,2)=vir3moli123(2,2)
     &                  +dtapr3b_ky*ep3moli123*y(i)
                    vir3moli123(3,2)=vir3moli123(3,2)
     &                  +dtapr3b_ky*ep3moli123*z(i)
                    vir3moli123(1,3)=vir3moli123(1,3)
     &                  +dtapr3b_kz*ep3moli123*x(i)
                    vir3moli123(2,3)=vir3moli123(2,3)
     &                  +dtapr3b_kz*ep3moli123*y(i)
                    vir3moli123(3,3)=vir3moli123(3,3)
     &                  +dtapr3b_kz*ep3moli123*z(i)
                    end if
                  end do
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
                      uindt(j,i) = uindt(j,i) + uind3moli123(j,l1)
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
