c
      subroutine All2bAll3bNew2SmoothInnerloop1_2cut_wskin_fij(
     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr,depmat3bt)
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use fmat
c      use cell
      use iounit
      implicit none
      real*8  ep2moli12
      real*8  dep2moli12(3,30)
      integer i,ii,j,l1,i1,i2,i3,k,l3,l4
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
      real*8 deptempmat(3,npole,npole),depmat2moli12(3,npole,npole)
      real*8 depmat3moli123(3,npole,npole)
      real*8 depmat3bt(3,npole,npole),r
      real*8 x2O,y2O,z2O,x2H1,y2H1,z2H1,x2H2,y2H2,z2H2
      real*8 x1O,y1O,z1O,x1H1,y1H1,z1H1,x1H2,y1H2,z1H2
      real*8 hbond1,hbond2,hbond3,hbond4,rlowhbond,xdonO,ydonO,zdonO
      real*8 xdonH,ydonH,zdonH,xacc,yacc,zacc,rdonOdonH2,rdonOdonH
      real*8 rdonOacc2,hbangle
      integer hcount,molnum(30)
      real*8 xnodonH,ynodonH,znodonH,xHacc,yHacc,zHacc,norm
      real*8 normal1_x,normal1_y,normal1_z
      real*8 normal2_x,normal2_y,normal2_z
      real*8 normal3_x,normal3_y,normal3_z,phi1,chi
      real*8 normal4_x,normal4_y,normal4_z
      real*8 xHacc2,yHacc2,zHacc2,theta2,phi2
      real*8 xbisect,ybisect,zbisect
      real*8 x3O,y3O,z3O,x3H1,y3H1,z3H1,x3H2,y3H2,z3H2
      real*8 rlowhbond_2,xdonO_2,ydonO_2,zdonO_2
      real*8 xdonH_2,ydonH_2,zdonH_2,xacc_2,yacc_2,zacc_2
      real*8 hbangle_2
      real*8 xnodonH_2,ynodonH_2,znodonH_2,xHacc_2,yHacc_2,zHacc_2
      real*8 xHacc2_2,yHacc2_2,zHacc2_2
      real*8 rlowhbond_3,xdonO_3,ydonO_3,zdonO_3
      real*8 xdonH_3,ydonH_3,zdonH_3,xacc_3,yacc_3,zacc_3
      real*8 hbangle_3
      real*8 xnodonH_3,ynodonH_3,znodonH_3,xHacc_3,yHacc_3,zHacc_3
      real*8 xHacc2_3,yHacc2_3,zHacc2_3
      logical printforcematbody3
      real*8 phi1_2,chi_2,theta2_2,phi2_2,phi1_3,chi_3,theta2_3,phi2_3

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
                   do k = 1, npole
                      do j = 1, 3
                        depmat3bt(j,k,i)=0.0d0
                      end do
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(j,i)=0.0d0
                   end do
                end do

      !printforcematbody=.false.
c       print*,"printforcematbody=",printforcematbody
c        print*,"Fij NOCut3b NoCut2b"
c        print*,"start last",start,last,moli1rmndr
         printforcematbody3=.true.

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start,last,x,y,z,name,
!$OMP& nmol,depmat3bt,printforcematbody,iout,
!$OMP& printforcematbody3)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,depmat3bt)
!$OMP& schedule(guided)
      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
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

            hcount=0
            do l1 = 1, np1
              i = pnum(l1)
              if(name(i).eq.'O') then
                x1O = x(i)
                y1O = y(i)
                z1O = z(i)
              else if((name(i).eq.'H').and.(hcount.eq.0)) then
                x1H1 = x(i)
                y1H1 = y(i)
                z1H1 = z(i)
                hcount=hcount+1
              else if((name(i).eq.'H').and.(hcount.eq.1)) then
                x1H2 = x(i)
                y1H2 = y(i)
                z1H2 = z(i)
                hcount=hcount+1
              end if
            end do

            hcount=0
            do l1 = np1+1, np2
              i = pnum(l1)
              if(name(i).eq.'O') then
                x2O = x(i)
                y2O = y(i)
                z2O = z(i)
              else if((name(i).eq.'H').and.(hcount.eq.0)) then
                x2H1 = x(i)
                y2H1 = y(i)
                z2H1 = z(i)
                hcount=hcount+1
              else if((name(i).eq.'H').and.(hcount.eq.1)) then
                x2H2 = x(i)
                y2H2 = y(i)
                z2H2 = z(i)
                hcount=hcount+1
              end if
            end do

          xr = x1O - x2H1
          yr = y1O - y2H1
          zr = z1O - z2H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond1= r
            rlowhbond = hbond1
            xdonO= x2O
            ydonO= y2O
            zdonO= z2O
            xdonH= x2H1
            ydonH= y2H1
            zdonH= z2H1
            xnodonH = x2H2
            ynodonH = y2H2
            znodonH = z2H2

            xacc = x1O
            yacc = y1O
            zacc = z1O
            xHacc = x1H1
            yHacc = y1H1
            zHacc = z1H1

            xHacc2 = x1H2
            yHacc2 = y1H2
            zHacc2 = z1H2

 
          xr = x1O - x2H2
          yr = y1O - y2H2
          zr = z1O - z2H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond2= r
            if(hbond2.lt.rlowhbond) then
              rlowhbond = hbond2
              xdonO= x2O
              ydonO= y2O
              zdonO= z2O
              xdonH= x2H2
              ydonH= y2H2
              zdonH= z2H2
              xnodonH=x2H1
              ynodonH=y2H1
              znodonH=z2H1

              xacc = x1O
              yacc = y1O
              zacc = z1O
              xHacc = x1H1
              yHacc = y1H1
              zHacc = z1H1
              
            xHacc2 = x1H2
            yHacc2 = y1H2
            zHacc2 = z1H2       
            end if

          xr = x2O - x1H1
          yr = y2O - y1H1
          zr = z2O - z1H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond3= r
            if(hbond3.lt.rlowhbond) then
              rlowhbond = hbond3
              xdonO=x1O
              ydonO=y1O
              zdonO=z1O
              xdonH=x1H1
              ydonH=y1H1
              zdonH=z1H1
              xnodonH=x1H2
              ynodonH=y1H2
              znodonH=z1H2              

              xacc=x2O
              yacc=y2O
              zacc=z2O
              xHacc=x2H1
              yHacc=y2H1
              zHacc=z2H1

              xHacc2=x2H2
              yHacc2=y2H2
              zHacc2=z2H2 
            end if

          xr = x2O - x1H2
          yr = y2O - y1H2
          zr = z2O - z1H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond4= r
            if(hbond4.lt.rlowhbond) then
              rlowhbond = hbond4
              xdonO=x1O
              ydonO=y1O
              zdonO=z1O
              xdonH=x1H2
              ydonH=y1H2
              zdonH=z1H2
              xnodonH=x1H1
              ynodonH=y1H1
              znodonH=z1H1

              xacc=x2O
              yacc=y2O
              zacc=z2O
              xHacc=x2H1
              yHacc=y2H1
              zHacc=z2H1

              xHacc2=x2H2
              yHacc2=y2H2
              zHacc2=z2H2
            end if

            xr = xdonO-xdonH
            yr = ydonO-ydonH
            zr = zdonO-zdonH
            r2 = xr*xr + yr*yr + zr*zr
            rdonOdonH2=r2
            rdonOdonH=sqrt(rdonOdonH2)

            xr = xdonO-xacc
            yr = ydonO-yacc
            zr = zdonO-zacc
            r2 = xr*xr + yr*yr + zr*zr
            rdonOacc2=r2

c            rdonOacc2 = rdonOdonH2 + rlowhbond*rlowhbond - 2*rlowhbond*rdonOdonH*cos(angle)
c            2*rlowhbond*rdonOdonH*cos(hbangle)=rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2

c            cos(hbangle)= (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)/(2*rlowhbond*rdonOdonH)
            hbangle=acos( (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)
     &                  /(2*rlowhbond*rdonOdonH))

c FIRST CALCULATE PLANE OF HYDROGEN BOND

         xr1= xacc-xdonH
         yr1= yacc-ydonH
         zr1= zacc-zdonH

         xr2= xdonO-xdonH
         yr2= ydonO-ydonH
         zr2= zdonO-zdonH

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) + 
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal1_x=(yr1*zr2-yr2*zr1)/norm
         normal1_y=(xr2*zr1-xr1*zr2)/norm
         normal1_z=(xr1*yr2-xr2*yr1)/norm

c NEXT CALCULATE PLANE FORMED BY DONATED HYDROGEN, DONOR OXYGEN, AND NONDONATED HYDROGEN

         xr1= xdonH-xdonO
         yr1= ydonH-ydonO
         zr1= zdonH-zdonO

         xr2= xnodonH-xdonO
         yr2= ynodonH-ydonO
         zr2= znodonH-zdonO

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) + 
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal2_x=(yr1*zr2-yr2*zr1)/norm
         normal2_y=(xr2*zr1-xr1*zr2)/norm
         normal2_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY HYDROGEN COVALENTLY BOUND TO ACCEPTOR OXYGEN, THE ACCEPTOR OXYGEN, AND THE DONATED HYDROGEN

         xr1 = xHacc-xacc
         yr1 = yHacc-yacc
         zr1 = zHacc-zacc

         xr2 = xdonH-xacc
         yr2 = ydonH-yacc
         zr2 = zdonH-zacc

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal3_x=(yr1*zr2-yr2*zr1)/norm
         normal3_y=(xr2*zr1-xr1*zr2)/norm
         normal3_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY ATOMS ON ACCEPTOR

         xr1 = xHacc-xacc
         yr1 = yHacc-yacc
         zr1 = zHacc-zacc

         xr2 = xHacc2-xacc
         yr2 = yHacc2-yacc
         zr2 = zHacc2-zacc

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal4_x=(yr1*zr2-yr2*zr1)/norm
         normal4_y=(xr2*zr1-xr1*zr2)/norm
         normal4_z=(xr1*yr2-xr2*yr1)/norm

         phi1= acos(normal1_x*normal2_x + normal1_y*normal2_y 
     &          + normal1_z*normal2_z)

         chi= acos(normal1_x*normal3_x + normal1_y*normal3_y 
     &          + normal1_z*normal3_z)

         theta2=acos(normal1_x*normal4_x + normal1_y*normal4_y
     &          + normal1_z*normal4_z)

         phi2= acos(normal2_x*normal4_x + normal2_y*normal4_y
     &          + normal2_z*normal4_z)

         xbisect=xnodonH+(xnodonH-xdonH)/2.0d0
         ybisect=ynodonH+(ynodonH-ydonH)/2.0d0
         zbisect=znodonH+(znodonH-zdonH)/2.0d0


         call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,pnum,
     &           eptemp,deptemp,virtemp,deptempmat)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,l1)=deptemp(j,l1)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                end do
                do l3=1,npole3b
                   k=pnum(l3)
                   do j = 1, 3
                      depmat2moli12(j,k,i)=deptempmat(j,k,i)
                      depmat3bt(j,k,i)=depmat3bt(j,k,i)
     &                                 +depmat2moli12(j,k,i)
                   end do
                end do
              end do

          if(printforcematbody) then
            molnum(1)=moli1
            molnum(2)=moli1
            molnum(3)=moli1
            molnum(4)=moli2
            molnum(5)=moli2
            molnum(6)=moli2
            do l1=1,npole3b
               i=pnum(l1)
               do l3=1,npole3b
                  k=pnum(l3)
                  if(i.lt.k) then
               !   if(molnum(l1).ne.molnum(l3)) then

                      xr = x(k)-x(i)
                      yr = y(k)-y(i)
                      zr = z(k)-z(i)
                      r=sqrt(xr*xr + yr*yr + zr*zr)
c                write(iout,31) moli1,moli2,i,k,r,xr,depmat2moli12(1,k,i)
c                write(iout,21) moli1,moli2,i,k,r,yr,depmat2moli12(2,k,i)
c                write(iout,11) moli1,moli2,i,k,r,zr,depmat2moli12(3,k,i)
c   31 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5)
c   21 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5)
c   11 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5)

                write(iout,31) moli1,moli2,i,k,r,xr,depmat2moli12(1,k,i)
     &               ,rlowhbond,hbangle,phi1,chi,theta2,phi2,ep2moli12
                write(iout,21) moli1,moli2,i,k,r,yr,depmat2moli12(2,k,i)
     &               ,rlowhbond,hbangle,phi1,chi,theta2,phi2,ep2moli12
                write(iout,11) moli1,moli2,i,k,r,zr,depmat2moli12(3,k,i)
     &               ,rlowhbond,hbangle,phi1,chi,theta2,phi2,ep2moli12
   31 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5,1x,"hbr=",f7.3,1x,
     &       "hbang=",f7.3,1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,1x,"ep2=",f7.3)
   21 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5,1x,"hbr=",f7.3,1x,
     &       "hbang=",f7.3,1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,1x,"ep2=",f7.3)
   11 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5,1x,"hbr=",f7.3,1x,
     &       "hbang=",f7.3,1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,1x,"ep2=",f7.3)
                  end if
               end do
            end do
          end if


              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do


          do moli3=moli2+1,nmol
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

            hcount=0
            do l1 = np2+1, np3
              i = pnum(l1)
              if(name(i).eq.'O') then
                x3O = x(i)
                y3O = y(i)
                z3O = z(i)
              else if((name(i).eq.'H').and.(hcount.eq.0)) then
                x3H1 = x(i)
                y3H1 = y(i)
                z3H1 = z(i)
                hcount=hcount+1
              else if((name(i).eq.'H').and.(hcount.eq.1)) then
                x3H2 = x(i)
                y3H2 = y(i)
                z3H2 = z(i)
                hcount=hcount+1
              end if
            end do

           xr = x1O - x3H1
           yr = y1O - y3H1
           zr = z1O - z3H1

c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond1= r
            rlowhbond_2 = hbond1
            xdonO_2= x3O
            ydonO_2= y3O
            zdonO_2= z3O
            xdonH_2= x3H1
            ydonH_2= y3H1
            zdonH_2= z3H1
            xnodonH_2 = x3H2
            ynodonH_2 = y3H2
            znodonH_2 = z3H2

            xacc_2 = x1O
            yacc_2 = y1O
            zacc_2 = z1O
            xHacc_2 = x1H1
            yHacc_2 = y1H1
            zHacc_2 = z1H1

            xHacc2_2 = x1H2
            yHacc2_2 = y1H2
            zHacc2_2 = z1H2


            xr = x1O - x3H2
            yr = y1O - y3H2
            zr = z1O - z3H2

            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond2= r
            if(hbond2.lt.rlowhbond_2) then
              rlowhbond_2 = hbond2
            xdonO_2= x3O
            ydonO_2= y3O
            zdonO_2= z3O
            xdonH_2= x3H2
            ydonH_2= y3H2
            zdonH_2= z3H2
            xnodonH_2 = x3H1
            ynodonH_2 = y3H1
            znodonH_2 = z3H1

            xacc_2 = x1O
            yacc_2 = y1O
            zacc_2 = z1O
            xHacc_2 = x1H1
            yHacc_2 = y1H1
            zHacc_2 = z1H1

            xHacc2_2 = x1H2
            yHacc2_2 = y1H2
            zHacc2_2 = z1H2
            end if              

          xr = x3O - x1H1
          yr = y3O - y1H1
          zr = z3O - z1H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond3= r
            if(hbond3.lt.rlowhbond_2) then
              rlowhbond_2 = hbond3
              xdonO_2=x1O
              ydonO_2=y1O
              zdonO_2=z1O
              xdonH_2=x1H1
              ydonH_2=y1H1
              zdonH_2=z1H1
              xnodonH_2=x1H2
              ynodonH_2=y1H2
              znodonH_2=z1H2

              xacc_2=x3O
              yacc_2=y3O
              zacc_2=z3O
              xHacc_2=x3H1
              yHacc_2=y3H1
              zHacc_2=z3H1

              xHacc2_2=x3H2
              yHacc2_2=y3H2
              zHacc2_2=z3H2
            end if

           xr = x3O - x1H2
           yr = y3O - y1H2
           zr = z3O - z1H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond4= r
            if(hbond4.lt.rlowhbond_2) then
              rlowhbond_2 = hbond4
              xdonO_2=x1O
              ydonO_2=y1O
              zdonO_2=z1O
              xdonH_2=x1H2
              ydonH_2=y1H2
              zdonH_2=z1H2
              xnodonH_2=x1H1
              ynodonH_2=y1H1
              znodonH_2=z1H1

              xacc_2=x3O
              yacc_2=y3O
              zacc_2=z3O
              xHacc_2=x3H1
              yHacc_2=y3H1
              zHacc_2=z3H1

              xHacc2_2=x3H2
              yHacc2_2=y3H2
              zHacc2_2=z3H2
            end if

            xr = xdonO_2-xdonH_2
            yr = ydonO_2-ydonH_2
            zr = zdonO_2-zdonH_2
            r2 = xr*xr + yr*yr + zr*zr
            rdonOdonH2=r2
            rdonOdonH=sqrt(rdonOdonH2)

            xr = xdonO_2-xacc_2
            yr = ydonO_2-yacc_2
            zr = zdonO_2-zacc_2
            r2 = xr*xr + yr*yr + zr*zr
            rdonOacc2=r2

c            rdonOacc2 = rdonOdonH2 + rlowhbond*rlowhbond - 2*rlowhbond*rdonOdonH*cos(angle)
c            2*rlowhbond*rdonOdonH*cos(hbangle)=rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2

c            cos(hbangle)= (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)/(2*rlowhbond*rdonOdonH)
        hbangle_2=acos( (rdonOdonH2 + rlowhbond_2*rlowhbond_2-rdonOacc2)
     &                  /(2*rlowhbond_2*rdonOdonH))

c FIRST CALCULATE PLANE OF HYDROGEN BOND

         xr1= xacc_2-xdonH_2
         yr1= yacc_2-ydonH_2
         zr1= zacc_2-zdonH_2

         xr2= xdonO_2-xdonH_2
         yr2= ydonO_2-ydonH_2
         zr2= zdonO_2-zdonH_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal1_x=(yr1*zr2-yr2*zr1)/norm
         normal1_y=(xr2*zr1-xr1*zr2)/norm
         normal1_z=(xr1*yr2-xr2*yr1)/norm

c NEXT CALCULATE PLANE FORMED BY DONATED HYDROGEN, DONOR OXYGEN, AND NONDONATED HYDROGEN

         xr1= xdonH_2-xdonO_2
         yr1= ydonH_2-ydonO_2
         zr1= zdonH_2-zdonO_2

         xr2= xnodonH_2-xdonO_2
         yr2= ynodonH_2-ydonO_2
         zr2= znodonH_2-zdonO_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal2_x=(yr1*zr2-yr2*zr1)/norm
         normal2_y=(xr2*zr1-xr1*zr2)/norm
         normal2_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY HYDROGEN COVALENTLY BOUND TO ACCEPTOR OXYGEN, THE ACCEPTOR OXYGEN, AND THE DONATED HYDROGEN

         xr1 = xHacc_2-xacc_2
         yr1 = yHacc_2-yacc_2
         zr1 = zHacc_2-zacc_2

         xr2 = xdonH_2-xacc_2
         yr2 = ydonH_2-yacc_2
         zr2 = zdonH_2-zacc_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal3_x=(yr1*zr2-yr2*zr1)/norm
         normal3_y=(xr2*zr1-xr1*zr2)/norm
         normal3_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY ATOMS ON ACCEPTOR

         xr1 = xHacc_2-xacc_2
         yr1 = yHacc_2-yacc_2
         zr1 = zHacc_2-zacc_2

         xr2 = xHacc2_2-xacc_2
         yr2 = yHacc2_2-yacc_2
         zr2 = zHacc2_2-zacc_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal4_x=(yr1*zr2-yr2*zr1)/norm
         normal4_y=(xr2*zr1-xr1*zr2)/norm
         normal4_z=(xr1*yr2-xr2*yr1)/norm

         phi1_2= acos(normal1_x*normal2_x + normal1_y*normal2_y
     &          + normal1_z*normal2_z)

         chi_2= acos(normal1_x*normal3_x + normal1_y*normal3_y
     &          + normal1_z*normal3_z)

         theta2_2=acos(normal1_x*normal4_x + normal1_y*normal4_y
     &          + normal1_z*normal4_z)

         phi2_2= acos(normal2_x*normal4_x + normal2_y*normal4_y
     &          + normal2_z*normal4_z)



           xr = x2O - x3H1
           yr = y2O - y3H1
           zr = z2O - z3H1

c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond1= r
            rlowhbond_3 = hbond1
            xdonO_3= x3O
            ydonO_3= y3O
            zdonO_3= z3O
            xdonH_3= x3H1
            ydonH_3= y3H1
            zdonH_3= z3H1
            xnodonH_3 = x3H2
            ynodonH_3 = y3H2
            znodonH_3 = z3H2

            xacc_3 = x2O
            yacc_3 = y2O
            zacc_3 = z2O
            xHacc_3 = x2H1
            yHacc_3 = y2H1
            zHacc_3 = z2H1

            xHacc2_3 = x2H2
            yHacc2_3 = y2H2
            zHacc2_3 = z2H2


            xr = x2O - x3H2
            yr = y2O - y3H2
            zr = z2O - z3H2

            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond2= r
            if(hbond2.lt.rlowhbond_3) then
              rlowhbond_3 = hbond2
            xdonO_3= x3O
            ydonO_3= y3O
            zdonO_3= z3O
            xdonH_3= x3H2
            ydonH_3= y3H2
            zdonH_3= z3H2
            xnodonH_3 = x3H1
            ynodonH_3 = y3H1
            znodonH_3 = z3H1

            xacc_3 = x2O
            yacc_3 = y2O
            zacc_3 = z2O
            xHacc_3 = x2H1
            yHacc_3 = y2H1
            zHacc_3 = z2H1

            xHacc2_3 = x2H2
            yHacc2_3 = y2H2
            zHacc2_3 = z2H2
            end if

          xr = x3O - x2H1
          yr = y3O - y2H1
          zr = z3O - z2H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond3= r
            if(hbond3.lt.rlowhbond_3) then
              rlowhbond_3 = hbond3
              xdonO_3=x2O
              ydonO_3=y2O
              zdonO_3=z2O
              xdonH_3=x2H1
              ydonH_3=y2H1
              zdonH_3=z2H1
              xnodonH_3=x2H2
              ynodonH_3=y2H2
              znodonH_3=z2H2

              xacc_3=x3O
              yacc_3=y3O
              zacc_3=z3O
              xHacc_3=x3H1
              yHacc_3=y3H1
              zHacc_3=z3H1

              xHacc2_3=x3H2
              yHacc2_3=y3H2
              zHacc2_3=z3H2
            end if

           xr = x3O - x2H2
           yr = y3O - y2H2
           zr = z3O - z2H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond4= r
            if(hbond4.lt.rlowhbond_3) then
              rlowhbond_3 = hbond4
              xdonO_3=x2O
              ydonO_3=y2O
              zdonO_3=z2O
              xdonH_3=x2H2
              ydonH_3=y2H2
              zdonH_3=z2H2
              xnodonH_3=x2H1
              ynodonH_3=y2H1
              znodonH_3=z2H1

              xacc_3=x3O
              yacc_3=y3O
              zacc_3=z3O
              xHacc_3=x3H1
              yHacc_3=y3H1
              zHacc_3=z3H1

              xHacc2_3=x3H2
              yHacc2_3=y3H2
              zHacc2_3=z3H2
            end if

            xr = xdonO_3-xdonH_3
            yr = ydonO_3-ydonH_3
            zr = zdonO_3-zdonH_3
            r2 = xr*xr + yr*yr + zr*zr
            rdonOdonH2=r2
            rdonOdonH=sqrt(rdonOdonH2)

            xr = xdonO_3-xacc_3
            yr = ydonO_3-yacc_3
            zr = zdonO_3-zacc_3
            r2 = xr*xr + yr*yr + zr*zr
            rdonOacc2=r2

c            rdonOacc2 = rdonOdonH2 + rlowhbond*rlowhbond - 2*rlowhbond*rdonOdonH*cos(angle)
c            2*rlowhbond*rdonOdonH*cos(hbangle)=rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2

c            cos(hbangle)= (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)/(2*rlowhbond*rdonOdonH)
        hbangle_3=acos( (rdonOdonH2 + rlowhbond_3*rlowhbond_3-rdonOacc2)
     &                  /(2*rlowhbond_3*rdonOdonH))

c FIRST CALCULATE PLANE OF HYDROGEN BOND

         xr1= xacc_3-xdonH_3
         yr1= yacc_3-ydonH_3
         zr1= zacc_3-zdonH_3

         xr2= xdonO_3-xdonH_3
         yr2= ydonO_3-ydonH_3
         zr2= zdonO_3-zdonH_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal1_x=(yr1*zr2-yr2*zr1)/norm
         normal1_y=(xr2*zr1-xr1*zr2)/norm
         normal1_z=(xr1*yr2-xr2*yr1)/norm

c NEXT CALCULATE PLANE FORMED BY DONATED HYDROGEN, DONOR OXYGEN, AND NONDONATED HYDROGEN

         xr1= xdonH_3-xdonO_3
         yr1= ydonH_3-ydonO_3
         zr1= zdonH_3-zdonO_3

         xr2= xnodonH_3-xdonO_3
         yr2= ynodonH_3-ydonO_3
         zr2= znodonH_3-zdonO_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal2_x=(yr1*zr2-yr2*zr1)/norm
         normal2_y=(xr2*zr1-xr1*zr2)/norm
         normal2_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY HYDROGEN COVALENTLY BOUND TO ACCEPTOR OXYGEN, THE ACCEPTOR OXYGEN, AND THE DONATED HYDROGEN

         xr1 = xHacc_3-xacc_3
         yr1 = yHacc_3-yacc_3
         zr1 = zHacc_3-zacc_3

         xr2 = xdonH_3-xacc_3
         yr2 = ydonH_3-yacc_3
         zr2 = zdonH_3-zacc_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal3_x=(yr1*zr2-yr2*zr1)/norm
         normal3_y=(xr2*zr1-xr1*zr2)/norm
         normal3_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY ATOMS ON ACCEPTOR

         xr1 = xHacc_3-xacc_3
         yr1 = yHacc_3-yacc_3
         zr1 = zHacc_3-zacc_3

         xr2 = xHacc2_3-xacc_3
         yr2 = yHacc2_3-yacc_3
         zr2 = zHacc2_3-zacc_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal4_x=(yr1*zr2-yr2*zr1)/norm
         normal4_y=(xr2*zr1-xr1*zr2)/norm
         normal4_z=(xr1*yr2-xr2*yr1)/norm

         phi1_3= acos(normal1_x*normal2_x + normal1_y*normal2_y
     &          + normal1_z*normal2_z)

         chi_3= acos(normal1_x*normal3_x + normal1_y*normal3_y
     &          + normal1_z*normal3_z)

         theta2_3=acos(normal1_x*normal4_x + normal1_y*normal4_y
     &          + normal1_z*normal4_z)

         phi2_3= acos(normal2_x*normal4_x + normal2_y*normal4_y
     &          + normal2_z*normal4_z)
                        
               call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,
     &               pnum,eptemp,deptemp,virtemp,deptempmat)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,l1) = deptemp(1,l1)
                    dep3moli123(2,l1) = deptemp(2,l1)
                    dep3moli123(3,l1) = deptemp(3,l1)
                    do l3 = 1,npole3b
                       k=pnum(l3)
                       depmat3moli123(1,k,i)=deptempmat(1,k,i)
                       depmat3moli123(2,k,i)=deptempmat(2,k,i)
                       depmat3moli123(3,k,i)=deptempmat(3,k,i)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=virtemp(j,i)
                    end do
                  end do

c          if(printforcematbody) then
c            do l1=1,npole3b
c               i=pnum(l1)
c               do l3=1,npole3b
c                  k=pnum(l3)
c                      xr = x(k)-x(i)
c                      yr = y(k)-y(i)
c                      zr = z(k)-z(i)
c                      r=sqrt(xr*xr + yr*yr + zr*zr)
c       write(iout,61) moli1,moli2,moli3,i,k,r,xr,deptempmat(1,k,i)
c       write(iout,51) moli1,moli2,moli3,i,k,r,yr,deptempmat(2,k,i)
c       write(iout,41) moli1,moli2,moli3,i,k,r,zr,deptempmat(3,k,i)
c  61  format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,
c     &      1x,'r=',
c     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5)
c  51  format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,
c     &      1x,'r=',
c     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5)
c  41  format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,
c     &      1x,'r=',
c     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5)
c               end do
c            end do
c          end if
            


               npole3b=6
               ep3moli123=ep3moli123-ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3moli123(j,l1)=dep3moli123(j,l1)-dep2moli12(j,l1)
                  end do
                  do l3 =1,npole3b
                     k=pnum(l3)
                     do j =1,3
                     depmat3moli123(j,k,i)=depmat3moli123(j,k,i)
     &                -depmat2moli12(j,k,i)
                     end do
                  end do 
               end do
               do i=1,3
                 do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)-vir2moli12(j,i)
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
               call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,
     &               pnum,eptemp,deptemp,virtemp,deptempmat)
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       end do
                     end do
                     do l1=1,npole3b
                        i=pnum(l1)
                        do l3=1,npole3b
                           k=pnum(l3)
                           do j=1,3
                     depmat3moli123(j,k,i)=depmat3moli123(j,k,i)
     &                -deptempmat(j,k,i)
                           end do
                        end do
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
               call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,
     &                pnum,eptemp,deptemp,virtemp,deptempmat)
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

                     do l1=1,npole3b
                        i=pnum(l1)
                        do l3=1,npole3b
                           k=pnum(l3)
                           do j=1,3
                     depmat3moli123(j,k,i)=depmat3moli123(j,k,i)
     &                -deptempmat(j,k,i)                            
                           end do
                        end do
                     end do

                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                      end do
                     end do

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

                  ep3bt = ep3bt + ep3moli123
                  do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,l1)
                    end do
                    do l3 = 1,npole3b
                       k = pnum(l3)
                       do j = 1,3
                        depmat3bt(j,k,i)=depmat3bt(j,k,i)
     &                    + depmat3moli123(j,k,i)
                       end do
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do

          if(printforcematbody3) then
            molnum(1)=moli1
            molnum(2)=moli1
            molnum(3)=moli1
            molnum(4)=moli2
            molnum(5)=moli2
            molnum(6)=moli2
            molnum(7)=moli3
            molnum(8)=moli3
            molnum(9)=moli3
            do l1=1,npole3b
               i=pnum(l1)
               do l3=1,npole3b
                  k=pnum(l3)
                  if(i.lt.k) then
               !   if(molnum(l1).ne.molnum(l3)) then

                      xr = x(k)-x(i)
                      yr = y(k)-y(i)
                      zr = z(k)-z(i)
                      r=sqrt(xr*xr + yr*yr + zr*zr)
c                write(iout,31) moli1,moli2,i,k,r,xr,depmat2moli12(1,k,i)
c                write(iout,21) moli1,moli2,i,k,r,yr,depmat2moli12(2,k,i)
c                write(iout,11) moli1,moli2,i,k,r,zr,depmat2moli12(3,k,i)
c   31 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5)
c   21 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5)
c   11 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5)

         write(iout,61) moli1,moli2,moli3,i,k,r,xr,depmat3moli123(1,k,i)
     &   ,rlowhbond,hbangle,rlowhbond_2,hbangle_2,rlowhbond_3,hbangle_3,
     &   phi1,chi,theta2,phi2,
     &   phi1_2,chi_2,theta2_2,phi2_2,
     &   phi1_3,chi_3,theta2_3,phi2_3,
     &   ep3moli123
         write(iout,51) moli1,moli2,moli3,i,k,r,yr,depmat3moli123(2,k,i)
     &   ,rlowhbond,hbangle,rlowhbond_2,hbangle_2,rlowhbond_3,hbangle_3,
     &   phi1,chi,theta2,phi2,
     &   phi1_2,chi_2,theta2_2,phi2_2,
     &   phi1_3,chi_3,theta2_3,phi2_3,
     &   ep3moli123
         write(iout,41) moli1,moli2,moli3,i,k,r,zr,depmat3moli123(3,k,i)
     &   ,rlowhbond,hbangle,rlowhbond_2,hbangle_2,rlowhbond_3,hbangle_3,
     &   phi1,chi,theta2,phi2,
     &   phi1_2,chi_2,theta2_2,phi2_2,
     &   phi1_3,chi_3,theta2_3,phi2_3,
     &   ep3moli123
   61 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,1x
     &    ,'r=',f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5,1x,"hbr1=",f7.3,1x,
     &       "hbang1=",f7.3,1x,"hbr2=",f7.3,1x,"hbang2=",f7.3,
     &       1x,"hbr3=",f7.3,1x,"hbang3=",f7.3,
     &       1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,
     &       1x,"phi1_2=",f7.3,1x,"chi_2=",f7.3,
     &       1x,"theta2_2=",f7.3,1x,"phi2_2=",f7.3,
     &       1x,"phi1_3=",f7.3,1x,"chi_3=",f7.3,
     &       1x,"theta2_3=",f7.3,1x,"phi2_3=",f7.3,1x,"ep3=",f7.3)
   51 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,1x
     &    ,'r=',f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5,1x,"hbr1=",f7.3,1x,
     &       "hbang1=",f7.3,1x,"hbr2=",f7.3,1x,"hbang2=",f7.3,
     &       1x,"hbr3=",f7.3,1x,"hbang3=",f7.3,
     &       1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,
     &       1x,"phi1_2=",f7.3,1x,"chi_2=",f7.3,
     &       1x,"theta2_2=",f7.3,1x,"phi2_2=",f7.3,
     &       1x,"phi1_3=",f7.3,1x,"chi_3=",f7.3,
     &       1x,"theta2_3=",f7.3,1x,"phi2_3=",f7.3,1x,"ep3=",f7.3)
   41 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,1x
     &    ,'r=',f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5,1x,"hbr1=",f7.3,1x,
     &       "hbang1=",f7.3,1x,"hbr2=",f7.3,1x,"hbang2=",f7.3,
     &       1x,"hbr3=",f7.3,1x,"hbang3=",f7.3,
     &       1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,
     &       1x,"phi1_2=",f7.3,1x,"chi_2=",f7.3,
     &       1x,"theta2_2=",f7.3,1x,"phi2_2=",f7.3,
     &       1x,"phi1_3=",f7.3,1x,"chi_3=",f7.3,
     &       1x,"theta2_3=",f7.3,1x,"phi2_3=",f7.3,1x,"ep3=",f7.3)
                  end if
               end do
            end do
          end if

          end do
        end do 
      end do  
!$OMP END DO
!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt

c
      if(moli1rmndr.eq.0) then
        return
      else
        goto 32
      end if

   32 continue


c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,depmat3bt,
c!$OMP& nmol,moli1rmndr)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,depmat3bt)
c!$OMP& schedule(guided)
c      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
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

            hcount=0
            do l1 = 1, np1
              i = pnum(l1)
              if(name(i).eq.'O') then
                x1O = x(i)
                y1O = y(i)
                z1O = z(i)
              else if((name(i).eq.'H').and.(hcount.eq.0)) then
                x1H1 = x(i)
                y1H1 = y(i)
                z1H1 = z(i)
                hcount=hcount+1
              else if((name(i).eq.'H').and.(hcount.eq.1)) then
                x1H2 = x(i)
                y1H2 = y(i)
                z1H2 = z(i)
                hcount=hcount+1
              end if
            end do

            hcount=0
            do l1 = np1+1, np2
              i = pnum(l1)
              if(name(i).eq.'O') then
                x2O = x(i)
                y2O = y(i)
                z2O = z(i)
              else if((name(i).eq.'H').and.(hcount.eq.0)) then
                x2H1 = x(i)
                y2H1 = y(i)
                z2H1 = z(i)
                hcount=hcount+1
              else if((name(i).eq.'H').and.(hcount.eq.1)) then
                x2H2 = x(i)
                y2H2 = y(i)
                z2H2 = z(i)
                hcount=hcount+1
              end if
            end do

          xr = x1O - x2H1
          yr = y1O - y2H1
          zr = z1O - z2H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond1= r
            rlowhbond = hbond1
            xdonO= x2O
            ydonO= y2O
            zdonO= z2O
            xdonH= x2H1
            ydonH= y2H1
            zdonH= z2H1
            xnodonH = x2H2
            ynodonH = y2H2
            znodonH = z2H2

            xacc = x1O
            yacc = y1O
            zacc = z1O
            xHacc = x1H1
            yHacc = y1H1
            zHacc = z1H1

            xHacc2 = x1H2
            yHacc2 = y1H2
            zHacc2 = z1H2
            
          xr = x1O - x2H2
          yr = y1O - y2H2
          zr = z1O - z2H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond2= r
            if(hbond2.lt.rlowhbond) then
              rlowhbond = hbond2
              xdonO= x2O
              ydonO= y2O
              zdonO= z2O
              xdonH= x2H2
              ydonH= y2H2
              zdonH= z2H2
              xnodonH=x2H1
              ynodonH=y2H1
              znodonH=z2H1

              xacc = x1O
              yacc = y1O
              zacc = z1O
              xHacc = x1H1
              yHacc = y1H1
              zHacc = z1H1
            xHacc2 = x1H2
            yHacc2 = y1H2
            zHacc2 = z1H2
            end if

          xr = x2O - x1H1
          yr = y2O - y1H1
          zr = z2O - z1H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond3= r
            if(hbond3.lt.rlowhbond) then
              rlowhbond = hbond3
              xdonO=x1O
              ydonO=y1O
              zdonO=z1O
              xdonH=x1H1
              ydonH=y1H1
              zdonH=z1H1
              xnodonH=x1H2
              ynodonH=y1H2
              znodonH=z1H2

              xacc=x2O
              yacc=y2O
              zacc=z2O
              xHacc=x2H1
              yHacc=y2H1
              zHacc=z2H1

              xHacc2=x2H2
              yHacc2=y2H2
              zHacc2=z2H2
            end if

          xr = x2O - x1H2
          yr = y2O - y1H2
          zr = z2O - z1H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond4= r
            if(hbond4.lt.rlowhbond) then
              rlowhbond = hbond4
              xdonO=x1O
              ydonO=y1O
              zdonO=z1O
              xdonH=x1H2
              ydonH=y1H2
              zdonH=z1H2
              xnodonH=x1H1
              ynodonH=y1H1
              znodonH=z1H1

              xacc=x2O
              yacc=y2O
              zacc=z2O
              xHacc=x2H1
              yHacc=y2H1
              zHacc=z2H1

              xHacc2=x2H2
              yHacc2=y2H2
              zHacc2=z2H2
            end if

            xr = xdonO-xdonH
            yr = ydonO-ydonH
            zr = zdonO-zdonH
            r2 = xr*xr + yr*yr + zr*zr
            rdonOdonH2=r2
            rdonOdonH=sqrt(rdonOdonH2)

            xr = xdonO-xacc
            yr = ydonO-yacc
            zr = zdonO-zacc
            r2 = xr*xr + yr*yr + zr*zr
            rdonOacc2=r2

c            rdonOacc2 = rdonOdonH2 + rlowhbond*rlowhbond - 2*rlowhbond*rdonOdonH*cos(angle)
c            2*rlowhbond*rdonOdonH*cos(hbangle)=rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2

c            cos(hbangle)= (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)/(2*rlowhbond*rdonOdonH)
            hbangle=acos( (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)
     &                  /(2*rlowhbond*rdonOdonH))

c FIRST CALCULATE PLANE OF HYDROGEN BOND

         xr1= xacc-xdonH
         yr1= yacc-ydonH
         zr1= zacc-zdonH

         xr2= xdonO-xdonH
         yr2= ydonO-ydonH
         zr2= zdonO-zdonH

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal1_x=(yr1*zr2-yr2*zr1)/norm
         normal1_y=(xr2*zr1-xr1*zr2)/norm
         normal1_z=(xr1*yr2-xr2*yr1)/norm

c NEXT CALCULATE PLANE FORMED BY DONATED HYDROGEN, DONOR OXYGEN, AND NONDONATED HYDROGEN

         xr1= xdonH-xdonO
         yr1= ydonH-ydonO
         zr1= zdonH-zdonO

         xr2= xnodonH-xdonO
         yr2= ynodonH-ydonO
         zr2= znodonH-zdonO

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal2_x=(yr1*zr2-yr2*zr1)/norm
         normal2_y=(xr2*zr1-xr1*zr2)/norm
         normal2_z=(xr1*yr2-xr2*yr1)/norm

c LASTLY CALCULATE PLANE FORMED BY HYDROGEN COVALENTLY BOUND TO ACCEPTOR OXYGEN, THE ACCEPTOR OXYGEN, AND THE DONATED HYDROGEN

         xr1 = xHacc-xacc
         yr1 = yHacc-yacc
         zr1 = zHacc-zacc

         xr2 = xdonH-xacc
         yr2 = ydonH-yacc
         zr2 = zdonH-zacc

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal3_x=(yr1*zr2-yr2*zr1)/norm
         normal3_y=(xr2*zr1-xr1*zr2)/norm
         normal3_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY ATOMS ON ACCEPTOR

         xr1 = xHacc-xacc
         yr1 = yHacc-yacc
         zr1 = zHacc-zacc

         xr2 = xHacc2-xacc
         yr2 = yHacc2-yacc
         zr2 = zHacc2-zacc

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal4_x=(yr1*zr2-yr2*zr1)/norm
         normal4_y=(xr2*zr1-xr1*zr2)/norm
         normal4_z=(xr1*yr2-xr2*yr1)/norm

         phi1= acos(normal1_x*normal2_x + normal1_y*normal2_y
     &          + normal1_z*normal2_z)

         chi= acos(normal1_x*normal3_x + normal1_y*normal3_y
     &          + normal1_z*normal3_z)

         theta2=acos(normal1_x*normal4_x + normal1_y*normal4_y
     &          + normal1_z*normal4_z)

         phi2= acos(normal2_x*normal4_x + normal2_y*normal4_y
     &          + normal2_z*normal4_z)

         xbisect=xnodonH+(xnodonH-xdonH)/2.0d0
         ybisect=ynodonH+(ynodonH-ydonH)/2.0d0
         zbisect=znodonH+(znodonH-zdonH)/2.0d0

         call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,pnum,
     &           eptemp,deptemp,virtemp,deptempmat)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,l1)=deptemp(j,l1)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
                end do
                do l3=1,npole3b
                   k=pnum(l3)
                   do j = 1, 3
                      depmat2moli12(j,k,i)=deptempmat(j,k,i)
                      depmat3bt(j,k,i)=depmat3bt(j,k,i)
     &                                 +depmat2moli12(j,k,i)
                   end do
                end do
              end do

          if(printforcematbody) then
            molnum(1)=moli1rmndr
            molnum(2)=moli1rmndr
            molnum(3)=moli1rmndr
            molnum(4)=moli2
            molnum(5)=moli2
            molnum(6)=moli2
            do l1=1,npole3b
               i=pnum(l1)
               do l3=1,npole3b
                  k=pnum(l3)
                  if(i.lt.k) then
                !  if(molnum(l1).ne.molnum(l3)) then
                      xr = x(k)-x(i)
                      yr = y(k)-y(i)
                      zr = z(k)-z(i)
                      r=sqrt(xr*xr + yr*yr + zr*zr)
c           write(iout,91) moli1rmndr,moli2,i,k,r,xr,depmat2moli12(1,k,i)
c           write(iout,81) moli1rmndr,moli2,i,k,r,yr,depmat2moli12(2,k,i)
c           write(iout,71) moli1rmndr,moli2,i,k,r,zr,depmat2moli12(3,k,i)
c   91 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5)
c   81 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5)
c   71 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5)

           write(iout,91) moli1rmndr,moli2,i,k,r,xr,depmat2moli12(1,k,i)
     &          ,rlowhbond,hbangle,phi1,chi,theta2,phi2,ep2moli12
           write(iout,81) moli1rmndr,moli2,i,k,r,yr,depmat2moli12(2,k,i)
     &          ,rlowhbond,hbangle,phi1,chi,theta2,phi2,ep2moli12
           write(iout,71) moli1rmndr,moli2,i,k,r,zr,depmat2moli12(3,k,i)
     &         ,rlowhbond,hbangle,phi1,chi,theta2,phi2,ep2moli12
   91 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5,1x,"hbr=",f7.3,1x,
     &       "hbang=",f7.3,1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,1x,"ep2=",f7.3)
   81 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5,1x,"hbr=",f7.3,1x,
     &        "hbang=",f7.3,1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,1x,"ep2=",f7.3)
   71 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5,1x,"hbr=",f7.3,1x,
     &        "hbang=",f7.3,1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,1x,"ep2=",f7.3)
                  end if  
               end do
            end do
          end if

              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do


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

            hcount=0
            do l1 = np2+1, np3
              i = pnum(l1)
              if(name(i).eq.'O') then
                x3O = x(i)
                y3O = y(i)
                z3O = z(i)
              else if((name(i).eq.'H').and.(hcount.eq.0)) then
                x3H1 = x(i)
                y3H1 = y(i)
                z3H1 = z(i)
                hcount=hcount+1
              else if((name(i).eq.'H').and.(hcount.eq.1)) then
                x3H2 = x(i)
                y3H2 = y(i)
                z3H2 = z(i)
                hcount=hcount+1
              end if
            end do

           xr = x1O - x3H1
           yr = y1O - y3H1
           zr = z1O - z3H1

c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond1= r
            rlowhbond_2 = hbond1
            xdonO_2= x3O
            ydonO_2= y3O
            zdonO_2= z3O
            xdonH_2= x3H1
            ydonH_2= y3H1
            zdonH_2= z3H1
            xnodonH_2 = x3H2
            ynodonH_2 = y3H2
            znodonH_2 = z3H2

            xacc_2 = x1O
            yacc_2 = y1O
            zacc_2 = z1O
            xHacc_2 = x1H1
            yHacc_2 = y1H1
            zHacc_2 = z1H1

            xHacc2_2 = x1H2
            yHacc2_2 = y1H2
            zHacc2_2 = z1H2


            xr = x1O - x3H2
            yr = y1O - y3H2
            zr = z1O - z3H2

            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond2= r
            if(hbond2.lt.rlowhbond_2) then
              rlowhbond_2 = hbond2
            xdonO_2= x3O
            ydonO_2= y3O
            zdonO_2= z3O
            xdonH_2= x3H2
            ydonH_2= y3H2
            zdonH_2= z3H2
            xnodonH_2 = x3H1
            ynodonH_2 = y3H1
            znodonH_2 = z3H1

            xacc_2 = x1O
            yacc_2 = y1O
            zacc_2 = z1O
            xHacc_2 = x1H1
            yHacc_2 = y1H1
            zHacc_2 = z1H1

            xHacc2_2 = x1H2
            yHacc2_2 = y1H2
            zHacc2_2 = z1H2
            end if

          xr = x3O - x1H1
          yr = y3O - y1H1
          zr = z3O - z1H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond3= r
            if(hbond3.lt.rlowhbond_2) then
              rlowhbond_2 = hbond3
              xdonO_2=x1O
              ydonO_2=y1O
              zdonO_2=z1O
              xdonH_2=x1H1
              ydonH_2=y1H1
              zdonH_2=z1H1
              xnodonH_2=x1H2
              ynodonH_2=y1H2
              znodonH_2=z1H2

              xacc_2=x3O
              yacc_2=y3O
              zacc_2=z3O
              xHacc_2=x3H1
              yHacc_2=y3H1
              zHacc_2=z3H1

              xHacc2_2=x3H2
              yHacc2_2=y3H2
              zHacc2_2=z3H2
            end if

           xr = x3O - x1H2
           yr = y3O - y1H2
           zr = z3O - z1H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond4= r
            if(hbond4.lt.rlowhbond_2) then
              rlowhbond_2 = hbond4
              xdonO_2=x1O
              ydonO_2=y1O
              zdonO_2=z1O
              xdonH_2=x1H2
              ydonH_2=y1H2
              zdonH_2=z1H2
              xnodonH_2=x1H1
              ynodonH_2=y1H1
              znodonH_2=z1H1

              xacc_2=x3O
              yacc_2=y3O
              zacc_2=z3O
              xHacc_2=x3H1
              yHacc_2=y3H1
              zHacc_2=z3H1

              xHacc2_2=x3H2
              yHacc2_2=y3H2
              zHacc2_2=z3H2
            end if

            xr = xdonO_2-xdonH_2
            yr = ydonO_2-ydonH_2
            zr = zdonO_2-zdonH_2
            r2 = xr*xr + yr*yr + zr*zr
            rdonOdonH2=r2
            rdonOdonH=sqrt(rdonOdonH2)

            xr = xdonO_2-xacc_2
            yr = ydonO_2-yacc_2
            zr = zdonO_2-zacc_2
            r2 = xr*xr + yr*yr + zr*zr
            rdonOacc2=r2

c            rdonOacc2 = rdonOdonH2 + rlowhbond*rlowhbond - 2*rlowhbond*rdonOdonH*cos(angle)
c            2*rlowhbond*rdonOdonH*cos(hbangle)=rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2

c            cos(hbangle)= (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)/(2*rlowhbond*rdonOdonH)
        hbangle_2=acos( (rdonOdonH2 + rlowhbond_2*rlowhbond_2-rdonOacc2)
     &                  /(2*rlowhbond_2*rdonOdonH))

c FIRST CALCULATE PLANE OF HYDROGEN BOND

         xr1= xacc_2-xdonH_2
         yr1= yacc_2-ydonH_2
         zr1= zacc_2-zdonH_2

         xr2= xdonO_2-xdonH_2
         yr2= ydonO_2-ydonH_2
         zr2= zdonO_2-zdonH_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal1_x=(yr1*zr2-yr2*zr1)/norm
         normal1_y=(xr2*zr1-xr1*zr2)/norm
         normal1_z=(xr1*yr2-xr2*yr1)/norm

c NEXT CALCULATE PLANE FORMED BY DONATED HYDROGEN, DONOR OXYGEN, AND NONDONATED HYDROGEN

         xr1= xdonH_2-xdonO_2
         yr1= ydonH_2-ydonO_2
         zr1= zdonH_2-zdonO_2

         xr2= xnodonH_2-xdonO_2
         yr2= ynodonH_2-ydonO_2
         zr2= znodonH_2-zdonO_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal2_x=(yr1*zr2-yr2*zr1)/norm
         normal2_y=(xr2*zr1-xr1*zr2)/norm
         normal2_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY HYDROGEN COVALENTLY BOUND TO ACCEPTOR OXYGEN, THE ACCEPTOR OXYGEN, AND THE DONATED HYDROGEN

         xr1 = xHacc_2-xacc_2
         yr1 = yHacc_2-yacc_2
         zr1 = zHacc_2-zacc_2

         xr2 = xdonH_2-xacc_2
         yr2 = ydonH_2-yacc_2
         zr2 = zdonH_2-zacc_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal3_x=(yr1*zr2-yr2*zr1)/norm
         normal3_y=(xr2*zr1-xr1*zr2)/norm
         normal3_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY ATOMS ON ACCEPTOR

         xr1 = xHacc_2-xacc_2
         yr1 = yHacc_2-yacc_2
         zr1 = zHacc_2-zacc_2

         xr2 = xHacc2_2-xacc_2
         yr2 = yHacc2_2-yacc_2
         zr2 = zHacc2_2-zacc_2

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal4_x=(yr1*zr2-yr2*zr1)/norm
         normal4_y=(xr2*zr1-xr1*zr2)/norm
         normal4_z=(xr1*yr2-xr2*yr1)/norm

         phi1_2= acos(normal1_x*normal2_x + normal1_y*normal2_y
     &          + normal1_z*normal2_z)

         chi_2= acos(normal1_x*normal3_x + normal1_y*normal3_y
     &          + normal1_z*normal3_z)

         theta2_2=acos(normal1_x*normal4_x + normal1_y*normal4_y
     &          + normal1_z*normal4_z)

         phi2_2= acos(normal2_x*normal4_x + normal2_y*normal4_y
     &          + normal2_z*normal4_z)

           xr = x2O - x3H1
           yr = y2O - y3H1
           zr = z2O - z3H1

c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond1= r
            rlowhbond_3 = hbond1
            xdonO_3= x3O
            ydonO_3= y3O
            zdonO_3= z3O
            xdonH_3= x3H1
            ydonH_3= y3H1
            zdonH_3= z3H1
            xnodonH_3 = x3H2
            ynodonH_3 = y3H2
            znodonH_3 = z3H2

            xacc_3 = x2O
            yacc_3 = y2O
            zacc_3 = z2O
            xHacc_3 = x2H1
            yHacc_3 = y2H1
            zHacc_3 = z2H1

            xHacc2_3 = x2H2
            yHacc2_3 = y2H2
            zHacc2_3 = z2H2


            xr = x2O - x3H2
            yr = y2O - y3H2
            zr = z2O - z3H2

            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond2= r
            if(hbond2.lt.rlowhbond_3) then
              rlowhbond_3 = hbond2
            xdonO_3= x3O
            ydonO_3= y3O
            zdonO_3= z3O
            xdonH_3= x3H2
            ydonH_3= y3H2
            zdonH_3= z3H2
            xnodonH_3 = x3H1
            ynodonH_3 = y3H1
            znodonH_3 = z3H1

            xacc_3 = x2O
            yacc_3 = y2O
            zacc_3 = z2O
            xHacc_3 = x2H1
            yHacc_3 = y2H1
            zHacc_3 = z2H1

            xHacc2_3 = x2H2
            yHacc2_3 = y2H2
            zHacc2_3 = z2H2
            end if

          xr = x3O - x2H1
          yr = y3O - y2H1
          zr = z3O - z2H1
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond3= r
            if(hbond3.lt.rlowhbond_3) then
              rlowhbond_3 = hbond3
              xdonO_3=x2O
              ydonO_3=y2O
              zdonO_3=z2O
              xdonH_3=x2H1
              ydonH_3=y2H1
              zdonH_3=z2H1
              xnodonH_3=x2H2
              ynodonH_3=y2H2
              znodonH_3=z2H2

              xacc_3=x3O
              yacc_3=y3O
              zacc_3=z3O
              xHacc_3=x3H1
              yHacc_3=y3H1
              zHacc_3=z3H1

              xHacc2_3=x3H2
              yHacc2_3=y3H2
              zHacc2_3=z3H2
            end if

           xr = x3O - x2H2
           yr = y3O - y2H2
           zr = z3O - z2H2
c          call image(xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r=sqrt(r2)
            hbond4= r
            if(hbond4.lt.rlowhbond_3) then
              rlowhbond_3 = hbond4
              xdonO_3=x2O
              ydonO_3=y2O
              zdonO_3=z2O
              xdonH_3=x2H2
              ydonH_3=y2H2
              zdonH_3=z2H2
              xnodonH_3=x2H1
              ynodonH_3=y2H1
              znodonH_3=z2H1

              xacc_3=x3O
              yacc_3=y3O
              zacc_3=z3O
              xHacc_3=x3H1
              yHacc_3=y3H1
              zHacc_3=z3H1

              xHacc2_3=x3H2
              yHacc2_3=y3H2
              zHacc2_3=z3H2
            end if

            xr = xdonO_3-xdonH_3
            yr = ydonO_3-ydonH_3
            zr = zdonO_3-zdonH_3
            r2 = xr*xr + yr*yr + zr*zr
            rdonOdonH2=r2
            rdonOdonH=sqrt(rdonOdonH2)

            xr = xdonO_3-xacc_3
            yr = ydonO_3-yacc_3
            zr = zdonO_3-zacc_3
            r2 = xr*xr + yr*yr + zr*zr
            rdonOacc2=r2

c            rdonOacc2 = rdonOdonH2 + rlowhbond*rlowhbond - 2*rlowhbond*rdonOdonH*cos(angle)
c            2*rlowhbond*rdonOdonH*cos(hbangle)=rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2

c            cos(hbangle)= (rdonOdonH2 + rlowhbond*rlowhbond-rdonOacc2)/(2*rlowhbond*rdonOdonH)
        hbangle_3=acos( (rdonOdonH2 + rlowhbond_3*rlowhbond_3-rdonOacc2)
     &                  /(2*rlowhbond_3*rdonOdonH))

c FIRST CALCULATE PLANE OF HYDROGEN BOND

         xr1= xacc_3-xdonH_3
         yr1= yacc_3-ydonH_3
         zr1= zacc_3-zdonH_3

         xr2= xdonO_3-xdonH_3
         yr2= ydonO_3-ydonH_3
         zr2= zdonO_3-zdonH_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal1_x=(yr1*zr2-yr2*zr1)/norm
         normal1_y=(xr2*zr1-xr1*zr2)/norm
         normal1_z=(xr1*yr2-xr2*yr1)/norm

c NEXT CALCULATE PLANE FORMED BY DONATED HYDROGEN, DONOR OXYGEN, AND NONDONATED HYDROGEN

         xr1= xdonH_3-xdonO_3
         yr1= ydonH_3-ydonO_3
         zr1= zdonH_3-zdonO_3

         xr2= xnodonH_3-xdonO_3
         yr2= ynodonH_3-ydonO_3
         zr2= znodonH_3-zdonO_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal2_x=(yr1*zr2-yr2*zr1)/norm
         normal2_y=(xr2*zr1-xr1*zr2)/norm
         normal2_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY HYDROGEN COVALENTLY BOUND TO ACCEPTOR OXYGEN, THE ACCEPTOR OXYGEN, AND THE DONATED HYDROGEN

         xr1 = xHacc_3-xacc_3
         yr1 = yHacc_3-yacc_3
         zr1 = zHacc_3-zacc_3

         xr2 = xdonH_3-xacc_3
         yr2 = ydonH_3-yacc_3
         zr2 = zdonH_3-zacc_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal3_x=(yr1*zr2-yr2*zr1)/norm
         normal3_y=(xr2*zr1-xr1*zr2)/norm
         normal3_z=(xr1*yr2-xr2*yr1)/norm

c  CALCULATE PLANE FORMED BY ATOMS ON ACCEPTOR

         xr1 = xHacc_3-xacc_3
         yr1 = yHacc_3-yacc_3
         zr1 = zHacc_3-zacc_3

         xr2 = xHacc2_3-xacc_3
         yr2 = yHacc2_3-yacc_3
         zr2 = zHacc2_3-zacc_3

         norm=sqrt( (yr1*zr2-yr2*zr1)*(yr1*zr2-yr2*zr1) +
     &     (xr2*zr1-xr1*zr2)*(xr2*zr1-xr1*zr2) +
     &     (xr1*yr2-xr2*yr1)*(xr1*yr2-xr2*yr1))

         normal4_x=(yr1*zr2-yr2*zr1)/norm
         normal4_y=(xr2*zr1-xr1*zr2)/norm
         normal4_z=(xr1*yr2-xr2*yr1)/norm

         phi1_3= acos(normal1_x*normal2_x + normal1_y*normal2_y
     &          + normal1_z*normal2_z)

         chi_3= acos(normal1_x*normal3_x + normal1_y*normal3_y
     &          + normal1_z*normal3_z)

         theta2_3=acos(normal1_x*normal4_x + normal1_y*normal4_y
     &          + normal1_z*normal4_z)

         phi2_3= acos(normal2_x*normal4_x + normal2_y*normal4_y
     &          + normal2_z*normal4_z)
            
               call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,
     &              pnum,eptemp,deptemp,virtemp,deptempmat)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,l1) = deptemp(1,l1)
                    dep3moli123(2,l1) = deptemp(2,l1)
                    dep3moli123(3,l1) = deptemp(3,l1)
                    do l3 = 1,npole3b
                       k=pnum(l3)
                       depmat3moli123(1,k,i)=deptempmat(1,k,i)
                       depmat3moli123(2,k,i)=deptempmat(2,k,i)
                       depmat3moli123(3,k,i)=deptempmat(3,k,i)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=virtemp(j,i)
                    end do
                  end do
            
c          if(printforcematbody) then
c            do l1=1,npole3b
c               i=pnum(l1)
c               do l3=1,npole3b
c                  k=pnum(l3)
c                      xr = x(k)-x(i)
c                      yr = y(k)-y(i)
c                      zr = z(k)-z(i)
c                      r=sqrt(xr*xr + yr*yr + zr*zr)
c       write(iout,121) moli1rmndr,moli2,moli3,i,k,r,xr,deptempmat(1,k,i)
c       write(iout,111) moli1rmndr,moli2,moli3,i,k,r,yr,deptempmat(2,k,i)
c       write(iout,101) moli1rmndr,moli2,moli3,i,k,r,zr,deptempmat(3,k,i)
c  121 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,
c     &      1x,'r=',
c     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5)
c  111 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,
c     &      1x,'r=',
c     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5)
c  101 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,
c     &      1x,'r=',
c     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5)
c               end do
c            end do
c          end if


               npole3b=6
               ep3moli123=ep3moli123-ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3moli123(j,l1)=dep3moli123(j,l1)-dep2moli12(j,l1)
                  end do
                  do l3 =1,npole3b
                     k=pnum(l3)
                     do j =1,3
                     depmat3moli123(j,k,i)=depmat3moli123(j,k,i)
     &                -depmat2moli12(j,k,i)
                     end do
                  end do 
               end do
               do i=1,3
                 do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)-vir2moli12(j,i)
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
               call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,
     &                pnum,eptemp,deptemp,virtemp,deptempmat)
                     ep3moli123 = ep3moli123 - eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,l2)=dep3moli123(j,l2)-deptemp(j,l1)
                       end do
                     end do
                     do l1=1,npole3b
                        i=pnum(l1)
                        do l3=1,npole3b
                           k=pnum(l3)
                           do j=1,3
                     depmat3moli123(j,k,i)=depmat3moli123(j,k,i)
     &                -deptempmat(j,k,i)
                           end do
                        end do
                     end do                     
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       end do
                     end do

                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
               call empole1a_3b_Polar_orig_nopriordir_fij(npole3b,
     &                pnum,eptemp,deptemp,virtemp,deptempmat)
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

                     do l1=1,npole3b
                        i=pnum(l1)
                        do l3=1,npole3b
                           k=pnum(l3)
                           do j=1,3
                     depmat3moli123(j,k,i)=depmat3moli123(j,k,i)
     &                -deptempmat(j,k,i)                            
                           end do
                        end do
                     end do

                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                      end do
                     end do

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

                  ep3bt = ep3bt + ep3moli123
                  do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,l1)
                    end do
                    do l3 = 1,npole3b
                       k = pnum(l3)
                       do j = 1,3
                        depmat3bt(j,k,i)=depmat3bt(j,k,i)
     &                    + depmat3moli123(j,k,i)
                       end do
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do

          if(printforcematbody3) then
            molnum(1)=moli1rmndr
            molnum(2)=moli1rmndr
            molnum(3)=moli1rmndr
            molnum(4)=moli2
            molnum(5)=moli2
            molnum(6)=moli2
            molnum(7)=moli3
            molnum(8)=moli3
            molnum(9)=moli3
            do l1=1,npole3b
               i=pnum(l1)
               do l3=1,npole3b
                  k=pnum(l3)
                  if(i.lt.k) then
               !   if(molnum(l1).ne.molnum(l3)) then

                      xr = x(k)-x(i)
                      yr = y(k)-y(i)
                      zr = z(k)-z(i)
                      r=sqrt(xr*xr + yr*yr + zr*zr)
c                write(iout,31) moli1,moli2,i,k,r,xr,depmat2moli12(1,k,i)
c                write(iout,21) moli1,moli2,i,k,r,yr,depmat2moli12(2,k,i)
c                write(iout,11) moli1,moli2,i,k,r,zr,depmat2moli12(3,k,i)
c   31 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5)
c   21 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5)
c   11 format(/,'2Bod mol1mol2=',i3,1x,i3,1x,'i k=',i4,1x,i4,1x,'r=',
c     &       f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5)

        write(iout,121) moli1rmndr,moli2,moli3,i,k,r,xr,
     &    depmat3moli123(1,k,i)
     &  ,rlowhbond,hbangle,rlowhbond_2,hbangle_2,rlowhbond_3,hbangle_3,
     &   phi1,chi,theta2,phi2,
     &   phi1_2,chi_2,theta2_2,phi2_2,
     &   phi1_3,chi_3,theta2_3,phi2_3,
     &   ep3moli123
        write(iout,111) moli1rmndr,moli2,moli3,i,k,r,yr,
     &    depmat3moli123(2,k,i)
     &  ,rlowhbond,hbangle,rlowhbond_2,hbangle_2,rlowhbond_3,hbangle_3,
     &   phi1,chi,theta2,phi2,
     &   phi1_2,chi_2,theta2_2,phi2_2,
     &   phi1_3,chi_3,theta2_3,phi2_3,
     &   ep3moli123
        write(iout,101) moli1rmndr,moli2,moli3,i,k,r,zr,
     &    depmat3moli123(3,k,i)
     &  ,rlowhbond,hbangle,rlowhbond_2,hbangle_2,rlowhbond_3,hbangle_3,
     &   phi1,chi,theta2,phi2,
     &   phi1_2,chi_2,theta2_2,phi2_2,
     &   phi1_3,chi_3,theta2_3,phi2_3,
     &   ep3moli123
  121 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,1x
     &    ,'r=',f7.3,1x,'xr=',f7.3,1x,'fijx=',f10.5,1x,"hbr1=",f7.3,1x,
     &       "hbang1=",f7.3,1x,"hbr2=",f7.3,1x,"hbang2=",f7.3,
     &       1x,"hbr3=",f7.3,1x,"hbang3=",f7.3,
     &       1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,
     &       1x,"phi1_2=",f7.3,1x,"chi_2=",f7.3,
     &       1x,"theta2_2=",f7.3,1x,"phi2_2=",f7.3,
     &       1x,"phi1_3=",f7.3,1x,"chi_3=",f7.3,
     &       1x,"theta2_3=",f7.3,1x,"phi2_3=",f7.3,1x,"ep3=",f7.3)
  111 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,1x
     &    ,'r=',f7.3,1x,'yr=',f7.3,1x,'fijy=',f10.5,1x,"hbr1=",f7.3,1x,
     &       "hbang1=",f7.3,1x,"hbr2=",f7.3,1x,"hbang2=",f7.3,
     &       1x,"hbr3=",f7.3,1x,"hbang3=",f7.3,
     &       1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,
     &       1x,"phi1_2=",f7.3,1x,"chi_2=",f7.3,
     &       1x,"theta2_2=",f7.3,1x,"phi2_2=",f7.3,
     &       1x,"phi1_3=",f7.3,1x,"chi_3=",f7.3,
     &       1x,"theta2_3=",f7.3,1x,"phi2_3=",f7.3,1x,"ep3=",f7.3)
  101 format(/,'3Bod mol1mol2mol3=',i3,1x,i3,1x,i3,1x,'i k=',i4,1x,i4,1x
     &    ,'r=',f7.3,1x,'zr=',f7.3,1x,'fijz=',f10.5,1x,"hbr1=",f7.3,1x,
     &       "hbang1=",f7.3,1x,"hbr2=",f7.3,1x,"hbang2=",f7.3,
     &       1x,"hbr3=",f7.3,1x,"hbang3=",f7.3,
     &       1x,"phi1=",f7.3,1x,"chi=",f7.3,
     &       1x,"theta2=",f7.3,1x,"phi2=",f7.3,
     &       1x,"phi1_2=",f7.3,1x,"chi_2=",f7.3,
     &       1x,"theta2_2=",f7.3,1x,"phi2_2=",f7.3,
     &       1x,"phi1_3=",f7.3,1x,"chi_3=",f7.3,
     &       1x,"theta2_3=",f7.3,1x,"phi2_3=",f7.3,1x,"ep3=",f7.3)
                  end if
               end do
            end do
          end if

          end do
        end do 
c      end do  
c!$OMP END DO
c!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1 single",ep3bt
      return
      end

