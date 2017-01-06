
      subroutine LoadBalNew2SmoothInnerloop0_2cut_wskinOrig(
     &  ep3bt) 
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
      !real*8  dep2moli12(3,30)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp!,deptemp(3,30)
      !real*8 vir2moli12(3,3)
      !real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3,np1,np2,np3
      real*8 ep3bt!,dep3bt(3,npole),virep3bt(3,3)
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
      real*8 ep3moli123!,vir3moli123(3,3),dep3moli123(3,30) 
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
      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b
c      rtapr3b=7.0d0
c      Cobarcut3b=8.5d0
      SwitchR3b=Cobarcut3b-rtapr3b
c      print*,"xcell2=",xcell2

      ep3bt=0.0d0
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

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& start_polar,last_polar,
!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
!$OMP& Rcut2b,SwitchR2b,SwitchR3b,molnew,taskid,
!$OMP& Rcut2b2)
!$OMP DO reduction(+:ep3bt)
!$OMP& schedule(guided)
      do kouter=start_polar(taskid),last_polar(taskid)
        moli1=molnew(kouter)
c      do moli1=start_polar(taskid),last_polar(taskid)

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

c        print*,"Anint moli12 xr1t yr1t zr1t",moli1,moli2,xr1t,yr1t,zr1t
         
          r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
          

C
C   Pair cutoffs
C
          if(r1_2.le.Rcut2b2) then
          r1=sqrt(r1_2)
          call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,eptemp)!,
          ep2moli12=eptemp
           if (r1.gt.rtapr2b) then
            term=(r1-rtapr2b)/SwitchR2b
            tapr2b_2=term*term
            tapr2b_3=tapr2b_2*term
            tapr2b_4=tapr2b_3*term
            tapr2b_5=tapr2b_4*term
            tapr2b = 1.0d0 - 10.0d0*tapr2b_3 + 15.0d0*tapr2b_4
     &                - 6.0d0*tapr2b_5

! NEED TO SMOOTHEN VIRIAL !

             ep3bt=ep3bt+tapr2b*ep2moli12
           else
              ep3bt=ep3bt+ep2moli12
           end if 
          else
            ep2moli12=0.0d0
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

            r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
            r2_2=xr2*xr2 + yr2*yr2 + zr2*zr2
            r3_2=xr3*xr3 + yr3*yr3 + zr3*zr3
            r1=sqrt(r1_2)         
            r2=sqrt(r2_2)
            r3=sqrt(r3_2)
            
           
            if ((r1.lt.r3).and.(r2.lt.r3) ) then 
              shellsum=r1+r2
              if (shellsum.le.Cobarcut3b) then
                call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp)!,deptemp,virtemp)
                ep3moli123=eptemp
            
                if (shellsum.gt.rtapr3b) then
                  term3=(shellsum-rtapr3b)/SwitchR3b
                  tapr3b_2=term3*term3
                  tapr3b_3=tapr3b_2*term3
                  tapr3b_4=tapr3b_3*term3
                  tapr3b_5=tapr3b_4*term3
                  tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
     &                       - 6.0d0*tapr3b_5
                end if
              end if             
            else if ((r1.lt.r2).and.(r3.lt.r2)) then
              shellsum=r1+r3
              if (shellsum.le.Cobarcut3b) then
                call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp)!,deptemp,virtemp)
                ep3moli123=eptemp

                 if (shellsum.gt.rtapr3b) then
                   term3=(shellsum-rtapr3b)/SwitchR3b
                   tapr3b_2=term3*term3
                   tapr3b_3=tapr3b_2*term3
                   tapr3b_4=tapr3b_3*term3
                   tapr3b_5=tapr3b_4*term3
                   tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
     &                - 6.0d0*tapr3b_5
                 end if
              end if
            else if ((r2.lt.r1).and.(r3.lt.r1)) then
              shellsum=r2+r3
               if (shellsum.le.Cobarcut3b) then
                call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp)!,deptemp,virtemp)
                  ep3moli123=eptemp

                 if (shellsum.gt.rtapr3b) then
                   term3=(shellsum-rtapr3b)/SwitchR3b
                   tapr3b_2=term3*term3
                   tapr3b_3=tapr3b_2*term3
                   tapr3b_4=tapr3b_3*term3
                   tapr3b_5=tapr3b_4*term3
                   tapr3b = 1.0d0 - 10.0d0*tapr3b_3 + 15.0d0*tapr3b_4
     &                - 6.0d0*tapr3b_5
                 end if
               end if
            end if

C LEFT OFF HERE!!!!!!!

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
                 else
                  call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &            eptemp)!,deptemp,virtemp)
                  ep3moli123=ep3moli123-eptemp
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
                  call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &            eptemp)!,deptemp,virtemp)
                     ep3moli123 = ep3moli123 - eptemp

                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole0a_3b_Polar_orig_nopriordir(npole3b,pnum,
     &           eptemp)!,deptemp,virtemp)
                     ep3moli123 = ep3moli123 - eptemp

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
               else
                  ep3bt = ep3bt + ep3moli123
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

