c
      subroutine All2bAll3bNew2SmoothInnerloop1_empole1c_totfield(
     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
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

        print*,"Totfield DEBUG NOCut3b NoCut2b"
        print*,"task start last",taskid,start,last,moli1rmndr

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,npole,
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
           call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
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

           call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &             deptemp,virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              !do i = 1, npole
                !i = pnum(l1)
              do i=1,npole
                do j = 1, 3
                  dep2moli12(j,i)=deptemp(j,i)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,i)
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
             call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &             deptemp,virtemp)
             ep1moli2=eptemp
             ep3bt=ep3bt-ep1moli2
              do i = 1, npole
                !i = pnum(l1)
               do j = 1, 3
                 dep1moli2(j,i)=deptemp(j,i)
                dep3bt(j,i)=dep3bt(j,i)-dep1moli2(j,i)
               end do
              end do
              do i=1,3
                do j=1,3
                 vir1moli2(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)-vir1moli2(j,i)
                end do
              end do


             pnum(1)=imol(1,moli1)
             pnum(2)=imol(1,moli1)+1
             pnum(3)=imol(2,moli1)
             ep3bt=ep3bt-ep1moli1
              do i = 1, npole
                !i = pnum(l1)
c              do i=1,npole
                do j = 1, 3
                   dep3bt(j,i)=dep3bt(j,i)-dep1moli1(j,i)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)-vir1moli1(j,i)
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

            
                 call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                ep3bt=ep3bt+eptemp
                ep3moli123=eptemp

                do i =1,npole
                  !i = pnum(l1)
c                do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,i)
                end do
                do i=1,3
                  do j=1,3
                  virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
                  end do
                end do
            


               npole3b=6
               ep3bt=ep3bt-ep2moli12
               ep3moli123=ep3moli123-ep2moli12
               do i = 1, npole
                  !i = pnum(l1)
                  do j = 1, 3
                    dep3bt(j,i)=dep3bt(j,i)-dep2moli12(j,i)
                  end do
               end do
              do i=1,3
                do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)-vir2moli12(j,i)
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
                  call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &            deptemp,virtemp)
                  ep3bt=ep3bt-eptemp
                  ep3moli123=ep3moli123-eptemp
                  do i =1,npole
                    !i = pnum(l1)
               ! do i=1,npole
                    dep3bt(1,i) = dep3bt(1,i)-deptemp(1,i)
                    dep3bt(2,i) = dep3bt(2,i)-deptemp(2,i)
                    dep3bt(3,i) = dep3bt(3,i)-deptemp(3,i)
                  end do
                do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i)-virtemp(j,i)
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
                 call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                  ep3bt = ep3bt - eptemp
                  ep3moli123=ep3moli123-eptemp
                do i =1,npole
                  !i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)-deptemp(3,i)
                end do
                do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i)-virtemp(j,i)
                  end do
                end do

                 np1=3
                 npole3b=3
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 ep3bt=ep3bt+ep1moli1
                 ep3moli123=ep3moli123+ep1moli1
                do i =1,npole
                  !i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli1(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli1(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+dep1moli1(3,i)
                end do
                do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i)+vir1moli1(j,i)
                  end do
                end do

                 npole3b=3
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                ep3bt=ep3bt+ep1moli2
                   ep3moli123=ep3moli123+ep1moli2
                do i =1,npole
                  !i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli2(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+dep1moli2(3,i)
                end do
                do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i)+vir1moli2(j,i)
                  end do
                end do


                 npole3b=3
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 call empole1c_3b_Polar_totfield(
     &           npole3b,pnum,eptemp,deptemp,virtemp)
                ep3bt=ep3bt+eptemp
                   ep3moli123=ep3moli123+eptemp

                do i =1,npole
                  !i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,i)
                end do
                do i=1,3
                  do j=1,3
                   virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
                  end do
                end do

c               print*,"3body Cobar ep3body",shellsum,ep3moli123


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
          call empole1c_3b_Polar_totfield(npole3b,
     &           pnum,eptemp,deptemp,virtemp)
          ep1moli1=eptemp
        !  print*,"moli1 ep1moli1=",moli1rmndr,ep1moli1
          ep3bt=ep3bt+ep1moli1
          do i = 1, npole
            !i = pnum(l1)
          !do i=1,npole
            do j = 1, 3
              dep1moli1(j,i)=deptemp(j,i)
              dep3bt(j,i)=dep3bt(j,i)+dep1moli1(j,i)
            end do
          end do
              do i=1,3
                do j=1,3
                 vir1moli1(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)+vir1moli1(j,i)
                end do
              end do

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,moli1rmndr,npole,
c!$OMP& nmol,ep1moli1,dep1moli1,vir1moli1)
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
           call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,deptemp,
     &             virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do i = 1, npole
                !i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,i)=deptemp(j,i)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,i)
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
            call empole1c_3b_Polar_totfield(npole3b,
     &           pnum,eptemp,deptemp,virtemp)
          ep1moli2=eptemp
          ep3bt=ep3bt-ep1moli2
          do i = 1, npole
            !i = pnum(l1)
         ! do i=1,npole
            do j = 1, 3
              dep1moli2(j,i)=deptemp(j,i)
              dep3bt(j,i)=dep3bt(j,i)-dep1moli2(j,i)
            end do
          end do
              do i=1,3
                do j=1,3
                 vir1moli2(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)-vir1moli2(j,i)
                end do
              end do

          ep3bt=ep3bt-ep1moli1
          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
              do i = 1, npole
                !i = pnum(l1)
             ! do i=1,npole
                do j = 1, 3
                   dep3bt(j,i)=dep3bt(j,i)-dep1moli1(j,i)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)-vir1moli1(j,i)
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

                 call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &           deptemp,virtemp)
                ep3bt=ep3bt+eptemp
                ep3moli123=eptemp
                do i =1,npole
                  !i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,i)
                end do
              do i=1,3
                do j=1,3
                 virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
                end do
              end do
            

                  npole3b=6
                   ep3bt=ep3bt-ep2moli12
                  ep3moli123=ep3moli123-ep2moli12

                   do i = 1, npole
                    !i = pnum(l1)
                  ! do i=1,npole
                    do j = 1, 3
                      dep3bt(j,i)=dep3bt(j,i)-dep2moli12(j,i)
                    end do
                   end do
                   do i=1,3
                     do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)-vir2moli12(j,i)
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
                  call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &                 deptemp,virtemp)
                  ep3bt=ep3bt-eptemp
                  ep3moli123=ep3moli123-eptemp
                 do i =1,npole
                  !i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)-deptemp(3,i)
                 end do
                   do i=1,3
                     do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)-virtemp(j,i)
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
                 call empole1c_3b_Polar_totfield(npole3b,pnum,eptemp,
     &                         deptemp,virtemp)
                   ep3moli123 = ep3moli123 - eptemp
                ep3bt=ep3bt-eptemp
                do i =1,npole
                  !i = pnum(l1)
               ! do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)-deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)-deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)-deptemp(3,i)
                end do
                   do i=1,3
                     do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)-virtemp(j,i)
                     end do
                   end do

                 np1=3
                 npole3b=3
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                ep3bt=ep3bt+ep1moli1
                  ep3moli123=ep3moli123+ep1moli1
                do i =1,npole
                  !i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli1(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli1(2,i)
                 dep3bt(3,i) = dep3bt(3,i)+dep1moli1(3,i)
                end do
                   do i=1,3
                     do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir1moli1(j,i)
                     end do
                   end do

                 npole3b=3
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 ep3bt=ep3bt+ep1moli2
                  ep3moli123=ep3moli123+ep1moli2
                do i =1,npole
                  !i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+dep1moli2(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+dep1moli2(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+dep1moli2(3,i)
                end do
                   do i=1,3
                     do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir1moli2(j,i)
                     end do
                   end do


                 npole3b=3
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 call empole1c_3b_Polar_totfield(
     &           npole3b,pnum,eptemp,deptemp,virtemp)
                ep3bt=ep3bt+eptemp
                  ep3moli123=ep3moli123+eptemp
                do i =1,npole
                  !i = pnum(l1)
                !do i=1,npole
                  dep3bt(1,i) = dep3bt(1,i)+deptemp(1,i)
                  dep3bt(2,i) = dep3bt(2,i)+deptemp(2,i)
                  dep3bt(3,i) = dep3bt(3,i)+deptemp(3,i)
                end do
c               print*,"3body Cobar ep3body",shellsum,ep3moli123
                   do i=1,3
                     do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+virtemp(j,i)
                     end do
                   end do


          end do
        end do 
c      end do  
c!$OMP END DO
c!$OMP END PARALLEL
c      end if
      print*,"Ewald Innerloop= ep3bt=",ep3bt
c       print*,"End inner1 single",ep3bt
      return
      end

