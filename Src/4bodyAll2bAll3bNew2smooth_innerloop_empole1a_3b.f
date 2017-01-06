
      subroutine All2bAll3bAll4bNew2SmoothInnerloop1_2cut_wskin(
     &  start,last,ep3bt,virep3bt,dep3bt,moli1rmndr)
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      implicit none
      real*8 ep2moli12,ep2moli13,ep2moli23
      real*8 dep2moli12(3,npole),dep2moli13(3,npole)
      real*8 dep2moli23(3,npole)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp,deptemp(3,npole)
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
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,npole)
      real*8 ep4moli1234,dep4moli1234(3,npole) 
      real*8 Cobarcut3b,Rcut2b,SwitchR2b,SwitchR3b
      integer moli1rmndr,l2,moli4
      real*8 vir2moli23(3,3),vir2moli13(3,3),vir4moli1234(3,3)
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

        print*,"4BODY NOCut4b NOCut3b NoCut2b"
        print*,"start last",start,last

!$OMP PARALLEL default(private) shared(imol,ep3bt,
!$OMP& virep3bt,dep3bt,start,last,
!$OMP& nmol)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
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

           call empole1a_mod(npole3b,pnum,eptemp,deptemp,virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
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

            
                 call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,i) = deptemp(1,i)
                    dep3moli123(2,i) = deptemp(2,i)
                    dep3moli123(3,i) = deptemp(3,i)
                  end do
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=virtemp(j,i)
                    end do
                  end do
            


               npole3b=6
               ep3moli123=ep3moli123-ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3moli123(j,i)=dep3moli123(j,i)-dep2moli12(j,i)
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
                  call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     ep2moli23=eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli23(j,i)=deptemp(j,i)
                       end do
                     end do
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       vir2moli23(j,i)=virtemp(j,i)
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
                 call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     ep2moli13=eptemp
                     do l1 = 1, np1
                      i = pnum(l1)
                      do j = 1, 3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli13(j,i)=deptemp(j,i)
                      end do
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1,3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli13(j,i)=deptemp(j,i)
                       end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       vir2moli13(j,i)=virtemp(j,i)
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
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,i)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
            
             do moli4=moli3+1,nmol
                npole3b=12
                pnum(1)=imol(1,moli1)
                pnum(2)=imol(1,moli1)+1
                pnum(3)=imol(2,moli1)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli3)
                pnum(8)=imol(1,moli3)+1
                pnum(9)=imol(2,moli3)
                pnum(10)=imol(1,moli4)
                pnum(11)=imol(1,moli4)+1
                pnum(12)=imol(2,moli4)
                 
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ep4moli1234=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep4moli1234(1,i) = deptemp(1,i)
                    dep4moli1234(2,i) = deptemp(2,i)
                    dep4moli1234(3,i) = deptemp(3,i)
                  end do
                  do i=1,3
                     do j=1,3
                        vir4moli1234(j,i)=virtemp(j,i)
                     end do
                  end do
                
                npole3b=9
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ep4moli1234=ep4moli1234-eptemp
                do l1=1,npole3b
                   i=pnum(l1)
                   do j=1,3
                    dep4moli1234(j,i)=dep4moli1234(j,i)
     &                              -deptemp(j,i)
                   end do
                end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do 

                npole3b=9
                pnum(1)=imol(1,moli1)
                pnum(2)=imol(1,moli1)+1
                pnum(3)=imol(2,moli1)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3             
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                       -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do


                npole3b=9
                pnum(1)=imol(1,moli1)
                pnum(2)=imol(1,moli1)+1
                pnum(3)=imol(2,moli1)
                pnum(4)=imol(1,moli3)
                pnum(5)=imol(1,moli3)+1
                pnum(6)=imol(2,moli3)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                    -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do

                npole3b=9
                pnum(1)=imol(1,moli2)
                pnum(2)=imol(1,moli2)+1
                pnum(3)=imol(2,moli2)
                pnum(4)=imol(1,moli3)
                pnum(5)=imol(1,moli3)+1
                pnum(6)=imol(2,moli3)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                      -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do

                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli2)
                 pnum(5)=imol(1,moli2)+1
                 pnum(6)=imol(2,moli2)
                 ep4moli1234=ep4moli1234+ep2moli12
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                      +dep2moli12(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli12(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 ep4moli1234=ep4moli1234+ep2moli13
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +dep2moli13(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli13(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 ep4moli1234=ep4moli1234+ep2moli23
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +dep2moli23(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli23(j,i)
                    end do
                 end do
 

                 npole3b=6
                 pnum(1)=imol(1,moli1)
                 pnum(2)=imol(1,moli1)+1
                 pnum(3)=imol(2,moli1)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                    +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                npole3b=12
                pnum(1)=imol(1,moli1)
                pnum(2)=imol(1,moli1)+1
                pnum(3)=imol(2,moli1)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli3)
                pnum(8)=imol(1,moli3)+1
                pnum(9)=imol(2,moli3)
                pnum(10)=imol(1,moli4)
                pnum(11)=imol(1,moli4)+1
                pnum(12)=imol(2,moli4)
                ep3bt = ep3bt + ep4moli1234
                do l1 = 1, npole3b
                    i = pnum(l1)
                   do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep4moli1234(j,i)
                   end do
                end do
                do i=1,3
                   do j=1,3
                      virep3bt(j,i)=virep3bt(j,i)+vir4moli1234(j,i)
                   end do
                end do

             end do
          end do
        end do 
      end do  
!$OMP END DO
!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      if(moli1rmndr.eq.0) then
        return
      else
        goto 31
      end if

   31 continue


c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,
c!$OMP& nmol)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
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

           call empole1a_mod(npole3b,pnum,eptemp,deptemp,virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
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

            
                 call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,i) = deptemp(1,i)
                    dep3moli123(2,i) = deptemp(2,i)
                    dep3moli123(3,i) = deptemp(3,i)
                  end do
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=virtemp(j,i)
                    end do
                  end do
            


               npole3b=6
               ep3moli123=ep3moli123-ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3moli123(j,i)=dep3moli123(j,i)-dep2moli12(j,i)
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
                  call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     ep2moli23=eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli23(j,i)=deptemp(j,i)
                       end do
                     end do
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       vir2moli23(j,i)=virtemp(j,i)
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
                 call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     ep2moli13=eptemp
                     do l1 = 1, np1
                      i = pnum(l1)
                      do j = 1, 3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli13(j,i)=deptemp(j,i)
                      end do
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1,3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli13(j,i)=deptemp(j,i)
                       end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       vir2moli13(j,i)=virtemp(j,i)
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
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,i)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
            
             do moli4=moli3+1,nmol
                npole3b=12
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli3)
                pnum(8)=imol(1,moli3)+1
                pnum(9)=imol(2,moli3)
                pnum(10)=imol(1,moli4)
                pnum(11)=imol(1,moli4)+1
                pnum(12)=imol(2,moli4)
                 
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ep4moli1234=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep4moli1234(1,i) = deptemp(1,i)
                    dep4moli1234(2,i) = deptemp(2,i)
                    dep4moli1234(3,i) = deptemp(3,i)
                  end do
                  do i=1,3
                     do j=1,3
                        vir4moli1234(j,i)=virtemp(j,i)
                     end do
                  end do
                
                npole3b=9
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ep4moli1234=ep4moli1234-eptemp
                do l1=1,npole3b
                   i=pnum(l1)
                   do j=1,3
                    dep4moli1234(j,i)=dep4moli1234(j,i)
     &                              -deptemp(j,i)
                   end do
                end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do 

                npole3b=9
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3             
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                       -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do


                npole3b=9
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli3)
                pnum(5)=imol(1,moli3)+1
                pnum(6)=imol(2,moli3)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                    -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do

                npole3b=9
                pnum(1)=imol(1,moli2)
                pnum(2)=imol(1,moli2)+1
                pnum(3)=imol(2,moli2)
                pnum(4)=imol(1,moli3)
                pnum(5)=imol(1,moli3)+1
                pnum(6)=imol(2,moli3)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                      -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do

                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli2)
                 pnum(5)=imol(1,moli2)+1
                 pnum(6)=imol(2,moli2)
                 ep4moli1234=ep4moli1234+ep2moli12
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                      +dep2moli12(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli12(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 ep4moli1234=ep4moli1234+ep2moli13
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +dep2moli13(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli13(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 ep4moli1234=ep4moli1234+ep2moli23
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +dep2moli23(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli23(j,i)
                    end do
                 end do
 

                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                    +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                npole3b=12
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli3)
                pnum(8)=imol(1,moli3)+1
                pnum(9)=imol(2,moli3)
                pnum(10)=imol(1,moli4)
                pnum(11)=imol(1,moli4)+1
                pnum(12)=imol(2,moli4)
                ep3bt = ep3bt + ep4moli1234
                do l1 = 1, npole3b
                    i = pnum(l1)
                   do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep4moli1234(j,i)
                   end do
                end do
                do i=1,3
                   do j=1,3
                      virep3bt(j,i)=virep3bt(j,i)+vir4moli1234(j,i)
                   end do
                end do

             end do
          end do
        end do 
c      end do  
c!$OMP END DO
c!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      return
      end
