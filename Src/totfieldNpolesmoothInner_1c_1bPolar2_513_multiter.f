c
      subroutine totfieldNpolesmoothInner_1c_1bPolar2_513_multiter(
     &  start,last,start2,last2,ep3bt,virep3bt,dep3bt,moli1rmndr)
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use aprx
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
      real*8 ep3moli123,vir3moli123(3,3),dep3moli123(3,npole) 
      real*8 Cobarcut3b,Rcut2b,SwitchR2b,SwitchR3b,Rcut2b2
      integer moli1rmndr,l2
      real*8 term,term2,term3,term4
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      integer, allocatable :: embedlst(:,:)
      integer, allocatable :: nembedlst(:)
      real*8 ewaldcut2
      integer start2,last2,length
      !allocate (nembedlst(itermode*3))
      !allocate (embedlst(maxelst,itermode*3))

      rtapr2b=1.0d12
      Rcut2b=1.0d12
      rtapr2b=6.0d0
      Rcut2b=8.0d0

      Rcut2b2=Rcut2b*Rcut2b
      SwitchR2b=Rcut2b-rtapr2b

      rtapr3b=6.0d0
      Cobarcut3b=8.0d0
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

        print*,"TOTFIELDEWALD 1body"
      print*,"task start last moli1rmndr",taskid,start,last,moli1rmndr
      print*,"task start2 last2",taskid,start2,last2
      print*,"itermode=",itermode

c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,npole,
c!$OMP& nmol)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
c!$OMP& schedule(guided)
      do moli1=start,last,itermode
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)

          do k=0,itermode-1
           pnum(3*k+1)=imol(1,moli1+k)
           pnum(3*k+2)=imol(1,moli1+k)+1
           pnum(3*k+3)=imol(2,moli1+k)
          end do

          !pnum(1)=imol(1,moli1)
          !pnum(2)=imol(1,moli1)+1
          !pnum(3)=imol(2,moli1)

          !pnum(4)=imol(1,moli1+1)
          !pnum(5)=imol(1,moli1+1)+1
          !pnum(6)=imol(2,moli1+1)

          !pnum(7)=imol(1,moli1+2)
          !pnum(8)=imol(1,moli1+2)+1
          !pnum(9)=imol(2,moli1+2)

          !pnum(10)=imol(1,moli1+3)
          !pnum(11)=imol(1,moli1+3)+1
          !pnum(12)=imol(2,moli1+3)

          !pnum(13)=imol(1,moli1+4)
          !pnum(14)=imol(1,moli1+4)+1
          !pnum(15)=imol(2,moli1+4)

          !npole3b=15
          npole3b=3*itermode
          !do l1 = 1, npole3b
          !   i=pnum(l1)
          !   nembedlst(l1) = 0
          !   do k = 1,npole
          !     xr = x(k) - x(i)
          !     yr = y(k) - y(i)
          !     zr = z(k) - z(i)
          !     call imagen (xr,yr,zr)
          !     r2 = xr*xr + yr*yr + zr*zr
          !     if ((r2 .le. ewaldcut2).and.(i.ne.k)) then
          !      nembedlst(l1) = nembedlst(l1)+1
          !      embedlst(nembedlst(l1),l1) = k
          !     end if
          !   end do
          !end do

           !call empole1c_3b_Polar_totfield_513(embedlst,nembedlst,
     &     !     npole3b,pnum,ep1moli1,dep1moli1,vir1moli1)
           call empole1c_3b_Polar_totfield(npole3b,pnum,ep1moli1,
     &          dep1moli1,vir1moli1)
          ep3bt=ep3bt+ep1moli1
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
      end do  
c!$OMP END DO
c!$OMP END PARALLEL

      print*,"Ewald 1body Innerloop= ep3bt=",ep3bt

      length=last2-start2+1

      !moli1=start2

      do k=0,length-1
         pnum(3*k+1)=imol(1,start2+k)
         pnum(3*k+2)=imol(1,start2+k)+1
         pnum(3*k+3)=imol(2,start2+k)
      end do
      npole3b=3*length
          !do l1 = 1, npole3b
          !   i=pnum(l1)
          !   nembedlst(l1) = 0
          !   do k = 1,npole
          !     xr = x(k) - x(i)
          !     yr = y(k) - y(i)
          !     zr = z(k) - z(i)
          !     call imagen (xr,yr,zr)
          !     r2 = xr*xr + yr*yr + zr*zr
          !     if ((r2 .le. ewaldcut2).and.(i.ne.k)) then
          !      nembedlst(l1) = nembedlst(l1)+1
          !      embedlst(nembedlst(l1),l1) = k
          !     end if
          !   end do
          !end do

           !call empole1c_3b_Polar_totfield_513(embedlst,nembedlst,
     &     !     npole3b,pnum,ep1moli1,dep1moli1,vir1moli1)
           call empole1c_3b_Polar_totfield(npole3b,pnum,ep1moli1,
     &          dep1moli1,vir1moli1)
          ep3bt=ep3bt+ep1moli1
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
          !do l1 = 1, npole3b
          !   i=pnum(l1)
          !   nembedlst(l1) = 0
          !   do k = 1,npole
          !     xr = x(k) - x(i)
          !     yr = y(k) - y(i)
          !     zr = z(k) - z(i)
          !     call imagen (xr,yr,zr)
          !     r2 = xr*xr + yr*yr + zr*zr
          !     if ((r2 .le. ewaldcut2).and.(i.ne.k)) then
          !      nembedlst(l1) = nembedlst(l1)+1
          !      embedlst(nembedlst(l1),l1) = k
          !     end if
          !   end do
          !end do

           !call empole1c_3b_Polar_totfield_513(embedlst,nembedlst,
     &     !     npole3b,pnum,ep1moli1,dep1moli1,vir1moli1)
           call empole1c_3b_Polar_totfield(npole3b,pnum,ep1moli1,
     &          dep1moli1,vir1moli1)

          ep3bt=ep3bt+ep1moli1
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

      print*,"Ewald 1body Innerloop= ep3bt=",ep3bt
      !deallocate(nembedlst)
      !deallocate(embedlst)
      return
      end

