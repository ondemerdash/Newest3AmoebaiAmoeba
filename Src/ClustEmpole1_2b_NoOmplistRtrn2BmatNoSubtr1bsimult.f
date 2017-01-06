
      subroutine ClustEmpole1_2b_NoOmplistRtrn2BmatNoSubtr1bsimult
     &  (ep3bt,virep3bt,dep3bt,ntript,em3bt,dem3bt,ep2btmat,dep2btmat,
     &   virep2btmat) 
      use sizes
      use atoms
      use atomid
      use molcul
      use mpole
      use neigh3b
      use mpidat
      use cobar
      use neigh2clust
      use cell
      use deriv1bmat
      implicit none
      real*8  ep2moli12,em2moli12
c      real*8  dep2moli12(3,36)
      real*8, allocatable :: dep2moli12(:,:)
      real*8, allocatable :: dem2moli12(:,:)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp,emtemp
      real*8, allocatable :: deptemp(:,:)
      real*8, allocatable :: demtemp(:,:)
      real*8 vir2moli12(3,3)
      real*8 virtemp(3,3)
      integer npole3b,moli1,moli2,moli3,np1,np2,np3
      integer, allocatable :: pnum(:)
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      real*8 em3bt,dem3bt(3,npole)
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
      integer clust1,clust2,clust3,npole3b12,k3
      real*8 ep1moli1,ep1moli2 
      real*8, allocatable :: dep1moli1(:,:) 
      real*8, allocatable :: dep1moli2(:,:)
      real*8 vir1moli1(3,3),vir1moli2(3,3)
      real*8 ep2moli12nosubtr 
      real*8, allocatable :: dep2moli12nosubtr(:,:)
      real*8 vir2moli12nosubtr(3,3)
      integer cnt,taskid_offset,ind
      real*8 ep2btmat(num2bsave)
      real*8 dep2btmat(3,2*3*maxsizeclust,num2bsave)
      real*8 virep2btmat(3,3,num2bsave)

      Rcut2b=cut2b_input
c      Rcut2b=7.0d0
      Rcut2b2=Rcut2b*Rcut2b

c        allocate(dep2moli12(3,3*3*maxsizeclust))
c        allocate(dem2moli12(3,3*3*maxsizeclust))
c        allocate(deptemp(3,3*3*maxsizeclust))

        allocate(dep2moli12(3,npole))
        allocate(dem2moli12(3,npole))

        allocate(pnum(2*3*maxsizeclust))
        taskid_offset=taskid-numtasks_polar1


c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start_polar,last_polar,
c!$OMP& nmollst,mollst,name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
c!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,molnew,taskid,
c!$OMP& Rcut2b2,sizeclust,clust,clust_cm,ep1bmat,dep1bmat,virep1bmat,
c!$OMP& start_mollst_polar,last_mollst_polar,num_mollst_chunk)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript)
c!$OMP& schedule(guided)
c      do kouter=start_polar(taskid),last_polar(taskid)
c        moli1=molnew(kouter)
c      do clust1=start_polar(taskid),last_polar(taskid)
      cnt=0
      do kouter=start_polar(taskid_offset),last_polar(taskid_offset)
        clust1=molnew(kouter)

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


c        do k1=1,nmollst(clust1)
c        do k1=start_mollst_polar(taskid),last_mollst_polar(taskid)
        do k3=1,num_mollst_chunk(clust1,taskid_offset)
         do k1=start_mollst_polar(k3,clust1,taskid_offset),
     &                        last_mollst_polar(k3,clust1,taskid_offset)
           clust2=mollst(k1,clust1)
           xr1=clust_cm(1,clust2)-clust_cm(1,clust1)
           yr1=clust_cm(2,clust2)-clust_cm(2,clust1)
           zr1=clust_cm(3,clust2)-clust_cm(3,clust1)
           call image(xr1,yr1,zr1)
           r1_2=xr1*xr1 + yr1*yr1 + zr1*zr1
           r1=sqrt(r1_2)
           cnt=cnt+1

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
c          cnt=cnt+1
          if(mlistclust) then
            !if(useboxclust) then
          !call empole1d_3b_PolarPerm2btestboxclust(npole3b,pnum,
     &    !     ep2moli12,dep2moli12,vir2moli12,cnt,em2moli12,dem2moli12)
            !else
          call empole1d_3b_PolarPerm2btest(npole3b,pnum,
     &         ep2moli12,dep2moli12,vir2moli12,cnt,em2moli12,dem2moli12)
            !end if  
          else
          call empole1c_3b_PolarPermtest(npole3b,pnum,
     &         ep2moli12,dep2moli12,vir2moli12,cnt,em2moli12,dem2moli12)
          end if

          npole3b=np1+np2

           
! NEED TO SMOOTHEN VIRIAL !

              ep3bt=ep3bt+ep2moli12
              em3bt=em3bt+em2moli12
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
c                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
c                  dem3bt(j,i) = dem3bt(j,i)+dem2moli12(j,l1)
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,i)
                  dem3bt(j,i) = dem3bt(j,i)+dem2moli12(j,i)

                end do
              end do
              do i=1,3
                do j=1,3
                  virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do

              if(save2bfor3b(k1,clust1)) then
                ind=ind2b_counter(clust2,clust1)
                ep2btmat(ind)=ep2moli12
                do l1 = 1, npole3b
                   i = pnum(l1)
                   do j = 1, 3
c                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,l1)
c                  dem3bt(j,i) = dem3bt(j,i)+dem2moli12(j,l1)
                     dep2btmat(j,l1,ind) = dep2moli12(j,i)
                   end do
                end do
                do i=1,3
                  do j=1,3
                  virep2btmat(j,i,ind)=vir2moli12(j,i)
                  end do
                end do
              end if

          end if
         end do
        end do 
      end do  
c!$OMP END DO
c!$OMP END PARALLEL
      !print*,"Energy,taskid,loadbalsmoothInner_1a_3bPolar",ep3bt,taskid
        deallocate(dep2moli12)
c        deallocate(deptemp)
        deallocate(pnum)
        deallocate(dem2moli12)
c        deallocate(dep3moli123)
c        deallocate(dep1moli1)
c        deallocate(dep1moli2)
c        deallocate(dep2moli12nosubtr)
      return
      end

