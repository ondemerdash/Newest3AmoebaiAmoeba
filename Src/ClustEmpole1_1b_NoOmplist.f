
      subroutine ClustEmpole1_1b_NoOmplist 
     &  (ep3bt,virep3bt,dep3bt,ntript,start,last,moli1rmndr,
     &  ep1btmat,dep1btmat,virep1btmat,em3bt,dem3bt,em1btmat,
     &  dem1btmat) 
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
      use pcg
      use polpot
      implicit none
      real*8  ep2moli12
      real*8, allocatable :: dep2moli12(:,:)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,start,last,ntript
      real*8 eptemp,emtemp
      real*8, allocatable :: deptemp(:,:)
      real*8, allocatable :: demtemp(:,:)
      real*8 vir2moli12(3,3)
      real*8 virtemp(3,3),viremtemp(3,3)
      integer npole3b,moli1,moli2,moli3,np1,np2,np3
      integer, allocatable :: pnum(:)
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      real*8 em3bt,dem3bt(3,npole),virem3bt(3,3)
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
      real*8 ep3moli123,vir3moli123(3,3)
c      real*8 dep3moli123(3,36) 
      real*8, allocatable :: dep3moli123(:,:)
      real*8 Rcut2b,Rcut2b2,SwitchR2b,SwitchR3b
      integer moli1rmndr,l2,kouter
      real*8 x2min,y2min,z2min,x3min,y3min,z3min
      real*8 xr3test,yr3test,zr3test,xr1t,yr1t,zr1t
      real*8 virtapr3b_xx,virtapr3b_xy,virtapr3b_xz
      real*8 virtapr3b_yy,virtapr3b_yz,virtapr3b_zz
      integer clust1,clust2,clust3,npole3b12
      real*8 ep1moli1,ep1moli2 
c      real*8 dep1moli1(3,36),dep1moli2(3,36)
      real*8, allocatable :: dep1moli1(:,:)
      real*8, allocatable :: dep1moli2(:,:)
      real*8 vir1moli1(3,3),vir1moli2(3,3)
      real*8 ep2moli12nosubtr
c      real*8 dep2moli12nosubtr(3,36)
      real*8, allocatable :: dep2moli12nosubtr(:,:)
      real*8 vir2moli12nosubtr(3,3)
      real*8 ep1btmat(clustcount) 
      real*8 dep1btmat(3,3*maxsizeclust,clustcount)
      real*8 virep1btmat(3,3,clustcount)
      integer cnt
      real*8 em1btmat(clustcount)
      real*8 dem1btmat(3,3*maxsizeclust,clustcount)
      real*8 virem1btmat(3,3,clustcount)

c        allocate(deptemp(3,3*maxsizeclust))
c        allocate(demtemp(3,3*maxsizeclust))
c        allocate(pnum(3*maxsizeclust))

        allocate(deptemp(3,npole))
        allocate(demtemp(3,npole))
        allocate(pnum(3*maxsizeclust))


c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,
c!$OMP& name,x,y,z,rtapr2b,rtapr3b,Cobarcut3b,
c!$OMP& Rcut2b,SwitchR2b,SwitchR3b,ntript,taskid,
c!$OMP& Rcut2b2,sizeclust,clust,clust_cm,clustcount,
c!$OMP& ep1btmat,dep1btmat,virep1btmat)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt,ntript,
c!$OMP& ep1btmat,dep1btmat,virep1btmat)
c!$OMP& schedule(guided)
      cnt=0
      do clust1=start,last

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
           cnt=cnt+1
c          call empole1d_3b_PolarPerm1b(npole3b,pnum,
c     &         eptemp,deptemp,virtemp,cnt,emtemp,demtemp)
          if(mlistclust) then
           if(useboxclust) then
          !call empole1d_3b_PolarPerm1btestboxclust(npole3b,pnum,
     &    !     eptemp,deptemp,virtemp,cnt,emtemp,demtemp) 
           else
          call empole1d_3b_PolarPerm1btest(npole3b,pnum,
     &         eptemp,deptemp,virtemp,cnt,emtemp,demtemp)
           end if
          else
          call empole1c_3b_PolarPermtest(npole3b,pnum,
     &         eptemp,deptemp,virtemp,cnt,emtemp,demtemp)
          end if
          ep3bt = ep3bt + eptemp
          em3bt = em3bt + emtemp
          ep1btmat(clust1)=eptemp
          em1btmat(clust1)=emtemp
          do l1 = 1, npole3b
             i = pnum(l1)
             do j = 1, 3
c              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
c              dep1btmat(j,l1,clust1)=deptemp(j,l1)
c
c              dem3bt(j,i) = dem3bt(j,i)+demtemp(j,l1)
c              dem1btmat(j,l1,clust1)=demtemp(j,l1)
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              dep1btmat(j,l1,clust1)=deptemp(j,i)

              dem3bt(j,i) = dem3bt(j,i)+demtemp(j,i)
              dem1btmat(j,l1,clust1)=demtemp(j,i)

             end do            
          end do
          do i=1,3
            do j=1,3
              virep3bt(j,i) = virep3bt(j,i)+virtemp(j,i)
              virep1btmat(j,i,clust1)=virtemp(j,i)
            end do
          end do


      end do  
c!$OMP END DO
c!$OMP END PARALLEL
      !print*,"Energy,taskid,loadbalsmoothInner_1a_3bPolar",ep3bt,taskid
c       print*,"Before rmndr In Clst1Btskid=",taskid,"ep3bt=",ep3bt

      if(moli1rmndr.eq.0) then
c        deallocate(dep2moli12)
        deallocate(deptemp)
        deallocate(demtemp)
        deallocate(pnum)
c        deallocate(dep3moli123)
c        deallocate(dep1moli1)
c        deallocate(dep1moli2)
c        deallocate(dep2moli12nosubtr)
        return
      else
        goto 31
      end if

   31 continue
           npole3b=0
           np1=0
           do i=1,sizeclust(moli1rmndr)
             moli1=clust(i,moli1rmndr)
             np1=np1+3
             pnum(3*(i-1)+1)=imol(1,moli1)
             pnum(3*(i-1)+2)=imol(1,moli1)+1
             pnum(3*(i-1)+3)=imol(2,moli1)
           end do
           npole3b=np1
           cnt=cnt+1
          if(mlistclust) then
           if(useboxclust) then
          !call empole1d_3b_PolarPerm1btestboxclust(npole3b,pnum,
     &    !     eptemp,deptemp,virtemp,cnt,emtemp,demtemp)
           else
          call empole1d_3b_PolarPerm1btest(npole3b,pnum,
     &         eptemp,deptemp,virtemp,cnt,emtemp,demtemp)
           end if
          else
          call empole1c_3b_PolarPermtest(npole3b,pnum,
     &         eptemp,deptemp,virtemp,cnt,emtemp,demtemp)
          end if
          ep3bt = ep3bt + eptemp
          em3bt = em3bt + emtemp
          ep1btmat(moli1rmndr)=eptemp
          em1btmat(moli1rmndr)=emtemp
          do l1 = 1, npole3b
             i = pnum(l1)
             do j = 1, 3
c              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
c              dep1btmat(j,l1,clust1)=deptemp(j,l1)
c
c              dem3bt(j,i) = dem3bt(j,i)+demtemp(j,l1)
c              dem1btmat(j,l1,clust1)=demtemp(j,l1)
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              dep1btmat(j,l1,moli1rmndr)=deptemp(j,i)

              dem3bt(j,i) = dem3bt(j,i)+demtemp(j,i)
              dem1btmat(j,l1,moli1rmndr)=demtemp(j,i)

             end do
          end do
          do i=1,3
            do j=1,3
              virep3bt(j,i) = virep3bt(j,i)+virtemp(j,i)
              virep1btmat(j,i,moli1rmndr)=virtemp(j,i)
            end do
          end do

c        deallocate(dep2moli12)
        deallocate(deptemp)
        deallocate(demtemp)
        deallocate(pnum)
c        deallocate(dep3moli123)
c        deallocate(dep1moli1)
c        deallocate(dep1moli2)
c        deallocate(dep2moli12nosubtr)
c       print*,"In Clst1Btskid=",taskid,"ep3bt=",ep3bt

        return
        end
