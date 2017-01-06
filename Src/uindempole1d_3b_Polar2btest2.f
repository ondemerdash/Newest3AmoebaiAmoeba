      subroutine uindempole1d_3b_Polar2btest2(npole3b,pnum,
     &   eptemp,deptemp,virtemp,cnt,uind)!,
c     &   emtemp,demtemp)
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polpot
      use polar, only: polarity, thole, pdamp
      use cho
      use pcg
      use neigh2clust
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
      integer npole3b
      integer pnum(*),l1,l3
      real*8 eptemp,deptemp(3,*),virtemp(3,3)
c      real*8 emtemp,demtemp(3,*),viremtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      integer cnt
      real*8 t1,t2,t3,t4,t5,t6,t7,t8

            t1=mpi_Wtime()

      eptemp = 0.0d0
c      emtemp = 0.0d0
c
c   Zero out temporary gradient of polarization energy
c

c      do l1 = 1, npole3b
      do i = 1, npole
        do j = 1, 3
          deptemp(j,i) = 0.0d0
c          demtemp(j,i) = 0.0d0
        end do
      end do

c
c     set the energy unit conversion factor
c
      f = electric / dielec

c
c   Zero out temporary virial
c
      do i=1,3
         do j=1,3
           virtemp(j,i)=0.0d0
c           viremtemp(j,i)=0.0d0
         end do
      end do

c      if(usecholesky) then
c       call induce0a_3b_cholesky_nopriordir(npole3b,pnum,
c     & uind,uinp) 
c      else if (usepcgscf) then
       call induce0a_3b_scf_vacomp2bclust(npole3b,pnum,uind,uinp,cnt)
c       print*,"Done w induce0a_3b_scf_vacomp2bclust"
c      else
c       call induce0a_3b_PolelecOnly_orig_nopriordir(npole3b,pnum,
c     & uind,uinp)
c      end if


      call emrecip1_3b_Polartest2(npole3b,pnum,uind,uinp,
     & eptemp,deptemp,virtemp)!,emtemp,demtemp)
c       call emrecip1_3b_PolarPermtest(npole3b,pnum,uind,uinp,
c     & eptemp,deptemp,virtemp,emtemp,demtemp) 
      !print*,"Done w emrecip1_3b_Polar_new in 2bclust"

      call ereal1d_3b_2bclustPolartest2(npole3b,pnum,uind,uinp,cnt,
     & eptemp,deptemp,virtemp)!,emtemp,demtemp)

c      print*,"Done w ereal1d_3b_3bclust in 3bclust"

      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      !do i = 1, npole
      do l1 =1, npole3b
         i=pnum(l1)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,l1)
         uiy = uind(2,l1)
         uiz = uind(3,l1)
c         cii = ci*ci
c         dii = dix*dix + diy*diy + diz*diz
c         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
c     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
c         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
c         emtemp = emtemp + e
         eptemp = eptemp + ei
      end do

      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      !do i = 1, npole
      do l1 = 1, npole3b
         i=pnum(l1)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,l1)+uinp(1,l1))
         uiy = 0.5d0 * (uind(2,l1)+uinp(2,l1))
         uiz = 0.5d0 * (uind(3,l1)+uinp(3,l1))
         trqi(1) = term * (diy*uiz-diz*uiy)
         trqi(2) = term * (diz*uix-dix*uiz)
         trqi(3) = term * (dix*uiy-diy*uix)
         !call torque (i,trq,trqi,frcx,frcy,frcz)
c        if(i.eq.0) then
c       print*,"Tilt!empole1d_PolarPerm3b b4 SelfTorq",taskid,cnt,npole3b
c        else 
c       print*,"In empole1d_PolarPerm3b b4 SelfTorque",taskid,cnt,npole3b
c        end if
c         call torque3_PolarPerm (npole3b,pnum,
c     &    i,trq,trqi,frcx,frcy,frcz,demtemp,deptemp)
         call torque3_Polar (
     &    i,trq,trqi,frcx,frcy,frcz,deptemp)
      end do
            t2=mpi_Wtime()
       !print*,"tsk=",taskid,"TotalCost Within empole1d_3b_3btest",t2-t1
c      print*,"Completed all steps of empole1d 3bclust"
      return
      end
c
