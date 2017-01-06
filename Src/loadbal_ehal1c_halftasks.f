c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1c  --  buffered 14-7 vdw derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1c" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine loadbal_ehal1c_half(evtmp,devtmp,virevtmp)
      use sizes
      use atomid
      use atoms
      use bound
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mutant
      !use neigh
      use parvdwneigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      use cell
      use mpidat
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rho,rho6,rho7
      real*8 tau,tau7,scal
      real*8 s1,s2,t1,t2
      real*8 dt1drho,dt2drho
      real*8 dtau,gtau
      real*8 taper,dtaper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 rik6,rik7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 evt,eintert
      real*8 virt(3,3)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: devt(:,:)
      logical proceed,usei
      logical muti,mutk
      character*6 mode
      real*8 evtmp,devtmp(3,*),virevtmp(3,3)
      integer iter3
      integer taskid_offset
c
c
c     zero out the van der Waals energy and first derivatives
c
      !ev = 0.0d0
      !do i = 1, n
      !   dev(1,i) = 0.0d0
      !   dev(2,i) = 0.0d0
      !   dev(3,i) = 0.0d0
      !end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
      allocate (devt(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
      ! print*,"taskid",taskid,"off2=",off2
      ! print*,"taskid",taskid,"cut2=",cut2

c      print*,"Ehal1c xcell ycell zcell",xcell,ycell,zcell
c      print *,"Ehal1c xcell2 ycell2 zcell2",xcell2,ycell2,zcell2
      !print*,"taskid=",taskid,"start_vdw2",start_vdw2(taskid)
      !print*,"taskid=",taskid,"last_vdw2",last_vdw2(taskid) 
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer global to local copies for OpenMP calculation
c
      evt = 0.0d0 
      do i = 1, n
         devt(1,i) = 0.0d0 
         devt(2,i) = 0.0d0
         devt(3,i) = 0.0d0
      end do
      do i = 1, 3
         virt(1,i) = 0.0d0 
         virt(2,i) = 0.0d0
         virt(3,i) = 0.0d0
      end do
      taskid_offset=taskid-numtasks_emreal
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nvdw,ivdw,ired,kred,
!$OMP& jvdw,xred,yred,zred,use,nvlst_recv,vlst_recv,n12,n13,n14,n15,
!$OMP& i12,i13,i14,i15,v2scale,v3scale,v4scale,v5scale,
!$OMP& use_group,off2,radmin,epsilon,radmin4,epsilon4,ghal,dhal,
!$OMP& cut2,vlambda,scalpha,scexp,mut,c0,c1,c2,c3,c4,c5,
!$OMP& start_vdw2,last_vdw2,taskid_offset)
!$OMP& firstprivate(vscale,iv14) shared(evt,devt,virt)
!$OMP DO reduction(+:evt,devt,virt) schedule(guided)
c
c     find van der Waals energy and derivatives via neighbor list
c
      !do ii = 1, nvdw
      do ii =start_vdw2(taskid_offset),last_vdw2(taskid_offset)
         iter3=ii-start_vdw2(taskid_offset)+1
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         !do kk = 1, nvlst(ii)
         !   k = ivdw(vlst(kk,ii))
         do kk = 1, nvlst_recv(iter3)
            k = ivdw(vlst_recv(kk,iter3))
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
c
c     get the energy and gradient, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     rho6 = rho**6
                     rho7 = rho6 * rho
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     s1 = 1.0d0 / (scal+(rho+dhal)**7)
                     s2 = 1.0d0 / (scal+rho7+ghal)
                     t1 = (1.0d0+dhal)**7 * s1
                     t2 = (1.0d0+ghal) * s2
                     dt1drho = -7.0d0*(rho+dhal)**6 * t1 * s1
                     dt2drho = -7.0d0*rho6 * t2 * s2
                     e = eps * t1 * (t2-2.0d0)
                     de = eps * (dt1drho*(t2-2.0d0)+t1*dt2drho) / rv
                  else
                     rv7 = rv**7
                     rik6 = rik2**3
                     rik7 = rik6 * rik
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     tau7 = tau**7
                     dtau = tau / (dhal+1.0d0)
                     gtau = eps*tau7*rik6*(ghal+1.0d0)*(rv7/rho)**2
                     e = eps*tau7*rv7*((ghal+1.0d0)*rv7/rho-2.0d0)
                     de = -7.0d0 * (dtau*e+gtau)
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  evt = evt + e
                  if (i .eq. iv) then
                     devt(1,i) = devt(1,i) + dedx
                     devt(2,i) = devt(2,i) + dedy
                     devt(3,i) = devt(3,i) + dedz
                  else
                     devt(1,i) = devt(1,i) + dedx*redi
                     devt(2,i) = devt(2,i) + dedy*redi
                     devt(3,i) = devt(3,i) + dedz*redi
                     devt(1,iv) = devt(1,iv) + dedx*rediv
                     devt(2,iv) = devt(2,iv) + dedy*rediv
                     devt(3,iv) = devt(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     devt(1,k) = devt(1,k) - dedx
                     devt(2,k) = devt(2,k) - dedy
                     devt(3,k) = devt(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     devt(1,k) = devt(1,k) - dedx*redk
                     devt(2,k) = devt(2,k) - dedy*redk
                     devt(3,k) = devt(3,k) - dedz*redk
                     devt(1,kv) = devt(1,kv) - dedx*redkv
                     devt(2,kv) = devt(2,kv) - dedy*redkv
                     devt(3,kv) = devt(3,kv) - dedz*redkv
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  virt(1,1) = virt(1,1) + vxx
                  virt(2,1) = virt(2,1) + vyx
                  virt(3,1) = virt(3,1) + vzx
                  virt(1,2) = virt(1,2) + vyx
                  virt(2,2) = virt(2,2) + vyy
                  virt(3,2) = virt(3,2) + vzy
                  virt(1,3) = virt(1,3) + vzx
                  virt(2,3) = virt(2,3) + vzy
                  virt(3,3) = virt(3,3) + vzz
c
c     increment the total intermolecular energy
c
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      evtmp = evt
      do i = 1, n
         devtmp(1,i) = devt(1,i)
         devtmp(2,i) = devt(2,i)
         devtmp(3,i) = devt(3,i)
      end do
      do i = 1, 3
         virevtmp(1,i) = virt(1,i)
         virevtmp(2,i) = virt(2,i)
         virevtmp(3,i) = virt(3,i)
      end do
c      print*,"End of ehal1c ev",ev
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      deallocate (devt)
      return
      end
