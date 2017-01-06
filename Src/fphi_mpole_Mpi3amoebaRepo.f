
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_mpole  --  multipole potential from grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_mpole" extracts the permanent multipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_mpole_Mpi3amoebaRepo (fphi_totfield)
      use sizes
      use mpole
      use pme
      use mpidat
      use chunks
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fphi_totfield(20,*)
      integer offset, start, finish
      integer new_iter
      integer lstart,lend
      integer rmndr

      call splitloop(lstart, lend, npole)
      !offset=int(npole/numtasks)
      !rmndr=mod(npole,numtasks)
      !if(taskid.ne.numtasks-1) then
      !  lstart=offset*taskid+1
      !  lend=lstart+offset-1
      !else
      !  lstart=offset*taskid+1
      !  lend=lstart+offset-1+rmndr
      !end if

      do i=1,npole
         do j=1,20
         fphi_totfield(j,i) = 0.0d0
         end do
      end do 

      !print*, "offset, start, finish, taskid", offset, lstart,
     & !    lend, taskid
       print*,"taskid",taskid,"lstart=",lstart,"lend=",lend
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,igrid,bsorder,
!$OMP& nfft3,thetai3,nfft2,thetai2,nfft1,thetai1,qgrid,
!$OMP& fphi_totfield,start,finish,lstart,lend)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
c      do isite = 1, npole
      do isite = lstart, lend
      !do new_iter=start, finish
c         iatm = ipole(isite)
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)
               t0 = 0.0d0
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq = qgrid(1,i,j,k)
                  t0 = t0 + tq*thetai1(1,it1,iatm)
                  t1 = t1 + tq*thetai1(2,it1,iatm)
                  t2 = t2 + tq*thetai1(3,it1,iatm)
                  t3 = t3 + tq*thetai1(4,it1,iatm)
               end do
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fphi_totfield(1,isite) = tuv000
         fphi_totfield(2,isite) = tuv100
         fphi_totfield(3,isite) = tuv010
         fphi_totfield(4,isite) = tuv001
         fphi_totfield(5,isite) = tuv200
         fphi_totfield(6,isite) = tuv020
         fphi_totfield(7,isite) = tuv002
         fphi_totfield(8,isite) = tuv110
         fphi_totfield(9,isite) = tuv101
         fphi_totfield(10,isite) = tuv011
         fphi_totfield(11,isite) = tuv300
         fphi_totfield(12,isite) = tuv030
         fphi_totfield(13,isite) = tuv003
         fphi_totfield(14,isite) = tuv210
         fphi_totfield(15,isite) = tuv201
         fphi_totfield(16,isite) = tuv120
         fphi_totfield(17,isite) = tuv021
         fphi_totfield(18,isite) = tuv102
         fphi_totfield(19,isite) = tuv012
         fphi_totfield(20,isite) = tuv111
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      !print*,"successful completion of MPI fphi_mpole!"
      return
      end
