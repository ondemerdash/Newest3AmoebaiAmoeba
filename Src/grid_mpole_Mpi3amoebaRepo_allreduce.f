c
      subroutine grid_mpole_Mpi3amoebaRepo_allreduce(fmp)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
      use mpidat
      implicit none
      include 'mpif.h'
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 v2,u2,t2
      real*8 term0,term1,term2
      real*8 fmp(10,*)
      integer offset, start, finish
      integer new_iter
      real*8 qgrid_temp(2,nfft1,nfft2,nfft3)
      integer ierr
c
c
c     zero out the particle mesh Ewald charge grid
c

      !print*, "hi from grid_mpole from id", taskid

      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do

      qgrid_temp = qgrid

      offset = npole*nchunk/numtasks
      start = offset*taskid
      if(taskid .eq. (numtasks-1)) then
         finish = npole*nchunk-1
      else
         finish=offset*(taskid+1)-1
      end if

      !print*, "offset, start, finish, taskid", offset, start,
     & !    finish, taskid
      
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,v2,u0,u1,u2,term0,term1,term2,t0,t1,t2,
!$OMP& new_iter)
!$OMP DO
c
c     put the permanent multipole moments onto the grid
c
c      do ichk = 1, nchunk

c      do new_iter=0, npole*nchunk 
      do new_iter=start, finish
         
         isite = mod(new_iter,npole) + 1
         ichk = new_iter/npole + 1

         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
c         do isite = 1, npole
            iatm = ipole(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  v2 = thetai3(3,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
                     u2 = thetai2(3,m,iatm)
                     term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
     &                     + fmp(4,isite)*u0*v1 + fmp(6,isite)*u2*v0
     &                     + fmp(7,isite)*u0*v2 + fmp(10,isite)*u1*v1
                     term1 = fmp(2,isite)*u0*v0 + fmp(8,isite)*u1*v0
     &                          + fmp(9,isite)*u0*v1
                     term2 = fmp(5,isite) * u0 * v0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
                        t2 = thetai1(3,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + 
     &                   term0*t0 + term1*t1 + term2*t2
                     end do
                  end do
               end do
            end if
         end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
         
         qgrid_temp=qgrid
c         qgrid = 0 
         
         !print*, "hi before allreduce from id", taskid
         call mpi_allreduce(qgrid_temp, qgrid, 2*nfft1*nfft2*nfft3, 
     &        MPI_REAL8,MPI_sum, mpi_comm_world, ierr)    
         !print*, "hi after allreduce from id", taskid
         
      return
      end

