      subroutine fftsetup3b (qgrid3b,planf3b,planb3b,iprime3b,
     &  ffttable3b)
c      subroutine fftsetup3b
      use sizes
      use fft
      use openmp
c      use pme, only: nfft1, nfft2, nfft3
      use pme
      implicit none
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      integer*8 planf3b
      integer*8 planb3b
!$    integer ifront,iback
!$    integer error,iguess
      integer iprime3b(maxprime,3)
      real*8 ffttable3b(maxtable,3)
c
c
c     initialization of Fast Fourier transform using FFTW;
c     comment "dfftw_init_threads" and "dfftw_plan_with_threads"
c     calls if serial FFTW is used in place of OpenMP FFTW
c

!$    if (ffttyp .eq. 'FFTW') then
!$       ifront = -1
!$       iback = 1
!$       error = 0
!$       iguess = 0
!$       call dfftw_init_threads (error)
c         print*,"After dfftw_init_threads"
!$       call dfftw_plan_with_nthreads (nthread)
c          print*,"After dfftw_plan_with_nthreads"
!$       call dfftw_plan_dft_3d (planf3b,nfft1,nfft2,nfft3,qgrid3b,
!$   &                              qgrid3b,ifront,iguess)
c          print*,"After dfftw_plan_dft_3d 1"
!$       call dfftw_plan_dft_3d (planb3b,nfft1,nfft2,nfft3,qgrid3b,
!$   &                              qgrid3b,iback,iguess)
c          print*,"After dfftw_plan_dft_3d 2"
!$    else

c
c     initialization of Fast Fourier transform using FFTPACK
c

         call cffti (nfft1,ffttable3b(1,1),iprime3b(1,1))
         call cffti (nfft2,ffttable3b(1,2),iprime3b(1,2))
         call cffti (nfft3,ffttable3b(1,3),iprime3b(1,3))
!$    end if
      return
      end

      subroutine fftfront3b(planf3b,qgrid3b,iprime3b,
     &  ffttable3b)
      use sizes
      use fft
c      use pme, only: nfft1,nfft2,nfft3
      use pme
      implicit none
      integer i,j,k
      real*8, allocatable :: work(:,:)
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      integer*8 planf3b
      integer iprime3b(maxprime,3)
      real*8 ffttable3b(maxtable,3)
      
c
c
c     perform a single 3-D forward transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       call dfftw_execute_dft (planf3b,qgrid3b,qgrid3b)
!$    else
c
c     perform three 1-D forward transforms using FFTPACK
c
         allocate (work(2,max(nfft1,nfft2,nfft3)))
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  work(1,i) = qgrid3b(1,i,j,k)
                  work(2,i) = qgrid3b(2,i,j,k)
               end do
               call cfftf (nfft1,work,ffttable3b(1,1),iprime3b(1,1))
               do i = 1, nfft1
                  qgrid3b(1,i,j,k) = work(1,i)
                  qgrid3b(2,i,j,k) = work(2,i)
               end do
            end do
         end do
         do k = 1, nfft3
            do i = 1, nfft1
               do j = 1, nfft2
                  work(1,j) = qgrid3b(1,i,j,k)
                  work(2,j) = qgrid3b(2,i,j,k)
               end do
               call cfftf (nfft2,work,ffttable3b(1,2),iprime3b(1,2))
               do j = 1, nfft2
                  qgrid3b(1,i,j,k) = work(1,j)
                  qgrid3b(2,i,j,k) = work(2,j)
               end do
            end do
         end do
         do i = 1, nfft1
            do j = 1, nfft2
               do k = 1, nfft3
                  work(1,k) = qgrid3b(1,i,j,k)
                  work(2,k) = qgrid3b(2,i,j,k)
               end do
               call cfftf (nfft3,work,ffttable3b(1,3),iprime3b(1,3))
               do k = 1, nfft3
                  qgrid3b(1,i,j,k) = work(1,k)
                  qgrid3b(2,i,j,k) = work(2,k)
               end do
            end do
         end do
         deallocate (work)
!$    end if
      return
      end

      subroutine fftback3b(planb3b,qgrid3b,iprime3b,ffttable3b)
      use sizes
      use fft
c      use pme, only: nfft1,nfft2,nfft3
      use pme
      implicit none
      integer i,j,k
      real*8, allocatable :: work(:,:)
      real*8 qgrid3b(2,nfft1,nfft2,nfft3)
      integer*8 planb3b
      integer iprime3b(maxprime,3)
      real*8 ffttable3b(maxtable,3)

c
c
c     perform a single 3-D backward transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       call dfftw_execute_dft (planb3b,qgrid3b,qgrid3b)
!$    else
c
c     perform three 1-D backward transforms using FFTPACK
c
         allocate (work(2,max(nfft1,nfft2,nfft3)))
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  work(1,i) = qgrid3b(1,i,j,k)
                  work(2,i) = qgrid3b(2,i,j,k)
               end do
               call cfftb (nfft1,work,ffttable3b(1,1),iprime3b(1,1))
               do i = 1, nfft1
                  qgrid3b(1,i,j,k) = work(1,i)
                  qgrid3b(2,i,j,k) = work(2,i)
               end do
            end do
         end do
         do k = 1, nfft3
            do i = 1, nfft1
               do j = 1, nfft2
                  work(1,j) = qgrid3b(1,i,j,k)
                  work(2,j) = qgrid3b(2,i,j,k)
               end do
               call cfftb (nfft2,work,ffttable3b(1,2),iprime3b(1,2))
               do j = 1, nfft2
                  qgrid3b(1,i,j,k) = work(1,j)
                  qgrid3b(2,i,j,k) = work(2,j)
               end do
            end do
         end do
         do i = 1, nfft1
            do j = 1, nfft2
               do k = 1, nfft3
                  work(1,k) = qgrid3b(1,i,j,k)
                  work(2,k) = qgrid3b(2,i,j,k)
               end do
               call cfftb (nfft3,work,ffttable3b(1,3),iprime3b(1,3))
               do k = 1, nfft3
                  qgrid3b(1,i,j,k) = work(1,k)
                  qgrid3b(2,i,j,k) = work(2,k)
               end do
            end do
         end do
         deallocate (work)
!$    end if
      return
      end

