      subroutine grad_covemrecip_vandr_noreduce_totfield
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      use totfield
      use math
      use ewald
      use dEtensor
      use sizes
      use polar
      use aprx
      use deriv3b
      use chgpot   
      use paremneigh
      use neigh
      use atoms
      use vdw
      use pme
      use chunks
      implicit none
      include 'mpif.h'
      real*8 emrealt,viremrealt(3,3)
      real*8, allocatable :: demrealt(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6
      real*8, allocatable :: fieldnpolet(:,:)
      real*8, allocatable :: fieldpnpolet(:,:)
      real*8 ucell(3),term,f
      integer k,maxsize_master,maxoffset_master,offset_cum
      integer stat(MPI_STATUS_SIZE),sz,maxlocal_allthread
      integer maxlocal_allthread_tau
      real*8 evtmp,virevtmp(3,3)
      real*8, allocatable :: devtmp(:,:)
      real*8, allocatable :: dep3bt1(:,:)
      real*8 virep3bt1(3,3)
      real*8, allocatable :: uindt(:,:)
      real*8, allocatable :: uinpt(:,:)
      logical first
      save first
      data first  / .true. /


      f = electric / dielec
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"
             !t1=mpi_Wtime()
      !print*,"numtasks_emreal2=",numtasks_emreal2

      if(taskid.eq.master) then
         call gradient_covalent2

      else if(taskid.eq.master1) then
      em = 0.0d0

      if (first) then
       first = .false.
       if (.not. allocated(dem))  allocate (dem(3,n))
      end if

      do i = 1, n
         do j=1,3
            dem(j,i) =0.0d0
         end do
      end do
      do i = 1, 3
         do j=1,3
            viremrecip(j,i) =0.0d0
         end do
      end do

        if(embedtyp.ne.'O') then
        call emrecip1_3b_Perm2_totfield
        else
        call emrecip1_3b_Perm2
        end if
c      call mpi_isend(em,1,mpi_real8,master,1,mpi_comm_world,
c     &      reqs(1),ierr)
c      call mpi_isend(dem,3*n,mpi_real8,master,2,mpi_comm_world,
c     &      reqs(2),ierr)
c      call mpi_isend(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
c     &      reqs(3),ierr)

c       call mpi_send(em,1,mpi_real8,master,1,mpi_comm_world,
c     &      ierr)
c      call mpi_send(dem,3*n,mpi_real8,master,2,mpi_comm_world,
c     &      ierr)
c      call mpi_send(viremrecip,3*3,mpi_real8,master,3,mpi_comm_world,
c     &      ierr)
 
c      if(taskid.lt.numtasks_emreal2) then
      !else if(taskid.eq.master2) then
       !call dfield0b_totfield
      else if(taskid.eq.2) then

        call ereal1d_3b_Perm2

          if((taskid.eq.master2).and.(embedtyp.ne.'O')) then
            call dfield0b_totfield
          end if
       call mpi_send(emreal_tmp,1,mpi_real8,master,1,mpi_comm_world,
     &      ierr)
      call mpi_send(demreal_tmp,3*n,mpi_real8,master,2,mpi_comm_world,
     &      ierr)
      call mpi_send(viremreal_tmp,3*3,mpi_real8,master,3,mpi_comm_world,
     &      ierr)

       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      else if(taskid.eq.3) then
        ev = 0.0d0
        if (.not. allocated(dev))  allocate (dev(3,n))

        do i = 1, n
          do j=1,3
             dev(j,i) = 0.0d0
          end do
        end do
        do i = 1, 3
          do j=1,3
            vir(j,i) =0.0d0
            virev(j,i)=0.0d0
          end do
        end do

         call ehal1

      do i = 1, 3
         do j=1,3
            virev(j,i)=vir(j,i)
         end do
      end do
       call mpi_send(ev,1,mpi_real8,master,4,mpi_comm_world,
     &      ierr)
      call mpi_send(dev,3*n,mpi_real8,master,5,mpi_comm_world,
     &      ierr)
      call mpi_send(virev,3*3,mpi_real8,master,6,mpi_comm_world,
     &      ierr)
         
      end if


c      call mpi_barrier(mpi_comm_world,ierr)

       if(embedtyp.ne.'O') then
         call mpi_bcast(fieldnpole,3*n,mpi_real8,master2,
     &      mpi_comm_world,ierr)
         call mpi_bcast(fieldpnpole,3*n,mpi_real8,master2,
     &      mpi_comm_world,ierr)
       end if
      return
      end

