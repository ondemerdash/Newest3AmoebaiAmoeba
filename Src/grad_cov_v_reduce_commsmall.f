      subroutine grad_cov_v_reduce_commsmall 
      use mpole
      use mpidat
      use mpidat2
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
      use mpidat4
      implicit none
      include 'mpif.h'
      real*8 emrealt,viremrealt(3,3)
      real*8, allocatable :: demrealt(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6
      real*8 t11,t21,t31,t41,t51,t61
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

      if((taskid.eq.master).or.((taskid.ge.numtasks_emreal+2)
     &    .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) ) then

        allocate (devtmp(3,nvdw))
        evtmp=0.0d0
        do i=1,n
         do j=1,3
            devtmp(j,i)=0.0d0
         end do
        end do
        do i=1,3
         do j=1,3
           virevtmp(j,i)=0.0d0
         end do
        end do
      end if
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"
             !t1=mpi_Wtime()
      !print*,"numtasks_emreal2=",numtasks_emreal2

      if(taskid.eq.master) then
         call gradient_covalent2

      else if((taskid.ge.numtasks_emreal+2)
     &    .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) then

            call loadbal_ehal1c_half_nolistsend(
     &        evtmp,devtmp,virevtmp)
      end if

      if((taskid.eq.master).or.((taskid.ge.numtasks_emreal+2)
     &  .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) ) then
       call mpi_ireduce(evtmp,ev,1,mpi_real8,
     &      mpi_sum,master,vdw_comm,reqs16,ierr)
       call mpi_ireduce(devtmp,dev,nvdw*3,mpi_real8,
     &      mpi_sum,master,vdw_comm,reqs17,ierr)
       call mpi_ireduce(virevtmp,virev,3*3,mpi_real8,
     &      mpi_sum,master,vdw_comm,reqs18,ierr)
      end if
c       call mpi_reduce(evtmp,ev,1,mpi_real8,
c     &      mpi_sum,master,mpi_comm_world,ierr)
c       call mpi_reduce(devtmp,dev,nvdw*3,mpi_real8,
c     &      mpi_sum,master,mpi_comm_world,ierr)
c       call mpi_reduce(virevtmp,virev,3*3,mpi_real8,
c     &      mpi_sum,master,mpi_comm_world,ierr)
      if((taskid.eq.master).or.((taskid.ge.numtasks_emreal+2)
     &   .and.(taskid.lt.(numtasks_vdw2+numtasks_emreal+2))) ) then

        deallocate(devtmp)
      end if

      return
      end

