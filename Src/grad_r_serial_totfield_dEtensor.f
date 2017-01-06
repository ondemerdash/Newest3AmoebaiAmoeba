      subroutine gradient_emreal_serial_totfield_dEtensor
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      use totfield
      use math
      use ewald
      implicit none
      include 'mpif.h'
      real*8 emrealt,viremrealt(3,3)
      real*8, allocatable :: demrealt(:,:)
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6
      real*8, allocatable :: fieldnpolet(:,:)
      real*8, allocatable :: fieldpnpolet(:,:)
      real*8 ucell(3),term

      allocate (demrealt(3,npole))
      allocate (fieldnpolet(3,npole))
      allocate (fieldpnpolet(3,npole))
      
      emrealt=0.0d0
      do i=1,npole
         do j=1,3
            demrealt(j,i)=0.0d0
            fieldnpolet(j,i)=0.0d0
            fieldpnpolet(j,i) = 0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           viremrealt(j,i)=0.0d0
        end do
      end do
c      print*,"Before call to LoadBalInnerloop_ereal1d_3b_Perm_sentlist2"

      !if(taskid.lt.numtasks_emreal2) then
      !call LoadBal_ereal1d_3b_Perm_sentlist2_totfield_dEtensorOmp(
     &! emrealt,viremrealt,demrealt,fieldnpolet,fieldpnpolet)

      call LoadBal_ereal1d_3b_Perm_totfield_dEtensorOmp

       !print*,"PermElec Real time w/o MPIReduce,taskid",t4-t3,taskid
      !end if

       !call mpi_reduce(emrealt,emreal_tmp,1,mpi_real8,
     & !     mpi_sum,master,emreal_comm,ierr)
       emreal_tmp=emrealt
       !call mpi_reduce(demrealt,demreal_tmp,npole*3,mpi_real8,
     & !     mpi_sum,master,emreal_comm,ierr)
       !call mpi_reduce(viremrealt,viremreal_tmp,3*3,mpi_real8,
     & !     mpi_sum,master,emreal_comm,ierr)
       do i=1,3
          do j=1,3
             viremreal_tmp(j,i)=viremrealt(j,i)
          end do
       end do
       !call mpi_reduce(fieldnpolet,fieldnpole,npole*3,mpi_real8,
     & !     mpi_sum,master1,emreal_comm,ierr)
       !call mpi_reduce(fieldpnpolet,fieldpnpole,npole*3,mpi_real8,
     & !     mpi_sum,master1,emreal_comm,ierr)
       !print*,"PermElec Real time MPIReduceOnly,taskid",t2-t1,taskid

c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      !if(taskid.eq.master1) then
c
c     get the self-energy portion of the electrostatic field
c
       term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
       do i = 1, npole
         do j = 1, 3
            fieldnpole(j,i)=fieldnpolet(j,i)
            fieldpnpole(j,i)=fieldpnpolet(j,i)
            demreal_tmp(j,i)=demrealt(j,i)
            !fieldnpole(j,i) = fieldnpole(j,i) +fieldnpole_rcp(j,i)
            !fieldpnpole(j,i) = fieldpnpole(j,i) +fieldnpole_rcp(j,i)
            !fieldnpole(j,i) = fieldnpole(j,i) + term*rpole(j+1,i)
            !fieldpnpole(j,i) = fieldpnpole(j,i) + term*rpole(j+1,i)
         end do
       end do

       !if (boundary .eq. 'VACUUM') then
       !  do i = 1, 3
       !     ucell(i) = 0.0d0
       !  end do
       !  do i = 1, npole
       !     ii = ipole(i)
       !     ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
       !     ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
       !     ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
       !  end do
       !  term = (4.0d0/3.0d0) * pi/volbox
       !  do i = 1, npole
       !     do j = 1, 3
       !        fieldnpole(j,i) = fieldnpole(j,i) - term*ucell(j)
       !        fieldpnpole(j,i) = fieldpnpole(j,i) - term*ucell(j)
       !     end do
       !  end do
       !end if

      !end if

      deallocate(demrealt)
      deallocate(fieldnpolet)
      deallocate(fieldpnpolet)
      return
      end

