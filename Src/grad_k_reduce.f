      subroutine gradient_emrecip_reduce
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      implicit none
      include 'mpif.h'
      real*8 emrecipt,viremrecipt(3,3)
      real*8, allocatable :: demrecipt(:,:)
      integer startem,lastem,atomind,i,j,ierr

      allocate (demrecipt(3,npole))

      emrecipt=0.0d0
      do i=1,npole
         do j=1,3
            demrecipt(j,i)=0.0d0
         end do 
      end do  
      do i=1,3
        do j=1,3
           viremrecipt(j,i)=0.0d0
        end do  
      end do 

c       call Innerloop_ereal1d_3b_Perm_fulllist2(atomind,startem,
c     &      lastem,emrealt,viremrealt,demrealt)

       if(taskid.eq.master1) then

         call emrecip1_3b_Perm2_reduce(emrecipt,demrecipt,viremrecipt)
       end if
       call mpi_reduce(emrecipt,em,1,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(demrecipt,dem,npole*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
       call mpi_reduce(viremrecipt,viremrecip,3*3,mpi_real8,
     &      mpi_sum,master,mpi_comm_world,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(demrecipt)
      return
      end

      subroutine gradient_emrecip_reduce_emrecipcomm
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      implicit none
      include 'mpif.h'
      real*8 emrecipt,viremrecipt(3,3)
      real*8, allocatable :: demrecipt(:,:)
      integer startem,lastem,atomind,i,j,ierr

      allocate (demrecipt(3,npole))

      emrecipt=0.0d0
      do i=1,npole
         do j=1,3
            demrecipt(j,i)=0.0d0
         end do
      end do
      do i=1,3
        do j=1,3
           viremrecipt(j,i)=0.0d0
        end do
      end do

c       call Innerloop_ereal1d_3b_Perm_fulllist2(atomind,startem,
c     &      lastem,emrealt,viremrealt,demrealt)

       if(taskid.eq.master1) then

         call emrecip1_3b_Perm2_reduce(emrecipt,demrecipt,viremrecipt)
       end if
       call mpi_reduce(emrecipt,em,1,mpi_real8,
     &      mpi_sum,master,emrecip_comm,ierr)
       call mpi_reduce(demrecipt,dem,npole*3,mpi_real8,
     &      mpi_sum,master,emrecip_comm,ierr)
       call mpi_reduce(viremrecipt,viremrecip,3*3,mpi_real8,
     &      mpi_sum,master,emrecip_comm,ierr)
c       print*,"In gradient_polar After mpi_reduce ",ep3bt_tot

      deallocate(demrecipt)
      return
      end

