      subroutine field_emreal_reduce_load_balanced_totfield_uindcalc
      use mpole
      use mpidat
      use deriv
      use energi
      use virial
      use totfield
      use math
      use ewald
      use polar
      use aprx
      implicit none
      include 'mpif.h'
      integer startem,lastem,atomind,i,j,ierr
      real*8 t1,t2,t3,t4,t5,t6
      real*8, allocatable :: fieldnpolet(:,:)
      real*8, allocatable :: fieldpnpolet(:,:)
      real*8 ucell(3),term

      if(uzepmedirpolz) then    
         allocate (fieldnpolet(3,npole))
         allocate (fieldpnpolet(3,npole))
      
         do i=1,npole
            do j=1,3
               fieldnpolet(j,i)=0.0d0
               fieldpnpolet(j,i) = 0.0d0
            end do
         end do

         if(taskid.lt.numtasks_emreal2) then
         call udirect2b_uinddir_mpi(fieldnpolet,fieldpnpolet) 
         end if


          call mpi_reduce(fieldnpolet,fieldnpole,npole*3,mpi_real8,
     &      mpi_sum,master1,mpi_comm_world,ierr)
          call mpi_reduce(fieldpnpolet,fieldpnpole,npole*3,mpi_real8,
     &      mpi_sum,master1,mpi_comm_world,ierr)
      else
         call dfield0b_totfield
      end if
       call mpi_barrier(mpi_comm_world,ierr)

      if(taskid.eq.master1) then
c
c     get the self-energy portion of the electrostatic field
c
       term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
        if(uzepmedirpolz) then
         do i = 1, npole
           do j = 1, 3
            fieldnpole(j,i) = fieldnpole(j,i) +fieldnpole_rcp(j,i)
            fieldpnpole(j,i) = fieldpnpole(j,i) +fieldnpole_rcp(j,i)
            fieldnpole(j,i) = fieldnpole(j,i) + term*rpole(j+1,i)
            fieldpnpole(j,i) = fieldpnpole(j,i) + term*rpole(j+1,i)

            uind(j,i)=polarity(i)*fieldnpole(j,i)
            uinp(j,i)=polarity(i)*fieldpnpole(j,i)
           end do
         end do
        else
         do i = 1, npole
           do j = 1, 3
            uind(j,i)=polarity(i)*fieldnpole(j,i)
            uinp(j,i)=polarity(i)*fieldpnpole(j,i)
           end do
         end do

        end if  
      end if

      deallocate(fieldnpolet)
      deallocate(fieldpnpolet)
      return
      end

