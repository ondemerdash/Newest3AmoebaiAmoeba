           allocate(taud1tmp_rcv(9,m1))
           allocate(taud2tmp_rcv(9,m2))
           allocate(taup1tmp_rcv(9,m1))
           allocate(taup2tmp_rcv(9,m2))

           call mpi_recv(taud1tmp_rcv,9*m1,mpi_real8,k,8,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(taud2tmp_rcv,9*m2,mpi_real8,k,9,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(taup1tmp_rcv,9*m1,mpi_real8,k,10,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(taup2tmp_rcv,9*m2,mpi_real8,k,11,
     &     mpi_comm_world,stat,ierr)

           do i=1,m1
              k2=ntpair_tau1offset+i
              do j=1,9
               taud1(j,k2)=taud1tmp_rcv(j,i)
               taup1(j,k2)=taup1tmp_rcv(j,i)
              end do
           end do

           do i=1,m2
              k2=ntpair_tau2offset+i
              do j=1,9
               taud2(j,k2)=taud2tmp_rcv(j,i)
               taup2(j,k2)=taup2tmp_rcv(j,i)
              end do
           end do

           deallocate(taud1tmp_rcv)
           deallocate(taud2tmp_rcv)
           deallocate(taup1tmp_rcv)
           deallocate(taup2tmp_rcv)


