           allocate(elsttau1_index_rcv(4,largest_nelst,sz))
           allocate(elsttau2_index_rcv(4,largest_nelst,sz))

           call mpi_recv(elsttau1_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,26,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau2_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,28,mpi_comm_world,stat,ierr)

           do i=start_emreal2(k-2),last_emreal2(k-2)
              iter=i-start_emreal2(k-2)+1
              do k2=1,nelst2(i)
                 elstgrad(k2,i)=elstgrad_rcv(k2,iter)+ntpair_dEoffset
                do k3=1,4
                  elsttau1_index(k3,k2,i)=elsttau1_index_rcv(k3,k2,iter)
                  elsttau2_index(k3,k2,i)=elsttau2_index_rcv(k3,k2,iter)

                  if(elsttau1_index_rcv(k3,k2,iter).ne.0) then
                    elsttau1(k3,k2,i)=elsttau1_rcv(k3,k2,iter)
     &              +ntpair_tau1offset
                  else
                    elsttau1(k3,k2,i)=elsttau1_rcv(k3,k2,iter)
                  end if

                  if(elsttau2_index_rcv(k3,k2,iter).ne.0) then
                    elsttau2(k3,k2,i)=elsttau2_rcv(k3,k2,iter)
     &              +ntpair_tau2offset
                  else
                    elsttau2(k3,k2,i)=elsttau2_rcv(k3,k2,iter)
                  end if
                end do
              end do
           end do

           deallocate(elstgrad_rcv)
           deallocate(elsttau1_index_rcv)
           deallocate(elsttau1_rcv)
           deallocate(elsttau2_index_rcv)
           deallocate(elsttau2_rcv)

           ntpair_dEoffset=ntpair_dEoffset+m
           ntpair_tau1offset=ntpair_tau1offset+m1
           ntpair_tau2offset=ntpair_tau2offset+m2

