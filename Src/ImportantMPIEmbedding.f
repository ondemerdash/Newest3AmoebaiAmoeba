         largest_nelst2=0

         do i=1,npole
           if(nelst2(i).gt.largest_nelst2) then
             largest_nelst2=nelst2(i)
           end if
         end do

      if(taskid.eq.numtasks-1) then


         ntpair_dEoffset=0
         ntpair_tau1offset=0
         ntpair_tau2offset=0

         allocate(ntpair_dE_array(0:numtasks_emreal2-1))
         allocate(ntpair_tau1_array(0:numtasks_emreal2-1))
         allocate(ntpair_tau2_array(0:numtasks_emreal2-1))

         do k=2,numtasks_emreal2+1
           call mpi_recv(ntpair_dEtmp_rcv,1,mpi_integer,k,1,
     &     mpi_comm_world,stat,ierr)
           ntpair_dEoffset=ntpair_dEoffset+ntpair_dEtmp_rcv
           ntpair_dE_array(k-2)=ntpair_dEtmp_rcv

           call mpi_recv(ntpair_tau1tmp_rcv,1,mpi_integer,k,6,
     &     mpi_comm_world,stat,ierr)
           ntpair_tau1offset=ntpair_tau1offset+ntpair_tau1tmp_rcv
           ntpair_tau1_array(k-2)=ntpair_tau1tmp_rcv

           call mpi_recv(ntpair_tau2tmp_rcv,1,mpi_integer,k,7,
     &     mpi_comm_world,stat,ierr)
           ntpair_tau2offset=ntpair_tau2offset+ntpair_tau2tmp_rcv
           ntpair_tau2_array(k-2)=ntpair_tau2tmp_rcv
         end do

         ntpair_dE=ntpair_dEoffset
         ntpair_tau1=ntpair_tau1offset
         ntpair_tau2=ntpair_tau2offset

         if(allocated(elst_grad)) deallocate(elst_grad)
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau1)) deallocate(elsttau1)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)
         if(allocated(elsttau2)) deallocate(elsttau2)

         if(allocated(dEd1)) deallocate(dEd1)
         if(allocated(dEd2)) deallocate(dEd2)
         if(allocated(dEp1)) deallocate(dEp1)
         if(allocated(dEp2)) deallocate(dEp2)

         if(allocated(taud1)) deallocate(taud1)
         if(allocated(taud2)) deallocate(taud2)
         if(allocated(taup1)) deallocate(taup1)
         if(allocated(taup2)) deallocate(taup2)

         if(allocated(frcztau1dtot)) deallocate(frcztau1dtot)
         if(allocated(frcytau1dtot)) deallocate(frcytau1dtot)
         if(allocated(frcxtau1dtot)) deallocate(frcxtau1dtot)
         if(allocated(frcztau1ptot)) deallocate(frcztau1ptot)
         if(allocated(frcytau1ptot)) deallocate(frcytau1ptot)
         if(allocated(frcxtau1ptot)) deallocate(frcxtau1ptot)
         if(allocated(frcztau2dtot)) deallocate(frcztau2dtot)
         if(allocated(frcytau2dtot)) deallocate(frcytau2dtot)
         if(allocated(frcxtau2dtot)) deallocate(frcxtau2dtot)
         if(allocated(frcztau2ptot)) deallocate(frcztau2ptot)
         if(allocated(frcytau2ptot)) deallocate(frcytau2ptot)
         if(allocated(frcxtau2ptot)) deallocate(frcxtau2ptot)

         allocate(elst_grad(largest_nelst2,npole))
         allocate(elsttau1_index(4,largest_nelst2,npole))
         allocate(elsttau1(4,largest_nelst2,npole))
         allocate(elsttau2_index(4,largest_nelst2,npole))
         allocate(elsttau2(4,largest_nelst2,npole))

         allocate(dEd1(9,ntpair_dE))
         allocate(dEd2(9,ntpair_dE))
         allocate(dEp1(9,ntpair_dE))
         allocate(dEp2(9,ntpair_dE))

         allocate(taud1(9,ntpair_tau1))
         allocate(taup1(9,ntpair_tau1))
         allocate(taud2(9,ntpair_tau2))
         allocate(taup2(9,ntpair_tau2))

         allocate(frcztau1dtot(9,ntpair_dE))
         allocate(frcytau1dtot(9,ntpair_dE))
         allocate(frcxtau1dtot(9,ntpair_dE))
         allocate(frcztau1ptot(9,ntpair_dE))
         allocate(frcytau1ptot(9,ntpair_dE))
         allocate(frcxtau1ptot(9,ntpair_dE))
         allocate(frcztau2dtot(9,ntpair_dE))
         allocate(frcytau2dtot(9,ntpair_dE))
         allocate(frcxtau2dtot(9,ntpair_dE))
         allocate(frcztau2ptot(9,ntpair_dE))
         allocate(frcytau2ptot(9,ntpair_dE))
         allocate(frcxtau2ptot(9,ntpair_dE))

         ntpair_dEoffset=0
         ntpair_tau1offset=0
         ntpair_tau2offset=0

         do k=2,numtasks_emreal2+1
           m=ntpair_dE_array(k-2)
           allocate(dEd1tmp_rcv(9,m))
           allocate(dEd2tmp_rcv(9,m))
           allocate(dEp1tmp_rcv(9,m))
           allocate(dEp2tmp_rcv(9,m))

           allocate(frcztau1d_rcv(9,m))
           allocate(frcytau1d_rcv(9,m))
           allocate(frcxtau1d_rcv(9,m))
           allocate(frcztau1p_rcv(9,m))
           allocate(frcytau1p_rcv(9,m))
           allocate(frcxtau1p_rcv(9,m))

           allocate(frcztau2d_rcv(9,m))
           allocate(frcytau2d_rcv(9,m))
           allocate(frcxtau2d_rcv(9,m))
           allocate(frcztau2p_rcv(9,m))
           allocate(frcytau2p_rcv(9,m))
           allocate(frcxtau2p_rcv(9,m))
           
NEED TO CORRECT THE TAGS
           call mpi_recv(dEd1tmp_rcv,9*m,mpi_real8,k,2,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEd2tmp_rcv,9*m,mpi_real8,k,3,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEp1tmp_rcv,9*m,mpi_real8,k,4,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(dEp2tmp_rcv,9*m,mpi_real8,k,5,
     &     mpi_comm_world,stat,ierr)

           call mpi_recv(frcztau1d_rcv,9*m,mpi_real8,k,12,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau1d_rcv,9*m,mpi_real8,k,13,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau1d_rcv,9*m,mpi_real8,k,14,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcztau1p_rcv,9*m,mpi_real8,k,15,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau1p_rcv,9*m,mpi_real8,k,16,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau1p_rcv,9*m,mpi_real8,k,17,
     &     mpi_comm_world,stat,ierr)

           call mpi_recv(frcztau2d_rcv,9*m,mpi_real8,k,18,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau2d_rcv,9*m,mpi_real8,k,19,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau2d_rcv,9*m,mpi_real8,k,20,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcztau2p_rcv,9*m,mpi_real8,k,21,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcytau2p_rcv,9*m,mpi_real8,k,22,
     &     mpi_comm_world,stat,ierr)
           call mpi_recv(frcxtau2p_rcv,9*m,mpi_real8,k,23,
     &     mpi_comm_world,stat,ierr)

           do i=1,m
              k2=ntpair_dEoffset+i
              do j=1,9
                dEd1(j,k2)=dEd1tmp_rcv(j,i)
                dEd2(j,k2)=dEd2tmp_rcv(j,i)
                dEp1(j,k2)=dEp1tmp_rcv(j,i)
                dEp2(j,k2)=dEp2tmp_rcv(j,i)

                frcztau1dtot(j,k2)=frcztau1d_rcv(j,i)
                frcytau1dtot(j,k2)=frcytau1d_rcv(j,i)
                frcxtau1dtot(j,k2)=frcxtau1d_rcv(j,i)
                frcztau1ptot(j,k2)=frcztau1p_rcv(j,i)
                frcytau1ptot(j,k2)=frcytau1p_rcv(j,i)
                frcxtau1ptot(j,k2)=frcxtau1p_rcv(j,i)

                frcztau2dtot(j,k2)=frcztau2d_rcv(j,i)
                frcytau2dtot(j,k2)=frcytau2d_rcv(j,i)
                frcxtau2dtot(j,k2)=frcxtau2d_rcv(j,i)
                frcztau2ptot(j,k2)=frcztau2p_rcv(j,i)
                frcytau2ptot(j,k2)=frcytau2p_rcv(j,i)
                frcxtau2ptot(j,k2)=frcxtau2p_rcv(j,i)
              end do
           end do

           deallocate(dEd1tmp_rcv)
           deallocate(dEd2tmp_rcv)
           deallocate(dEp1tmp_rcv)
           deallocate(dEp2tmp_rcv)

           deallocate(frcztau1d_rcv)
           deallocate(frcytau1d_rcv)
           deallocate(frcxtau1d_rcv)
           deallocate(frcztau1p_rcv)
           deallocate(frcytau1p_rcv)
           deallocate(frcxtau1p_rcv)
           deallocate(frcztau2d_rcv)
           deallocate(frcytau2d_rcv)
           deallocate(frcxtau2d_rcv)
           deallocate(frcztau2p_rcv)
           deallocate(frcytau2p_rcv)
           deallocate(frcxtau2p_rcv)


           m1=ntpair_tau1_array(k-2)
           m2=ntpair_tau2_array(k-2)
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

           call mpi_recv(largest_nelst,1,mpi_integer,k,24,
     &     mpi_comm_world,stat,ierr)

           sz=last_emreal2(k-2)-start_emreal2(k-2)+1

           allocate(elstgrad_rcv(largest_nelst,sz))
           allocate(elsttau1_index_rcv(4,largest_nelst,sz))
           allocate(elsttau1_rcv(4,largest_nelst,sz))
           allocate(elsttau2_index_rcv(4,largest_nelst,sz))
           allocate(elsttau2_rcv(4,largest_nelst,sz))

           call mpi_recv(elstgrad_rcv,sz*largest_nelst,
     &     mpi_integer,k,25,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau1_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,26,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau1_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,27,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau2_index_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,28,mpi_comm_world,stat,ierr)
           call mpi_recv(elsttau2_rcv,4*sz*largest_nelst,
     &     mpi_integer,k,29,mpi_comm_world,stat,ierr)
 
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
         end do

      end if

      if(allocated(dEd1tmp)) deallocate(dEd1tmp)
      if(allocated(dEd2tmp)) deallocate(dEd2tmp)
      if(allocated(dEp1tmp)) deallocate(dEp1tmp)
      if(allocated(dEp2tmp)) deallocate(dEp2tmp)

      if(allocated(taud1tmp)) deallocate(taud1tmp)
      if(allocated(taud2tmp)) deallocate(taud2tmp)
      if(allocated(taup1tmp)) deallocate(taup1tmp)
      if(allocated(taup2tmp)) deallocate(taup2tmp)

      !if(allocated(tau1indextmp)) deallocate(tau1indextmp)
      !if(allocated(tau2indextmp)) deallocate(tau2indextmp)

      if(allocated(frcztau1d)) deallocate (frcztau1d)
      if(allocated(frcytau1d)) deallocate (frcytau1d)
      if(allocated(frcxtau1d)) deallocate (frcxtau1d)
      if(allocated(frcztau1p)) deallocate (frcztau1p)
      if(allocated(frcytau1p)) deallocate (frcytau1p)
      if(allocated(frcxtau1p)) deallocate (frcxtau1p)
      if(allocated(frcztau2d)) deallocate (frcztau2d)
      if(allocated(frcytau2d)) deallocate (frcytau2d)
      if(allocated(frcxtau2d)) deallocate (frcxtau2d)
      if(allocated(frcztau2p)) deallocate (frcztau2p)
      if(allocated(frcytau2p)) deallocate (frcytau2p)
      if(allocated(frcxtau2p)) deallocate (frcxtau2p)

      if(allocated(elstgrad_tmp)) deallocate(elstgrad_tmp)
      if(allocated(elsttau1_index_tmp)) deallocate(elsttau1_index_tmp)
      if(allocated(elsttau1_tmp)) deallocate(elsttau1_tmp)
      if(allocated(elsttau2_index_tmp)) deallocate(elsttau2_index_tmp)
      if(allocated(elsttau2_tmp)) deallocate(elsttau2_tmp)

         largest_nelst2=0

         do i=1,npole
           if(nelst2(i).gt.largest_nelst2) then
             largest_nelst2=nelst2(i)
           end if 
         end do

       call mpi_bcast(ntpair_dE,1,mpi_integer,numtasks-1,
     &   mpi_comm_world,ierr)
       call mpi_bcast(ntpair_tau1,1,mpi_integer,numtasks-1,
     &   mpi_comm_world,ierr)
       call mpi_bcast(ntpair_tau2,1,mpi_integer,numtasks-1,
     &   mpi_comm_world,ierr)

      if(taskid.ne.numtasks-1) then   
         if(allocated(elst_grad)) deallocate(elst_grad)
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau1)) deallocate(elsttau1)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)
         if(allocated(elsttau2)) deallocate(elsttau2)

         if(allocated(dEd1)) deallocate(dEd1)
         if(allocated(dEd2)) deallocate(dEd2)
         if(allocated(dEp1)) deallocate(dEp1)
         if(allocated(dEp2)) deallocate(dEp2)

         if(allocated(taud1)) deallocate(taud1)
         if(allocated(taud2)) deallocate(taud2)
         if(allocated(taup1)) deallocate(taup1)
         if(allocated(taup2)) deallocate(taup2)

         if(allocated(frcztau1dtot)) deallocate(frcztau1dtot)
         if(allocated(frcytau1dtot)) deallocate(frcytau1dtot)
         if(allocated(frcxtau1dtot)) deallocate(frcxtau1dtot)
         if(allocated(frcztau1ptot)) deallocate(frcztau1ptot)
         if(allocated(frcytau1ptot)) deallocate(frcytau1ptot)
         if(allocated(frcxtau1ptot)) deallocate(frcxtau1ptot)
         if(allocated(frcztau2dtot)) deallocate(frcztau2dtot)
         if(allocated(frcytau2dtot)) deallocate(frcytau2dtot)
         if(allocated(frcxtau2dtot)) deallocate(frcxtau2dtot)
         if(allocated(frcztau2ptot)) deallocate(frcztau2ptot)
         if(allocated(frcytau2ptot)) deallocate(frcytau2ptot)
         if(allocated(frcxtau2ptot)) deallocate(frcxtau2ptot)

         allocate(elst_grad(largest_nelst2,npole))
         allocate(elsttau1_index(4,largest_nelst2,npole))
         allocate(elsttau1(4,largest_nelst2,npole))
         allocate(elsttau2_index(4,largest_nelst2,npole))
         allocate(elsttau2(4,largest_nelst2,npole))

         allocate(dEd1(9,ntpair_dE))
         allocate(dEd2(9,ntpair_dE))
         allocate(dEp1(9,ntpair_dE))
         allocate(dEp2(9,ntpair_dE))

         allocate(taud1(9,ntpair_tau1))
         allocate(taup1(9,ntpair_tau1))
         allocate(taud2(9,ntpair_tau2))
         allocate(taup2(9,ntpair_tau2))

         allocate(frcztau1dtot(9,ntpair_dE))
         allocate(frcytau1dtot(9,ntpair_dE))
         allocate(frcxtau1dtot(9,ntpair_dE)) 
         allocate(frcztau1ptot(9,ntpair_dE)) 
         allocate(frcytau1ptot(9,ntpair_dE)) 
         allocate(frcxtau1ptot(9,ntpair_dE)) 
         allocate(frcztau2dtot(9,ntpair_dE)) 
         allocate(frcytau2dtot(9,ntpair_dE)) 
         allocate(frcxtau2dtot(9,ntpair_dE)) 
         allocate(frcztau2ptot(9,ntpair_dE)) 
         allocate(frcytau2ptot(9,ntpair_dE)) 
         allocate(frcxtau2ptot(9,ntpair_dE)) 
      end if


      call mpi_bcast(dEd1,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dEd2,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)      
      call mpi_bcast(dEp1,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(dEp2,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(taud1,9*ntpair_tau1,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taud2,9*ntpair_tau2,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taup1,9*ntpair_tau1,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(taup2,9*ntpair_tau2,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(frcztau1dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau1dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau1dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcztau1ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau1ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau1ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcztau2dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau2dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau2dtot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcztau2ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcytau2ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)
      call mpi_bcast(frcxtau2ptot,9*ntpair_dE,mpi_real8,numtasks-1,
     &   mpi_comm_world,ierr)

      call mpi_bcast(elst_grad,largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau1_index,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau1,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau2_index,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
      call mpi_bcast(elsttau2,4*largest_nelst2*npole,mpi_integer,
     &   numtasks-1,mpi_comm_world,ierr)
