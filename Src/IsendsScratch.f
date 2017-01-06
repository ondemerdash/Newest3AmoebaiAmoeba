        call mpi_isend(dEd1_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,2,mpi_comm_world,reqs22,ierr)
        call mpi_isend(dEd2_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,3,mpi_comm_world,reqs23,ierr)
        call mpi_isend(dEp1_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,4,mpi_comm_world,reqs24,ierr)
        call mpi_isend(dEp2_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,5,mpi_comm_world,reqs25,ierr)

        call mpi_isend(taud1_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,8,mpi_comm_world,reqs28,ierr)
        call mpi_isend(taud2_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,9,mpi_comm_world,reqs29,ierr)
       call mpi_isend(taup1_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,10,mpi_comm_world,reqs30,ierr)
       call mpi_isend(taup2_3tmp,9*4*largest_nelst*sz,mpi_real8,
     &       numtasks-1,11,mpi_comm_world,reqs31,ierr)


        call mpi_isend(frcztau1dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,12,mpi_comm_world,reqs32,ierr)
        call mpi_isend(frcytau1dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,13,mpi_comm_world,reqs33,ierr)
        call mpi_isend(frcxtau1dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,14,mpi_comm_world,reqs34,ierr)
        call mpi_isend(frcztau1ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,15,mpi_comm_world,reqs35,ierr)
        call mpi_isend(frcytau1ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,16,mpi_comm_world,reqs36,ierr)
        call mpi_isend(frcxtau1ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,17,mpi_comm_world,reqs37,ierr)

        call mpi_isend(frcztau2dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,18,mpi_comm_world,reqs38,ierr)
        call mpi_isend(frcytau2dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,19,mpi_comm_world,reqs39,ierr)
        call mpi_isend(frcxtau2dtot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,20,mpi_comm_world,reqs40,ierr)
        call mpi_isend(frcztau2ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,21,mpi_comm_world,reqs41,ierr)
        call mpi_isend(frcytau2ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,22,mpi_comm_world,reqs42,ierr)
        call mpi_isend(frcxtau2ptot_3tmp,9*largest_nelst*sz,mpi_real8,
     &       numtasks-1,23,mpi_comm_world,reqs43,ierr)

        call mpi_isend(largest_nelst,1,mpi_integer,
     &       numtasks-1,24,mpi_comm_world,reqs44,ierr)
        call mpi_isend(elsttau1_index_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,26,mpi_comm_world,reqs46,ierr)
        call mpi_isend(elsttau2_index_tmp,4*sz*largest_nelst,
     &   mpi_integer,numtasks-1,28,mpi_comm_world,reqs48,ierr)

