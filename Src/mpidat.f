      module mpidat
      implicit none 
      integer master,taskid,numtasks,master1,master2
      integer emreal_comm,numtasks_emreal,mpi_comm_emreal
      integer remainder_emreal,offset_emreal
      integer emrecip_comm,numtasks_emrecip
      integer numtasks_polar,numtasks_emreal2
      integer, allocatable :: start_emreal2(:)
      integer, allocatable :: last_emreal2(:)
      integer, allocatable :: maxsize_elst(:)
c      integer, allocatable :: offset_emreal2(:)
      integer, allocatable :: start_polar(:)
      integer, allocatable :: last_polar(:)
      integer, allocatable :: molnew(:)
      integer offset_emreal3
      integer, allocatable :: start_polar3(:)
      integer, allocatable :: last_polar3(:)
      integer, allocatable :: maxsize_mollst3(:)
      integer numtasks_polar3
      integer maxsize_elst2
      integer, allocatable :: start_vdw2(:)
      integer, allocatable :: last_vdw2(:)
      integer, allocatable :: maxsize_vlst(:)
      integer numtasks_vdw2,offset_vdw3
      integer numtasks_vdw
      integer reqs(6)
      integer numtasks_polz
      save 
      end
