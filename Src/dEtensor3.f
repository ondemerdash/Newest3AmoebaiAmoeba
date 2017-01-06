c
c
      module dEtensor3
      implicit none
      real*8, allocatable :: dEd1_3(:,:,:)
      real*8, allocatable :: dEd2_3(:,:,:)
      real*8, allocatable :: dEp1_3(:,:,:)
      real*8, allocatable :: dEp2_3(:,:,:)

      real*8, allocatable :: taud1_3(:,:,:,:)
      real*8, allocatable :: taud2_3(:,:,:,:)
      real*8, allocatable :: taup1_3(:,:,:,:)
      real*8, allocatable :: taup2_3(:,:,:,:)

      real*8, allocatable :: frcztau1dtot_3(:,:,:)
      real*8, allocatable :: frcytau1dtot_3(:,:,:)
      real*8, allocatable :: frcxtau1dtot_3(:,:,:)
      real*8, allocatable :: frcztau1ptot_3(:,:,:)
      real*8, allocatable :: frcytau1ptot_3(:,:,:)
      real*8, allocatable :: frcxtau1ptot_3(:,:,:)
      real*8, allocatable :: frcztau2dtot_3(:,:,:)
      real*8, allocatable :: frcytau2dtot_3(:,:,:)
      real*8, allocatable :: frcxtau2dtot_3(:,:,:)
      real*8, allocatable :: frcztau2ptot_3(:,:,:)
      real*8, allocatable :: frcytau2ptot_3(:,:,:)
      real*8, allocatable :: frcxtau2ptot_3(:,:,:)

      real*8, allocatable :: dEd1_3tmp(:,:,:)
      real*8, allocatable :: dEd2_3tmp(:,:,:)
      real*8, allocatable :: dEp1_3tmp(:,:,:)
      real*8, allocatable :: dEp2_3tmp(:,:,:)

      real*8, allocatable :: taud1_3tmp(:,:,:,:)
      real*8, allocatable :: taud2_3tmp(:,:,:,:)
      real*8, allocatable :: taup1_3tmp(:,:,:,:)
      real*8, allocatable :: taup2_3tmp(:,:,:,:)

      real*8, allocatable :: frcztau1dtot_3tmp(:,:,:)
      real*8, allocatable :: frcytau1dtot_3tmp(:,:,:)
      real*8, allocatable :: frcxtau1dtot_3tmp(:,:,:)
      real*8, allocatable :: frcztau1ptot_3tmp(:,:,:)
      real*8, allocatable :: frcytau1ptot_3tmp(:,:,:)
      real*8, allocatable :: frcxtau1ptot_3tmp(:,:,:)
      real*8, allocatable :: frcztau2dtot_3tmp(:,:,:)
      real*8, allocatable :: frcytau2dtot_3tmp(:,:,:)
      real*8, allocatable :: frcxtau2dtot_3tmp(:,:,:)
      real*8, allocatable :: frcztau2ptot_3tmp(:,:,:)
      real*8, allocatable :: frcytau2ptot_3tmp(:,:,:)
      real*8, allocatable :: frcxtau2ptot_3tmp(:,:,:)


      real*8, allocatable :: dEd1_3tmp_rcv(:,:,:)
      real*8, allocatable :: dEd2_3tmp_rcv(:,:,:)
      real*8, allocatable :: dEp1_3tmp_rcv(:,:,:)
      real*8, allocatable :: dEp2_3tmp_rcv(:,:,:)

      real*8, allocatable :: taud1_3tmp_rcv(:,:,:,:)
      real*8, allocatable :: taud2_3tmp_rcv(:,:,:,:)
      real*8, allocatable :: taup1_3tmp_rcv(:,:,:,:)
      real*8, allocatable :: taup2_3tmp_rcv(:,:,:,:)


      real*8, allocatable :: frcztau1d_3rcv(:,:,:)
      real*8, allocatable :: frcytau1d_3rcv(:,:,:)
      real*8, allocatable :: frcxtau1d_3rcv(:,:,:)
      real*8, allocatable :: frcztau1p_3rcv(:,:,:)
      real*8, allocatable :: frcytau1p_3rcv(:,:,:)
      real*8, allocatable :: frcxtau1p_3rcv(:,:,:)
      real*8, allocatable :: frcztau2d_3rcv(:,:,:)
      real*8, allocatable :: frcytau2d_3rcv(:,:,:)
      real*8, allocatable :: frcxtau2d_3rcv(:,:,:)
      real*8, allocatable :: frcztau2p_3rcv(:,:,:)
      real*8, allocatable :: frcytau2p_3rcv(:,:,:)
      real*8, allocatable :: frcxtau2p_3rcv(:,:,:)
   
      save 
      end
