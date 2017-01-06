c
c
      module dEtensor2
      implicit none

      real*8, allocatable :: taud1(:,:)
      real*8, allocatable :: taud2(:,:)
      real*8, allocatable :: taup1(:,:)
      real*8, allocatable :: taup2(:,:)
      integer ntpair_tau1,ntpair_tau2

      real*8, allocatable :: frcztau1dtot(:,:)
      real*8, allocatable :: frcytau1dtot(:,:)
      real*8, allocatable :: frcxtau1dtot(:,:)
      real*8, allocatable :: frcztau1ptot(:,:)
      real*8, allocatable :: frcytau1ptot(:,:)
      real*8, allocatable :: frcxtau1ptot(:,:)
      real*8, allocatable :: frcztau2dtot(:,:)
      real*8, allocatable :: frcytau2dtot(:,:)
      real*8, allocatable :: frcxtau2dtot(:,:)
      real*8, allocatable :: frcztau2ptot(:,:)
      real*8, allocatable :: frcytau2ptot(:,:)
      real*8, allocatable :: frcxtau2ptot(:,:)


      integer, allocatable :: elstgrad(:,:)
      integer, allocatable :: elsttau1_index(:,:,:)
      integer, allocatable :: elsttau1(:,:,:)
      integer, allocatable :: elsttau2_index(:,:,:)
      integer, allocatable :: elsttau2(:,:,:)

      integer, allocatable :: elstgrad_tmp(:,:)
      integer, allocatable :: elsttau1_index_tmp(:,:,:)
      integer, allocatable :: elsttau1_tmp(:,:,:)
      integer, allocatable :: elsttau2_index_tmp(:,:,:)
      integer, allocatable :: elsttau2_tmp(:,:,:)

      real*8, allocatable :: taud1tmp_rcv(:,:)
      real*8, allocatable :: taud2tmp_rcv(:,:)
      real*8, allocatable :: taup1tmp_rcv(:,:)
      real*8, allocatable :: taup2tmp_rcv(:,:)

      real*8, allocatable :: frcztau1d_rcv(:,:)
      real*8, allocatable :: frcytau1d_rcv(:,:)
      real*8, allocatable :: frcxtau1d_rcv(:,:)
      real*8, allocatable :: frcztau1p_rcv(:,:)
      real*8, allocatable :: frcytau1p_rcv(:,:)
      real*8, allocatable :: frcxtau1p_rcv(:,:)
      real*8, allocatable :: frcztau2d_rcv(:,:)
      real*8, allocatable :: frcytau2d_rcv(:,:)
      real*8, allocatable :: frcxtau2d_rcv(:,:)
      real*8, allocatable :: frcztau2p_rcv(:,:)
      real*8, allocatable :: frcytau2p_rcv(:,:)
      real*8, allocatable :: frcxtau2p_rcv(:,:)
     
      integer, allocatable :: elstgrad_rcv(:,:)
      integer, allocatable :: elsttau1_index_rcv(:,:,:)
      integer, allocatable :: elsttau1_rcv(:,:,:)
      integer, allocatable :: elsttau2_index_rcv(:,:,:)
      integer, allocatable :: elsttau2_rcv(:,:,:)

      integer ntpair_tau1tmp_rcv,ntpair_tau2tmp_rcv
      save
      end
 
