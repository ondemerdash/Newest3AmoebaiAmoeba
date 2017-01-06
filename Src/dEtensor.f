c
c
      module dEtensor
      implicit none
      real*8, allocatable :: dEd1(:,:)
      real*8, allocatable :: dEd2(:,:)
      real*8, allocatable :: dEp1(:,:)
      real*8, allocatable :: dEp2(:,:)
      integer, allocatable :: dEindex(:,:) 
      integer ntpair_dE

      real*8, allocatable :: dEd1tmp(:,:)
      real*8, allocatable :: dEd2tmp(:,:)
      real*8, allocatable :: dEp1tmp(:,:)
      real*8, allocatable :: dEp2tmp(:,:)
      integer, allocatable :: dEindextmp(:,:)
      integer ntpair_dEtmp

      real*8, allocatable :: dEd1tmp_rcv(:,:)
      real*8, allocatable :: dEd2tmp_rcv(:,:)
      real*8, allocatable :: dEp1tmp_rcv(:,:)
      real*8, allocatable :: dEp2tmp_rcv(:,:)
      integer, allocatable :: dEindextmp_rcv(:,:)
      integer ntpair_dEtmp_rcv

      integer, allocatable :: ntpair_start(:)
      integer, allocatable :: ntpair_last(:)
      integer, allocatable :: ntpair_start_rmndr(:)
      integer, allocatable :: ntpair_last_rmndr(:)

      real*8, allocatable :: taud1tmp(:,:)
      real*8, allocatable :: taud2tmp(:,:)
      real*8, allocatable :: taup1tmp(:,:)
      real*8, allocatable :: taup2tmp(:,:)
      integer, allocatable :: tau1indextmp(:,:)
      integer, allocatable :: tau2indextmp(:,:)
      integer ntpair_tau1tmp,ntpair_tau2tmp

      real*8, allocatable :: frcztau1d(:,:)
      real*8, allocatable :: frcytau1d(:,:)
      real*8, allocatable :: frcxtau1d(:,:)
      real*8, allocatable :: frcztau1p(:,:)
      real*8, allocatable :: frcytau1p(:,:)
      real*8, allocatable :: frcxtau1p(:,:)
      real*8, allocatable :: frcztau2d(:,:)
      real*8, allocatable :: frcytau2d(:,:)
      real*8, allocatable :: frcxtau2d(:,:)
      real*8, allocatable :: frcztau2p(:,:)
      real*8, allocatable :: frcytau2p(:,:)
      real*8, allocatable :: frcxtau2p(:,:)

      save 
      end
