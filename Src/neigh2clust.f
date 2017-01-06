c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  neigh.i  --  pairwise neighbor list indices and storage  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lbuffer     width of the neighbor list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     cbuf2       square of charge cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     vbufx       square of vdw cutoff plus twice the list buffer
c     cbufx       square of charge cutoff plus twice the list buffer
c     mbufx       square of multipole cutoff plus twice the list buffer
c     xvold       x-coordinate at last vdw neighbor list update
c     yvold       y-coordinate at last vdw neighbor list update
c     zvold       z-coordinate at last vdw neighbor list update
c     xcold       x-coordinate at last charge neighbor list update
c     ycold       y-coordinate at last charge neighbor list update
c     zcold       z-coordinate at last charge neighbor list update
c     xmold       x-coordinate at last multipole neighbor list update
c     ymold       y-coordinate at last multipole neighbor list update
c     zmold       z-coordinate at last multipole neighbor list update
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     dovlst      logical flag to rebuild vdw neighbor list
c     doclst      logical flag to rebuild charge neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c     domollst2bod
c     domollst3bod
      module neigh2clust
      implicit none
      integer, allocatable :: clust(:,:)
      integer, allocatable :: sizeclust(:)
      integer clustcount
      real*8, allocatable :: clust_cm(:,:)
      real*8, allocatable :: distmax(:)
      real*8, allocatable :: xmolold2(:)
      real*8, allocatable :: ymolold2(:)
      real*8, allocatable :: zmolold2(:)
      logical doclust,do2waterclustlist,doclust3,doclust4
      real*8 clustcut
      logical inputkmeansclust
      integer maxsizeclust
      integer, allocatable :: start_mollst_polar(:,:,:)
      integer, allocatable :: last_mollst_polar(:,:,:)
      integer, allocatable :: num_mollst_chunk(:,:)
      integer max_nmollst
      logical ompouterloop3b
      integer, allocatable :: elst1b(:,:,:)
      integer, allocatable :: nelst1b(:,:)
      logical, allocatable :: domlstclust1b(:)
      real*8, allocatable :: xmold1b(:,:)
      real*8, allocatable :: ymold1b(:,:)
      real*8, allocatable :: zmold1b(:,:)
      integer, allocatable :: ulst1b(:,:,:)
      integer, allocatable :: nulst1b(:,:)
      logical, allocatable :: doulstclust1b(:)
      real*8, allocatable :: xuold1b(:,:)
      real*8, allocatable :: yuold1b(:,:)
      real*8, allocatable :: zuold1b(:,:)
      integer, allocatable :: elst2b(:,:,:)
      integer, allocatable :: nelst2b(:,:)
      logical, allocatable :: domlstclust2b(:)
      real*8, allocatable :: xmold2b(:,:)
      real*8, allocatable :: ymold2b(:,:)
      real*8, allocatable :: zmold2b(:,:)
      integer, allocatable :: ulst2b(:,:,:)
      integer, allocatable :: nulst2b(:,:)
      logical, allocatable :: doulstclust2b(:)
      real*8, allocatable :: xuold2b(:,:)
      real*8, allocatable :: yuold2b(:,:)
      real*8, allocatable :: zuold2b(:,:)
      real*8 mbuf2clust,mbufxclust
      logical mlistclust,use_ewaldclust
      integer maxulst2
      logical all2bclust,useboxclust
      integer cnt3bwork
      logical, allocatable :: save2bfor3b(:,:)
      integer, allocatable :: counter2b_ind(:,:)
      integer, allocatable :: ind2b_counter(:,:)
      integer num2bsave,numtasks_polar3b,max_nmollst3b
      integer, allocatable :: start_mollst_polar3b(:,:,:)
      integer, allocatable :: last_mollst_polar3b(:,:,:)
      integer, allocatable :: num_mollst_chunk3b(:,:)
      integer, allocatable :: start_polar3b(:)
      integer, allocatable :: last_polar3b(:)
      integer, allocatable :: molnew3b(:)
      integer, allocatable :: elst3b(:,:,:)
      integer, allocatable :: nelst3b(:,:)
      logical, allocatable :: domlstclust3b(:)
      real*8, allocatable :: xmold3b(:,:)
      real*8, allocatable :: ymold3b(:,:)
      real*8, allocatable :: zmold3b(:,:)
      integer, allocatable :: ulst3b(:,:,:)
      integer, allocatable :: nulst3b(:,:)
      logical, allocatable :: doulstclust3b(:)
      real*8, allocatable :: xuold3b(:,:)
      real*8, allocatable :: yuold3b(:,:)
      real*8, allocatable :: zuold3b(:,:)
      save
      end
