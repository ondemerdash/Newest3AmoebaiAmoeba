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
      module neigh3b
      implicit none
      integer, allocatable :: nmollst(:)
      integer, allocatable :: mollst(:,:)
      integer, allocatable :: nmollst2(:)
      integer, allocatable :: mollst2(:,:)      
      integer, allocatable :: nmollst3(:)
      integer, allocatable :: mollst3(:,:)
      integer, allocatable :: mol3new(:)
      integer count3new
      integer, allocatable :: nmollst3_recv(:)
      integer, allocatable :: mollst3_recv(:,:)
      real*8 molbuf2,molbufx,molbufcobar
      real*8 molbufcobar_x 
      real*8 molbuf2l,molbuf2s
      real*8 molbuf2l_x,molbuf2s_x
      real*8, allocatable :: xmolold(:)
      real*8, allocatable :: ymolold(:)
      real*8, allocatable :: zmolold(:)
      real*8, allocatable :: xmolold3(:)
      real*8, allocatable :: ymolold3(:)
      real*8, allocatable :: zmolold3(:)
      real*8 rtapr2b_input,Rcut2b_input
      real*8 cut2b_input
      real*8 rtapr3b,R123cut3b
      logical domollst2bod,domollst3bod,listbcast
      logical listbcast3b
      save
      end
c      common /neigh3b/ domollst2bod,domollst3bod,nmollst,mollst,
c     &                 nmollst2,mollst2,xmolold3,ymolold3,zmolold3,
c     &                 xmolold,ymolold,zmolold,molbuf2,molbufx,
c     &                 molbuf_cobar,molbuf2l,molbuf2s,mol3new,
c     &                 count3new,nmollst3mod,mollst3mod,molbuf2l_x,
c     &                 molbuf2s_x,listbcast
