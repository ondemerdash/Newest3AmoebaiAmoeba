
c     ##############################################################
c     ##                                                          ##
c     ##  ewreg.i  --  exponential factors for regular Ewald sum  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxvec    maximum number of k-vectors per reciprocal axis
c
c     ejc       exponental factors for cosine along the j-axis
c     ejs       exponental factors for sine along the j-axis
c     ekc       exponental factors for cosine along the k-axis
c     eks       exponental factors for sine along the k-axis
c     elc       exponental factors for cosine along the l-axis
c     els       exponental factors for sine along the l-axis
c
c
      module ewreg3bpolz
      use sizes
      implicit none
      integer maxvec
      parameter (maxvec=50)
c      real*8 ejc(maxatm,0:maxvec),ejs(maxatm,0:maxvec)
c      real*8 ekc(maxatm,-maxvec:maxvec)
c      real*8 eks(maxatm,-maxvec:maxvec)
c      real*8 elc(maxatm,-maxvec:maxvec)
c      real*8 els(maxatm,-maxvec:maxvec)
      integer jmax3b,kmax3b,lmax3b
      save
      end
c
