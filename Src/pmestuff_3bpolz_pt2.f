c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2010 by T. Darden, D. Gohara & Jay W. Ponder  ##
c     ##                      All Rights Reserved                     ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below implement various B-spline and coordinate  ##
c     ##  manipulations for particle mesh Ewald summation; spatial  ##
c     ##  grid assignment by David Gohara; modified from original   ##
c     ##  PME code by Thomas Darden, NIEHS, Research Triangle, NC   ##
c     ##                                                            ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine bspline_fill  --  get PME B-spline coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "bspline_fill" finds B-spline coefficients and derivatives
c     for PME atomic sites along the fractional coordinate axes
c
c
      subroutine bspline_fill3b(igrid3b,thetai1_3b,thetai2_3b,
     &           thetai3_3b)
      use sizes
      use atoms
      use boxes
      use pme
      implicit none
      integer i,ifr
      real*8 xi,yi,zi
      real*8 w,fr,eps
      logical first
      save first
      data first  / .true. /
      real*8 thetai1_3b(4,bsorder,n)
      real*8 thetai2_3b(4,bsorder,n)
      real*8 thetai3_3b(4,bsorder,n)
      integer igrid3b(3,n)
      
c
c
c     perform dynamic allocation of some global arrays
c

c      if (first) then
c         first = .false.
c         if (.not. allocated(igrid))  allocate (igrid(3,n))
c      end if

c
c     offset used to shift sites off exact lattice bounds
c
      eps = 1.0d-8
c
c     get the B-spline coefficients for each atomic site
c
      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr = dble(nfft1) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
c         igrid(1,i) = ifr - bsorder
         igrid3b(1,i) = ifr - bsorder
c         call bsplgen (w,thetai1(1,1,i))
         call bsplgen (w,thetai1_3b(1,1,i))
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
c         igrid(2,i) = ifr - bsorder
c         call bsplgen (w,thetai2(1,1,i))
         igrid3b(2,i) = ifr - bsorder
         call bsplgen (w,thetai2_3b(1,1,i))
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr = dble(nfft3) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
c         igrid(3,i) = ifr - bsorder
c         call bsplgen (w,thetai3(1,1,i))
         igrid3b(3,i) = ifr - bsorder
         call bsplgen (w,thetai3_3b(1,1,i))         
      end do
c      print*,"Successful end of bspline_fill"
      return
      end
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine table_fill  --  spatial chunks for each site  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "table_fill" constructs an array which stores the spatial
c     regions of the particle mesh Ewald grid with contributions
c     from each electrostatic site
c
c
      subroutine table_fill3b(igrid3b,pmetable3b)
      use sizes
      use atoms
      use chunks
      use pme
      implicit none
      integer i,k
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      logical negx,negy,negz
      logical posx,posy,posz
      logical midx,midy,midz
      integer pmetable3b(n,nchunk)
      integer igrid3b(3,n)

c
c
c     zero out the PME table marking chunks per site
c
c      print*,"nchunk in table_fill",nchunk
c      print*,"grdoff in table_fill",grdoff
c      print*,"nlpts nrpts in table_fill",nlpts,nrpts
c      print*,"nfft123 in table_fill",nfft1,nfft2,nfft3
 
      do k = 1, nchunk
         do i = 1, n
c            pmetable(i,k) = 0
            pmetable3b(i,k) =0
         end do
      end do
c
c     loop over sites to find the spatial chunks for each
c
      do i = 1, n
         nearpt(1) = igrid3b(1,i) + grdoff
         nearpt(2) = igrid3b(2,i) + grdoff
         nearpt(3) = igrid3b(3,i) + grdoff
         if (nearpt(1) .lt. 1) then
            nearpt(1) = mod(nearpt(1),nfft1) + nfft1
         else if (nearpt(1) .gt. nfft1) then
            nearpt(1) = mod(nearpt(1),nfft1)
         end if
         if (nearpt(2) .lt. 1) then
            nearpt(2) = mod(nearpt(2),nfft2) + nfft2
         else if (nearpt(2) .gt. nfft2) then
            nearpt(2) = mod(nearpt(2),nfft2)
         end if
         if (nearpt(3) .lt. 1) then
            nearpt(3) = mod(nearpt(3),nfft3) + nfft3
         else if (nearpt(3) .gt. nfft3) then
            nearpt(3) = mod(nearpt(3),nfft3)
         end if
         abound(1) = nearpt(1) - nlpts
         abound(2) = nearpt(1) + nrpts
         abound(3) = nearpt(2) - nlpts
         abound(4) = nearpt(2) + nrpts
         abound(5) = nearpt(3) - nlpts
         abound(6) = nearpt(3) + nrpts
         cid(1) = (nearpt(1)-1)/ngrd1 + 1
         cid(2) = (nearpt(2)-1)/ngrd2 + 1
         cid(3) = (nearpt(3)-1)/ngrd3 + 1
         cbound(1) = (cid(1)-1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = (cid(2)-1)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = (cid(3)-1)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
c
c     set and store central chunk where the site is located
c
         k = (cid(3)-1)*nchk1*nchk2 + (cid(2)-1)*nchk1 + cid(1)
c         pmetable(i,k) = 1
         pmetable3b(i,k) =1
c
c     flags for atom bounds to left or right of central chunk
c
         negx = (abound(1) .lt. cbound(1))
         negy = (abound(3) .lt. cbound(3))
         negz = (abound(5) .lt. cbound(5))
         posx = (abound(2) .gt. cbound(2))
         posy = (abound(4) .gt. cbound(4))
         posz = (abound(6) .gt. cbound(6))
c
c     flags for atom bounds fully inside the central chunk
c
         midx = (.not.negx .and. .not.posx)
         midy = (.not.negy .and. .not.posy)
         midz = (.not.negz .and. .not.posz)
         if (midx .and. midy .and. midz)  goto 10
c
c     flags for atom bounds that overlap the central chunk
c
         midx = (.not.negx .or. .not.posx)
         midy = (.not.negy .or. .not.posy)
         midz = (.not.negz .or. .not.posz)
c
c     check for overlap with any of the neighboring chunks
c
         if (midx .and. midy .and. negz) 
     &      call setchunk3b (i,cid,0,0,-1,pmetable3b)
         if (midx .and. midy .and. posz)
     &      call setchunk3b (i,cid,0,0,1,pmetable3b)
         if (midx .and. negy .and. midz)
     &      call setchunk3b (i,cid,0,-1,0,pmetable3b)
         if (midx .and. posy .and. midz)
     &      call setchunk3b (i,cid,0,1,0,pmetable3b)
         if (negx .and. midy .and. midz)
     &       call setchunk3b (i,cid,-1,0,0,pmetable3b)
         if (posx .and. midy .and. midz)
     &       call setchunk3b (i,cid,1,0,0,pmetable3b)
         if (midx .and. negy .and. negz)
     &       call setchunk3b (i,cid,0,-1,-1,pmetable3b)
         if (midx .and. negy .and. posz)
     &       call setchunk3b (i,cid,0,-1,1,pmetable3b)
         if (midx .and. posy .and. negz)
     &       call setchunk3b (i,cid,0,1,-1,pmetable3b)
         if (midx .and. posy .and. posz)
     &       call setchunk3b (i,cid,0,1,1,pmetable3b)
         if (negx .and. midy .and. negz)
     &       call setchunk3b (i,cid,-1,0,-1,pmetable3b)
         if (negx .and. midy .and. posz)
     &       call setchunk3b (i,cid,-1,0,1,pmetable3b)
         if (posx .and. midy .and. negz)
     &       call setchunk3b (i,cid,1,0,-1,pmetable3b)
         if (posx .and. midy .and. posz)
     &        call setchunk3b (i,cid,1,0,1,pmetable3b)
         if (negx .and. negy .and. midz) 
     &        call setchunk3b (i,cid,-1,-1,0,pmetable3b)
         if (negx .and. posy .and. midz)
     &        call setchunk3b (i,cid,-1,1,0,pmetable3b)
         if (posx .and. negy .and. midz)
     &        call setchunk3b (i,cid,1,-1,0,pmetable3b)
         if (posx .and. posy .and. midz)
     &       call setchunk3b (i,cid,1,1,0,pmetable3b)
         if (negx .and. negy .and. negz)
     &       call setchunk3b (i,cid,-1,-1,-1,pmetable3b)
         if (negx .and. negy .and. posz)
     &       call setchunk3b (i,cid,-1,-1,1,pmetable3b)
         if (negx .and. posy .and. negz)
     &       call setchunk3b (i,cid,-1,1,-1,pmetable3b)
         if (posx .and. negy .and. negz)
     &       call setchunk3b (i,cid,1,-1,-1,pmetable3b)
         if (negx .and. posy .and. posz)
     &       call setchunk3b (i,cid,-1,1,1,pmetable3b)
         if (posx .and. negy .and. posz)
     &       call setchunk3b (i,cid,1,-1,1,pmetable3b)
         if (posx .and. posy .and. negz)
     &       call setchunk3b (i,cid,1,1,-1,pmetable3b)
         if (posx .and. posy .and. posz)
     &       call setchunk3b (i,cid,1,1,1,pmetable3b)
   10    continue
      end do
c      print*,"Successful end of table_fill"
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setchunk  --  site overlaps neighboring chunk  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setchunk" marks a chunk in the PME spatial table which is
c     overlapped by the B-splines for an electrostatic site
c
c
c      subroutine setchunk (i,cid,off1,off2,off3)
      subroutine setchunk3b (i,cid,off1,off2,off3,pmetable3b)
      use sizes
      use chunks
      use pme
      use atoms
      implicit none
      integer i,k
      integer off1,off2,off3
      integer cid(3),temp(3)
      integer pmetable3b(n,nchunk)

c
c
c     mark neighboring chunk overlapped by an electrostatic site
c
      temp(1) = cid(1) + off1
      if (temp(1) .lt. 1)  temp(1) = nchk1
      if (temp(1) .gt. nchk1)  temp(1) = 1
      temp(2) = cid(2) + off2
      if (temp(2) .lt. 1)  temp(2) = nchk2
      if (temp(2) .gt. nchk2)  temp(2) = 1
      temp(3) = cid(3) + off3
      if (temp(3) .lt. 1)  temp(3) = nchk3
      if (temp(3) .gt. nchk3)  temp(3) = 1
      k = (temp(3)-1)*nchk1*nchk2 + (temp(2)-1)*nchk1 + temp(1)
c      pmetable(i,k) = 1
      pmetable3b(i,k) =1
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c
c
