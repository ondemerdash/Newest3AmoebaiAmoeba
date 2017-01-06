c     ###############################################################
c     ##                                                           ##
c     ##  subroutine table_fill_omp  --  spatial chunks for each site  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "table_fill_omp" constructs an array which stores the spatial
c     regions of the particle mesh Ewald grid with contributions
c     from each electrostatic site
c
c
      subroutine table_fill_omp
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
c
c
c     zero out the PME table marking chunks per site
c
c      print*,"nchunk in table_fill_omp",nchunk
c      print*,"grdoff in table_fill_omp",grdoff
c      print*,"nlpts nrpts in table_fill_omp",nlpts,nrpts
c      print*,"nfft123 in table_fill_omp",nfft1,nfft2,nfft3

      do k = 1, nchunk
         do i = 1, n
            pmetable(i,k) = 0
         end do
      end do
c
c     loop over sites to find the spatial chunks for each
c
!$OMP PARALLEL default(shared) private(i,nearpt,abound,cid,cbound,
!$OMP& negx,negy,negz,posx,posy,posz,midx,midy,midz,k
!$OMP& )
!$OMP DO
      do i = 1, n
         nearpt(1) = igrid(1,i) + grdoff
         nearpt(2) = igrid(2,i) + grdoff
         nearpt(3) = igrid(3,i) + grdoff
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
         pmetable(i,k) = 1
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
         if (midx.and.midy .and. negz)  call setchunk_omp (i,cid,0,0,-1)
         if (midx.and.midy .and. posz)  call setchunk_omp (i,cid,0,0,1)
         if (midx.and.negy .and. midz)  call setchunk_omp (i,cid,0,-1,0)
         if (midx.and.posy .and. midz)  call setchunk_omp (i,cid,0,1,0)
         if (negx.and.midy .and. midz)  call setchunk_omp (i,cid,-1,0,0)
         if (posx.and.midy .and. midz)  call setchunk_omp (i,cid,1,0,0)
         if (midx.and.negy .and. negz) call setchunk_omp (i,cid,0,-1,-1)
         if (midx.and.negy .and. posz)  call setchunk_omp (i,cid,0,-1,1)
         if (midx.and.posy .and. negz)  call setchunk_omp (i,cid,0,1,-1)
         if (midx.and.posy .and. posz)  call setchunk_omp (i,cid,0,1,1)
         if (negx.and.midy .and. negz) call setchunk_omp (i,cid,-1,0,-1)
         if (negx.and.midy .and. posz)  call setchunk_omp (i,cid,-1,0,1)
         if (posx.and.midy .and. negz)  call setchunk_omp (i,cid,1,0,-1)
         if (posx.and.midy .and. posz)  call setchunk_omp (i,cid,1,0,1)
         if (negx.and.negy .and. midz) call setchunk_omp (i,cid,-1,-1,0)
         if (negx.and.posy .and. midz)  call setchunk_omp (i,cid,-1,1,0)
         if (posx.and.negy .and. midz)  call setchunk_omp (i,cid,1,-1,0)
         if (posx.and.posy .and. midz)  call setchunk_omp (i,cid,1,1,0)
         if (negx.and.negy.and. negz) call setchunk_omp (i,cid,-1,-1,-1)
         if (negx.and.negy .and. posz) call setchunk_omp (i,cid,-1,-1,1)
         if (negx.and.posy .and. negz) call setchunk_omp (i,cid,-1,1,-1)
         if (posx.and.negy .and. negz) call setchunk_omp (i,cid,1,-1,-1)
         if (negx.and.posy .and. posz)  call setchunk_omp (i,cid,-1,1,1)
         if (posx.and.negy .and. posz)  call setchunk_omp (i,cid,1,-1,1)
         if (posx.and.posy .and. negz)  call setchunk_omp (i,cid,1,1,-1)
         if (posx.and.posy .and. posz)  call setchunk_omp (i,cid,1,1,1)
   10    continue
      end do
!$OMP END DO
!$OMP END PARALLEL

      !print*,"Successful end of table_fill_omp"
      return
      end
c
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setchunk_omp  --  site overlaps neighboring chunk  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setchunk_omp" marks a chunk in the PME spatial table which is
c     overlapped by the B-splines for an electrostatic site
c
c
      subroutine setchunk_omp (i,cid,off1,off2,off3)
      use sizes
      use chunks
      use pme
      implicit none
      integer i,k
      integer off1,off2,off3
      integer cid(3),temp(3)
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
      pmetable(i,k) = 1
      return
      end

