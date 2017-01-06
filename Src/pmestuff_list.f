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
      subroutine table_fill_list
      use sizes
      use atoms
      use chunks
      use pme
      use chunklight
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
c      print*,"nchunk in table_fill",nchunk
c      print*,"grdoff in table_fill",grdoff
c      print*,"nlpts nrpts in table_fill",nlpts,nrpts
c      print*,"nfft123 in table_fill",nfft1,nfft2,nfft3
 
      do k = 1, nchunk
         numpolechunk(k)=0
         do i = 1, n
            pmetable(i,k) = 0
         end do
      end do
c
c     loop over sites to find the spatial chunks for each
c
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
         numpolechunk(k)=numpolechunk(k)+1
         chunksitelist(numpolechunk(k),k)=i
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
        if (midx .and. midy .and. negz) call setchunklist (i,cid,0,0,-1)
         if (midx .and. midy .and. posz) call setchunklist (i,cid,0,0,1)
        if (midx .and. negy .and. midz) call setchunklist (i,cid,0,-1,0)
         if (midx .and. posy .and. midz) call setchunklist (i,cid,0,1,0)
        if (negx .and. midy .and. midz) call setchunklist (i,cid,-1,0,0)
         if (posx .and. midy .and. midz) call setchunklist (i,cid,1,0,0)
       if (midx .and. negy .and. negz) call setchunklist (i,cid,0,-1,-1)
        if (midx .and. negy .and. posz) call setchunklist (i,cid,0,-1,1)
        if (midx .and. posy .and. negz) call setchunklist (i,cid,0,1,-1)
         if (midx .and. posy .and. posz) call setchunklist (i,cid,0,1,1)
       if (negx .and. midy .and. negz) call setchunklist (i,cid,-1,0,-1)
        if (negx .and. midy .and. posz) call setchunklist (i,cid,-1,0,1)
        if (posx .and. midy .and. negz) call setchunklist (i,cid,1,0,-1)
         if (posx .and. midy .and. posz) call setchunklist (i,cid,1,0,1)
       if (negx .and. negy .and. midz) call setchunklist (i,cid,-1,-1,0)
        if (negx .and. posy .and. midz) call setchunklist (i,cid,-1,1,0)
        if (posx .and. negy .and. midz) call setchunklist (i,cid,1,-1,0)
         if (posx .and. posy .and. midz) call setchunklist (i,cid,1,1,0)
      if (negx .and. negy .and. negz) call setchunklist (i,cid,-1,-1,-1)
       if (negx .and. negy .and. posz) call setchunklist (i,cid,-1,-1,1)
       if (negx .and. posy .and. negz) call setchunklist (i,cid,-1,1,-1)
       if (posx .and. negy .and. negz) call setchunklist (i,cid,1,-1,-1)
        if (negx .and. posy .and. posz) call setchunklist (i,cid,-1,1,1)
        if (posx .and. negy .and. posz) call setchunklist (i,cid,1,-1,1)
        if (posx .and. posy .and. negz) call setchunklist (i,cid,1,1,-1)
         if (posx .and. posy .and. posz) call setchunklist (i,cid,1,1,1)
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
      subroutine setchunklist (i,cid,off1,off2,off3)
      use sizes
      use chunks
      use pme
      use chunklight
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

       numpolechunk(k)=numpolechunk(k)+1
       chunksitelist(numpolechunk(k),k)=i
      
      pmetable(i,k) = 1
      return
      end
c
c
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine grid_mpole  --  put multipoles on PME grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "grid_mpole" places the fractional atomic multipoles onto
c     the particle mesh Ewald grid
c
c
      subroutine grid_mpole_list (fmp)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
      use chunklight
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 v2,u2,t2
      real*8 term0,term1,term2
      real*8 fmp(10,*)
      integer iter
c
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,v2,u0,u1,u2,term0,term1,term2,t0,t1,t2,iter)
!$OMP DO
c
c     put the permanent multipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
c         do isite = 1, npole
         do iter = 1,numpolechunk(ichk)
            isite=chunksitelist(iter,ichk)
            iatm = ipole(isite)
c            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  v2 = thetai3(3,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
                     u2 = thetai2(3,m,iatm)
                     term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
     &                     + fmp(4,isite)*u0*v1 + fmp(6,isite)*u2*v0
     &                     + fmp(7,isite)*u0*v2 + fmp(10,isite)*u1*v1
                     term1 = fmp(2,isite)*u0*v0 + fmp(8,isite)*u1*v0
     &                          + fmp(9,isite)*u0*v1
                     term2 = fmp(5,isite) * u0 * v0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
                        t2 = thetai1(3,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
     &                                      + term1*t1 + term2*t2
                     end do
                  end do
               end do
c            end if
         end do
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
