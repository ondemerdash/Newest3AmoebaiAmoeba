
      subroutine gradient_covalent(energy,derivs)
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use deriv3b
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical first
      save first
      data first  / .true. /


c      allocate (dep3b(3,npole))
c      ep3b2=0.0d0
c      ep3b3=0.0d0
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eopb = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      em = 0.0d0

      ep3b=0.0d0

      if (first) then
         first = .false.
c         if (.not. allocated(dep))  allocate (dep(3,n))
         if (.not. allocated(desum))  allocate (desum(3,n))
         if (.not. allocated(deb))  allocate (deb(3,n))
         if (.not. allocated(dea))  allocate (dea(3,n))
         if (.not. allocated(deba))  allocate (deba(3,n))
         if (.not. allocated(deub))  allocate (deub(3,n))
         if (.not. allocated(deopb))  allocate (deopb(3,n))
         if (.not. allocated(det))  allocate (det(3,n))
         if (.not. allocated(dept))  allocate (dept(3,n))
         if (.not. allocated(dett))  allocate (dett(3,n))
         if (.not. allocated(dev))  allocate (dev(3,n))
         if (.not. allocated(dem))  allocate (dem(3,n))
         if (.not. allocated(dep3b)) allocate (dep3b(3,n))
      end if

        do i=1,3
          do j=1,3
c            virep3b2(i,j)=0.0d0
c            virep3b3(i,j)=0.0d0
            virep3b(i,j)=0.0d0

            vir(i,j)=0.0d0
c            viremreal(i,j)=0.0d0
c            virev(i,j)=0.0d0
          end do
        end do

      do i = 1, n
         do j=1,3
c            dev_total(j,i)=0.0d0
c            dep(j,i) = 0.0d0
c            dep3b2(j,i) = 0.0d0
c            dep3b3(j,i) = 0.0d0
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep3b(j,i) = 0.0d0

c            demreal(j,i) =0.0d0
         end do
      end do

c       call nblist_perm_2bod

c
c     call the van der Waals energy component routines
c
      call ebond1
      call eangle1
      call estrbnd1
      call eurey1
      call eopbend1
      call etors1
      call epitors1
      call etortor1
c      call ehal1
c         if (vdwtyp .eq. 'BUFFERED-14-7') call ehal1c_new
c      print*,"After vdw"

c        call empole1c_3b_Perm_PME
c        call empole1d_3b_Perm_PME
c
c     call any miscellaneous energy component routines
c

c      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
c     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
c     &          + ep + er + es + elf + eg + ex
c      energy = esum
c      do i = 1, n
c         do j = 1, 3
c            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
c     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
c     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
c     &                      + det(j,i) + dept(j,i) + debt(j,i)
c     &                      + dett(j,i) + dev(j,i) + dec(j,i)
c     &                      + decd(j,i) + ded(j,i) + dem(j,i)
c     &                      + dep(j,i) + der(j,i) + des(j,i)
c     &                      + delf(j,i) + deg(j,i) + dex(j,i)
c         end do
c      end do

      esum =eb + ea + eba + eub + eopb + et + ept  + ett + ev + em 
      energy = esum
c      print*,"em ev",em,ev
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i) + deopb(j,i)
     &                   + det(j,i) + dept(j,i) + dett(j,i)
     &                   + dev(j,i) + dem(j,i) 
            derivs(j,i) = desum(j,i)
         end do
      end do
      return
      end
