      subroutine induce0a_3b_cholesky_nopriordir(npole3b,pnum,
     & uind,uinp)
      use sizes
      use atoms
      use polar, only: polarity, thole, pdamp
      use polpot
      implicit none
      integer i,j,k,m
      integer ii,kk
      integer npole3b,pnum(*)
      real*8 field(3,npole3b),off3b
      real*8 fieldp(3,npole3b)
      character*6 mode
      integer l1,l2,l3,k1,k2,i1,i2
      real*8 M_tot(3*npole3b,3*npole3b)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 M_tot_column((9*npole3b*npole3b +3*npole3b)/2)
      real*8 fdir_column(3*npole3b)
      integer counter
c        print*,"npole3b=",npole3b

        counter=1

        do l1=1,npole3b
           do j=1,3
              field(j,l1)=0.0d0
              fieldp(j,l1)=0.0d0
               uind(j,l1) = 0.0d0
               uinp(j,l1) = 0.0d0
            fdir_column(counter)=0.0d0
            counter=counter+1
           end do
        end do

         do l1 = 1, 3*npole3b
            do l3 = 1, 3*npole3b
               M_tot(l1,l3) = 0.0d0
            end do
         end do

      counter=1

      do i = 1, 3*npole3b
         do j = i, 3*npole3b
            M_tot_column(counter) = 0.0d0
            counter=counter+1
         end do
      end do

       call field_noewald_umutual_rl_3b (field,fieldp,M_tot,
     &   npole3b,pnum)

      if (poltyp .eq. 'MUTUAL') then
         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               i1 = 3*(l1-1)+j
               M_tot(i1,i1) = M_tot(i1,i1)+1.0d0/polarity(i)
            end do
c             print*,"polarity(i)",polarity(i)
         end do

c         call invert(3*npole3b,M_tot)
c
c         do l1 = 1, npole3b
c            i = pnum(l1)
c            do i1 = 1, 3
c               i2 = 3*(l1-1) + i1
c               do l3 = 1, npole3b
c                  k = pnum(l3)
c                  k2 = 3*(l3-1)
c                  uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*field(1,l3)+
c     &                                      M_tot(i2,k2+2)*field(2,l3)+
c     &                                      M_tot(i2,k2+3)*field(3,l3)
c                  uinp(i1,l1)=uinp(i1,l1)+ M_tot(i2,k2+1)*fieldp(1,l3)+
c     &                                     M_tot(i2,k2+2)*fieldp(2,l3)+
c     &                                     M_tot(i2,k2+3)*fieldp(3,l3)
c
c               end do
c            end do
c         end do
         counter=1

         do l1 = 1, npole3b
c            i = pnum(l1)
            do j = 1, 3
c               fdir_column(counter)=field(j,i)
               fdir_column(counter)=field(j,l1)
               counter=counter+1
            end do
         end do

         counter=1

         do i = 1, 3*npole3b
            do j = i, 3*npole3b
               M_tot_column(counter) = M_tot(i,j)
               counter=counter+1
            end do
         end do

         call cholesky(3*npole3b,M_tot_column,fdir_column)

         do l1 = 1, npole3b
c            i = pnum(l1)
            do j = 1, 3
               i1=3*(l1-1)+j
c               uind(j,i)=fdir_column(i1)
               uind(j,l1)=fdir_column(i1)
               uinp(j,l1)=uind(j,l1)
            end do
         end do


      else
         do l1 = 1, npole3b
            i = pnum(l1)
            do j=1,3
               uind(j,l1)=polarity(i)*field(j,l1)
               uinp(j,l1)=polarity(i)*fieldp(j,l1)
            end do
         end do
      end if

      return
      end
c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
