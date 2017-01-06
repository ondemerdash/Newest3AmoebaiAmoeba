
      subroutine induce0a_3b_PolelecOnly_totfieldnpolevac(npole3b,
     &  pnum,npole3b1,pnum1,pnum1ref,uind,uinp,field,fieldp)
      use sizes
      use atoms
      use polar, only: polarity, thole, pdamp
      use totfield
      use mpole
      use polpot
      implicit none
      integer i,j,k,m
      integer ii,kk
      integer npole3b,pnum(*),npole3b1,pnum1(*),pnum1ref(*)
      real*8 field(3,npole3b)
      real*8 fieldp(3,npole3b)
      character*6 mode
      integer l1,l2,l3,k1,k2,i1,i2
      real*8 M_tot(3*npole3b1,3*npole3b1)
      real*8 uind(3,npole3b1)
      real*8 uinp(3,npole3b1)
c      real*8 M_tot(3*npole,3*npole)
c      real*8 uind(3,npole)
c      real*8 uinp(3,npole)
      logical doi
c        print*,"npole3b=",npole3b
c        do l1=1,npole3b
c           do j=1,3
c              field(j,l1)=0.0d0
c              fieldp(j,l1)=0.0d0
c           end do
c        end do

         do l1 = 1, 3*npole3b1
            do l3 = 1, 3*npole3b1
               M_tot(l1,l3) = 0.0d0
            end do
         end do

c         do l1 = 1, 3*npole
c            do l3 = 1, 3*npole
c               M_tot(l1,l3) = 0.0d0
c            end do
c         end do

         do l1 = 1, npole3b1
c         do l1 = 1,npole
            do j = 1, 3
               uind(j,l1) = 0.0d0
               uinp(j,l1) = 0.0d0
            end do
         end do


c       call totfield_noewald_umutual_rl_3b (field,fieldp,M_tot,
c     &   npole3b,pnum)

c       call dfield0a_totfieldnpole3b (field,npole3b,pnum)

c       call dfield0a_totfieldnpole3bmod (field,npole3b,pnum)
       
c       call totfieldnpole_noewald_umutual_rl_3b (field,fieldp,
c     &   M_tot,npole3b,pnum)

       !print*,"poltyp=",poltyp
      if (poltyp .eq. 'MUTUAL') then
    
       !call totfieldnpole3b_noewald_umutual_rl_3b (field,fieldp,
     & !  M_tot,npole3b,pnum)

       call totfieldnpole3b_noewald_umutual_rl_3b_nofield(
     &   M_tot,npole3b1,pnum1)

         do l1 = 1, npole3b1
            i = pnum1(l1)
c         do i =1,npole
c            doi=.false.
c            do k1=1,npole3b
c              if(pnum(k1).eq.i) then
c                doi=.true.
c                l1=k1
c                goto 21
c              end if
c            end do
c   21    continue
c            if(doi) then
            do j = 1, 3
               i1 = 3*(l1-1)+j
c               i1 = 3*(i-1)+j
               M_tot(i1,i1) = M_tot(i1,i1)+1.0d0/polarity(i)
            end do
c            end if
         end do

         call invert(3*npole3b1,M_tot)
         

         do l1 = 1, npole3b1
          !  i = pnum1(l1)
            do i1 = 1, 3
               i2 = 3*(l1-1) + i1
               do l3 = 1, npole3b1
               !   k = pnum1(l3)
                  k=pnum1ref(l3)
                  k2 = 3*(l3-1)

                  uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*field(1,k)+
     &                                      M_tot(i2,k2+2)*field(2,k)+
     &                                      M_tot(i2,k2+3)*field(3,k)
c                  uind(i1,i)=uind(i1,i)+ M_tot(i2,k2+1)*field(1,l3)+
c     &                                      M_tot(i2,k2+2)*field(2,l3)+
c     &                                      M_tot(i2,k2+3)*field(3,l3)
c               uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*fieldnpole(1,k)+
c     &                                  M_tot(i2,k2+2)*fieldnpole(2,k)+
c     &                                  M_tot(i2,k2+3)*fieldnpole(3,k)
                  uinp(i1,l1)=uinp(i1,l1)+ M_tot(i2,k2+1)*fieldp(1,k)+
     &                                     M_tot(i2,k2+2)*fieldp(2,k)+
     &                                     M_tot(i2,k2+3)*fieldp(3,k)

c                  uinp(i1,i)=uinp(i1,i)+ M_tot(i2,k2+1)*fieldp(1,l3)+
c     &                                     M_tot(i2,k2+2)*fieldp(2,l3)+
c     &                                     M_tot(i2,k2+3)*fieldp(3,l3)
c               uinp(i1,l1)=uinp(i1,l1)+M_tot(i2,k2+1)*fieldpnpole(1,k)+
c     &                                 M_tot(i2,k2+2)*fieldpnpole(2,k)+
c     &                                 M_tot(i2,k2+3)*fieldpnpole(3,k)

               end do
            end do
         end do
      else
         do l1 = 1, npole3b
            i = pnum(l1)
            do j=1,3
               uind(j,l1)=polarity(i)*fieldnpole(j,i)
               uinp(j,l1)=polarity(i)*fieldpnpole(j,i)
            end do
         end do
      end if

      return
      end
c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
