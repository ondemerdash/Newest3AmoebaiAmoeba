      subroutine empole1a_3b_Polar_new_subtract_1bodymut(npole3b,pnum,
     &   eptemp,deptemp,virtemp,order)
      implicit none
      integer npole3b,pnum(*)
      real*8 eptemp,deptemp(3,npole3b),virtemp(3,3)
      real*8 eptemp1,deptemp1(3,npole3b),virtemp1(3,3)
      integer order,npole3b1,pnum1(npole3b),pnum1ref(npole3b)
      real*8 field(3,npole3b),fieldp(3,npole3b)
      integer i,l1,l2
 
      call empole1a_3b_Polar_rtrnfield(npole3b,pnum,eptemp,
     &   deptemp,virtemp,field,fieldp)
     
      if(order.eq.2) then
        npole3b1=3
        do l1=1,npole3b1
           pnum1(l1)=pnum(l1)          
           pnum1ref(l1)=l1
        end do
       
        call empole1a_3b_Polar_totfieldnpolevac(npole3b,pnum,
     &  npole3b1,pnum1,pnum1ref,eptemp1,
     &   deptemp1,virtemp1,field,fieldp)   
 
         eptemp=eptemp - eptemp1
         
         do l1=1,npole3b
            deptemp(1,l1) = deptemp(1,l1) - deptemp1(1,l1)
            deptemp(2,l1) = deptemp(2,l1) - deptemp1(2,l1)
            deptemp(3,l1) = deptemp(3,l1) - deptemp1(3,l1)
         end do

         do i=1,3
            virtemp(1,i)=virtemp(1,i) - virtemp1(1,i)
            virtemp(2,i)=virtemp(2,i) - virtemp1(2,i)
            virtemp(3,i)=virtemp(3,i) - virtemp1(3,i)
         end do

         do l1=1,npole3b1
               l2=l1+npole3b1
            pnum1(l1)=pnum(l2)
            pnum1ref(l1)=l2
         end do

        call empole1a_3b_Polar_totfieldnpolevac(npole3b,pnum,
     &  npole3b1,pnum1,pnum1ref,eptemp1,
     &   deptemp1,virtemp1,field,fieldp)

         eptemp=eptemp - eptemp1
         
         do l1=1,npole3b
            deptemp(1,l1) = deptemp(1,l1) - deptemp1(1,l1)
            deptemp(2,l1) = deptemp(2,l1) - deptemp1(2,l1)
            deptemp(3,l1) = deptemp(3,l1) - deptemp1(3,l1)
         end do

         do i=1,3
            virtemp(1,i)=virtemp(1,i) - virtemp1(1,i)
            virtemp(2,i)=virtemp(2,i) - virtemp1(2,i)
            virtemp(3,i)=virtemp(3,i) - virtemp1(3,i)
         end do

      else

        npole3b1=3
        do l1=1,npole3b1
           pnum1(l1)=pnum(l1)
           pnum1ref(l1)=l1
        end do

        call empole1a_3b_Polar_totfieldnpolevac(npole3b,pnum,
     &  npole3b1,pnum1,pnum1ref,eptemp1,
     &   deptemp1,virtemp1,field,fieldp)

         eptemp=eptemp - eptemp1

         do l1=1,npole3b
            deptemp(1,l1) = deptemp(1,l1) - deptemp1(1,l1)
            deptemp(2,l1) = deptemp(2,l1) - deptemp1(2,l1)
            deptemp(3,l1) = deptemp(3,l1) - deptemp1(3,l1)
         end do

         do i=1,3
            virtemp(1,i)=virtemp(1,i) - virtemp1(1,i)
            virtemp(2,i)=virtemp(2,i) - virtemp1(2,i)
            virtemp(3,i)=virtemp(3,i) - virtemp1(3,i)
         end do

         do l1=1,npole3b1
               l2=l1+npole3b1
            pnum1(l1)=pnum(l2)
            pnum1ref(l1)=l2
         end do

        call empole1a_3b_Polar_totfieldnpolevac(npole3b,pnum,
     &  npole3b1,pnum1,pnum1ref,eptemp1,
     &   deptemp1,virtemp1,field,fieldp)

         eptemp=eptemp - eptemp1

         do l1=1,npole3b
            deptemp(1,l1) = deptemp(1,l1) - deptemp1(1,l1)
            deptemp(2,l1) = deptemp(2,l1) - deptemp1(2,l1)
            deptemp(3,l1) = deptemp(3,l1) - deptemp1(3,l1)
         end do

         do i=1,3
            virtemp(1,i)=virtemp(1,i) - virtemp1(1,i)
            virtemp(2,i)=virtemp(2,i) - virtemp1(2,i)
            virtemp(3,i)=virtemp(3,i) - virtemp1(3,i)
         end do

         do l1=1,npole3b1
               l2=l1+npole3b1+npole3b1
            pnum1(l1)=pnum(l2)
            pnum1ref(l1)=l2
         end do

        call empole1a_3b_Polar_totfieldnpolevac(npole3b,pnum,
     &  npole3b1,pnum1,pnum1ref,eptemp1,
     &   deptemp1,virtemp1,field,fieldp)

         eptemp=eptemp - eptemp1

         do l1=1,npole3b
            deptemp(1,l1) = deptemp(1,l1) - deptemp1(1,l1)
            deptemp(2,l1) = deptemp(2,l1) - deptemp1(2,l1)
            deptemp(3,l1) = deptemp(3,l1) - deptemp1(3,l1)
         end do

         do i=1,3
            virtemp(1,i)=virtemp(1,i) - virtemp1(1,i)
            virtemp(2,i)=virtemp(2,i) - virtemp1(2,i)
            virtemp(3,i)=virtemp(3,i) - virtemp1(3,i)
         end do
      end if
      return
      end 
