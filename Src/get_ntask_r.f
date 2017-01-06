      subroutine get_numtasks_emreal
      use keys
      use mpidat
      implicit none
      integer i,next
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string

      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:12) .eq. 'NTASKEMREAL ') then
            read (string,*,err=10,end=10) numtasks_emreal 
         end if
   10    continue
      end do

      return 
      end  
