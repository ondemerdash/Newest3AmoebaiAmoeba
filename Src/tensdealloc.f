      subroutine tensdealloc
      use dEtensor2
      use dEtensor3
      implicit none

c         if(allocated(elstgrad)) deallocate(elstgrad)
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
c         if(allocated(elsttau1)) deallocate(elsttau1)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)
c         if(allocated(elsttau2)) deallocate(elsttau2)

         if(allocated(dEd1_3)) deallocate(dEd1_3)
         if(allocated(dEd2_3)) deallocate(dEd2_3)
         if(allocated(dEp1_3)) deallocate(dEp1_3)
         if(allocated(dEp2_3)) deallocate(dEp2_3)

         if(allocated(taud1_3)) deallocate(taud1_3)
         if(allocated(taud2_3)) deallocate(taud2_3)
         if(allocated(taup1_3)) deallocate(taup1_3)
         if(allocated(taup2_3)) deallocate(taup2_3)

         if(allocated(frcztau1dtot_3)) deallocate(frcztau1dtot_3)
         if(allocated(frcytau1dtot_3)) deallocate(frcytau1dtot_3)
         if(allocated(frcxtau1dtot_3)) deallocate(frcxtau1dtot_3)
         if(allocated(frcztau1ptot_3)) deallocate(frcztau1ptot_3)
         if(allocated(frcytau1ptot_3)) deallocate(frcytau1ptot_3)
         if(allocated(frcxtau1ptot_3)) deallocate(frcxtau1ptot_3)
         if(allocated(frcztau2dtot_3)) deallocate(frcztau2dtot_3)
         if(allocated(frcytau2dtot_3)) deallocate(frcytau2dtot_3)
         if(allocated(frcxtau2dtot_3)) deallocate(frcxtau2dtot_3)
         if(allocated(frcztau2ptot_3)) deallocate(frcztau2ptot_3)
         if(allocated(frcytau2ptot_3)) deallocate(frcytau2ptot_3)
         if(allocated(frcxtau2ptot_3)) deallocate(frcxtau2ptot_3)
      return
      end
