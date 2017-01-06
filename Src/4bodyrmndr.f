
c!$OMP PARALLEL default(private) shared(imol,ep3bt,
c!$OMP& virep3bt,dep3bt,start,last,
c!$OMP& nmol)
c!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
c!$OMP& schedule(guided)
c      do moli1=start,last
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
        do moli2=moli1rmndr+1,nmol
          np1=3
          np2=6
          np3=9
          pnum(1)=imol(1,moli1rmndr)
          pnum(2)=imol(1,moli1rmndr)+1
          pnum(3)=imol(2,moli1rmndr)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)

           call empole1a_mod(npole3b,pnum,eptemp,deptemp,virtemp)

              ep2moli12=eptemp
              ep3bt=ep3bt+ep2moli12

              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep2moli12(j,i)=deptemp(j,i)
                  
                  dep3bt(j,i) = dep3bt(j,i)+dep2moli12(j,i)
                end do
              end do
              do i=1,3
                do j=1,3
                 vir2moli12(j,i)=virtemp(j,i)
                 virep3bt(j,i)=virep3bt(j,i)+vir2moli12(j,i)
                end do
              end do


          do moli3=moli2+1,nmol
            npole3b=9
            pnum(1)=imol(1,moli1rmndr)
            pnum(2)=imol(1,moli1rmndr)+1
            pnum(3)=imol(2,moli1rmndr)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)

            
                 call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ntript=ntript+1
                  ep3moli123=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep3moli123(1,i) = deptemp(1,i)
                    dep3moli123(2,i) = deptemp(2,i)
                    dep3moli123(3,i) = deptemp(3,i)
                  end do
                  do i=1,3
                    do j=1,3
                       vir3moli123(j,i)=virtemp(j,i)
                    end do
                  end do
            


               npole3b=6
               ep3moli123=ep3moli123-ep2moli12
               do l1 = 1, npole3b
                  i = pnum(l1)
                  do j = 1, 3
                    dep3moli123(j,i)=dep3moli123(j,i)-dep2moli12(j,i)
                  end do
               end do
               do i=1,3
                 do j=1,3
                   vir3moli123(j,i)=vir3moli123(j,i)-vir2moli12(j,i)
                 end do
               end do

                  npole3b=6
                  np1 =3
                  pnum(1)=imol(1,moli2)
                  pnum(2)=imol(1,moli2)+1
                  pnum(3)=imol(2,moli2)
                  pnum(4)=imol(1,moli3)
                  pnum(5)=imol(1,moli3)+1
                  pnum(6)=imol(2,moli3)              
                  call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &            virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     ep2moli23=eptemp
                     do l1 = 1, npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1, 3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli23(j,i)=deptemp(j,i)
                       end do
                     end do
                     do i=1,3
                       do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       vir2moli23(j,i)=virtemp(j,i)
                       end do
                     end do

                 np1=3
                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                     ep3moli123 = ep3moli123 - eptemp
                     ep2moli13=eptemp
                     do l1 = 1, np1
                      i = pnum(l1)
                      do j = 1, 3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli13(j,i)=deptemp(j,i)
                      end do
                     end do
                     do l1=np1+1,npole3b
                       i = pnum(l1)
                       l2=l1+np1
                       do j = 1,3
                       dep3moli123(j,i)=dep3moli123(j,i)-deptemp(j,i)
                       dep2moli13(j,i)=deptemp(j,i)
                       end do
                     end do
                     do i=1,3
                      do j=1,3
                       vir3moli123(j,i)=vir3moli123(j,i)-virtemp(j,i)
                       vir2moli13(j,i)=virtemp(j,i)
                      end do
                     end do

               npole3b=9
               pnum(1)=imol(1,moli1rmndr)
               pnum(2)=imol(1,moli1rmndr)+1
               pnum(3)=imol(2,moli1rmndr)
               pnum(4)=imol(1,moli2)
               pnum(5)=imol(1,moli2)+1
               pnum(6)=imol(2,moli2)
               pnum(7)=imol(1,moli3)
               pnum(8)=imol(1,moli3)+1
               pnum(9)=imol(2,moli3)

                  ep3bt = ep3bt + ep3moli123
                  do l1 = 1, npole3b
                    i = pnum(l1)
                    do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep3moli123(j,i)
                    end do
                  end do
                  do i=1,3
                    do j=1,3
                    virep3bt(j,i)=virep3bt(j,i)+vir3moli123(j,i)
                    end do
                  end do
            
             do moli4=moli3+1,nmol
                npole3b=12
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli3)
                pnum(8)=imol(1,moli3)+1
                pnum(9)=imol(2,moli3)
                pnum(10)=imol(1,moli4)
                pnum(11)=imol(1,moli4)+1
                pnum(12)=imol(2,moli4)
                 
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                  ep4moli1234=eptemp
                  do l1 = 1,npole3b
                    i = pnum(l1)
                    dep4moli1234(1,i) = deptemp(1,i)
                    dep4moli1234(2,i) = deptemp(2,i)
                    dep4moli1234(3,i) = deptemp(3,i)
                  end do
                  do i=1,3
                     do j=1,3
                        vir4moli1234(j,i)=virtemp(j,i)
                     end do
                  end do
                
                npole3b=9
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                ep4moli1234=ep4moli1234-eptemp
                do l1=1,npole3b
                   i=pnum(l1)
                   do j=1,3
                    dep4moli1234(j,i)=dep4moli1234(j,i)
     &                              -deptemp(j,i)
                   end do
                end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do 

                npole3b=9
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3             
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                       -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do


                npole3b=9
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli3)
                pnum(5)=imol(1,moli3)+1
                pnum(6)=imol(2,moli3)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                    -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do

                npole3b=9
                pnum(1)=imol(1,moli2)
                pnum(2)=imol(1,moli2)+1
                pnum(3)=imol(2,moli2)
                pnum(4)=imol(1,moli3)
                pnum(5)=imol(1,moli3)+1
                pnum(6)=imol(2,moli3)
                pnum(7)=imol(1,moli4)
                pnum(8)=imol(1,moli4)+1
                pnum(9)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234-eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                      -deptemp(j,i)
                    end do
                 end do
                do i=1,3
                   do j=1,3
                    vir4moli1234(j,i)=vir4moli1234(j,i)-virtemp(j,i)
                   end do
                end do

                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli2)
                 pnum(5)=imol(1,moli2)+1
                 pnum(6)=imol(2,moli2)
                 ep4moli1234=ep4moli1234+ep2moli12
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                      +dep2moli12(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli12(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 ep4moli1234=ep4moli1234+ep2moli13
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +dep2moli13(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli13(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 pnum(4)=imol(1,moli3)
                 pnum(5)=imol(1,moli3)+1
                 pnum(6)=imol(2,moli3)
                 ep4moli1234=ep4moli1234+ep2moli23
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +dep2moli23(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +vir2moli23(j,i)
                    end do
                 end do
 

                 npole3b=6
                 pnum(1)=imol(1,moli1rmndr)
                 pnum(2)=imol(1,moli1rmndr)+1
                 pnum(3)=imol(2,moli1rmndr)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli2)
                 pnum(2)=imol(1,moli2)+1
                 pnum(3)=imol(2,moli2)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                     +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                 npole3b=6
                 pnum(1)=imol(1,moli3)
                 pnum(2)=imol(1,moli3)+1
                 pnum(3)=imol(2,moli3)
                 pnum(4)=imol(1,moli4)
                 pnum(5)=imol(1,moli4)+1
                 pnum(6)=imol(2,moli4)
                call empole1a_mod(npole3b,pnum,eptemp,deptemp,
     &           virtemp)
                 ep4moli1234=ep4moli1234+eptemp
                 do l1=1,npole3b
                    i=pnum(l1)
                    do j=1,3
                     dep4moli1234(j,i)=dep4moli1234(j,i)
     &                    +deptemp(j,i)
                    end do
                 end do
                 do i=1,3
                    do j=1,3
                      vir4moli1234(j,i)=vir4moli1234(j,i)
     &                  +virtemp(j,i)
                    end do
                 end do

                npole3b=12
                pnum(1)=imol(1,moli1rmndr)
                pnum(2)=imol(1,moli1rmndr)+1
                pnum(3)=imol(2,moli1rmndr)
                pnum(4)=imol(1,moli2)
                pnum(5)=imol(1,moli2)+1
                pnum(6)=imol(2,moli2)
                pnum(7)=imol(1,moli3)
                pnum(8)=imol(1,moli3)+1
                pnum(9)=imol(2,moli3)
                pnum(10)=imol(1,moli4)
                pnum(11)=imol(1,moli4)+1
                pnum(12)=imol(2,moli4)
                ep3bt = ep3bt + ep4moli1234
                do l1 = 1, npole3b
                    i = pnum(l1)
                   do j = 1, 3
                      dep3bt(j,i) = dep3bt(j,i)+dep4moli1234(j,i)
                   end do
                end do
                do i=1,3
                   do j=1,3
                      virep3bt(j,i)=virep3bt(j,i)+vir4moli1234(j,i)
                   end do
                end do

             end do
          end do
        end do 
      end do  
c!$OMP END DO
c!$OMP END PARALLEL

c      print*,"Innerloop= ep3bt=",ep3bt
c       print*,"End inner1",ep3bt
      return
      end

