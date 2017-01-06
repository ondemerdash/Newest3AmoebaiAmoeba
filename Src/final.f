c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine final  --  final actions before program exit  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "final" performs any final program actions such as deallocation
c     of global memory, prints a status message, and then pauses if
c     necessary to avoid closing the execution window
c
c
      subroutine final
      use sizes
      use align
      use paremneigh
      use parvdwneigh
      use analyz
      use angang
      use angbnd
      use angtor
      use atmlst
      use bitor
      use bndstr
      use charge
      use chunks
      use couple
      use deriv
      use dipole
      use disgeo
      use domega
      use faces
      use fracs
      use freeze
      use group
      use hessn
      use hpmf
      use improp
      use imptor
      use inform
      use iounit
      use light
      use merck
      use molcul
      use moldyn
      use mpole
      use mutant
      use neigh
      use nonpol
      use opbend
      use opdist
      use orbits
      use paths
      use pdb
      use piorbs
      use pistuf
      use pitors
      use pme
      use polar
      use polgrp
      use qmstuf
      use refer
      use restrn
      use rgddyn
      use rigid
      use ring
      use socket
      use solute
      use stodyn
      use strbnd
      use strtor
      use syntrn
      use tarray
      use tors
      use tortor
      use uprior
      use urey
      use usage
      use usolve
      use vdw
      use vibs
      use warp
      use neigh3b
      use deriv3b
      use ewaldneigh
      use mpidat
      use totfield
      use dEtensor
      use dEtensor2
      use dEtensor3
      use chunklight
      use neigh2
      use deriv1bmat
      use neigh2clust
      use deriv2bmat
      use uprior1b
      use uprior2b
      use boxes1bclust
      use boxes2bclust
      use cell1bclust
      use cell2bclust
      use uind1bmatrix
      use uind2bmat
      implicit none

c     deallocation of arrays for cluster method
         if(allocated(xbox2b))
     &       deallocate(xbox2b)
         if(allocated(ybox2b))
     &        deallocate(ybox2b)
         if(allocated(zbox2b))
     &        deallocate(zbox2b)
         if(allocated(xcell2b))
     &       deallocate(xcell2b)
         if(allocated(ycell2b))
     &       deallocate(ycell2b)
         if(allocated(zcell2b))
     &       deallocate(zcell2b)
         if(allocated(xbox2_2b))
     &       deallocate(xbox2_2b)
         if(allocated(ybox2_2b))
     &        deallocate(ybox2_2b)
         if(allocated(zbox2_2b))
     &        deallocate(zbox2_2b)
         if(allocated(xcell2_2b))
     &       deallocate(xcell2_2b)
         if(allocated(ycell2_2b))
     &       deallocate(ycell2_2b)
         if(allocated(zcell2_2b))
     &       deallocate(zcell2_2b)
         if(allocated(ncell2b))
     &       deallocate(ncell2b)
         if(allocated(icell2b))
     &       deallocate(icell2b)
         if(allocated(volbox2b))
     &       deallocate(volbox2b)
         if(allocated(lvec2b))
     &       deallocate(lvec2b)
         if(allocated(recip2b))
     &       deallocate(recip2b)

c     deallocation of arrays for cluster method
         if(allocated(xbox1b))
     &       deallocate(xbox1b)
         if(allocated(ybox1b))
     &        deallocate(ybox1b)
         if(allocated(zbox1b))
     &        deallocate(zbox1b)
         if(allocated(xcell1b))
     &       deallocate(xcell1b)
         if(allocated(ycell1b))
     &       deallocate(ycell1b)
         if(allocated(zcell1b))
     &       deallocate(zcell1b)
         if(allocated(xbox2_1b))
     &       deallocate(xbox2_1b)
         if(allocated(ybox2_1b))
     &        deallocate(ybox2_1b)
         if(allocated(zbox2_1b))
     &        deallocate(zbox2_1b)
         if(allocated(xcell2_1b))
     &       deallocate(xcell2_1b)
         if(allocated(ycell2_1b))
     &       deallocate(ycell2_1b)
         if(allocated(zcell2_1b))
     &       deallocate(zcell2_1b)
         if(allocated(ncell1b))
     &       deallocate(ncell1b)
         if(allocated(icell1b))
     &       deallocate(icell1b)
         if(allocated(volbox1b))
     &       deallocate(volbox1b)
         if(allocated(lvec1b))
     &       deallocate(lvec1b)
         if(allocated(recip1b))
     &       deallocate(recip1b)


        if(allocated(uind1bmat)) deallocate(uind1bmat)
        if(allocated(uind2bmatmod)) deallocate(uind2bmatmod)
        if(allocated(save2bfor3b)) deallocate(save2bfor3b)
        if (allocated(ind2b_counter)) deallocate(ind2b_counter)
        if(allocated(ep2bmatmod)) deallocate(ep2bmatmod)
        if(allocated(dep2bmatmod)) deallocate(dep2bmatmod)
        if(allocated(virep2bmatmod)) deallocate(virep2bmatmod)

        if(allocated(ep1bmat)) deallocate(ep1bmat)
        if(allocated(dep1bmat)) deallocate(dep1bmat)
        if(allocated(em1bmat)) deallocate(em1bmat)
        if(allocated(dem1bmat)) deallocate(dem1bmat)
        if(allocated(virep1bmat)) deallocate(virep1bmat)
        if(allocated(ep2bmat)) deallocate(ep2bmat)
        if(allocated(dep2bmat)) deallocate(dep2bmat)
        if(allocated(em2bmat)) deallocate(em2bmat)
        if(allocated(dem2bmat)) deallocate(dem2bmat)
        if(allocated(virep2bmat)) deallocate(virep2bmat)

        if (allocated(clust)) deallocate (clust)
        if (allocated(sizeclust)) deallocate (sizeclust)
        if (allocated(clust_cm)) deallocate (clust_cm)
        if (allocated(distmax)) deallocate (distmax)

         if(allocated(elst1b))
     &       deallocate(elst1b)
         if(allocated(nelst1b))
     &        deallocate(nelst1b)
         if(allocated(domlstclust1b))
     &        deallocate(domlstclust1b)

         if(allocated(xmold1b))
     &      deallocate(xmold1b)
         if(allocated(ymold1b))
     &      deallocate(ymold1b)
         if(allocated(zmold1b))
     &      deallocate(zmold1b)


         if(allocated(ulst1b))
     &       deallocate(ulst1b)
         if(allocated(nulst1b))
     &        deallocate(nulst1b)
         if(allocated(doulstclust1b))
     &        deallocate(doulstclust1b)

         if(allocated(xuold1b))
     &      deallocate(xuold1b)
         if(allocated(yuold1b))
     &      deallocate(yuold1b)
         if(allocated(zuold1b))
     &      deallocate(zuold1b)


         if(allocated(elst2b))
     &       deallocate(elst2b)
         if(allocated(nelst2b))
     &        deallocate(nelst2b)
         if(allocated(domlstclust2b))
     &        deallocate(domlstclust2b)

         if(allocated(xmold2b))
     &      deallocate(xmold2b)
         if(allocated(ymold2b))
     &      deallocate(ymold2b)
         if(allocated(zmold2b))
     &      deallocate(zmold2b)


         if(allocated(ulst2b))
     &       deallocate(ulst2b)
         if(allocated(nulst2b))
     &        deallocate(nulst2b)
         if(allocated(doulstclust2b))
     &        deallocate(doulstclust2b)

         if(allocated(xuold2b))
     &      deallocate(xuold2b)
         if(allocated(yuold2b))
     &      deallocate(yuold2b)
         if(allocated(zuold2b))
     &      deallocate(zuold2b)

c
c
c     deallocation of global arrays from module align
c

      if (allocated(ifit))  deallocate (ifit)
      if (allocated(wfit))  deallocate (wfit)
c
c     deallocation of global arrays from module analyz
c
      if (allocated(aesum))  deallocate (aesum)
      if (allocated(aeb))  deallocate (aeb)
      if (allocated(aea))  deallocate (aea)
      if (allocated(aeba))  deallocate (aeba)
      if (allocated(aeub))  deallocate (aeub)
      if (allocated(aeaa))  deallocate (aeaa)
      if (allocated(aeopb))  deallocate (aeopb)
      if (allocated(aeopd))  deallocate (aeopd)
      if (allocated(aeid))  deallocate (aeid)
      if (allocated(aeit))  deallocate (aeit)
      if (allocated(aet))  deallocate (aet)
      if (allocated(aept))  deallocate (aept)
      if (allocated(aebt))  deallocate (aebt)
      if (allocated(aeat))  deallocate (aeat)
      if (allocated(aett))  deallocate (aett)
      if (allocated(aev))  deallocate (aev)
      if (allocated(aec))  deallocate (aec)
      if (allocated(aecd))  deallocate (aecd)
      if (allocated(aed))  deallocate (aed)
      if (allocated(aem))  deallocate (aem)
      if (allocated(aep))  deallocate (aep)
      if (allocated(aer))  deallocate (aer)
      if (allocated(aes))  deallocate (aes)
      if (allocated(aelf))  deallocate (aelf)
      if (allocated(aeg))  deallocate (aeg)
      if (allocated(aex))  deallocate (aex)
c
c     deallocation of global arrays from module angang
c
      if (allocated(iaa))  deallocate (iaa)
      if (allocated(kaa))  deallocate (kaa)
c
c     deallocation of global arrays from module angbnd
c
      if (allocated(iang))  deallocate (iang)
      if (allocated(ak))  deallocate (ak)
      if (allocated(anat))  deallocate (anat)
      if (allocated(afld))  deallocate (afld)
c
c     deallocation of global arrays from module angtor
c
      if (allocated(iat))  deallocate (iat)
      if (allocated(kant))  deallocate (kant)
c
c     deallocation of global arrays from module atmlst
c
      if (allocated(bndlist))  deallocate (bndlist)
      if (allocated(anglist))  deallocate (anglist)
c
c     deallocation of global arrays from module bitor
c
      if (allocated(ibitor))  deallocate (ibitor)
c
c     deallocation of global arrays from module bndstr
c
      if (allocated(ibnd))  deallocate (ibnd)
      if (allocated(bk))  deallocate (bk)
      if (allocated(bl))  deallocate (bl)
c
c     deallocation of global arrays from module charge
c
      if (allocated(iion))  deallocate (iion)
      if (allocated(jion))  deallocate (jion)
      if (allocated(kion))  deallocate (kion)
      if (allocated(chglist))  deallocate (chglist)
      if (allocated(pchg))  deallocate (pchg)
c
c     deallocation of global arrays from module chunks
c
      if (allocated(pmetable))  deallocate (pmetable)

      if (allocated(numpolechunk)) deallocate (numpolechunk)
      if (allocated(chunksitelist))
     &       deallocate (chunksitelist)

c
c     deallocation of global arrays from module couple
c
      if (allocated(n13))  deallocate (n13)
      if (allocated(n14))  deallocate (n14)
      if (allocated(n15))  deallocate (n15)
      if (allocated(i13))  deallocate (i13)
      if (allocated(i14))  deallocate (i14)
      if (allocated(i15))  deallocate (i15)
c
c     deallocation of global arrays from module deriv
c
      if (allocated(desum))  deallocate (desum)
      if (allocated(deb))  deallocate (deb)
      if (allocated(dea))  deallocate (dea)
      if (allocated(deba))  deallocate (deba)
      if (allocated(deub))  deallocate (deub)
      if (allocated(deaa))  deallocate (deaa)
      if (allocated(deopb))  deallocate (deopb)
      if (allocated(deopd))  deallocate (deopd)
      if (allocated(deid))  deallocate (deid)
      if (allocated(deit))  deallocate (deit)
      if (allocated(det))  deallocate (det)
      if (allocated(dept))  deallocate (dept)
      if (allocated(debt))  deallocate (debt)
      if (allocated(deat))  deallocate (deat)
      if (allocated(dett))  deallocate (dett)
      if (allocated(dev))  deallocate (dev)
      if (allocated(dec))  deallocate (dec)
      if (allocated(decd))  deallocate (decd)
      if (allocated(ded))  deallocate (ded)
      if (allocated(dem))  deallocate (dem)
      if (allocated(demreal))  deallocate (demreal)
      if (allocated(demreal_tmp))  deallocate (demreal_tmp)
      if (allocated(dep))  deallocate (dep)
      if (allocated(der))  deallocate (der)
      if (allocated(des))  deallocate (des)
      if (allocated(delf))  deallocate (delf)
      if (allocated(deg))  deallocate (deg)
      if (allocated(dex))  deallocate (dex)
      if (allocated(dep3b))  deallocate (dep3b)

      if (allocated(dep2b)) deallocate (dep2b)
      if (allocated(uind3b)) deallocate (uind3b)
      if (allocated(fieldnpole)) deallocate (fieldnpole)
      if (allocated(fieldpnpole)) deallocate (fieldpnpole)
      if (allocated(fphi_totfield)) deallocate(fphi_totfield)
      if (allocated(fieldnpole_rcp)) deallocate (fieldnpole_rcp)
        if(allocated(dEd1)) deallocate (dEd1)
        if(allocated(dEd2)) deallocate (dEd2)
        if(allocated(dEp1)) deallocate (dEp1)
        if(allocated(dEp2)) deallocate (dEp2)
        if(allocated(dEindex)) deallocate (dEindex)

        if(allocated(dEindextmp)) deallocate (dEindextmp)
        if(allocated(dEd1tmp)) deallocate (dEd1tmp)
        if(allocated(dEd2tmp)) deallocate (dEd2tmp)
        if(allocated(dEp1tmp)) deallocate (dEp1tmp)
        if(allocated(dEp2tmp)) deallocate (dEp2tmp)
        if(allocated(tau1indextmp)) deallocate(tau1indextmp)
        if(allocated(tau2indextmp)) deallocate(tau2indextmp)
        if(allocated(taud1tmp)) deallocate (taud1tmp)
        if(allocated(taud2tmp)) deallocate (taud2tmp)
        if(allocated(taup1tmp)) deallocate (taup1tmp)
        if(allocated(taup2tmp)) deallocate (taup2tmp)

      if(allocated(frcztau1d)) deallocate (frcztau1d)
      if(allocated(frcytau1d)) deallocate (frcytau1d)
      if(allocated(frcxtau1d)) deallocate (frcxtau1d)
      if(allocated(frcztau1p)) deallocate (frcztau1p)
      if(allocated(frcytau1p)) deallocate (frcytau1p)
      if(allocated(frcxtau1p)) deallocate (frcxtau1p)
      if(allocated(frcztau2d)) deallocate (frcztau2d)
      if(allocated(frcytau2d)) deallocate (frcytau2d)
      if(allocated(frcxtau2d)) deallocate (frcxtau2d)
      if(allocated(frcztau2p)) deallocate (frcztau2p)
      if(allocated(frcytau2p)) deallocate (frcytau2p)
      if(allocated(frcxtau2p)) deallocate (frcxtau2p)

        if(allocated(ntpair_start)) deallocate(ntpair_start)
        if(allocated(ntpair_last)) deallocate(ntpair_last)
        if(allocated(ntpair_start_rmndr)) deallocate(ntpair_start_rmndr)
        if(allocated(ntpair_last_rmndr)) deallocate(ntpair_last_rmndr)

         if(allocated(elstgrad)) deallocate(elstgrad)
         if(allocated(elsttau1_index)) deallocate(elsttau1_index)
         if(allocated(elsttau1)) deallocate(elsttau1)
         if(allocated(elsttau2_index)) deallocate(elsttau2_index)
         if(allocated(elsttau2)) deallocate(elsttau2)


         if(allocated(taud1)) deallocate(taud1)
         if(allocated(taud2)) deallocate(taud2)
         if(allocated(taup1)) deallocate(taup1)
         if(allocated(taup2)) deallocate(taup2)

         if(allocated(frcztau1dtot)) deallocate(frcztau1dtot)
         if(allocated(frcytau1dtot)) deallocate(frcytau1dtot)
         if(allocated(frcxtau1dtot)) deallocate(frcxtau1dtot)
         if(allocated(frcztau1ptot)) deallocate(frcztau1ptot)
         if(allocated(frcytau1ptot)) deallocate(frcytau1ptot)
         if(allocated(frcxtau1ptot)) deallocate(frcxtau1ptot)
         if(allocated(frcztau2dtot)) deallocate(frcztau2dtot)
         if(allocated(frcytau2dtot)) deallocate(frcytau2dtot)
         if(allocated(frcxtau2dtot)) deallocate(frcxtau2dtot)
         if(allocated(frcztau2ptot)) deallocate(frcztau2ptot)
         if(allocated(frcytau2ptot)) deallocate(frcytau2ptot)
         if(allocated(frcxtau2ptot)) deallocate(frcxtau2ptot)

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

c
c     deallocation of global arrays from module dipole
c
      if (allocated(idpl))  deallocate (idpl)
      if (allocated(bdpl))  deallocate (bdpl)
      if (allocated(sdpl))  deallocate (sdpl)
c
c     deallocation of global arrays from module disgeo
c
      if (allocated(dbnd))  deallocate (dbnd)
      if (allocated(georad))  deallocate (georad)
c
c     deallocation of global arrays from module domega
c
      if (allocated(tesum))  deallocate (tesum)
      if (allocated(teb))  deallocate (teb)
      if (allocated(tea))  deallocate (tea)
      if (allocated(teba))  deallocate (teba)
      if (allocated(teub))  deallocate (teub)
      if (allocated(teaa))  deallocate (teaa)
      if (allocated(teopb))  deallocate (teopb)
      if (allocated(teopd))  deallocate (teopd)
      if (allocated(teid))  deallocate (teid)
      if (allocated(teit))  deallocate (teit)
      if (allocated(tet))  deallocate (tet)
      if (allocated(tept))  deallocate (tept)
      if (allocated(tebt))  deallocate (tebt)
      if (allocated(teat))  deallocate (teat)
      if (allocated(tett))  deallocate (tett)
      if (allocated(tev))  deallocate (tev)
      if (allocated(tec))  deallocate (tec)
      if (allocated(tecd))  deallocate (tecd)
      if (allocated(ted))  deallocate (ted)
      if (allocated(tem))  deallocate (tem)
      if (allocated(tep))  deallocate (tep)
      if (allocated(ter))  deallocate (ter)
      if (allocated(tes))  deallocate (tes)
      if (allocated(telf))  deallocate (telf)
      if (allocated(teg))  deallocate (teg)
      if (allocated(tex))  deallocate (tex)
c
c     deallocation of global arrays from module faces
c
      if (allocated(ar))  deallocate (ar)
      if (allocated(axyz))  deallocate (axyz)
      if (allocated(skip))  deallocate (skip)
      if (allocated(nosurf))  deallocate (nosurf)
      if (allocated(afree))  deallocate (afree)
      if (allocated(abur))  deallocate (abur)
      if (allocated(cls))  deallocate (cls)
      if (allocated(clst))  deallocate (clst)
      if (allocated(acls))  deallocate (acls)
      if (allocated(ttfe))  deallocate (ttfe)
      if (allocated(ttle))  deallocate (ttle)
      if (allocated(enext))  deallocate (enext)
      if (allocated(tta))  deallocate (tta)
      if (allocated(ttbur))  deallocate (ttbur)
      if (allocated(ttfree))  deallocate (ttfree)
      if (allocated(tfe))  deallocate (tfe)
      if (allocated(ta))  deallocate (ta)
      if (allocated(tr))  deallocate (tr)
      if (allocated(t))  deallocate (t)
      if (allocated(tax))  deallocate (tax)
      if (allocated(tfree))  deallocate (tfree)
      if (allocated(pa))  deallocate (pa)
      if (allocated(p))  deallocate (p)
      if (allocated(va))  deallocate (va)
      if (allocated(vp))  deallocate (vp)
      if (allocated(vxyz))  deallocate (vxyz)
      if (allocated(env))  deallocate (env)
      if (allocated(fnen))  deallocate (fnen)
      if (allocated(ca))  deallocate (ca)
      if (allocated(ct))  deallocate (ct)
      if (allocated(cr))  deallocate (cr)
      if (allocated(c))  deallocate (c)
      if (allocated(epc))  deallocate (epc)
      if (allocated(epv))  deallocate (epv)
      if (allocated(afe))  deallocate (afe)
      if (allocated(ale))  deallocate (ale)
      if (allocated(epnext))  deallocate (epnext)
      if (allocated(fsen))  deallocate (fsen)
      if (allocated(fsep))  deallocate (fsep)
      if (allocated(cynep))  deallocate (cynep)
      if (allocated(cyep))  deallocate (cyep)
      if (allocated(fpa))  deallocate (fpa)
      if (allocated(fpncy))  deallocate (fpncy)
      if (allocated(fpcy))  deallocate (fpcy)
c
c     deallocation of global arrays from module fracs
c
      if (allocated(xfrac))  deallocate (xfrac)
      if (allocated(yfrac))  deallocate (yfrac)
      if (allocated(zfrac))  deallocate (zfrac)
c
c     deallocation of global arrays from module freeze
c
      if (allocated(iratx))  deallocate (iratx)
      if (allocated(kratx))  deallocate (kratx)
      if (allocated(irat))  deallocate (irat)
      if (allocated(krat))  deallocate (krat)
      if (allocated(ratimage))  deallocate (ratimage)
c
c     deallocation of global arrays from module group
c
      if (allocated(kgrp))  deallocate (kgrp)
      if (allocated(grplist))  deallocate (grplist)
      if (allocated(igrp))  deallocate (igrp)
      if (allocated(grpmass))  deallocate (grpmass)
      if (allocated(wgrp))  deallocate (wgrp)
c
c     deallocation of global arrays from module hessn
c
      if (allocated(hessx))  deallocate (hessx)
      if (allocated(hessy))  deallocate (hessy)
      if (allocated(hessz))  deallocate (hessz)
c
c     deallocation of global arrays from module hpmf
c
      if (allocated(ipmf))  deallocate (ipmf)
      if (allocated(rpmf))  deallocate (rpmf)
      if (allocated(acsa))  deallocate (acsa)
c
c     deallocation of global arrays from module improp
c
      if (allocated(iiprop))  deallocate (iiprop)
      if (allocated(kprop))  deallocate (kprop)
      if (allocated(vprop))  deallocate (vprop)
c
c     deallocation of global arrays from module imptor
c
      if (allocated(iitors))  deallocate (iitors)
      if (allocated(itors1))  deallocate (itors1)
      if (allocated(itors2))  deallocate (itors2)
      if (allocated(itors3))  deallocate (itors3)
c
c     deallocation of global arrays from module light
c
      if (allocated(kbx))  deallocate (kbx)
      if (allocated(kby))  deallocate (kby)
      if (allocated(kbz))  deallocate (kbz)
      if (allocated(kex))  deallocate (kex)
      if (allocated(key))  deallocate (key)
      if (allocated(kez))  deallocate (kez)
      if (allocated(locx))  deallocate (locx)
      if (allocated(locy))  deallocate (locy)
      if (allocated(locz))  deallocate (locz)
      if (allocated(rgx))  deallocate (rgx)
      if (allocated(rgy))  deallocate (rgy)
      if (allocated(rgz))  deallocate (rgz)
c
c     deallocation of global arrays from module merck
c
      if (allocated(mmff_ka))  deallocate (mmff_ka)
      if (allocated(mmff_ka1))  deallocate (mmff_ka1)
      if (allocated(mmff_ka2))  deallocate (mmff_ka2)
      if (allocated(mmff_ka3))  deallocate (mmff_ka3)
      if (allocated(mmff_ka4))  deallocate (mmff_ka4)
      if (allocated(mmff_ka5))  deallocate (mmff_ka5)
      if (allocated(mmff_ka6))  deallocate (mmff_ka6)
      if (allocated(mmff_ka7))  deallocate (mmff_ka7)
      if (allocated(mmff_ka8))  deallocate (mmff_ka8)
      if (allocated(mmff_ang0))  deallocate (mmff_ang0)
      if (allocated(mmff_ang1))  deallocate (mmff_ang1)
      if (allocated(mmff_ang2))  deallocate (mmff_ang2)
      if (allocated(mmff_ang3))  deallocate (mmff_ang3)
      if (allocated(mmff_ang4))  deallocate (mmff_ang4)
      if (allocated(mmff_ang5))  deallocate (mmff_ang5)
      if (allocated(mmff_ang6))  deallocate (mmff_ang6)
      if (allocated(mmff_ang7))  deallocate (mmff_ang7)
      if (allocated(mmff_ang8))  deallocate (mmff_ang8)
      if (allocated(stbn_abc))  deallocate (stbn_abc)
      if (allocated(stbn_cba))  deallocate (stbn_cba)
      if (allocated(stbn_abc1))  deallocate (stbn_abc1)
      if (allocated(stbn_cba1))  deallocate (stbn_cba1)
      if (allocated(stbn_abc2))  deallocate (stbn_abc2)
      if (allocated(stbn_cba2))  deallocate (stbn_cba2)
      if (allocated(stbn_abc3))  deallocate (stbn_abc3)
      if (allocated(stbn_cba3))  deallocate (stbn_cba3)
      if (allocated(stbn_abc4))  deallocate (stbn_abc4)
      if (allocated(stbn_cba4))  deallocate (stbn_cba4)
      if (allocated(stbn_abc5))  deallocate (stbn_abc5)
      if (allocated(stbn_cba5))  deallocate (stbn_cba5)
      if (allocated(stbn_abc6))  deallocate (stbn_abc6)
      if (allocated(stbn_cba6))  deallocate (stbn_cba6)
      if (allocated(stbn_abc7))  deallocate (stbn_abc7)
      if (allocated(stbn_cba7))  deallocate (stbn_cba7)
      if (allocated(stbn_abc8))  deallocate (stbn_abc8)
      if (allocated(stbn_cba8))  deallocate (stbn_cba8)
      if (allocated(stbn_abc9))  deallocate (stbn_abc9)
      if (allocated(stbn_cba9))  deallocate (stbn_cba9)
      if (allocated(stbn_abc10))  deallocate (stbn_abc10)
      if (allocated(stbn_cba10))  deallocate (stbn_cba10)
      if (allocated(stbn_abc11))  deallocate (stbn_abc11)
      if (allocated(stbn_cba11))  deallocate (stbn_cba11)
c
c     deallocation of global arrays from module molcul
c
      if (allocated(imol))  deallocate (imol)
      if (allocated(kmol))  deallocate (kmol)
      if (allocated(molcule))  deallocate (molcule)
      if (allocated(molmass))  deallocate (molmass)
c
c     deallocation of global arrays from module moldyn
c
      if (allocated(v))  deallocate (v)
      if (allocated(a))  deallocate (a)
      if (allocated(aalt))  deallocate (aalt)
c      if (allocated(vnopolz)) deallocate (vnopolz)
c      if (allocated(anopolz)) deallocate (anopolz)
c
c     deallocation of global arrays from module mpole
c
      if (allocated(ipole))  deallocate (ipole)
      if (allocated(polsiz))  deallocate (polsiz)
      if (allocated(pollist))  deallocate (pollist)
      if (allocated(zaxis))  deallocate (zaxis)
      if (allocated(xaxis))  deallocate (xaxis)
      if (allocated(yaxis))  deallocate (yaxis)
      if (allocated(pole))  deallocate (pole)
      if (allocated(rpole))  deallocate (rpole)
      if (allocated(polaxe))  deallocate (polaxe)
c
c     deallocation of global arrays from module mutant
c
      if (allocated(imut))  deallocate (imut)
      if (allocated(type0))  deallocate (type0)
      if (allocated(class0))  deallocate (class0)
      if (allocated(type1))  deallocate (type1)
      if (allocated(class1))  deallocate (class1)
      if (allocated(mut))  deallocate (mut)
c
c     deallocation of global arrays from module neigh
c
      if (allocated(xvold))  deallocate (xvold)
      if (allocated(yvold))  deallocate (yvold)
      if (allocated(zvold))  deallocate (zvold)
      if (allocated(xcold))  deallocate (xcold)
      if (allocated(ycold))  deallocate (ycold)
      if (allocated(zcold))  deallocate (zcold)
      if (allocated(xmold))  deallocate (xmold)
      if (allocated(ymold))  deallocate (ymold)
      if (allocated(zmold))  deallocate (zmold)
      if (allocated(nvlst))  deallocate (nvlst)
      if (allocated(vlst))  deallocate (vlst)
      if (allocated(nelst))  deallocate (nelst)
      if (allocated(elst))  deallocate (elst)
      if (allocated(nulst))  deallocate (nulst)
      if (allocated(ulst))  deallocate (ulst)

      if (allocated(xmold2))  deallocate (xmold2)
      if (allocated(ymold2))  deallocate (ymold2)
      if (allocated(zmold2))  deallocate (zmold2)
      if (allocated(nelst2))  deallocate (nelst2)
      if (allocated(elst2))  deallocate (elst2)


      if (allocated(elst_recv)) deallocate(elst_recv)
      if (allocated(nelst_recv)) deallocate(nelst_recv)
      if (allocated(elst_recv_single)) deallocate(elst_recv_single)

      if (allocated(vlst_recv)) deallocate(vlst_recv)
      if (allocated(nvlst_recv)) deallocate(nvlst_recv)
       
      if (allocated(mollst))  deallocate (mollst)
      if (allocated(nmollst))  deallocate (nmollst)
      if (allocated(mollst2))  deallocate (mollst2)
      if (allocated(nmollst2))  deallocate (nmollst2)
      if (allocated(xmolold))  deallocate (xmolold)
      if (allocated(ymolold))  deallocate (ymolold)
      if (allocated(zmolold))  deallocate (zmolold)

      if (allocated(njlst)) deallocate(njlst)
      if (allocated(jlst)) deallocate(jlst)
      if (allocated(starterecip)) deallocate(starterecip)
      if (allocated(enderecip)) deallocate(enderecip)
      if (allocated(jiter)) deallocate(jiter)

      if (allocated(starterecip2)) deallocate(starterecip2)
      if (allocated(enderecip2)) deallocate(enderecip2)
      if (allocated(jiter2)) deallocate(jiter2)
      if (allocated(doremainder)) deallocate(doremainder)

      if (allocated(jlst_recv)) deallocate(jlst_recv)
      if (allocated(jlst2_recv)) deallocate(jlst2_recv)

      if (allocated(start_polar)) deallocate(start_polar)
      if (allocated(last_polar))  deallocate(last_polar)
      if (allocated(molnew)) deallocate(molnew)

      if (allocated(start_emreal2)) deallocate(start_emreal2)
      if (allocated(last_emreal2))  deallocate(last_emreal2)
      if (allocated(maxsize_elst)) deallocate(maxsize_elst)
    
      if (allocated(start_vdw2)) deallocate(start_vdw2)
      if (allocated(last_vdw2)) deallocate(last_vdw2)
      if (allocated(maxsize_vlst)) deallocate(maxsize_vlst)
      
      
      if (allocated(start_polar3)) deallocate(start_polar3)
      if (allocated(last_polar3)) deallocate(last_polar3)
      if (allocated(maxsize_mollst3)) deallocate(maxsize_mollst3)
      if (allocated(xmolold3)) deallocate(xmolold3) 
      if (allocated(ymolold3)) deallocate(ymolold3)
      if (allocated(zmolold3)) deallocate(zmolold3)
      if (allocated(mollst3)) deallocate(mollst3)
      if (allocated(nmollst3)) deallocate(nmollst3)
      if (allocated(mollst3_recv)) deallocate(mollst3_recv)
      if (allocated(nmollst3_recv)) deallocate(nmollst3_recv)


c      if (allocated(offset_emreal2))  deallocate(offset_emreal2)

c
c     deallocation of global arrays from module nonpol
c
      if (allocated(rcav))  deallocate (rcav)
      if (allocated(rdisp))  deallocate (rdisp)
      if (allocated(cdisp))  deallocate (cdisp)
c
c     deallocation of global arrays from module opbend
c
      if (allocated(iopb))  deallocate (iopb)
      if (allocated(opbk))  deallocate (opbk)
c
c     deallocation of global arrays from module opdist
c
      if (allocated(iopd))  deallocate (iopd)
      if (allocated(opdk))  deallocate (opdk)
c
c     deallocation of global arrays from module orbits
c
      if (allocated(qorb))  deallocate (qorb)
      if (allocated(worb))  deallocate (worb)
      if (allocated(emorb))  deallocate (emorb)
c
c     deallocation of global arrays from module paths
c
      if (allocated(pc0))  deallocate (pc0)
      if (allocated(pc1))  deallocate (pc1)
      if (allocated(pvect))  deallocate (pvect)
      if (allocated(pstep))  deallocate (pstep)
      if (allocated(pzet))  deallocate (pzet)
      if (allocated(gc))  deallocate (gc)
c
c     deallocation of global arrays from module pdb
c
      if (allocated(resnum))  deallocate (resnum)
      if (allocated(resatm))  deallocate (resatm)
      if (allocated(npdb12))  deallocate (npdb12)
      if (allocated(ipdb12))  deallocate (ipdb12)
      if (allocated(pdblist))  deallocate (pdblist)
      if (allocated(xpdb))  deallocate (xpdb)
      if (allocated(ypdb))  deallocate (ypdb)
      if (allocated(zpdb))  deallocate (zpdb)
      if (allocated(pdbres))  deallocate (pdbres)
      if (allocated(pdbatm))  deallocate (pdbatm)
      if (allocated(pdbtyp))  deallocate (pdbtyp)
c
c     deallocation of global arrays from module piorbs
c
      if (allocated(iorbit))  deallocate (iorbit)
      if (allocated(iconj))  deallocate (iconj)
      if (allocated(kconj))  deallocate (kconj)
      if (allocated(piperp))  deallocate (piperp)
      if (allocated(ibpi))  deallocate (ibpi)
      if (allocated(itpi))  deallocate (itpi)
      if (allocated(pbpl))  deallocate (pbpl)
      if (allocated(pnpl))  deallocate (pnpl)
      if (allocated(listpi))  deallocate (listpi)
c
c     deallocation of global arrays from module pistuf
c
      if (allocated(bkpi))  deallocate (bkpi)
      if (allocated(blpi))  deallocate (blpi)
      if (allocated(kslope))  deallocate (kslope)
      if (allocated(lslope))  deallocate (lslope)
      if (allocated(torsp2))  deallocate (torsp2)
c
c     deallocation of global arrays from module pitors
c
      if (allocated(ipit))  deallocate (ipit)
      if (allocated(kpit))  deallocate (kpit)
c
c     deallocation of global arrays from module pme
c
      if (allocated(igrid))  deallocate (igrid)
      if (allocated(thetai1))  deallocate (thetai1)
      if (allocated(thetai2))  deallocate (thetai2)
      if (allocated(thetai3))  deallocate (thetai3)
      if (allocated(qgrid))  deallocate (qgrid)
c     if (allocated(qfac))  deallocate (qfac)
c
c     deallocation of global arrays from module polar
c
      if (allocated(polarity))  deallocate (polarity)
      if (allocated(thole))  deallocate (thole)
      if (allocated(pdamp))  deallocate (pdamp)
      if (allocated(uind))  deallocate (uind)
      if (allocated(uinp))  deallocate (uinp)
      if (allocated(uinds))  deallocate (uinds)
      if (allocated(uinps))  deallocate (uinps)
c
c     deallocation of global arrays from module polgrp
c
      if (allocated(ip11))  deallocate (ip11)
      if (allocated(ip12))  deallocate (ip12)
      if (allocated(ip13))  deallocate (ip13)
      if (allocated(ip14))  deallocate (ip14)
c
c     deallocation of global arrays from module qmstuf
c
      if (allocated(gx))  deallocate (gx)
      if (allocated(gy))  deallocate (gy)
      if (allocated(gz))  deallocate (gz)
      if (allocated(gfreq))  deallocate (gfreq)
      if (allocated(gforce))  deallocate (gforce)
      if (allocated(gh))  deallocate (gh)
c
c     deallocation of global arrays from module refer
c
      if (allocated(reftyp))  deallocate (reftyp)
      if (allocated(n12ref))  deallocate (n12ref)
      if (allocated(i12ref))  deallocate (i12ref)
      if (allocated(xref))  deallocate (xref)
      if (allocated(yref))  deallocate (yref)
      if (allocated(zref))  deallocate (zref)
      if (allocated(refnam))  deallocate (refnam)
c
c     deallocation of global arrays from module restrn
c
      if (allocated(ipfix))  deallocate (ipfix)
      if (allocated(kpfix))  deallocate (kpfix)
      if (allocated(idfix))  deallocate (idfix)
      if (allocated(iafix))  deallocate (iafix)
      if (allocated(itfix))  deallocate (itfix)
      if (allocated(igfix))  deallocate (igfix)
      if (allocated(ichir))  deallocate (ichir)
      if (allocated(xpfix))  deallocate (xpfix)
      if (allocated(ypfix))  deallocate (ypfix)
      if (allocated(zpfix))  deallocate (zpfix)
      if (allocated(pfix))  deallocate (pfix)
      if (allocated(dfix))  deallocate (dfix)
      if (allocated(afix))  deallocate (afix)
      if (allocated(tfix))  deallocate (tfix)
      if (allocated(gfix))  deallocate (gfix)
      if (allocated(chir))  deallocate (chir)
c
c     deallocation of global arrays from module rgddyn
c
      if (allocated(xcmo))  deallocate (xcmo)
      if (allocated(ycmo))  deallocate (ycmo)
      if (allocated(zcmo))  deallocate (zcmo)
      if (allocated(vcm))  deallocate (vcm)
      if (allocated(wcm))  deallocate (wcm)
      if (allocated(lm))  deallocate (lm)
      if (allocated(vc))  deallocate (vc)
      if (allocated(wc))  deallocate (wc)
      if (allocated(linear))  deallocate (linear)
c
c     deallocation of global arrays from module rigid
c
      if (allocated(xrb))  deallocate (xrb)
      if (allocated(yrb))  deallocate (yrb)
      if (allocated(zrb))  deallocate (zrb)
      if (allocated(rbc))  deallocate (rbc)
c
c     deallocation of global arrays from module ring
c
      if (allocated(iring3))  deallocate (iring3)
      if (allocated(iring4))  deallocate (iring4)
      if (allocated(iring5))  deallocate (iring5)
      if (allocated(iring6))  deallocate (iring6)
c
c     deallocation of global arrays from module solute
c
      if (allocated(rsolv))  deallocate (rsolv)
      if (allocated(asolv))  deallocate (asolv)
      if (allocated(rborn))  deallocate (rborn)
      if (allocated(drb))  deallocate (drb)
      if (allocated(drbp))  deallocate (drbp)
      if (allocated(drobc))  deallocate (drobc)
      if (allocated(gpol))  deallocate (gpol)
      if (allocated(shct))  deallocate (shct)
      if (allocated(aobc))  deallocate (aobc)
      if (allocated(bobc))  deallocate (bobc)
      if (allocated(gobc))  deallocate (gobc)
      if (allocated(vsolv))  deallocate (vsolv)
      if (allocated(wace))  deallocate (wace)
      if (allocated(s2ace))  deallocate (s2ace)
      if (allocated(uace))  deallocate (uace)
c
c     deallocation of global arrays from module stodyn
c
      if (allocated(fgamma))  deallocate (fgamma)
c
c     deallocation of global arrays from module strbnd
c
      if (allocated(isb))  deallocate (isb)
      if (allocated(sbk))  deallocate (sbk)
c
c     deallocation of global arrays from module strtor
c
      if (allocated(ist))  deallocate (ist)
      if (allocated(kst))  deallocate (kst)
c
c     deallocation of global arrays from module syntrn
c
      if (allocated(xmin1))  deallocate (xmin1)
      if (allocated(xmin2))  deallocate (xmin2)
      if (allocated(xm))  deallocate (xm)
c
c     deallocation of global arrays from module tarray
c
      if (allocated(tindex))  deallocate (tindex)
      if (allocated(tdipdip))  deallocate (tdipdip)
c
c     deallocation of global arrays from module tors
c
      if (allocated(itors))  deallocate (itors)
      if (allocated(tors1))  deallocate (tors1)
      if (allocated(tors2))  deallocate (tors2)
      if (allocated(tors3))  deallocate (tors3)
      if (allocated(tors4))  deallocate (tors4)
      if (allocated(tors5))  deallocate (tors5)
      if (allocated(tors6))  deallocate (tors6)
c
c     deallocation of global arrays from module tortor
c
      if (allocated(itt))  deallocate (itt)
c
c     deallocation of global arrays from module uprior
c
      if (allocated(udalt))  deallocate (udalt)
      if (allocated(upalt))  deallocate (upalt)
      if (allocated(usalt))  deallocate (usalt)
      if (allocated(upsalt))  deallocate (upsalt)

      if (allocated(udalt1b))  deallocate (udalt1b)
      if (allocated(upalt1b))  deallocate (upalt1b)
      if (allocated(nualt1b)) deallocate (nualt1b)
      if (allocated(bpred1b)) deallocate (bpred1b)
      if (allocated(bpredp1b)) deallocate (bpredp1b)

      if (allocated(udalt2b))  deallocate (udalt2b)
      if (allocated(upalt2b))  deallocate (upalt2b)
      if (allocated(nualt2b)) deallocate (nualt2b)
      if (allocated(bpred2b)) deallocate (bpred2b)
      if (allocated(bpredp2b)) deallocate (bpredp2b)

c
c     deallocation of global arrays from module urey
c
      if (allocated(iury))  deallocate (iury)
      if (allocated(uk))  deallocate (uk)
      if (allocated(ul))  deallocate (ul)
c
c     deallocation of global arrays from module usage
c
      if (allocated(iuse))  deallocate (iuse)
      if (allocated(use))  deallocate (use)
c
c     deallocation of global arrays from module usolve
c
      if (allocated(mindex))  deallocate (mindex)
      if (allocated(minv))  deallocate (minv)
c
c     deallocation of global arrays from module vdw
c
      if (allocated(ivdw))  deallocate (ivdw)
      if (allocated(jvdw))  deallocate (jvdw)
      if (allocated(ired))  deallocate (ired)
      if (allocated(kred))  deallocate (kred)
      if (allocated(radmin))  deallocate (radmin)
      if (allocated(epsilon))  deallocate (epsilon)
      if (allocated(radmin4))  deallocate (radmin4)
      if (allocated(epsilon4))  deallocate (epsilon4)
      if (allocated(radhbnd))  deallocate (radhbnd)
      if (allocated(epshbnd))  deallocate (epshbnd)
c
c     deallocation of global arrays from module vibs
c
      if (allocated(phi))  deallocate (phi)
      if (allocated(phik))  deallocate (phik)
      if (allocated(pwork))  deallocate (pwork)
c
c     deallocation of global arrays from module warp
c
      if (allocated(m2))  deallocate (m2)
c
c     free memory used by the APBS Poisson-Boltzmann solver
c
      if (solvtyp .eq. 'PB') then
         call apbsfinal
      end if
c
c     close any open socket used for external communication
c
      if (use_socket) then
         call sktkill
      end if
c
c     print a final status message before exiting TINKER
c
      if (debug) then
         write (iout,10)
   10    format (/,' TINKER is Exiting following Normal Termination',
     &              ' of the Program',/)
      end if
c
c     may need a pause to avoid closing the execution window
c
      if (holdup) then
         read (input,20)
   20    format ()
      end if
      return
      end
