##
###################################################################
##                                                               ##
##  Makefile for Building the Tinker Molecular Modeling Package  ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all the TINKER executables
##   2. make rename           Move the executables to BINDIR
##   3. make create_links     Create soft links in LINKDIR
##   4. make remove_links     Remove soft links from LINKDIR
##   6. make listing          Concatenate source to tinker.txt
##   5. make clean            Delete objects and executables
##
##  Original version of this file is due to Peter Happersberger
##  and Jochen Buehler of the University of Konstanz, Jan. 1998,
##  Modifications by Reece Hart & Jay Ponder, Washington University
##
##  Note: Building TINKER requires libraries for the FFTW Fast
##  Fourier Transform package and (optionally) the APBS Poisson-
##  Boltzmann solver. These libraries are assumed to already be
##  present, and are not built by this Makefile.
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##  TINKERDIR    TINKER Distribution Directory
##  BINDIR       Hard Copies of TINKER Executables
##  LINKDIR      Linked Copies of TINKER Executables

TINKERDIR = $(HOME)/3-amoeba
BINDIR = $(TINKERDIR)/bin
LINKDIR = /usr/local/bin

####################################################################
##  Known Machine Types; Uncomment One of the Following Sections  ##
##  May Need Editing to Match Your Desired OS & Compiler Version  ##
####################################################################

##  Machine:  Generic Linux or Macintosh
##  CPU Type: Intel x86 Compatible
##  Compiler: GNU gfortran

#F77 = gfortran
#LIBS =
#F77FLAGS = -c
#OPTFLAGS = -O3 -ffast-math
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Generic Linux or Macintosh
##  CPU Type: Intel x86 Compatible
##  Compiler: GNU gfortran
##  Parallel: OpenMP

#F77 = gfortran
#LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c
#OPTFLAGS = -O3 -ffast-math -fopenmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 14.0

#F77 = ifort
#LIBS =
#F77FLAGS = -c -xHost
#OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = -static-libgcc -static-intel

##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 14.0
##  Parallel: OpenMP

#F77 = ifort
#LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c -xHost
#OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive -openmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc -static-intel

# COMPILING ON HOPPER
F77 = ftn
LIBS = -L/opt/cray/fftw/3.3.4.1/x86_64/lib -lfftw3_threads -lfftw3
F77FLAGS = -c
OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive -openmp -dynamic \
           -mcmodel=large -check bounds -traceback
LIBFLAGS = -crusv
RANLIB = ranlib -c
LINKFLAGS = $(OPTFLAGS)


# COMPILING ON EDISON
F77 = ftn
LIBS = -L/opt/fftw/3.3.0.4/sandybridge/lib -lfftw3_threads -lfftw3
F77FLAGS = -c
OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive -openmp -dynamic \
           -mcmodel=large #-check bounds -traceback #-check all -traceback
LIBFLAGS = -crusv
RANLIB = ranlib -c
LINKFLAGS = $(OPTFLAGS)

# Find out the system hostname.

OLDALBAUGHDIR = /home/aalbaugh/tinker/tinker6.3.03/fftw/lib
NEWALBAUGHDIR = /home/aalbaugh/tinker/tinker7.1.2/fftw/lib
# New on Armada OpenMPI w/ OpenMP !!Gfortran!!
#F77 = mpif77
#LIBS = -L$(NEWALBAUGHDIR) -lfftw3_threads -lfftw3
#F77FLAGS = -c -fno-align-commons -fstrict-aliasing
#OPTFLAGS = -O3 -ffast-math -march=native -mno-avx -fopenmp \
           -mcmodel=large #-fbounds-check -fbacktrace # -fcheck=all
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc


HOSTNAME:=$(shell hostname)

ifeq ($(HOSTNAME),ultra)

  # Compiling on ultra
  F77 = mpif77
  LIBS = -L$(TINKERDIR)/thirdparty/fftw/lib -lfftw3_threads -lfftw3
  F77FLAGS = -c -xHost
  OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive -openmp
  LIBFLAGS = -crusv
  RANLIB = ranlib #-c
  LINKFLAGS = $(OPTFLAGS) -static-libgcc -static-intel

else ifeq ($(HOSTNAME),indy0)

  # Compiling options for indy
  F77 = mpif77
  LIBS = -L$(TINKERDIR)/thirdparty/fftw/lib -lfftw3_threads -lfftw3
  F77FLAGS = -c
  OPTFLAGS = -O3 -ffast-math -fopenmp
  LIBFLAGS = -crusv
  RANLIB = ranlib 
  LINKFLAGS = $(OPTFLAGS) -static-libgcc

else ifeq ($(HOSTNAME),phi.hydra)

  # Compiling options for indy
  F77 = mpiifort
  LIBS = -L$(TINKERDIR)/thirdparty/fftw/lib -lfftw3_threads -lfftw3
  F77FLAGS = -c
  OPTFLAGS = -O3 -g -fopenmp
  LIBFLAGS = -crusv
  RANLIB = ranlib 
  LINKFLAGS = $(OPTFLAGS) -static-intel

endif


# New on Armada OpenMPI w/ OpenMP !!Gfortran!! Non-static linking
#F77 = mpif77
#LIBS = -L/home/aalbaugh/tinker/tinker6.3.03/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c -fno-align-commons -fstrict-aliasing
#OPTFLAGS = -O3 -ffast-math -march=native -mno-avx -fopenmp \
           -mcmodel=large 
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS)


# New on Armada OpenMPI w/ OpenMP !!Gfortran!! Vanilla Compilation Flags!!
#F77 = mpif77
#LIBS = -L/home/aalbaugh/tinker/tinker6.3.03/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c 
#OPTFLAGS = -O3 -fopenmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.9 (Mavericks)
##  Compiler: Intel Fortran for Mac 14.0

#F77 = ifort
#LIBS =
#F77FLAGS = -c -axSSSE3
#OPTFLAGS = -O3 -no-ipo -no-prec-div -mdynamic-no-pic
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.6 \
             -Wl,-stack_size,0x10000000

##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.9 (Mavericks)
##  Compiler: Intel Fortran for Mac 14.0
##  Parallel: OpenMP

#F77 = ifort
#LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c -axSSSE3
#OPTFLAGS = -O3 -no-ipo -no-prec-div -mdynamic-no-pic -openmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.6 \
            -Wl,-stack_size,0x10000000

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = \
       dynamic_iAMOEBA_tensor.o \
       dynamic_3mAMOEBA.o \
       dynamic_3M_AMOEBA.o \
       dynamic_2M_AMOEBA.o \
       mdsave_writeDyns.o \
       prtdyn_writeDyns.o \
       uind2bmat.o \
       init2bmatmodPolar_uind.o \
       uind1bmatrix.o \
       init1bmat_uind.o \
       uindclust_gradient_Polar1b2b3bsimult_ireducemat_commsmall2.o \
       uindempole1d_3b_Polar1btest2.o \
       uindClustEmpole1plz_1b_NoOmplist.o \
       uindClustEmpole1plz_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o \
       uindClustEmpole1plz_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o \
       uindempole1d_3b_Polar2btest2.o \
       uindempole1d_3b_Polar3btest2.o \
       pol_loadbal_clust_splitlist_offset_2.o \
       pol_loadbal_clust_splitlist_offset3b_2.o \
       ClustEmpole1plz_2b_NoOmplistRtrn2BNoSubtr1bsimult.o \
       clust_gradient_Polar1b2bsimult_ireducemat_commsmall2.o \
       clust_gradient_Polar1b2b3bsimult_ireducemat_commsmall2.o \
       empole1d_3b_Polar1btest2.o \
       empole1d_3b_Polar2btest2.o \
       empole1d_3b_Polar3btest2.o \
       ClustEmpole1plz_1b_NoOmplist.o \
       ClustEmpole1plz_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o \
       ClustEmpole1plz_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o \
       clust_gradient_Polar1b2b3bsimult_reducemat_commsmall2.o \
       grad_covemrecip_vandr_reduce_totfield_commsmall.o \
       clust_gradient_PolarPerm1b2b3bsimult_reducemat_commsmall2.o \
       empole1d_3b_Polar3btest.o \
       ulspred3b.o \
       uprior3b.o \
       initpolarpred3b_offset.o \
       clust3blists_offset.o \
       mlist_3bclust.o \
       ulist_3bclust.o \
       pol_loadbal_clust_splitlist_offset3b_simple.o \
       init2bmatmodPolar.o \
       ClustEmpole1_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o \
       clust_gradient_PolarPerm1b2b3bsimult_reducemat_commsmall.o \
       ClustEmpole1_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o \
       empole1d_3b_PolarPerm3btest.o \
       induce0a_3b_scf_vacomp3bclustlist.o \
       dfield0c_3bPolar_3bclust.o \
       grad_cov_v_reduce_commsmall.o \
       mpidat4.o \
       clust_gradient_PolarPerm1b2bsimult_reducemat_commsmall.o \
       ClustEmpole1_2b_NoOmplistRtrn2BNoSubtr1bsimult.o \
       clust_gradient_covvdwPolarPerm1b2b_ireducemat_commsmall.o \
       ClustEmpole1_2b_NoOmplistRtrn2BNoSubtr1b.o \
       clust_gradient_covvdwPolarPerm1b2b_ireducemat.o \
       nblist_par_simultclust.o \
       loadbal_ehal1c_halftasks_nolistsend_simultclust.o \
       clust_gradient_covvdwPolarPerm1b2b.o \
       cell1bclust.o \
       boxes1bclust.o \
       cell2bclust.o \
       boxes2bclust.o \
       ulspred1b.o \
       initpolarpred1b.o \
       uprior1b.o \
       uprior2b.o \
       initpolarpred2b_offset.o \
       initpolarpred2b.o \
       ulspred2b.o \
       clust_gradient_PolarPerm1b2bsimultmod.o \
       ClustEmpole1_1b_NoOmplistsimult.o \
       ClustEmpole1_2b_NoOmplistsimult.o \
       clust_gradient_PolarPerm1b2bsimult.o \
       init1b2bmatPolarPerm.o \
       deriv2bmat.o \
       pol_loadbal_clust_splitlist_offset.o \
       clust2blists_offset.o \
       clust_gradient_PolarPerm1b2ball.o \
       ClustEmpole1_2b_NoOmpAll2b.o \
       empole1c_3b_PolarPermtest.o \
       empole1d_3b_PolarPerm2btest.o \
       empole1d_3b_PolarPerm1btest.o \
       empole1c_3b_PolarPerm.o \
       ClustEmpole1_2b_NoOmp.o \
       ClustEmpole1_1b_NoOmp.o \
       grad_cov_v_reduce.o \
       init1bmatPolarPerm.o \
       torque3_PolarPerm.o \
       clust_gradient_PolarPerm1b2b.o \
       ClustEmpole1_1b_NoOmplist.o \
       ClustEmpole1_2b_NoOmplist.o \
       empole1d_3b_PolarPerm1b.o \
       empole1d_3b_PolarPerm2b.o \
       dfield0c_3bPolar_2bclust.o \
       dfield0c_3bPolar_1bclust.o \
       ufield0c_3bPolar.o \
       empole1d_3bPolar_origversion_nopriordirectcalc2_omp2b.o \
       empole1d_3bPolar_origversion_nopriordirectcalc2_omp1b.o \
       ereal1d_3b_1bclust.o \
       ereal1d_3b_2bclust.o \
       emrecip1_3b_Polar_new.o \
       pmestuff_3bnew.o \
       mlist_2bclust.o \
       mlist_1bclust.o \
       ulist_2bclust.o \
       induce0a_3b_scf_vacomp1bclustlist.o \
       induce0a_3b_scf_vacomp2bclustlist.o \
       empole1b_3bPolar_origversion_nopriordirectcalc2_omp2b.o \
       empole1b_3bPolar_origversion_nopriordirectcalc2_omp1b.o \
       ulist_1bclust.o \
       clust2blists.o \
       clust1blists.o \
       clustloadbalsmoothInner_1a_2bPolar_vac_precalc1bNoOmplist.o \
       clustnolistsmoothInner_1a_1bPolar_vac_rtrn1bNoOmplist.o \
       clustnolistsmoothInner_1a_1bPolar_vac_rtrn1bNoOmp.o \
       clustloadbalsmoothInner_1a_2bPolar_vac_precalc1bNoOmp.o \
       induce0a_3b_scf_vacomp.o \
       empole1a_3bPolar_origversion_nopriordirectcalc2_omp.o \
       pol_loadbal_clust_splitlist.o \
       pol_loadbal_clust.o \
       clustnolistsmoothInner_1a_1bPolar_vac_rtrn1b.o \
       clustloadbalsmoothInner_1a_2bPolar_vac_precalc1b.o \
       clust_smooth_gradpolar1b2b.o \
       init1bmat.o \
       deriv1bmat.o \
       nblist5kmeansclust.o \
       clustloadbalsmoothInner_1a_2bPolar_vac.o \
       clustnolistsmoothInner_1a_1bPolar_vac.o \
       pcg.o \
       basefileclust.o \
       getclust.o \
       readclust.o \
       induce0a_3b_scf_vac.o \
       induce0a_3b_cholesky.o \
       cho.o \
       pmestuffNoQPoleNoDPole.o \
       emrecip1_3b_Perm2NoQPoleNoDPole.o \
       fmat.o \
       totfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2fij.o \
       empole1a_3bPolar_totfieldnpolefij.o \
       empole1a1_3b_Polar_totfieldnpole_modfij.o \
       empole1a2_3b_Polar_totfieldnpole_modfij.o \
       torque_3b_new_npolefij.o \
       totfieldNpoleloadbalsmoothInner_1a_3bPolar_w1bodyfij.o \
       empole1a_3bPolar_origversion_nopriordirectcalc2_fijscale2.o \
       empole1a_3bPolar_origversion_nopriordirectcalc2_fijscale.o \
       Newempole1c_3bPolar_totfield_dEtensor_fewertens1.o \
       totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi.o \
       totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi1.o \
       totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi1.o \
       empole1a_3bPolar_origversion_nopriordirectcalc2_fij.o \
       torque_3b_new_fij.o \
       All2bAll3bNew2smooth_innerloop_empole1a_3b_fij.o \
       4bodyAll2bAll3bNew2smooth_innerloop_empole1a_3b.o \
       empole1a_mod.o \
       sulf.o \
       sAll2bNew2smooth_innerloop_empole1a_3b.o \
       sAll2bAll3bNew2smooth_innerloop_empole1a_3b.o \
       sAll2bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
       sAll2bAll3bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
       All2bAll3bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
       All2bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
       All2bNew2smooth_innerloop_empole1a_3b.o \
       sAll2bNew2smooth_innerloop_empole1a_3b.o \
       s2btotfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2.o \
       s2btotfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o \
       s2btotfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o \
       stotfieldloadbalsmooth_findiff_test.o \
       ewpol_loadbal7_clust.o \
       loadbalsmoothInner_1a_2bPolar_vac.o \
       empole1a2_3b_Polar_totfieldnpole3bmod_dol1_2.o \
       totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi1.o \
       totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi1.o \
       Newempole1c_3bPolar_totfield_dEtensor_fewertens1.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens1.o \
       thg5_small_vir_nolistsend_Nn1_noewtens_notenssends_fewertens.o \
       torque3_dEtensor_72015_taunew4_pt1_vir_notenssends_fewertens.o \
       torque3_dEtensor_72015_taunew4_pt2_vir_notenssends_fewertens.o \
       empole1a1_3b_Polar_totfieldnpole3b.o \
       empole1a_3bPolar_totfieldnpole3b_rtrndeptemp2.o \
       thg5_small_vir_nolistsend_Nn_noewtens_fewertens.o \
       torque3_dEtensor_72015_taunew3_pt1_vir_fewertens.o \
       torque3_dEtensor_72015_taunew3_pt2_vir_fewertens.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_fewertens.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens2.o \
       totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi.o \
       tensdealloc.o \
       Newempole1c_3bPolar_totfield_dEtensor_fewertens.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens.o \
       thg5_small_vir_nolistsend_Nn_noewtens_notenssends_fewertens.o \
       torque3_dEtensor_72015_taunew3_pt1_vir_notenssends_fewertens.o \
       torque3_dEtensor_72015_taunew3_pt2_vir_notenssends_fewertens.o \
       dEtensor3.o \
       thg5_small_vir_nolistsend_noewtens2.o \
       thg5_small_vir_nolistsend_Nn_notenssends.o \
       torque3_dEtensor_72015_taunew3_pt1_vir_notenssends.o \
       torque3_dEtensor_72015_taunew3_pt2_vir_notenssends.o \
       thg5_small_vir_nolistsend_Nn_noewtens_notenssends.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends.o \
       mpidat3.o \
       dEtensor2.o \
       torque3_dEtensor_72015_taunew3_pt1_vir.o \
       torque3_dEtensor_72015_taunew3_pt2_vir.o \
       nblist_short.o \
       limits2.o \
       neigh2.o \
       switch2.o \
       sizes2.o \
       Newempole1c_3bPolar_totfield_dEtensor.o \
       thg5_small_vir_nolistsend_Nn_noewtens.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3.o \
       cutoffs_ewald_short.o \
       ewpol_loadbal7.o \
       nblist2_hbond.o \
       nblist2_test.o \
       nblist4_com.o \
       nblist3_com.o \
       clustnolistsmoothInner_1c_3bPolar_vac.o \
       nblist4.o \
       nblist3.o \
       mpidat2.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2_no3bmut.o \
       moments_onlyDipole_fixedcharge.o \
       uindInner_1a_3bPolar_vac.o \
       pmonte_par2.o \
       empole0a_3bPolar_origversion_nopriordirectcalc2.o \
       eng_r_reduce.o \
       loadbalsmoothInner_0a_3bPolar_vac.o \
       smooth_engpolar.o \
       empole_Perm.o \
       eng_covalent_vdw.o \
       verlet_par_pmonte.o \
       pmonte_par.o \
       uindempole1a_3bPolar_origversion_nopriordirectcalc2.o \
       uindloadbalsmoothInner_1a_3bPolar_vac.o \
       uindsmooth_gradpolar.o \
       clust_smooth_gradpolar.o \
       clustnolistsmoothInner_1a_3bPolar_vac.o \
       clustloadbalsmoothInner_1a_3bPolar_vac.o \
       nblist2.o \
       nblist5.o \
       neigh2clust.o \
       pol_loadbal7.o \
       pol_loadbal6.o \
       thg5_small_vir_nolistsend_noewtens.o \
       empole1c_3bPolar_totfield_dEtensor.o \
       virialnamoeba.o \
       dfield0a_ufield0a_uscale0a.o \
       empole1a_PolarStripped.o \
       induce0a_nobranch.o \
       torque_pol.o \
       ewpol_loadbal_halftasks.o \
       dEtensor_test_nolistsend_separatereduce_funcdecomp_simultpolz.o \
       loadbalsmoothInner_1a_3bPolar1bmut_halftasks.o \
       loadbalsmoothInner_1a_3bPolar_halftasks.o \
       loadbalsmoothInner_1a_3bPolar_vac.o \
       loadbalsmoothInner_1a_3bPolar1bmut.o \
       dfield0b_totfield.o \
       LoadBalInner_r_Perm_nolistsend.o \
       grad_covemrecip_vandr_reduce_totfield.o \
       empole1a_3bPolar_origversion_nopriordirectcalc_rtrnfield.o \
       induce0a_3b_PolelecOnly_totfieldnpolevac.o \
       empole1a2_3b_Polar_totfieldnpolevac.o \
       empole1a_3bPolar_new_subtract_1bodymut.o \
       empole1a_3bPolar_totfieldnpolevac.o \
       empole1a_3bPolar_origversion_nopriordirectcalc2.o \
       emrecip1_1_3b_Polar_totfield_parl_513.o \
       empole1c_3bPolar_totfield_901.o \
       emrecip1_3b_Polar_usemlist_mutual_virialonly.o \
       emrecip1_3b_Polar_usemlist_mutual.o \
       induce0c_field_amat.o \
       ereal1c1_3b_totfieldnpole.o \
       totfieldnpole3b_noewald_umutual_rl_3b_nofield.o \
       empole1a_3bPolar_totfieldnpole.o \
       empole1a1_3b_Polar_totfieldnpole_mod.o \
       empole1a2_3b_Polar_totfieldnpole_mod.o \
       induce0a_3b_PolelecOnly_totfieldnpole.o \
       empole1a1_3bPolar_totfield_dEtensor1bmut.o \
       emrecip1_3b_Polar_usemlist_direct_virialonly.o \
       b2test.o \
       emrecip1_3b_Polar_usemlist_direct_test2.o \
       emrecip1_3b_Polar_usemlist_direct_test.o \
       empole1d_3b_PolarPerm_selfeng_bcast2.o \
       grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2.o \
       bcast_cov.o \
       grad_r_allreduce_totfield_dEtensor_vir.o \
       thg6_small_vir_nolistsend.o \
       grad_vandr_ireduce_totfield_dEtensor_vir.o \
       grad_vandr_ireduce_totfield_dEtensor_vir2.o \
       loadbalsmoothInner_1a_2bPolar.o \
       nblist3b_1list_par_newesttink7.o \
       chkpole_omp.o \
       rotpole_omp.o \
       grad_vandr_reduce_totfield_dEtensor_vir2.o \
       nblist_par.o \
       loadbal_ehal1c_halftasks_nolistsend.o \
       ereal1d_Perm_savetensor.o \
       kewald_eastman_list.o \
       chunklight.o \
       pmestuff_list.o \
       v_loadbal_halftasks.o \
       loadbal_ehal1c_halftasks.o \
       grad_vandr_reduce_totfield_dEtensor_vir.o \
       empole1c_3bPolar_totfield_dEtensor1b_nouselist_mpivir.o \
       torque3_dEtensor_72015_taunew2_pt1_vir.o \
       torque3_dEtensor_72015_taunew2_pt2_vir.o \
       thg5_vir.o \
       thg5_small_vir.o \
       grad_r_reduce_totfield_dEtensor_vir.o \
       thg6_save.o \
       thg6_save_small.o \
       torque2_omp.o \
       table_fill_setchunk_omp.o \
       bspline_fill_bsplgen_omp.o \
       LoadBalInner_thg4omp_mpi_72015torque2.o \
       torque3_dEtensor_72015_taunew2_pt1.o \
       torque3_dEtensor_72015_taunew2_pt2.o \
       thg5.o \
       thg5_small.o \
       emrecip1_3b_Perm2_totfield_mpi_gridmpole.o \
       emrecip1_3b_Perm2_totfield_mpi_fphimpole.o \
       dEtensorewtotfield_dyn_rMPIOMP_kMPIgridmpoleBcast_rloadbal_vMPIOMP.o \
       dEtensorewtotfield_dyn_rMPIOMP_kMPIfphimpoleBcast_rloadbal_vMPIOMP.o \
       nose_tuckerman06_par.o \
       tuckerman06.o \
       grid_mpole_Mpi3amoebaRepo_allreduce.o \
       emrecip1_3b_Perm2_totfield_mpi2.o \
       grid_mpole_Mpi3amoebaRepo.o \
       gradient_emrecip_totfield_mpi.o \
       emrecip1_3b_Perm2_totfield_mpi.o \
       fphi_mpole_Mpi3amoebaRepo.o \
       kewald_nopmealloc_eastman.o \
       mech_pme4_eastman.o \
       kewald_eastman.o \
       allocvdw_master2.o \
       bcast_vdw.o \
       v_loadbal.o \
       parvdwneigh.o \
       loadbal_ehal1c.o \
       gradient_covalent2.o \
       grad_v_reduce.o \
       field_r_reduce_uindcalc.o \
       udirect2b_uinddir_mpi.o \
       LoadBalInnerloop_ereal1d_sentlist2.o \
       udirect2b_uinddir.o \
       empole1c_3bPolar_totfield_dEtensor1b_nouselist_mpi.o \
       empole1c_3bPolar_totfield_dEtensor1b_nouselist.o \
       totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi.o \
       totfield_smooth_gradpolar_dEtensor.o \
       LoadBalInner_thg4omp_mpi.o \
       grad_r_reduce_totfield_dEtensor.o \
       grad_r_serial_totfield_dEtensor.o \
       ewtotfieldsmooth2_gradient_polar_serial_dEtensor.o \
       totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensor.o \
       dEtensor.o \
       LoadBalInner_r_Perm_totfield_bcastlist.o \
       grad_r_reduce_totfield_bcast.o \
       empole1c_3bPolar_totfield_usemlist.o \
       emrecip1_3b_Polar_usemlist_direct.o \
       empole1d_3b_PolarPerm_selfeng_bcast.o \
       grad_r_reduce_totfield_uindcalc.o \
       ereal1d2_3b_totfieldnpole_mlist.o \
       ereal1c2_3b_totfieldnpole_omp.o \
       totfieldNpolesmoothInner_1c_1bPolar2_513_alliter.o \
       totfieldNpolesmoothInner_1c_1bPolar2_513_multiter.o \
       totfieldNpolesmoothInner_1c_1bPolar2_513.o \
       ewtotfieldsmallsmooth2_gradient_polar_serial.o \
       grad_r_serial_totfield.o \
       sendmlist_load_balanced3_sizelist_serial.o \
       pmestuff_3bpolz_513.o \
       induce_totfield.o \
       pmestuff_3bpolz_totfield_513.o \
       pmestuff_3bpolz_totfield_pt2_513.o \
       emrecip1_1and2_3b_Polar_totfield_parl_513.o \
       kewald3b_totfield.o \
       bspline_table_recipsum.o \
       totfieldNpoleloadbalsmoothInner_1c_2bPolar2.o \
       aprx.o \
       totfieldNpoleloadbalsmoothInner_1c_2bPolar.o \
       nblist3b_1list.o \
       alloc3b_vdwlist_1list.o \
       totfieldNpolesmoothInner_1c_3bPolar.o \
       totfieldNpoleloadbalsmoothInner_1c_3bPolar.o \
       ereal1d2_3b_totfieldnpole.o \
       r_loadbal2.o \
       gradient_emrecip_totfield.o \
       gradient_emrecip_totfield_send.o \
       emrecip1_3b_Perm2_totfield.o \
       LoadBalInner_r_Perm_totfield.o \
       grad_r_reduce_totfield.o \
       pmestuff_3bpolz_totfield_pt2.o \
       pmestuff_3bpolz_totfield.o \
       dfield0c_totfield.o \
       dfield0a_totfield.o \
       ereal1c1_3b_totfieldnpole.o \
       ereal1c2_3b_totfieldnpole.o \
       empole1c_3bPolar_totfield.o \
       All2bAll3bNew2smooth_innerloop_empole1c_3b_totfield.o \
       All2bAll3bNew2smooth_innerloop_empole1c_3b.o \
       induce.o \
       empole1c_3bPolar.o \
       pmestuff_3bpolz_pt2.o \
       torque3_Polar.o \
       torque_3b_new_npole.o \
       totfield.o \
       totfield_smooth_gradpolar.o \
       fft3d_3b.o \
       kewald3b.o \
       pmestuff_3bpolz.o \
       emrecip1_3b_Polar.o \
       All2bAll3bNew2smooth_innerloop_empole1a_3b.o \
       ewreg3bpolz.o \
       kewaldreg3b.o \
       ewpol_loadbal.o \
       action.o \
       mdstat_propertiesDaveDeTo.o \
       alloc3b_uind.o \
       moments_onlyDipole.o \
       get_ntask_r.o \
       r_loadbal.o \
       emreal_load_balance.o \
       sendmlist_load_balanced2.o \
       sendmlist_load_balanced.o \
       LoadBalInner_r_Perm.o \
       loadbalsmoothInner_1a_3bPolar.o \
       LoadBalcom_smooth2b_innerloop_empole1a_3b.o \
       pol_loadbal.o \
       empole1a_Perm.o \
       torque_Perm.o \
       nolistcom_smooth2b_innerloop_empole1a_3b.o \
       ptest_parallel.o \
       emrecip1PermReduce.o \
       mech_pme4.o \
       kewald_nopmealloc.o \
       grad_k_reduce.o \
       paremneigh.o \
       sendmlist.o \
       Inner_r_Perm.o \
       ewaldneigh.o \
       parewreg.o \
       kewaldreg.o \
       ewreg.o \
       smooth2b_gradient_polar.o \
       smooth_gradpolar.o \
       grad_r_reduce.o \
       mech_pme3.o \
       PermElec2.o \
       Modgrad_3b.o \
       neigh3b.o \
       mechanic_pme.o \
       mech_par.o \
       deriv3b.o \
       PermElec.o \
       nblist3b.o \
       par_nblist_2b.o \
       mpidat.o \
       empole1a_3bPolar.o \
       gradient_3b.o \
       nose_par.o \
       verlet_par.o \
       beeman_pt1_pt2_master.o \
       mdinit_par.o \
       gradslow_serial.o \
       cobar.o \
       active.o \
       alchemy.o \
       align.o \
       analysis.o \
       analyz.o \
       analyze.o \
       angang.o \
       angbnd.o \
       angles.o \
       angpot.o \
       angtor.o \
       anneal.o \
       archive.o \
       argue.o \
       ascii.o \
       atmlst.o \
       atomid.o \
       atoms.o \
       attach.o \
       bar.o \
       basefile.o \
       bath.o \
       beeman.o \
       bicubic.o \
       bitor.o \
       bitors.o \
       bndpot.o \
       bndstr.o \
       bonds.o \
       born.o \
       bound.o \
       bounds.o \
       boxes.o \
       bussi.o \
       calendar.o \
       cell.o \
       center.o \
       charge.o \
       chgpot.o \
       chkpole.o \
       chkring.o \
       chkxyz.o \
       cholesky.o \
       chrono.o \
       chunks.o \
       clock.o \
       cluster.o \
       column.o \
       command.o \
       connect.o \
       connolly.o \
       control.o \
       correlate.o \
       couple.o \
       crystal.o \
       cspline.o \
       cutoffs.o \
       cutoffs_vdw_nolist.o \
       deflate.o \
       delete.o \
       deriv.o \
       diagq.o \
       diffeq.o \
       diffuse.o \
       dipole.o \
       disgeo.o \
       distgeom.o \
       dma.o \
       document.o \
       domega.o \
       dynamic.o \
       eangang.o \
       eangang1.o \
       eangang2.o \
       eangang3.o \
       eangle.o \
       eangle1.o \
       eangle2.o \
       eangle3.o \
       eangtor.o \
       eangtor1.o \
       eangtor2.o \
       eangtor3.o \
       ebond.o \
       ebond1.o \
       ebond2.o \
       ebond3.o \
       ebuck.o \
       ebuck1.o \
       ebuck2.o \
       ebuck3.o \
       echarge.o \
       echarge1.o \
       echarge2.o \
       echarge3.o \
       echgdpl.o \
       echgdpl1.o \
       echgdpl2.o \
       echgdpl3.o \
       edipole.o \
       edipole1.o \
       edipole2.o \
       edipole3.o \
       egauss.o \
       egauss1.o \
       egauss2.o \
       egauss3.o \
       egeom.o \
       egeom1.o \
       egeom2.o \
       egeom3.o \
       ehal.o \
       ehal1.o \
       ehal2.o \
       ehal3.o \
       eimprop.o \
       eimprop1.o \
       eimprop2.o \
       eimprop3.o \
       eimptor.o \
       eimptor1.o \
       eimptor2.o \
       eimptor3.o \
       elj.o \
       elj1.o \
       elj2.o \
       elj3.o \
       embed.o \
       emetal.o \
       emetal1.o \
       emetal2.o \
       emetal3.o \
       emm3hb.o \
       emm3hb1.o \
       emm3hb2.o \
       emm3hb3.o \
       empole.o \
       empole1.o \
       empole2.o \
       empole3.o \
       energi.o \
       energy.o \
       eopbend.o \
       eopbend1.o \
       eopbend2.o \
       eopbend3.o \
       eopdist.o \
       eopdist1.o \
       eopdist2.o \
       eopdist3.o \
       epitors.o \
       epitors1.o \
       epitors2.o \
       epitors3.o \
       erf.o \
       erxnfld.o \
       erxnfld1.o \
       erxnfld2.o \
       erxnfld3.o \
       esolv.o \
       esolv1.o \
       esolv2.o \
       esolv3.o \
       estrbnd.o \
       estrbnd1.o \
       estrbnd2.o \
       estrbnd3.o \
       estrtor.o \
       estrtor1.o \
       estrtor2.o \
       estrtor3.o \
       etors.o \
       etors1.o \
       etors2.o \
       etors3.o \
       etortor.o \
       etortor1.o \
       etortor2.o \
       etortor3.o \
       eurey.o \
       eurey1.o \
       eurey2.o \
       eurey3.o \
       evcorr.o \
       ewald.o \
       extra.o \
       extra1.o \
       extra2.o \
       extra3.o \
       faces.o \
       fatal.o \
       fft.o \
       fft3d.o \
       fftpack.o \
       field.o \
       fields.o \
       files.o \
       final.o \
       flatten.o \
       fracs.o \
       freeunit.o \
       freeze.o \
       gda.o \
       geometry.o \
       getint.o \
       getkey.o \
       getmol.o \
       getmol2.o \
       getnumb.o \
       getpdb.o \
       getprm.o \
       getref.o \
       getstring.o \
       gettext.o \
       getword.o \
       getxyz.o \
       ghmcstep.o \
       gkstuf.o \
       gradient.o \
       gradrgd.o \
       gradrot.o \
       group.o \
       groups.o \
       grpline.o \
       gyrate.o \
       hescut.o \
       hessian.o \
       hessn.o \
       hessrgd.o \
       hessrot.o \
       hpmf.o \
       hybrid.o \
       image.o \
       impose.o \
       improp.o \
       imptor.o \
       inertia.o \
       inform.o \
       initatom.o \
       initial.o \
       initprm.o \
       initres.o \
       initrot.o \
       insert.o \
       intedit.o \
       inter.o \
       intxyz.o \
       invbeta.o \
       invert.o \
       iounit.o \
       jacobi.o \
       kanang.o \
       kangang.o \
       kangle.o \
       kangs.o \
       kangtor.o \
       kantor.o \
       katom.o \
       katoms.o \
       kbond.o \
       kbonds.o \
       kcharge.o \
       kchrge.o \
       kdipol.o \
       kdipole.o \
       kewald.o \
       kextra.o \
       keys.o \
       kgeom.o \
       khbond.o \
       kimprop.o \
       kimptor.o \
       kinetic.o \
       kiprop.o \
       kitors.o \
       kmetal.o \
       kmpole.o \
       kmulti.o \
       kopbend.o \
       kopbnd.o \
       kopdist.o \
       kopdst.o \
       korbit.o \
       korbs.o \
       kpitor.o \
       kpitors.o \
       kpolar.o \
       kpolr.o \
       ksolv.o \
       kstbnd.o \
       kstrbnd.o \
       kstrtor.o \
       ksttor.o \
       ktors.o \
       ktorsn.o \
       ktortor.o \
       ktrtor.o \
       kurey.o \
       kurybr.o \
       kvdw.o \
       kvdwpr.o \
       kvdws.o \
       lattice.o \
       lbfgs.o \
       light.o \
       lights.o \
       limits.o \
       linmin.o \
       makeint.o \
       makeref.o \
       makexyz.o \
       math.o \
       maxwell.o \
       mdinit.o \
       mdrest.o \
       mdsave.o \
       mdstat.o \
       mdstuf.o \
       mechanic.o \
       merck.o \
       merge.o \
       minima.o \
       minimize.o \
       minirot.o \
       minrigid.o \
       molcul.o \
       moldyn.o \
       molecule.o \
       molxyz.o \
       moment.o \
       moments.o \
       monte.o \
       mplpot.o \
       mpole.o \
       mutant.o \
       mutate.o \
       nblist.o \
       neigh.o \
       newton.o \
       newtrot.o \
       nextarg.o \
       nexttext.o \
       nonpol.o \
       nose.o \
       nspline.o \
       nucleic.o \
       nucleo.o \
       number.o \
       numeral.o \
       numgrad.o \
       ocvm.o \
       omega.o \
       opbend.o \
       opdist.o \
       openend.o \
       openmp.o \
       optimize.o \
       optirot.o \
       optrigid.o \
       optsave.o \
       orbital.o \
       orbits.o \
       orient.o \
       orthog.o \
       output.o \
       overlap.o \
       params.o \
       path.o \
       paths.o \
       pbstuf.o \
       pdb.o \
       pdbxyz.o \
       phipsi.o \
       picalc.o \
       piorbs.o \
       pistuf.o \
       pitors.o \
       pme.o \
       pmestuff.o \
       pmpb.o \
       polar.o \
       polarize.o \
       poledit.o \
       polgrp.o \
       polpot.o \
       polymer.o \
       potent.o \
       potential.o \
       potfit.o \
       precis.o \
       precise.o \
       pressure.o \
       prmedit.o \
       prmkey.o \
       promo.o \
       protein.o \
       prtdyn.o \
       prterr.o \
       prtint.o \
       prtmol2.o \
       prtpdb.o \
       prtprm.o \
       prtseq.o \
       prtxyz.o \
       pss.o \
       pssrigid.o \
       pssrot.o \
       ptable.o \
       qmstuf.o \
       quatfit.o \
       radial.o \
       random.o \
       rattle.o \
       readdyn.o \
       readgau.o \
       readint.o \
       readmol.o \
       readmol2.o \
       readpdb.o \
       readprm.o \
       readseq.o \
       readxyz.o \
       refer.o \
       replica.o \
       resdue.o \
       respa.o \
       restrn.o \
       rgddyn.o \
       rgdstep.o \
       rigid.o \
       ring.o \
       rings.o \
       rmsfit.o \
       rotbnd.o \
       rotlist.o \
       rotpole.o \
       rxnfld.o \
       rxnpot.o \
       saddle.o \
       scales.o \
       scan.o \
       sdstep.o \
       search.o \
       sequen.o \
       server.o \
       shakeup.o \
       shunt.o \
       sigmoid.o \
       sizes.o \
       sizes.o \
       sktstuff.o \
       sniffer.o \
       socket.o \
       solute.o \
       sort.o \
       spacefill.o \
       spectrum.o \
       square.o \
       stodyn.o \
       strbnd.o \
       strtor.o \
       suffix.o \
       superpose.o \
       surface.o \
       surfatom.o \
       switch.o \
       sybylxyz.o \
       syntrn.o \
       tarray.o \
       temper.o \
       testgrad.o \
       testhess.o \
       testpair.o \
       testpol.o \
       testrot.o \
       timer.o \
       timerot.o \
       titles.o \
       tncg.o \
       torphase.o \
       torpot.o \
       torque.o \
       tors.o \
       torsfit.o \
       torsions.o \
       tortor.o \
       tree.o \
       trimtext.o \
       unitcell.o \
       units.o \
       uprior.o \
       urey.o \
       urypot.o \
       usage.o \
       usolve.o \
       valence.o \
       valfit.o \
       vdw.o \
       vdwpot.o \
       verlet.o \
       version.o \
       vibbig.o \
       vibrate.o \
       vibrot.o \
       vibs.o \
       virial.o \
       volume.o \
       warp.o \
       xtalfit.o \
       xtalmin.o \
       xtals.o \
       xyzatm.o \
       xyzedit.o \
       xyzint.o \
       xyzpdb.o \
       xyzsybyl.o \
       zatom.o \
       zclose.o \
       zcoord.o

EXEFILES = \
           dynamic_iAMOEBA_tensor.x \
           dynamic_3mAMOEBA.x \
           dynamic_3M_AMOEBA.x \
           dynamic_2M_AMOEBA.x \
           alchemy.x \
           analyze.x \
           anneal.x \
           archive.x \
           bar.x \
           correlate.x \
           crystal.x \
           diffuse.x \
           distgeom.x \
           document.x \
           dynamic.x \
           gda.x \
           intedit.x \
           intxyz.x \
           minimize.x \
           minirot.x \
           minrigid.x \
           molxyz.x \
           monte.x \
           newton.x \
           newtrot.x \
           nucleic.x \
           optimize.x \
           optirot.x \
           optrigid.x \
           path.x \
           pdbxyz.x \
           polarize.x \
           poledit.x \
           potential.x \
           prmedit.x \
           protein.x \
           pss.x \
           pssrigid.x \
           pssrot.x \
           radial.x \
           saddle.x \
           scan.x \
           sniffer.x \
           spacefill.x \
           spectrum.x \
           superpose.x \
           sybylxyz.x \
           testgrad.x \
           testhess.x \
           testpair.x \
           testpol.x \
           testrot.x \
           timer.x \
           timerot.x \
           torsfit.x \
           valence.x \
           vibbig.x \
           vibrate.x \
           vibrot.x \
           xtalfit.x \
           xtalmin.x \
           xyzedit.x \
           xyzint.x \
           xyzpdb.x \
           xyzsybyl.x

%.o: %.f
	${F77} ${F77FLAGS} ${OPTFLAGS} $< -o $@ 

%.o: %.c
	${CC} ${CFLAGS} ${INCLUDEDIR} ${OPTFLAGS} $<

%.x: %.o libtinker.a
	${F77} ${LINKFLAGS} -o $@ $^ ${LIBS}; strip $@

all:	${EXEFILES}

clean:
	rm -f *.o *.mod *.a *.x

listing:
	cat *.i *.f *.c > tinker.txt

rename:
	mv  alchemy.x    $(BINDIR)/alchemy
	mv  analyze.x    $(BINDIR)/analyze
	mv  anneal.x     $(BINDIR)/anneal
	mv  archive.x    $(BINDIR)/archive
	mv  bar.x        $(BINDIR)/bar
	mv  correlate.x  $(BINDIR)/correlate
	mv  crystal.x    $(BINDIR)/crystal
	mv  diffuse.x    $(BINDIR)/diffuse
	mv  distgeom.x   $(BINDIR)/distgeom
	mv  document.x   $(BINDIR)/document
	mv  dynamic.x    $(BINDIR)/dynamic
	mv  gda.x        $(BINDIR)/gda
	mv  intedit.x    $(BINDIR)/intedit
	mv  intxyz.x     $(BINDIR)/intxyz
	mv  minimize.x   $(BINDIR)/minimize
	mv  minirot.x    $(BINDIR)/minirot
	mv  minrigid.x   $(BINDIR)/minrigid
	mv  molxyz.x     $(BINDIR)/molxyz
	mv  monte.x      $(BINDIR)/monte
	mv  newton.x     $(BINDIR)/newton
	mv  newtrot.x    $(BINDIR)/newtrot
	mv  nucleic.x    $(BINDIR)/nucleic
	mv  optimize.x   $(BINDIR)/optimize
	mv  optirot.x    $(BINDIR)/optirot
	mv  optrigid.x   $(BINDIR)/optrigid
	mv  path.x       $(BINDIR)/path
	mv  pdbxyz.x     $(BINDIR)/pdbxyz
	mv  polarize.x   $(BINDIR)/polarize
	mv  poledit.x    $(BINDIR)/poledit
	mv  potential.x  $(BINDIR)/potential
	mv  prmedit.x    $(BINDIR)/prmedit
	mv  protein.x    $(BINDIR)/protein
	mv  pss.x        $(BINDIR)/pss
	mv  pssrigid.x   $(BINDIR)/pssrigid
	mv  pssrot.x     $(BINDIR)/pssrot
	mv  radial.x     $(BINDIR)/radial
	mv  saddle.x     $(BINDIR)/saddle
	mv  scan.x       $(BINDIR)/scan
	mv  sniffer.x    $(BINDIR)/sniffer
	mv  spacefill.x  $(BINDIR)/spacefill
	mv  spectrum.x   $(BINDIR)/spectrum
	mv  superpose.x  $(BINDIR)/superpose
	mv  sybylxyz.x   $(BINDIR)/sybylxyz
	mv  testgrad.x   $(BINDIR)/testgrad
	mv  testhess.x   $(BINDIR)/testhess
	mv  testpair.x   $(BINDIR)/testpair
	mv  testpol.x    $(BINDIR)/testpol
	mv  testrot.x    $(BINDIR)/testrot
	mv  timer.x      $(BINDIR)/timer
	mv  timerot.x    $(BINDIR)/timerot
	mv  torsfit.x    $(BINDIR)/torsfit
	mv  valence.x    $(BINDIR)/valence
	mv  vibbig.x     $(BINDIR)/vibbig
	mv  vibrate.x    $(BINDIR)/vibrate
	mv  vibrot.x     $(BINDIR)/vibrot
	mv  xtalfit.x    $(BINDIR)/xtalfit
	mv  xtalmin.x    $(BINDIR)/xtalmin
	mv  xyzedit.x    $(BINDIR)/xyzedit
	mv  xyzint.x     $(BINDIR)/xyzint
	mv  xyzpdb.x     $(BINDIR)/xyzpdb
	mv  xyzsybyl.x   $(BINDIR)/xyzsybyl

remove_links:
	rm -f $(LINKDIR)/alchemy
	rm -f $(LINKDIR)/analyze
	rm -f $(LINKDIR)/anneal
	rm -f $(LINKDIR)/archive
	rm -f $(LINKDIR)/bar
	rm -f $(LINKDIR)/correlate
	rm -f $(LINKDIR)/crystal
	rm -f $(LINKDIR)/diffuse
	rm -f $(LINKDIR)/distgeom
	rm -f $(LINKDIR)/document
	rm -f $(LINKDIR)/dynamic
	rm -f $(LINKDIR)/gda
	rm -f $(LINKDIR)/intedit
	rm -f $(LINKDIR)/intxyz
	rm -f $(LINKDIR)/minimize
	rm -f $(LINKDIR)/minirot
	rm -f $(LINKDIR)/minrigid
	rm -f $(LINKDIR)/molxyz
	rm -f $(LINKDIR)/monte
	rm -f $(LINKDIR)/newton
	rm -f $(LINKDIR)/newtrot
	rm -f $(LINKDIR)/nucleic
	rm -f $(LINKDIR)/optimize
	rm -f $(LINKDIR)/optirot
	rm -f $(LINKDIR)/optrigid
	rm -f $(LINKDIR)/path
	rm -f $(LINKDIR)/pdbxyz
	rm -f $(LINKDIR)/polarize
	rm -f $(LINKDIR)/poledit
	rm -f $(LINKDIR)/potential
	rm -f $(LINKDIR)/prmedit
	rm -f $(LINKDIR)/protein
	rm -f $(LINKDIR)/pss
	rm -f $(LINKDIR)/pssrigid
	rm -f $(LINKDIR)/pssrot
	rm -f $(LINKDIR)/radial
	rm -f $(LINKDIR)/saddle
	rm -f $(LINKDIR)/scan
	rm -f $(LINKDIR)/sniffer
	rm -f $(LINKDIR)/spacefill
	rm -f $(LINKDIR)/spectrum
	rm -f $(LINKDIR)/superpose
	rm -f $(LINKDIR)/sybylxyz
	rm -f $(LINKDIR)/testgrad
	rm -f $(LINKDIR)/testhess
	rm -f $(LINKDIR)/testpair
	rm -f $(LINKDIR)/testpol
	rm -f $(LINKDIR)/testrot
	rm -f $(LINKDIR)/timer
	rm -f $(LINKDIR)/timerot
	rm -f $(LINKDIR)/torsfit
	rm -f $(LINKDIR)/valence
	rm -f $(LINKDIR)/vibbig
	rm -f $(LINKDIR)/vibrate
	rm -f $(LINKDIR)/vibrot
	rm -f $(LINKDIR)/xtalfit
	rm -f $(LINKDIR)/xtalmin
	rm -f $(LINKDIR)/xyzedit
	rm -f $(LINKDIR)/xyzint
	rm -f $(LINKDIR)/xyzpdb
	rm -f $(LINKDIR)/xyzsybyl

create_links:
	ln -s $(BINDIR)/alchemy    $(LINKDIR)/alchemy
	ln -s $(BINDIR)/analyze    $(LINKDIR)/analyze
	ln -s $(BINDIR)/anneal     $(LINKDIR)/anneal
	ln -s $(BINDIR)/archive    $(LINKDIR)/archive
	ln -s $(BINDIR)/bar        $(LINKDIR)/bar
	ln -s $(BINDIR)/correlate  $(LINKDIR)/correlate
	ln -s $(BINDIR)/crystal    $(LINKDIR)/crystal
	ln -s $(BINDIR)/diffuse    $(LINKDIR)/diffuse
	ln -s $(BINDIR)/distgeom   $(LINKDIR)/distgeom
	ln -s $(BINDIR)/document   $(LINKDIR)/document
	ln -s $(BINDIR)/dynamic    $(LINKDIR)/dynamic
	ln -s $(BINDIR)/gda        $(LINKDIR)/gda
	ln -s $(BINDIR)/intedit    $(LINKDIR)/intedit
	ln -s $(BINDIR)/intxyz     $(LINKDIR)/intxyz
	ln -s $(BINDIR)/minimize   $(LINKDIR)/minimize
	ln -s $(BINDIR)/minirot    $(LINKDIR)/minirot
	ln -s $(BINDIR)/minrigid   $(LINKDIR)/minrigid
	ln -s $(BINDIR)/molxyz     $(LINKDIR)/molxyz
	ln -s $(BINDIR)/monte      $(LINKDIR)/monte
	ln -s $(BINDIR)/newton     $(LINKDIR)/newton
	ln -s $(BINDIR)/newtrot    $(LINKDIR)/newtrot
	ln -s $(BINDIR)/nucleic    $(LINKDIR)/nucleic
	ln -s $(BINDIR)/optimize   $(LINKDIR)/optimize
	ln -s $(BINDIR)/optirot    $(LINKDIR)/optirot
	ln -s $(BINDIR)/optrigid   $(LINKDIR)/optrigid
	ln -s $(BINDIR)/path       $(LINKDIR)/path
	ln -s $(BINDIR)/pdbxyz     $(LINKDIR)/pdbxyz
	ln -s $(BINDIR)/polarize   $(LINKDIR)/polarize
	ln -s $(BINDIR)/poledit    $(LINKDIR)/poledit
	ln -s $(BINDIR)/potential  $(LINKDIR)/potential
	ln -s $(BINDIR)/prmedit    $(LINKDIR)/prmedit
	ln -s $(BINDIR)/protein    $(LINKDIR)/protein
	ln -s $(BINDIR)/pss        $(LINKDIR)/pss
	ln -s $(BINDIR)/pssrigid   $(LINKDIR)/pssrigid
	ln -s $(BINDIR)/pssrot     $(LINKDIR)/pssrot
	ln -s $(BINDIR)/radial     $(LINKDIR)/radial
	ln -s $(BINDIR)/saddle     $(LINKDIR)/saddle
	ln -s $(BINDIR)/scan       $(LINKDIR)/scan
	ln -s $(BINDIR)/sniffer    $(LINKDIR)/sniffer
	ln -s $(BINDIR)/spacefill  $(LINKDIR)/spacefill
	ln -s $(BINDIR)/spectrum   $(LINKDIR)/spectrum
	ln -s $(BINDIR)/superpose  $(LINKDIR)/superpose
	ln -s $(BINDIR)/sybylxyz   $(LINKDIR)/sybylxyz
	ln -s $(BINDIR)/testgrad   $(LINKDIR)/testgrad
	ln -s $(BINDIR)/testhess   $(LINKDIR)/testhess
	ln -s $(BINDIR)/testpair   $(LINKDIR)/testpair
	ln -s $(BINDIR)/testpol    $(LINKDIR)/testpol
	ln -s $(BINDIR)/testrot    $(LINKDIR)/testrot
	ln -s $(BINDIR)/timer      $(LINKDIR)/timer
	ln -s $(BINDIR)/timerot    $(LINKDIR)/timerot
	ln -s $(BINDIR)/torsfit    $(LINKDIR)/torsfit
	ln -s $(BINDIR)/valence    $(LINKDIR)/valence
	ln -s $(BINDIR)/vibbig     $(LINKDIR)/vibbig
	ln -s $(BINDIR)/vibrate    $(LINKDIR)/vibrate
	ln -s $(BINDIR)/vibrot     $(LINKDIR)/vibrot
	ln -s $(BINDIR)/xtalfit    $(LINKDIR)/xtalfit
	ln -s $(BINDIR)/xtalmin    $(LINKDIR)/xtalmin
	ln -s $(BINDIR)/xyzedit    $(LINKDIR)/xyzedit
	ln -s $(BINDIR)/xyzint     $(LINKDIR)/xyzint
	ln -s $(BINDIR)/xyzpdb     $(LINKDIR)/xyzpdb
	ln -s $(BINDIR)/xyzsybyl   $(LINKDIR)/xyzsybyl

libtinker.a: ${OBJS}
	ar ${LIBFLAGS} libtinker.a \
        mdsave_writeDyns.o \
        prtdyn_writeDyns.o \
        uind2bmat.o \
        init2bmatmodPolar_uind.o \
        uind1bmatrix.o \
        init1bmat_uind.o \
        uindclust_gradient_Polar1b2b3bsimult_ireducemat_commsmall2.o \
        uindempole1d_3b_Polar1btest2.o \
        uindClustEmpole1plz_1b_NoOmplist.o \
        uindClustEmpole1plz_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o \
        uindClustEmpole1plz_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o \
        uindempole1d_3b_Polar2btest2.o \
        uindempole1d_3b_Polar3btest2.o \
        pol_loadbal_clust_splitlist_offset_2.o \
        pol_loadbal_clust_splitlist_offset3b_2.o \
        ClustEmpole1plz_2b_NoOmplistRtrn2BNoSubtr1bsimult.o \
        clust_gradient_Polar1b2bsimult_ireducemat_commsmall2.o \
        clust_gradient_Polar1b2b3bsimult_ireducemat_commsmall2.o \
        empole1d_3b_Polar1btest2.o \
        empole1d_3b_Polar2btest2.o \
        empole1d_3b_Polar3btest2.o \
        ClustEmpole1plz_1b_NoOmplist.o \
        ClustEmpole1plz_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o \
        ClustEmpole1plz_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o \
        clust_gradient_Polar1b2b3bsimult_reducemat_commsmall2.o \
        grad_covemrecip_vandr_reduce_totfield_commsmall.o \
        clust_gradient_PolarPerm1b2b3bsimult_reducemat_commsmall2.o \
        empole1d_3b_Polar3btest.o \
        ulspred3b.o \
        uprior3b.o \
        initpolarpred3b_offset.o \
        clust3blists_offset.o \
        mlist_3bclust.o \
        ulist_3bclust.o \
        pol_loadbal_clust_splitlist_offset3b_simple.o \
        init2bmatmodPolar.o \
        ClustEmpole1_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o \
        clust_gradient_PolarPerm1b2b3bsimult_reducemat_commsmall.o \
        dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust3TinkEmpSimult1b2b3b.o \
        ClustEmpole1_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o \
        empole1d_3b_PolarPerm3btest.o \
        induce0a_3b_scf_vacomp3bclustlist.o \
        dfield0c_3bPolar_3bclust.o \
        grad_cov_v_reduce_commsmall.o \
        mpidat4.o \
        clust_gradient_PolarPerm1b2bsimult_reducemat_commsmall.o \
        ClustEmpole1_2b_NoOmplistRtrn2BNoSubtr1bsimult.o \
        clust_gradient_covvdwPolarPerm1b2b_ireducemat_commsmall.o \
        ClustEmpole1_2b_NoOmplistRtrn2BNoSubtr1b.o \
        clust_gradient_covvdwPolarPerm1b2b_ireducemat.o \
        nblist_par_simultclust.o \
        loadbal_ehal1c_halftasks_nolistsend_simultclust.o \
        clust_gradient_covvdwPolarPerm1b2b.o \
        cell1bclust.o \
        boxes1bclust.o \
        cell2bclust.o \
        boxes2bclust.o \
        ulspred1b.o \
        initpolarpred1b.o \
        uprior1b.o \
        uprior2b.o \
        initpolarpred2b_offset.o \
        initpolarpred2b.o \
        ulspred2b.o \
        clust_gradient_PolarPerm1b2bsimultmod.o \
        ClustEmpole1_1b_NoOmplistsimult.o \
        ClustEmpole1_2b_NoOmplistsimult.o \
        clust_gradient_PolarPerm1b2bsimult.o \
        init1b2bmatPolarPerm.o \
        deriv2bmat.o \
        pol_loadbal_clust_splitlist_offset.o \
        clust2blists_offset.o \
        clust_gradient_PolarPerm1b2ball.o \
        ClustEmpole1_2b_NoOmpAll2b.o \
        empole1c_3b_PolarPermtest.o \
        empole1d_3b_PolarPerm2btest.o \
        empole1d_3b_PolarPerm1btest.o \
        empole1c_3b_PolarPerm.o \
        ClustEmpole1_2b_NoOmp.o \
        ClustEmpole1_1b_NoOmp.o \
        grad_cov_v_reduce.o \
        init1bmatPolarPerm.o \
        torque3_PolarPerm.o \
        clust_gradient_PolarPerm1b2b.o \
        ClustEmpole1_1b_NoOmplist.o \
        ClustEmpole1_2b_NoOmplist.o \
        empole1d_3b_PolarPerm1b.o \
        empole1d_3b_PolarPerm2b.o \
        dfield0c_3bPolar_2bclust.o \
        dfield0c_3bPolar_1bclust.o \
        ufield0c_3bPolar.o \
        empole1d_3bPolar_origversion_nopriordirectcalc2_omp2b.o \
        empole1d_3bPolar_origversion_nopriordirectcalc2_omp1b.o \
        ereal1d_3b_1bclust.o \
        ereal1d_3b_2bclust.o \
        emrecip1_3b_Polar_new.o \
        pmestuff_3bnew.o \
        mlist_2bclust.o \
        mlist_1bclust.o \
        ulist_2bclust.o \
        induce0a_3b_scf_vacomp1bclustlist.o \
        induce0a_3b_scf_vacomp2bclustlist.o \
        empole1b_3bPolar_origversion_nopriordirectcalc2_omp2b.o \
        empole1b_3bPolar_origversion_nopriordirectcalc2_omp1b.o \
        ulist_1bclust.o \
        clust2blists.o \
        clust1blists.o \
        clustloadbalsmoothInner_1a_2bPolar_vac_precalc1bNoOmplist.o \
        clustnolistsmoothInner_1a_1bPolar_vac_rtrn1bNoOmplist.o \
        clustnolistsmoothInner_1a_1bPolar_vac_rtrn1bNoOmp.o \
        clustloadbalsmoothInner_1a_2bPolar_vac_precalc1bNoOmp.o \
        induce0a_3b_scf_vacomp.o \
        empole1a_3bPolar_origversion_nopriordirectcalc2_omp.o \
        pol_loadbal_clust_splitlist.o \
        pol_loadbal_clust.o \
        clustnolistsmoothInner_1a_1bPolar_vac_rtrn1b.o \
        clustloadbalsmoothInner_1a_2bPolar_vac_precalc1b.o \
        clust_smooth_gradpolar1b2b.o \
        init1bmat.o \
        deriv1bmat.o \
        nblist5kmeansclust.o \
        clustloadbalsmoothInner_1a_2bPolar_vac.o \
        clustnolistsmoothInner_1a_1bPolar_vac.o \
        pcg.o \
        basefileclust.o \
        getclust.o \
        readclust.o \
        induce0a_3b_scf_vac.o \
        induce0a_3b_cholesky.o \
        cho.o \
        fmat.o \
        totfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2fij.o \
        empole1a_3bPolar_totfieldnpolefij.o \
        empole1a1_3b_Polar_totfieldnpole_modfij.o \
        empole1a2_3b_Polar_totfieldnpole_modfij.o \
        torque_3b_new_npolefij.o \
        totfieldNpoleloadbalsmoothInner_1a_3bPolar_w1bodyfij.o \
        empole1a_3bPolar_origversion_nopriordirectcalc2_fijscale2.o \
        empole1a_3bPolar_origversion_nopriordirectcalc2_fijscale.o \
        Newempole1c_3bPolar_totfield_dEtensor_fewertens1.o \
        totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi.o \
        totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi1.o \
        totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi1.o \
        empole1a_3bPolar_origversion_nopriordirectcalc2_fij.o \
        torque_3b_new_fij.o \
        All2bAll3bNew2smooth_innerloop_empole1a_3b_fij.o \
        4bodyAll2bAll3bNew2smooth_innerloop_empole1a_3b.o \
        empole1a_mod.o \
        sAll2bNew2smooth_innerloop_empole1a_3b.o \
        sAll2bAll3bNew2smooth_innerloop_empole1a_3b.o \
        sAll2bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
        sAll2bAll3bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
        All2bAll3bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
        All2bNew2smooth_innerloop_empole1a_3b_subtrdir.o \
        All2bNew2smooth_innerloop_empole1a_3b.o \
        sAll2bNew2smooth_innerloop_empole1a_3b.o \
        s2btotfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2.o \
        s2btotfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o \
        s2btotfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o \
        stotfieldloadbalsmooth_findiff_test.o \
        sulf.o \
        pmestuffNoQPoleNoDPole.o \
        emrecip1_3b_Perm2NoQPoleNoDPole.o \
        ewpol_loadbal7_clust.o \
        loadbalsmoothInner_1a_2bPolar_vac.o \
        empole1a2_3b_Polar_totfieldnpole3bmod_dol1_2.o \
        totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi1.o \
        totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi1.o \
        Newempole1c_3bPolar_totfield_dEtensor_fewertens1.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens1.o \
        thg5_small_vir_nolistsend_Nn1_noewtens_notenssends_fewertens.o \
        torque3_dEtensor_72015_taunew4_pt1_vir_notenssends_fewertens.o \
        torque3_dEtensor_72015_taunew4_pt2_vir_notenssends_fewertens.o \
        empole1a1_3b_Polar_totfieldnpole3b.o \
        empole1a_3bPolar_totfieldnpole3b_rtrndeptemp2.o \
        thg5_small_vir_nolistsend_Nn_noewtens_fewertens.o \
        torque3_dEtensor_72015_taunew3_pt1_vir_fewertens.o \
        torque3_dEtensor_72015_taunew3_pt2_vir_fewertens.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_fewertens.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens2.o \
        tensdealloc.o \
        totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi.o \
        Newempole1c_3bPolar_totfield_dEtensor_fewertens.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens.o \
        thg5_small_vir_nolistsend_Nn_noewtens_notenssends_fewertens.o \
        torque3_dEtensor_72015_taunew3_pt1_vir_notenssends_fewertens.o \
        torque3_dEtensor_72015_taunew3_pt2_vir_notenssends_fewertens.o \
        dEtensor3.o \
        thg5_small_vir_nolistsend_noewtens2.o \
        thg5_small_vir_nolistsend_Nn_notenssends.o \
        torque3_dEtensor_72015_taunew3_pt1_vir_notenssends.o \
        torque3_dEtensor_72015_taunew3_pt2_vir_notenssends.o \
        thg5_small_vir_nolistsend_Nn_noewtens_notenssends.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends.o \
        mpidat3.o \
        dEtensor2.o \
        torque3_dEtensor_72015_taunew3_pt1_vir.o \
        torque3_dEtensor_72015_taunew3_pt2_vir.o \
        nblist_short.o \
        limits2.o \
        neigh2.o \
        switch2.o \
        sizes2.o \
        Newempole1c_3bPolar_totfield_dEtensor.o \
        thg5_small_vir_nolistsend_Nn_noewtens.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3.o \
        cutoffs_ewald_short.o \
        ewpol_loadbal7.o \
        nblist2_hbond.o \
        nblist2_test.o \
        clustnolistsmoothInner_1c_3bPolar_vac.o \
        nblist4_com.o \
        nblist3_com.o \
        nblist4.o \
        nblist3.o \
        mpidat2.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2_no3bmut.o \
        moments_onlyDipole_fixedcharge.o \
        uindInner_1a_3bPolar_vac.o \
        pmonte_par2.o \
        empole0a_3bPolar_origversion_nopriordirectcalc2.o \
        eng_r_reduce.o \
        loadbalsmoothInner_0a_3bPolar_vac.o \
        smooth_engpolar.o \
        empole_Perm.o \
        eng_covalent_vdw.o \
        verlet_par_pmonte.o \
        pmonte_par.o \
        uindempole1a_3bPolar_origversion_nopriordirectcalc2.o \
        uindloadbalsmoothInner_1a_3bPolar_vac.o \
        uindsmooth_gradpolar.o \
        clust_smooth_gradpolar.o \
        clustnolistsmoothInner_1a_3bPolar_vac.o \
        clustloadbalsmoothInner_1a_3bPolar_vac.o \
        nblist2.o \
        nblist5.o \
        neigh2clust.o \
        pol_loadbal7.o \
        pol_loadbal6.o \
        thg5_small_vir_nolistsend_noewtens.o \
        empole1c_3bPolar_totfield_dEtensor.o \
        virialnamoeba.o \
        dfield0a_ufield0a_uscale0a.o \
        empole1a_PolarStripped.o \
        induce0a_nobranch.o \
        torque_pol.o \
        ewpol_loadbal_halftasks.o \
        loadbalsmoothInner_1a_3bPolar1bmut_halftasks.o \
        loadbalsmoothInner_1a_3bPolar_halftasks.o \
        loadbalsmoothInner_1a_3bPolar_vac.o \
        loadbalsmoothInner_1a_3bPolar1bmut.o \
        dfield0b_totfield.o \
        LoadBalInner_r_Perm_nolistsend.o \
        grad_covemrecip_vandr_reduce_totfield.o \
        empole1a_3bPolar_origversion_nopriordirectcalc_rtrnfield.o \
        induce0a_3b_PolelecOnly_totfieldnpolevac.o \
        empole1a2_3b_Polar_totfieldnpolevac.o \
        empole1a2_3b_Polar_totfieldnpole_mod.o \
        empole1a_3bPolar_new_subtract_1bodymut.o \
        empole1a_3bPolar_totfieldnpolevac.o \
        empole1a_3bPolar_totfieldnpole.o \
        empole1a_3bPolar_origversion_nopriordirectcalc2.o \
        emrecip1_1_3b_Polar_totfield_parl_513.o \
        empole1c_3bPolar_totfield_901.o \
        emrecip1_3b_Polar_usemlist_mutual_virialonly.o \
        emrecip1_3b_Polar_usemlist_mutual.o \
        induce0c_field_amat.o \
        ereal1c1_3b_totfieldnpole.o \
        totfieldnpole3b_noewald_umutual_rl_3b_nofield.o \
        empole1a1_3b_Polar_totfieldnpole_mod.o \
        induce0a_3b_PolelecOnly_totfieldnpole.o \
        empole1a1_3bPolar_totfield_dEtensor1bmut.o \
        emrecip1_3b_Polar_usemlist_direct_virialonly.o \
        b2test.o \
        emrecip1_3b_Polar_usemlist_direct_test2.o \
        emrecip1_3b_Polar_usemlist_direct_test.o \
        empole1d_3b_PolarPerm_selfeng_bcast2.o \
        grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2.o \
        bcast_cov.o \
        grad_r_allreduce_totfield_dEtensor_vir.o \
        thg6_small_vir_nolistsend.o \
        grad_vandr_ireduce_totfield_dEtensor_vir.o \
        grad_vandr_ireduce_totfield_dEtensor_vir2.o \
        loadbalsmoothInner_1a_2bPolar.o \
        nblist3b_1list_par_newesttink7.o \
        chkpole_omp.o \
        rotpole_omp.o \
        grad_vandr_reduce_totfield_dEtensor_vir2.o \
        nblist_par.o \
        loadbal_ehal1c_halftasks_nolistsend.o \
        ereal1d_Perm_savetensor.o \
        kewald_eastman_list.o \
        chunklight.o \
        pmestuff_list.o \
        v_loadbal_halftasks.o \
        loadbal_ehal1c_halftasks.o \
        grad_vandr_reduce_totfield_dEtensor_vir.o \
        empole1c_3bPolar_totfield_dEtensor1b_nouselist_mpivir.o \
        torque3_dEtensor_72015_taunew2_pt1_vir.o \
        torque3_dEtensor_72015_taunew2_pt2_vir.o \
        thg5_vir.o \
        thg5_small_vir.o \
        grad_r_reduce_totfield_dEtensor_vir.o \
        thg6_save.o \
        thg6_save_small.o \
        torque2_omp.o \
        table_fill_setchunk_omp.o \
        bspline_fill_bsplgen_omp.o \
        LoadBalInner_thg4omp_mpi_72015torque2.o \
        torque3_dEtensor_72015_taunew2_pt1.o \
        torque3_dEtensor_72015_taunew2_pt2.o \
        thg5.o \
        thg5_small.o \
        emrecip1_3b_Perm2_totfield_mpi_gridmpole.o \
        emrecip1_3b_Perm2_totfield_mpi_fphimpole.o \
        nose_tuckerman06_par.o \
        tuckerman06.o \
        emrecip1_3b_Perm2_totfield_mpi2.o \
        grid_mpole_Mpi3amoebaRepo.o \
        gradient_emrecip_totfield_mpi.o \
        emrecip1_3b_Perm2_totfield_mpi.o \
        fphi_mpole_Mpi3amoebaRepo.o \
        kewald_nopmealloc_eastman.o \
        mech_pme4_eastman.o \
        kewald_eastman.o \
        gradient_emrecip_totfield_send.o \
        allocvdw_master2.o \
        bcast_vdw.o \
        v_loadbal.o \
        parvdwneigh.o \
        loadbal_ehal1c.o \
        gradient_covalent2.o \
        grad_v_reduce.o \
        field_r_reduce_uindcalc.o \
        udirect2b_uinddir_mpi.o \
        LoadBalInnerloop_ereal1d_sentlist2.o \
        udirect2b_uinddir.o \
        empole1c_3bPolar_totfield_dEtensor1b_nouselist_mpi.o \
        empole1c_3bPolar_totfield_dEtensor1b_nouselist.o \
        totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi.o \
        totfield_smooth_gradpolar_dEtensor.o \
        LoadBalInner_thg4omp_mpi.o \
        grad_r_reduce_totfield_dEtensor.o \
        grad_r_serial_totfield_dEtensor.o \
        ewtotfieldsmooth2_gradient_polar_serial_dEtensor.o \
        totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensor.o \
        dEtensor.o \
        empole1c_3bPolar_totfield_dEtensor.o \
        LoadBalInner_r_Perm_totfield_bcastlist.o \
        grad_r_reduce_totfield_bcast.o \
        empole1c_3bPolar_totfield_usemlist.o \
        emrecip1_3b_Polar_usemlist_direct.o \
        empole1d_3b_PolarPerm_selfeng_bcast.o \
        grad_r_reduce_totfield_uindcalc.o \
        ereal1d2_3b_totfieldnpole_mlist.o \
        ereal1c2_3b_totfieldnpole_omp.o \
        totfieldNpolesmoothInner_1c_1bPolar2_513_alliter.o \
        totfieldNpolesmoothInner_1c_1bPolar2_513_multiter.o \
        totfieldNpolesmoothInner_1c_1bPolar2_513.o \
        ewtotfieldsmallsmooth2_gradient_polar_serial.o \
        grad_r_serial_totfield.o \
        sendmlist_load_balanced3_sizelist_serial.o \
        pmestuff_3bpolz_513.o \
        induce_totfield.o \
        pmestuff_3bpolz_totfield_513.o \
        pmestuff_3bpolz_totfield_pt2_513.o \
        emrecip1_1and2_3b_Polar_totfield_parl_513.o \
        kewald3b_totfield.o \
        bspline_table_recipsum.o \
        aprx.o \
        totfieldNpoleloadbalsmoothInner_1c_2bPolar2.o \
        totfieldNpoleloadbalsmoothInner_1c_2bPolar.o \
        nblist3b_1list.o \
        alloc3b_vdwlist_1list.o \
        totfieldNpolesmoothInner_1c_3bPolar.o \
        totfieldNpoleloadbalsmoothInner_1c_3bPolar.o \
        ereal1d2_3b_totfieldnpole.o \
        r_loadbal2.o \
        gradient_emrecip_totfield.o \
        emrecip1_3b_Perm2_totfield.o \
        LoadBalInner_r_Perm_totfield.o \
        grad_r_reduce_totfield.o \
        pmestuff_3bpolz_totfield_pt2.o \
        pmestuff_3bpolz_totfield.o \
        dfield0c_totfield.o \
        dfield0a_totfield.o \
        ereal1c1_3b_totfieldnpole.o \
        ereal1c2_3b_totfieldnpole.o \
        empole1c_3bPolar_totfield.o \
        All2bAll3bNew2smooth_innerloop_empole1c_3b_totfield.o \
        All2bAll3bNew2smooth_innerloop_empole1c_3b.o \
        induce.o \
        empole1c_3bPolar.o \
        pmestuff_3bpolz_pt2.o \
        torque3_Polar.o \
        torque_3b_new_npole.o \
        totfield.o \
        totfield_smooth_gradpolar.o \
        fft3d_3b.o \
        kewald3b.o \
        pmestuff_3bpolz.o \
        emrecip1_3b_Polar.o \
        All2bAll3bNew2smooth_innerloop_empole1a_3b.o \
        ewreg3bpolz.o \
        kewaldreg3b.o \
        ewpol_loadbal.o \
        mdstat_propertiesDaveDeTo.o \
        alloc3b_uind.o \
        moments_onlyDipole.o \
        get_ntask_r.o \
        LoadBalcom_smooth2b_innerloop_empole1a_3b.o \
        r_loadbal.o \
        emreal_load_balance.o \
        sendmlist_load_balanced2.o \
        sendmlist_load_balanced.o \
        LoadBalInner_r_Perm.o \
        loadbalsmoothInner_1a_3bPolar.o \
        pol_loadbal.o \
        empole1a_Perm.o \
        torque_Perm.o \
        nolistcom_smooth2b_innerloop_empole1a_3b.o \
        ptest_parallel.o \
        emrecip1PermReduce.o \
        mech_pme4.o \
        kewald_nopmealloc.o \
        grad_k_reduce.o \
        paremneigh.o \
        sendmlist.o \
        Inner_r_Perm.o \
        ewaldneigh.o \
        parewreg.o \
        kewaldreg.o \
        ewreg.o \
        smooth_gradpolar.o \
        smooth2b_gradient_polar.o \
        grad_r_reduce.o \
        PermElec2.o \
        Modgrad_3b.o \
        neigh3b.o \
        mechanic_pme.o \
        mech_pme3.o \
        mech_par.o \
        deriv3b.o \
        PermElec.o \
        nblist3b.o \
        par_nblist_2b.o \
        mpidat.o \
        empole1a_3bPolar.o \
        gradient_3b.o \
        nose_par.o \
        verlet_par.o \
        beeman_pt1_pt2_master.o \
        mdinit_par.o \
        cobar.o \
        gradslow_serial.o \
        action.o \
        active.o \
        align.o \
        analysis.o \
        analyz.o \
        angang.o \
        angbnd.o \
        angles.o \
        angpot.o \
        angtor.o \
        argue.o \
        ascii.o \
        atmlst.o \
        atomid.o \
        atoms.o \
        attach.o \
        basefile.o \
        bath.o \
        beeman.o \
        bicubic.o \
        bitor.o \
        bitors.o \
        bndpot.o \
        bndstr.o \
        bonds.o \
        born.o \
        bound.o \
        bounds.o \
        boxes.o \
        bussi.o \
        calendar.o \
        cell.o \
        center.o \
        charge.o \
        chgpot.o \
        chkpole.o \
        chkring.o \
        chkxyz.o \
        cholesky.o \
        chrono.o \
        chunks.o \
        clock.o \
        cluster.o \
        column.o \
        command.o \
        connect.o \
        connolly.o \
        control.o \
        couple.o \
        cspline.o \
        cutoffs.o \
        cutoffs_vdw_nolist.o \
        deflate.o \
        delete.o \
        deriv.o \
        diagq.o \
        diffeq.o \
        dipole.o \
        disgeo.o \
        dma.o \
        domega.o \
        eangang.o \
        eangang1.o \
        eangang2.o \
        eangang3.o \
        eangle.o \
        eangle1.o \
        eangle2.o \
        eangle3.o \
        eangtor.o \
        eangtor1.o \
        eangtor2.o \
        eangtor3.o \
        ebond.o \
        ebond1.o \
        ebond2.o \
        ebond3.o \
        ebuck.o \
        ebuck1.o \
        ebuck2.o \
        ebuck3.o \
        echarge.o \
        echarge1.o \
        echarge2.o \
        echarge3.o \
        echgdpl.o \
        echgdpl1.o \
        echgdpl2.o \
        echgdpl3.o \
        edipole.o \
        edipole1.o \
        edipole2.o \
        edipole3.o \
        egauss.o \
        egauss1.o \
        egauss2.o \
        egauss3.o \
        egeom.o \
        egeom1.o \
        egeom2.o \
        egeom3.o \
        ehal.o \
        ehal1.o \
        ehal2.o \
        ehal3.o \
        eimprop.o \
        eimprop1.o \
        eimprop2.o \
        eimprop3.o \
        eimptor.o \
        eimptor1.o \
        eimptor2.o \
        eimptor3.o \
        elj.o \
        elj1.o \
        elj2.o \
        elj3.o \
        embed.o \
        emetal.o \
        emetal1.o \
        emetal2.o \
        emetal3.o \
        emm3hb.o \
        emm3hb1.o \
        emm3hb2.o \
        emm3hb3.o \
        empole.o \
        empole1.o \
        empole2.o \
        empole3.o \
        energi.o \
        energy.o \
        eopbend.o \
        eopbend1.o \
        eopbend2.o \
        eopbend3.o \
        eopdist.o \
        eopdist1.o \
        eopdist2.o \
        eopdist3.o \
        epitors.o \
        epitors1.o \
        epitors2.o \
        epitors3.o \
        erf.o \
        erxnfld.o \
        erxnfld1.o \
        erxnfld2.o \
        erxnfld3.o \
        esolv.o \
        esolv1.o \
        esolv2.o \
        esolv3.o \
        estrbnd.o \
        estrbnd1.o \
        estrbnd2.o \
        estrbnd3.o \
        estrtor.o \
        estrtor1.o \
        estrtor2.o \
        estrtor3.o \
        etors.o \
        etors1.o \
        etors2.o \
        etors3.o \
        etortor.o \
        etortor1.o \
        etortor2.o \
        etortor3.o \
        eurey.o \
        eurey1.o \
        eurey2.o \
        eurey3.o \
        evcorr.o \
        ewald.o \
        extra.o \
        extra1.o \
        extra2.o \
        extra3.o \
        faces.o \
        fatal.o \
        fft.o \
        fft3d.o \
        fftpack.o \
        field.o \
        fields.o \
        files.o \
        final.o \
        flatten.o \
        fracs.o \
        freeunit.o \
        freeze.o \
        geometry.o \
        getint.o \
        getkey.o \
        getmol.o \
        getmol2.o \
        getnumb.o \
        getpdb.o \
        getprm.o \
        getref.o \
        getstring.o \
        gettext.o \
        getword.o \
        getxyz.o \
        ghmcstep.o \
        gkstuf.o \
        gradient.o \
        gradrgd.o \
        gradrot.o \
        group.o \
        groups.o \
        grpline.o \
        gyrate.o \
        hescut.o \
        hessian.o \
        hessn.o \
        hessrgd.o \
        hessrot.o \
        hpmf.o \
        hybrid.o \
        image.o \
        impose.o \
        improp.o \
        imptor.o \
        inertia.o \
        inform.o \
        initatom.o \
        initial.o \
        initprm.o \
        initres.o \
        initrot.o \
        insert.o \
        inter.o \
        invbeta.o \
        invert.o \
        iounit.o \
        jacobi.o \
        kanang.o \
        kangang.o \
        kangle.o \
        kangs.o \
        kangtor.o \
        kantor.o \
        katom.o \
        katoms.o \
        kbond.o \
        kbonds.o \
        kcharge.o \
        kchrge.o \
        kdipol.o \
        kdipole.o \
        kewald.o \
        kextra.o \
        keys.o \
        kgeom.o \
        khbond.o \
        kimprop.o \
        kimptor.o \
        kinetic.o \
        kiprop.o \
        kitors.o \
        kmetal.o \
        kmpole.o \
        kmulti.o \
        kopbend.o \
        kopbnd.o \
        kopdist.o \
        kopdst.o \
        korbit.o \
        korbs.o \
        kpitor.o \
        kpitors.o \
        kpolar.o \
        kpolr.o \
        ksolv.o \
        kstbnd.o \
        kstrbnd.o \
        kstrtor.o \
        ksttor.o \
        ktors.o \
        ktorsn.o \
        ktortor.o \
        ktrtor.o \
        kurey.o \
        kurybr.o \
        kvdw.o \
        kvdwpr.o \
        kvdws.o \
        lattice.o \
        lbfgs.o \
        light.o \
        lights.o \
        limits.o \
        linmin.o \
        makeint.o \
        makeref.o \
        makexyz.o \
        math.o \
        maxwell.o \
        mdinit.o \
        mdrest.o \
        mdsave.o \
        mdstat.o \
        mdstuf.o \
        mechanic.o \
        merck.o \
        merge.o \
        minima.o \
        molcul.o \
        moldyn.o \
        molecule.o \
        moment.o \
        moments.o \
        mplpot.o \
        mpole.o \
        mutant.o \
        mutate.o \
        nblist.o \
        neigh.o \
        nextarg.o \
        nexttext.o \
        nonpol.o \
        nose.o \
        nspline.o \
        nucleo.o \
        number.o \
        numeral.o \
        numgrad.o \
        ocvm.o \
        omega.o \
        opbend.o \
        opdist.o \
        openend.o \
        openmp.o \
        optsave.o \
        orbital.o \
        orbits.o \
        orient.o \
        orthog.o \
        output.o \
        overlap.o \
        params.o \
        paths.o \
        pbstuf.o \
        pdb.o \
        phipsi.o \
        picalc.o \
        piorbs.o \
        pistuf.o \
        pitors.o \
        pme.o \
        pmestuff.o \
        pmpb.o \
        polar.o \
        polgrp.o \
        polpot.o \
        polymer.o \
        potent.o \
        potfit.o \
        precis.o \
        precise.o \
        pressure.o \
        prmkey.o \
        promo.o \
        prtdyn.o \
        prterr.o \
        prtint.o \
        prtmol2.o \
        prtpdb.o \
        prtprm.o \
        prtseq.o \
        prtxyz.o \
        ptable.o \
        qmstuf.o \
        quatfit.o \
        random.o \
        rattle.o \
        readdyn.o \
        readgau.o \
        readint.o \
        readmol.o \
        readmol2.o \
        readpdb.o \
        readprm.o \
        readseq.o \
        readxyz.o \
        refer.o \
        replica.o \
        resdue.o \
        respa.o \
        restrn.o \
        rgddyn.o \
        rgdstep.o \
        rigid.o \
        ring.o \
        rings.o \
        rmsfit.o \
        rotbnd.o \
        rotlist.o \
        rotpole.o \
        rxnfld.o \
        rxnpot.o \
        scales.o \
        sdstep.o \
        search.o \
        sequen.o \
        server.o \
        shakeup.o \
        shunt.o \
        sigmoid.o \
        sizes.o \
        sktstuff.o \
        socket.o \
        solute.o \
        sort.o \
        square.o \
        stodyn.o \
        strbnd.o \
        strtor.o \
        suffix.o \
        surface.o \
        surfatom.o \
        switch.o \
        syntrn.o \
        tarray.o \
        temper.o \
        titles.o \
        tncg.o \
        torphase.o \
        torpot.o \
        torque.o \
        tors.o \
        torsions.o \
        tortor.o \
        tree.o \
        trimtext.o \
        unitcell.o \
        units.o \
        uprior.o \
        urey.o \
        urypot.o \
        usage.o \
        usolve.o \
        valfit.o \
        vdw.o \
        vdwpot.o \
        verlet.o \
        version.o \
        vibs.o \
        virial.o \
        volume.o \
        warp.o \
        xtals.o \
        xyzatm.o \
        zatom.o \
        zclose.o \
        zcoord.o
	${RANLIB} libtinker.a

###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################
dynamic_iAMOEBA_tensor.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o mpidat2.o 
dynamic_3mAMOEBA.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o mpidat2.o mpidat3.o
dynamic_3M_AMOEBA.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o mpidat3.o
dynamic_2M_AMOEBA.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o mpidat3.o
totfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2fij.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
empole1a_3bPolar_totfieldnpolefij.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a1_3b_Polar_totfieldnpole_modfij.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole_modfij.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
torque_3b_new_npolefij.o: atoms.o deriv.o mpole.o sizes.o
totfieldNpoleloadbalsmoothInner_1a_3bPolar_w1bodyfij.o:sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o fmat.o
empole1a_3bPolar_origversion_nopriordirectcalc2_fijscale2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_3bPolar_origversion_nopriordirectcalc2_fijscale.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_3bPolar_origversion_nopriordirectcalc2_fij.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
torque_3b_new_fij.o: atoms.o mpole.o sizes.o
All2bAll3bNew2smooth_innerloop_empole1a_3b_fij.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
4bodyAll2bAll3bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
sulf.o:
fmat.o:
thzuinddyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust4Simult1b2b3b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o mpidat3.o uind1bmatrix.o uind2bmat.o
mdsave_writeDyns.o: atomid.o atoms.o bound.o boxes.o files.o group.o inform.o iounit.o mdstuf.o moldyn.o mpole.o output.o polar.o potent.o rgddyn.o sizes.o socket.o titles.o units.o
prtdyn_writeDyns.o: atoms.o boxes.o files.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o titles.o
uind2bmat.o:
uind1bmatrix.o:
init2bmatmodPolar_uind.o: uind2bmat.o neigh2clust.o
init1bmat_uind.o: uind1bmatrix.o neigh2clust.o
uindclust_gradient_Polar1b2b3bsimult_ireducemat_commsmall2.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o mpidat4.o uind1bmatrix.o uind2bmat.o cobar.o
uindempole1d_3b_Polar1btest2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
uindClustEmpole1plz_1b_NoOmplist.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
uindClustEmpole1plz_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
uindClustEmpole1plz_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
uindempole1d_3b_Polar2btest2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
uindempole1d_3b_Polar3btest2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
uinddyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust4Simult1b2b3b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o mpidat3.o uind1bmatrix.o uind2bmat.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust5Simult1b2b3b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o mpidat3.o
pol_loadbal_clust_splitlist_offset_2.o: mpidat.o neigh3b.o molcul.o neigh2clust.o
pol_loadbal_clust_splitlist_offset3b_2.o: mpidat.o neigh3b.o molcul.o neigh2clust.o
ClustEmpole1plz_2b_NoOmplistRtrn2BNoSubtr1bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o cobar.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust4Simult1b2b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o mpidat3.o
clust_gradient_Polar1b2bsimult_ireducemat_commsmall2.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o mpidat4.o
clust_gradient_Polar1b2b3bsimult_ireducemat_commsmall2.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o mpidat4.o
empole1d_3b_Polar1btest2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1d_3b_Polar2btest2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1d_3b_Polar3btest2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
ClustEmpole1plz_1b_NoOmplist.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
ClustEmpole1plz_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
ClustEmpole1plz_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
clust_gradient_Polar1b2b3bsimult_reducemat_commsmall2.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
grad_covemrecip_vandr_reduce_totfield_commsmall.o: paremneigh.o  mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o mpidat4.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust4Simult1b2b3b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o
clust_gradient_PolarPerm1b2b3bsimult_reducemat_commsmall2.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o mpidat4.o
empole1d_3b_Polar3btest.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
ulspred3b.o: sizes.o mpole.o uprior.o uprior3b.o 
grad_cov_v_reduce_commsmall.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o dEtensor3.o paremneigh.o mpidat4.o
mpidat4.o:
clust_gradient_PolarPerm1b2bsimult_reducemat_commsmall.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
ClustEmpole1_2b_NoOmplistRtrn2BNoSubtr1bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
clust_gradient_covvdwPolarPerm1b2b_ireducemat_commsmall.o: 
ClustEmpole1_2b_NoOmplistRtrn2BNoSubtr1b.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
clust_gradient_covvdwPolarPerm1b2b_ireducemat.o: 
nblist_par_simultclust.o:
loadbal_ehal1c_halftasks_nolistsend_simultclust.o: parvdwneigh.o
clust_gradient_covvdwPolarPerm1b2b.o:
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust3TinkEmpSimult1bVdw.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o
boxes1bclust.o:
unitcell1bclust.o:
lattice1bclust.o: math.o cell1bclust.o boxes1bclust.o
boxes2bclust.o:
cell1bclust.o:
cell2bclust.o:
clust1blatticeboxes.o: cell1bclust.o boxes1bclust.o
lattice2bclust.o: math.o cell2bclust.o boxes2bclust.o
clust2blatticeboxes.o: cell2bclust.o boxes2bclust.o
replicabox1bclust.o:
mlist_1bclustboxclust.o:
lightnimagenbox1bclust.o: cell1bclust.o boxes1bclust.o
ulist_1bclustboxclust.o:
lightnimagenbox2bclust.o: cell2bclust.o boxes2bclust.o
mlist_2bclustboxclust.o:
ulist_2bclustboxclust.o:
replicabox2bclust.o:
clust2blistsboxclust.o: 
clust1blistsboxclust.o:
bspline_fill_omp1bclust.o:
bspline_fill_omp2bclust.o:
switchbox1bclust.o:
image1b.o: cell1bclust.o
udirect2b_3b_vacomp1btestbox1bclust.o:
ufield0c_3bPolarbox1bclust.o:
uscale0b_3b_vacomp1btestbox1bclust.o:
uscale0b_3b_vacomp2btestbox2bclust.o:
udirect1_3b_newbox1bclust.o: boxes1bclust.o
udirect1_3b_newbox2bclust.o: boxes2bclust.o
ClustScratchpmestuff_3bpolztest.o:
ClustScratchpmestuff.o:
image2b.o: cell2bclust.o
switchbox2bclust.o:
udirect2b_3b_vacomp2btestbox2bclust.o:
ufield0c_3bPolarbox2bclust.o:
kewald1bclust.o: boxes1bclust.o
kewald2bclust.o: boxes2bclust.o
empole1d_3b_PolarPerm1btestboxclust.o: boxes1bclust.o
empole1d_3b_PolarPerm2btestboxclust.o: boxes2bclust.o
initpolarpred3b_offset.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o uprior.o uprior3b.o
uprior3b.o:
clust3blists_offset.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o
mlist_3bclust.o: neigh2clust.o atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o neigh2.o sizes2.o
ulist_3bclust.o: neigh2clust.o atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
pol_loadbal_clust_splitlist_offset3b_simple.o: mpidat.o neigh3b.o molcul.o neigh2clust.o
init2bmatmodPolar.o: deriv2bmat.o neigh2clust.o
ClustEmpole1_2b_NoOmplistRtrn2BmatNoSubtr1bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o cobar.o
clust_gradient_PolarPerm1b2b3bsimult_reducemat_commsmall.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o mpidat3.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust3TinkEmpSimult1b2b3b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o
ClustEmpole1_3b_NoOmplistRtrn3BNoSubtr1b2bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o
empole1d_3b_PolarPerm3btest.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
induce0a_3b_scf_vacomp3bclustlist.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o uprior.o usolve.o mpidat.o fft.o openmp.o
dfield0c_3bPolar_3bclust.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o tarray.o openmp.o
ulspred1b.o: sizes.o mpole.o uprior.o uprior1b.o
initpolarpred1b.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o uprior.o uprior1b.o
uprior1b.o:
uprior2b.o:
initpolarpred2b_offset.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o uprior.o uprior2b.o
initpolarpred2b.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o uprior.o uprior2b.o
ulspred2b.o: sizes.o mpole.o uprior.o uprior2b.o
clust_gradient_PolarPerm1b2bsimultmod.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o neigh3b.o 
ClustEmpole1_1b_NoOmplistsimult.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o pcg.o
ClustEmpole1_2b_NoOmplistsimult.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o pcg.o
clust_gradient_PolarPerm1b2bsimult.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o deriv2bmat.o molcul.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust3TinkEmpSimult1b2b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o deriv2bmat.o polpot.o uprior.o mpidat4.o
init1b2bmatPolarPerm.o: deriv1bmat.o neigh2clust.o deriv2bmat.o
deriv2bmat.o:
pol_loadbal_clust_splitlist_offset.o: mpidat.o neigh3b.o molcul.o neigh2clust.o aprx.o
clust2blists_offset.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o 
clust_gradient_PolarPerm1b2ball.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o
ClustEmpole1_2b_NoOmpAll2b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
empole1c_3b_PolarPermtest.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1d_3b_PolarPerm2btest.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1d_3b_PolarPerm1btest.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1c_3b_PolarPerm.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
ClustEmpole1_2b_NoOmp.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
ClustEmpole1_1b_NoOmp.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
grad_cov_v_reduce.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o dEtensor3.o paremneigh.o
init1bmatPolarPerm.o: deriv1bmat.o neigh2clust.o
torque3_PolarPerm.o: atoms.o mpole.o sizes.o
clust_gradient_PolarPerm1b2b.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o deriv.o energi.o
ClustEmpole1_1b_NoOmplist.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
ClustEmpole1_2b_NoOmplist.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
empole1d_3b_PolarPerm1b.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1d_3b_PolarPerm2b.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust3TinkEmp.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o mpidat2.o usolve.o tarray.o polpot.o
dfield0c_3bPolar_2bclust.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o tarray.o openmp.o
dfield0c_3bPolar_1bclust.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o
ufield0c_3bPolar.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o
empole1d_3bPolar_origversion_nopriordirectcalc2_omp2b.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1d_3bPolar_origversion_nopriordirectcalc2_omp1b.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
ereal1d_3b_1bclust.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
ereal1d_3b_2bclust.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
emrecip1_3b_Polar_new.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o fft.o
pmestuff_3bnew.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
mlist_2bclust.o: neigh2clust.o atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
mlist_1bclust.o: neigh2clust.o atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
ulist_2bclust.o: neigh2clust.o atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
induce0a_3b_scf_vacomp1bclustlist.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o uprior.o usolve.o mpidat.o
induce0a_3b_scf_vacomp2bclustlist.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o neigh2clust.o uprior.o usolve.o
empole1b_3bPolar_origversion_nopriordirectcalc2_omp2b.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
empole1b_3bPolar_origversion_nopriordirectcalc2_omp1b.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o neigh2clust.o
ulist_1bclust.o: neigh2clust.o atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
clust2blists.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o
clust1blists.o: neigh2clust.o mpidat.o limits.o molcul.o sizes.o neigh.o neigh3b.o
clustloadbalsmoothInner_1a_2bPolar_vac_precalc1bNoOmplist.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
clustnolistsmoothInner_1a_1bPolar_vac_rtrn1bNoOmplist.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
clustnolistsmoothInner_1a_1bPolar_vac_rtrn1bNoOmp.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
clustloadbalsmoothInner_1a_2bPolar_vac_precalc1bNoOmp.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
induce0a_3b_scf_vacomp.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_3bPolar_origversion_nopriordirectcalc2_omp.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o
pol_loadbal_clust_splitlist.o: mpidat.o neigh3b.o molcul.o neigh2clust.o
pol_loadbal_clust.o: mpidat.o neigh3b.o molcul.o neigh2clust.o
clustnolistsmoothInner_1a_1bPolar_vac_rtrn1b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
clustloadbalsmoothInner_1a_2bPolar_vac_precalc1b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o deriv1bmat.o
clust_smooth_gradpolar1b2b.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o
init1bmat.o: deriv1bmat.o neigh2clust.o
deriv1bmat.o:
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust3.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o
nblist5kmeansclust.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o aprx.o cobar.o
clustloadbalsmoothInner_1a_2bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o
clustnolistsmoothInner_1a_1bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o
pcg.o:
basefileclust.o: ascii.o files.o
getclust.o: inform.o iounit.o output.o
readclust.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o neigh2clust.o molcul.o mpidat.o
induce0a_3b_scf_vac.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
induce0a_3b_cholesky.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
cho.o:
ewpol_loadbal7_clust.o: mpidat.o neigh3b.o molcul.o neigh2clust.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_nopolz.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o mpidat2.o
loadbalsmoothInner_1a_2bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o 
totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi1.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi1.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
Newempole1c_3bPolar_totfield_dEtensor_fewertens1.o: sizes2.o dEtensor2.o sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o dEtensor3.o limits2.o neigh2.o 
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens1.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o dEtensor3.o
thg5_small_vir_nolistsend_Nn1_noewtens_notenssends_fewertens.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o totfield.o dEtensor3.o
torque3_dEtensor_72015_taunew4_pt1_vir_notenssends_fewertens.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o dEtensor3.o
torque3_dEtensor_72015_taunew4_pt2_vir_notenssends_fewertens.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o dEtensor3.o
thg5_small_vir_nolistsend_Nn_noewtens_fewertens.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o totfield.o dEtensor3.o
torque3_dEtensor_72015_taunew3_pt1_vir_fewertens.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o dEtensor3.o
torque3_dEtensor_72015_taunew3_pt2_vir_fewertens.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o dEtensor3.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_fewertens.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o dEtensor3.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens2.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o dEtensor3.o
dEtensor_test_nolistsend_separatereduce_funcdecomp_NN2b.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o mpidat2.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o
nblist2_hbond.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
totfieldNpoleloadbalsmoothInner_1c_2bPolar_dEtensorMpi.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
tensdealloc.o: dEtensor2.o dEtensor3.o
Newempole1c_3bPolar_totfield_dEtensor_fewertens.o: sizes2.o dEtensor2.o sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o dEtensor3.o limits2.o neigh2.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends_fewertens.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o dEtensor3.o
thg5_small_vir_nolistsend_Nn_noewtens_notenssends_fewertens.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o totfield.o dEtensor3.o
torque3_dEtensor_72015_taunew3_pt1_vir_notenssends_fewertens.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o dEtensor3.o
torque3_dEtensor_72015_taunew3_pt2_vir_notenssends_fewertens.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o dEtensor3.o
dEtensor3.o:
thg5_small_vir_nolistsend_noewtens2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o 
thg5_small_vir_nolistsend_Nn_notenssends.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o totfield.o
torque3_dEtensor_72015_taunew3_pt1_vir_notenssends.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o
torque3_dEtensor_72015_taunew3_pt2_vir_notenssends.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o
thg5_small_vir_nolistsend_Nn_noewtens_notenssends.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o totfield.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3_notenssends.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o
mpidat3.o:
dEtensor2.o:
torque3_dEtensor_72015_taunew3_pt1_vir.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o
torque3_dEtensor_72015_taunew3_pt2_vir.o: atoms.o mpole.o sizes.o dEtensor.o dEtensor2.o
nblist_short.o: atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh2.o potent.o sizes.o vdw.o sizes2.o neigh.o
limits2.o:
neigh2.o:
switch2.o: limits.o nonpol.o shunt.o sizes.o
sizes2.o: 
Newempole1c_3bPolar_totfield_dEtensor.o: sizes2.o dEtensor2.o sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o
thg5_small_vir_nolistsend_Nn_noewtens.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o sizes2.o dEtensor2.o  
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir3.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o dEtensor2.o neigh2.o sizes2.o mpidat3.o
cutoffs_ewald_short.o: sizes2.o neigh2.o atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o sizes.o tarray.o neigh3b.o cobar.o totfield.o aprx.o neigh2clust.o
ewpol_loadbal7.o: mpidat.o neigh3b.o molcul.o
dEtensor_test_nolistsend_separatereduce_funcdecomp_NN.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o mpidat2.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o
nblist2_hbond.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
nblist2_test.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust2com.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o
clustnolistsmoothInner_1c_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o
nblist4.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
nblist3.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
nblist4_com.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
nblist3_com.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2_no3bmut.o: atoms.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2.o: atoms.o
grad_covemrecip_vandr_reduce_totfield.o: atoms.o
grad_vandr_ireduce_totfield_dEtensor_vir2.o: atoms.o
grad_vandr_reduce_totfield_dEtensor_vir2.o: atoms.o
grad_vandr_ireduce_totfield_dEtensor_vir.o: atoms.o
dEtensor_test_nolistsend_separatereduce_funcdecomp_no3bmut.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o mpidat2.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2_no3bmut.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o mpidat2.o
moments_onlyDipole_fixedcharge.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o charge.o dipole.o moment.o rigid.o
runmoments.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o sizes.o solute.o units.o usage.o deriv3b.o
uinddyn_rMPIOMP_kBcast_rloadbal_nocut2.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
uinddyn_rMPIOMP_kBcast_rloadbal_nocut.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
uindInner_1a_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
pmonte_par2.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o sizes.o units.o usage.o virial.o
empole0a_3bPolar_origversion_nopriordirectcalc2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
eng_r_reduce.o: mpole.o mpidat.o deriv.o energi.o virial.o
loadbalsmoothInner_0a_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
smooth_engpolar.o: mpole.o mpidat.o deriv3b.o
empole_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o mpidat.o paremneigh.o
eng_covalent_vdw.o: bound.o group.o sizes.o atoms.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o vdwpot.o virial.o deriv3b.o mpole.o mpidat.o couple.o molcul.o usage.o boxes.o polar.o polgrp.o neigh3b.o parewreg.o
verlet_par_pmonte.o: sizes.o atomid.o atoms.o freeze.o moldyn.o units.o usage.o
pmonte_par.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o sizes.o units.o usage.o virial.o
dyn_rMPIOMP_kBcast_rloadbal_pmonte_par.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o polar.o
uindempole1a_3bPolar_origversion_nopriordirectcalc2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
uindloadbalsmoothInner_1a_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
uindsmooth_gradpolar.o: mpole.o mpidat.o deriv3b.o
uinddyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
clust_smooth_gradpolar.o: mpole.o mpidat.o deriv3b.o neigh2clust.o aprx.o deriv1bmat.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp_clust2.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o neigh2clust.o deriv1bmat.o
clustnolistsmoothInner_1a_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o
clustloadbalsmoothInner_1a_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o neigh2clust.o
nblist2.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
nblist5.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o neigh2clust.o
neigh2clust.o:
pol_loadbal7.o: mpidat.o neigh3b.o molcul.o
pol_loadbal6.o: mpidat.o neigh3b.o molcul.o
thg5_small_vir_nolistsend_noewtens.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
dyn_rMPIOMP_kBcast_rloadbal_Namoeba2.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o polar.o virialnamoeba.o
virialnamoeba.o: 
dfield0a_ufield0a_uscale0a.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o
dyn_rMPIOMP_kBcast_rloadbal_Namoeba.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o polar.o virialnamoeba.o
empole1a_PolarStripped.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o virialnamoeba.o
induce0a_nobranch.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o
torque_pol.o: atoms.o deriv.o mpole.o sizes.o
dyn_rMPIOMP_kBcast_rloadbal_nolistsend_funcdecomp.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o mpidat2.o mpidat3.o
ewpol_loadbal_halftasks.o: mpidat.o neigh3b.o molcul.o
dEtensor_test_nolistsend_separatereduce_funcdecomp_simultpolz.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o moldyn_nopolz_vel.o
loadbalsmoothInner_1a_3bPolar1bmut_halftasks.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
loadbalsmoothInner_1a_3bPolar_halftasks.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
kinetic_nopolz_vel.o: atomid.o atoms.o bath.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o usage.o moldyn_nopolz_vel.o
hoover_nopolz_vel.o: sizes.o atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o units.o usage.o virial.o mpidat.o boxes.o cell.o moldyn_nopolz_vel.o
nose_par_nopolz_vel.o: sizes.o atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o units.o usage.o virial.o mpidat.o boxes.o cell.o moldyn_nopolz_vel.o
dEtensor_test_nolistsend_separatereduce_funcdecomp_noEast_nopolz_vel.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o moldyn_nopolz_vel.o
moldyn_nopolz_vel.o:
loadbalsmoothInner_1a_3bPolar_vac.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o pcg.o cho.o
loadbalsmoothInner_1a_3bPolar1bmut.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o 
dEtensor_test_nolistsend_separatereduce_funcdecomp_noEast.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
LoadBalInner_r_Perm_nolistsend.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o shunt.o
grad_covemrecip_vandr_reduce_totfield.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o mpidat2.o
empole1a_3bPolar_origversion_nopriordirectcalc_rtrnfield.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
induce0a_3b_PolelecOnly_totfieldnpolevac.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a2_3b_Polar_totfieldnpolevac.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_3bPolar_new_subtract_1bodymut.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_3bPolar_totfieldnpolevac.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_3bPolar_origversion_nopriordirectcalc2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o cho.o pcg.o
empole1c_3bPolar_totfield_901.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o
emrecip1_3b_Polar_usemlist_mutual_virialonly.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
emrecip1_3b_Polar_usemlist_mutual.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
totfieldnpole3b_noewald_umutual_rl_3b_nofield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a1_3bPolar_totfield_dEtensor1bmut.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o molcul.o aprx.o deriv3b.o
emrecip1_3b_Polar_usemlist_direct_virialonly.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o 
emrecip1_3b_Polar_usemlist_direct_test2.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
emrecip1_3b_Polar_usemlist_direct_test.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
grad_covemrecip_vandr_reduce_totfield_dEtensor_vir2.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o mpidat2.o
dEtensor_test_nolistsend_separatereduce_funcdecomp.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o mpidat2.o
grad_r_allreduce_totfield_dEtensor_vir.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o
AllreducedEtensorewtotfield_dyn_rMPIOMP_kBcast_rloadbal_vMPIOMP.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o bndpot.o bndstr.o angbnd.o angpot.o kurybr.o urey.o
bcast_cov.o: mpidat.o sizes.o atoms.o bndpot.o bndstr.o angbnd.o angpot.o kurybr.o urey.o limits.o potent.o strbnd.o urypot.o
thg6_small_vir_nolistsend.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
grad_vandr_ireduce_totfield_dEtensor_vir.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o
grad_vandr_ireduce_totfield_dEtensor_vir2.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o
loadbalsmoothInner_1a_2bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o
nblist3b_1list_par_newesttink7.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o mpidat.o
chkpole_omp.o: atoms.o mpole.o sizes.o
rotpole_omp.o: atoms.o mpole.o sizes.o
dEtensor_test_nolistsend_separatereduce.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
grad_vandr_reduce_totfield_dEtensor_vir2.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o
nblist_par.o: atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o mpidat.o
loadbal_ehal1c_halftasks_nolistsend.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o parvdwneigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o mpidat.o
ereal1d_Perm_savetensor.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
dEtensor_test_nolistsend.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
kewald_eastman_list.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o chunklight.o
chunklight.o: atoms.o
pmestuff_list.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o chunklight.o
v_loadbal_halftasks.o: sizes.o mpidat.o neigh.o parvdwneigh.o vdw.o
loadbal_ehal1c_halftasks.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o parvdwneigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o mpidat.o
grad_vandr_reduce_totfield_dEtensor_vir.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o atoms.o
dEtensor_test.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o 
empole1c_3bPolar_totfield_dEtensor1b_nouselist_mpivir.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o molcul.o aprx.o deriv3b.o 
torque3_dEtensor_72015_taunew2_pt1_vir.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_taunew2_pt2_vir.o: atoms.o mpole.o sizes.o
thg5_vir.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
thg5_small_vir.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
grad_r_reduce_totfield_dEtensor_vir.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o neigh.o
thg5.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o 
thg6_save.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
thg6_save_small.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
torque2_omp.o: sizes.o atoms.o mpole.o
table_fill_setchunk_omp.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
bspline_fill_bsplgen_omp.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
torque3_dEtensor_72015_taunew2_pt1.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_taunew2_pt2.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_taunew_pt1.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_taunew_pt2.o: atoms.o mpole.o sizes.o
LoadBalInner_thg4omp_mpi_72015torque2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
thg6.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
grad_r_reduce_thg5.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o
thg5.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
thg5_small.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
emrecip1_3b_Perm2_totfield_mpi_gridmpole.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o totfield.o mpidat.o
emrecip1_3b_Perm2_totfield_mpi_fphimpole.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o totfield.o mpidat.o
dEtensorewtotfield_dyn_rMPIOMP_kMPIgridmpoleBcast_rloadbal_vMPIOMP.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
dEtensorewtotfield_dyn_rMPIOMP_kMPIfphimpoleBcast_rloadbal_vMPIOMP.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
torque3_dEtensor_72015_tau_pt1.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_tau_pt2.o: atoms.o mpole.o sizes.o
nose_tuckerman06_par.o: sizes.o atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o units.o usage.o virial.o mpidat.o boxes.o cell.o tuckerman06.o
tuckerman06.o: 
torque3_dEtensor_72015_dE2only_pt2.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_dE1only_pt1.o: atoms.o mpole.o sizes.o
LoadBalInner_thg4omp_mpi_72015torque.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
torque3_dEtensor_72015_pt1.o: atoms.o mpole.o sizes.o
torque3_dEtensor_72015_pt2.o: atoms.o mpole.o sizes.o
grid_mpole_Mpi3amoebaRepo_allreduce.o: sizes.o atoms.o chunks.o mpole.o pme.o mpidat.o
emrecip1_3b_Perm2_totfield_mpi2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o totfield.o mpidat.o
fphi_mpole_Mpi3amoebaRepo.o: sizes.o mpole.o pme.o mpidat.o
grid_mpole_Mpi3amoebaRepo.o: sizes.o atoms.o chunks.o mpole.o pme.o mpidat.o
gradient_emrecip_totfield_mpi.o: sizes.o atoms.o energi.o deriv.o virial.o
emrecip1_3b_Perm2_totfield_mpi.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o totfield.o mpidat.o
fphi_mpole_Mpi3amoebaRepo.o: sizes.o mpole.o pme.o mpidat.o
dEtensorewtotfield_dyn_rMPIOMP_kMPIBcast_rloadbal_vMPIOMP.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
dEtensorewtotfield_dyn_rMPIOMP_kBcast_rloadbal_vMPIOMP.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o vdw.o vdwpot.o mutant.o
kewald_nopmealloc_eastman.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
mech_pme4_eastman.o: inform.o iounit.o limits.o potent.o vdwpot.o sizes.o atoms.o
kewald_eastman.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
grad_v_reduce_master2_send.o: mpidat.o deriv.o energi.o virial.o vdw.o
allocvdw_master2.o: sizes.o atoms.o bound.o deriv.o energi.o iounit.o limits.o potent.o virial.o
v_loadbal.o: sizes.o mpidat.o neigh.o parvdwneigh.o vdw.o
dir1bewtotfield_dyn_rMPIOMP_kBcast_rloadbal_vMPIOMPm2.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o vdw.o vdwpot.o mutant.o
parvdwneigh.o:
loadbal_ehal1c.o:  atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o parvdwneigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o mpidat.o 
gradient_covalent2.o: sizes.o atoms.o bound.o deriv.o energi.o iounit.o limits.o potent.o virial.o
grad_v_reduce.o: mpidat.o deriv.o energi.o virial.o vdw.o
dir1bewtotfield_dyn_rMPIOMP_kBcast_rloadbal_vMPIOMP.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o vdw.o vdwpot.o mutant.o 
field_r_reduce_uindcalc.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o
udirect2b_uinddir_mpi.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o chunks.o fft.o totfield.o paremneigh.o mpidat.o
LoadBalInnerloop_ereal1d_sentlist2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
udirect2b_uinddir.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o chunks.o fft.o totfield.o
empole1c_3bPolar_totfield_dEtensor1b_nouselist_mpi.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o molcul.o aprx.o deriv3b.o
get_ntpair_dE_task_nouselist.o: mpidat.o dEtensor.o neigh.o molcul.o mpole.o
empole1c_3bPolar_totfield_dEtensor1b_nouselist.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o molcul.o aprx.o deriv3b.o
get_ntpair_dE_task.o: mpidat.o dEtensor.o neigh.o molcul.o mpole.o
empole1c_3bPolar_totfield_dEtensor1b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o molcul.o
torque3_dEtensor_pt1.o: atoms.o mpole.o sizes.o
torque3_dEtensor_pt2.o: atoms.o mpole.o sizes.o
torque3_dEtensorMod2.o: atoms.o mpole.o sizes.o 
dEtensorewtotfield_dyn_rOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o 
totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensorMpi.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
totfield_smooth_gradpolar_dEtensor.o: mpole.o mpidat.o deriv3b.o totfield.o molcul.o aprx.o mpidat2.o
LoadBalInner_thg4omp_mpi.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o
LoadBalInner_thg4omp.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o openmp.o energi.o totfield.o 
grad_r_reduce_totfield_dEtensor.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o dEtensor.o aprx.o deriv3b.o chgpot.o
dEtensorewtotfield_dyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o aprx.o
grad_r_serial_totfield_dEtensor.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o
ewtotfieldsmooth2_gradient_polar_serial_dEtensor.o: mpole.o mpidat.o deriv3b.o totfield.o aprx.o
dEtensorewtotfield_dyn_serial.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o dEtensor.o
totfieldNpoleloadbalsmoothInner_1c_3bPolar_dEtensor.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
dEtensor.o:
LoadBalInner_thg4omp.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o dEtensor.o shunt.o
empole1c_3bPolar_totfield_dEtensor.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o dEtensor.o
empole1c_3bPolar_totfield_usemlist.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o neigh.o
emrecip1_3b_Polar_usemlist_direct.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o 
empole1d_3b_PolarPerm_selfeng_bcast.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
empole1d_3b_PolarPerm_selfeng_bcast2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
grad_r_reduce_totfield_uindcalc.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o polar.o
dir1bewtotfield_dyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o
ereal1d2_3b_totfieldnpole_mlist.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o neigh.o shunt.o mpidat.o paremneigh.o
r_loadbal2.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
ereal1c2_3b_totfieldnpole_omp.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o neigh.o shunt.o
totfieldNpolesmoothInner_1c_1bPolar2_513_alliter.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
totfieldNpolesmoothInner_1c_1bPolar2_513_multiter.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
totfieldNpolesmoothInner_1c_1bPolar2_513.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o aprx.o
totfieldNpoleloadbalsmoothInner_1c_3bPolar2_513.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o
ewtotfield_dyn_serial.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o boxes.o cell.o limits.o
ewtotfieldsmallsmooth2_gradient_polar_serial.o: mpole.o mpidat.o deriv3b.o totfield.o aprx.o
grad_r_serial_totfield.o: mpole.o mpidat.o deriv.o energi.o virial.o totfield.o math.o ewald.o
sendmlist_load_balanced3_sizelist_serial.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o 
pmestuff_3bpolz_513.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
induce_totfield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o ewald.o fft.o
ewtotfield_dyn_serial.o: sizes.o atoms.o bath.o bndstr.o bound.o boxes.o cell.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o limits.o
ewtotfieldsmallsmooth2_gradient_polar_serial.o: mpole.o mpidat.o deriv3b.o totfield.o aprx.o
grad_r_serial_totfield.o: math.o mpole.o mpidat.o deriv.o energi.o ewald.o virial.o totfield.o
sendmlist_load_balanced3_sizelist_serial.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o 
pmestuff_3bpolz_513.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
induce_totfield.o: sizes.o bound.o boxes.o ewald.o fft.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o uprior.o
pmestuff_3bpolz_totfield_513.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmestuff_3bpolz_totfield_pt2_513.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o 
emrecip1_1and2_3b_Polar_totfield_parl_513.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
empole1c_3bPolar_totfield_513.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o
kewald3b_totfield.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
empole1c_3bPolar_totfield_513offset.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o 
totfieldNpoleloadbalsmoothInner_1c_2bPolar2_513.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o
bspline_table_recipsum.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o pme.o chunks.o atoms.o fft.o 
empole1c_3bPolar_totfield_511offset.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o
ereal1d2_3b_totfieldnpole_511offset.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o neigh.o shunt.o 
empole1c_3bPolar_totfield_511.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o
totfieldNpoleloadbalsmoothInner_1c_2bPolar2_511.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o limits.o
aprx.o:
totfieldNpoleloadbalsmoothInner_1c_2bPolar2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
totfieldNpoleloadbalsmoothInner_1c_2bPolar.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
nblist3b_1list.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o
alloc3b_vdwlist_1list.o: sizes.o atoms.o deriv3b.o limits.o aprx.o
ewtotfield_dyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o pme.o chunks.o aprx.o
totfieldNpolesmoothInner_1c_3bPolar.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
totfieldNpoleloadbalsmoothInner_1c_3bPolar.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
ereal1d2_3b_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o neigh.o shunt.o 
r_loadbal2.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
LoadBalInner_r_Perm_totfield_bcastlist.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o
grad_r_reduce_totfield_bcast.o: mpole.o mpidat.o deriv.o energi.o virial.o totfield.o
bcastmpolelistparalleltotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
paralleltotfieldloadbalsmooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
emrecip1_1and2_3b_Polar_totfield_parl.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o totfield.o
paralleltotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
gradient_emrecip_totfield.o: sizes.o atoms.o energi.o deriv.o virial.o
gradient_emrecip_totfield_send.o: sizes.o atoms.o energi.o deriv.o virial.o
emrecip1_3b_Perm2_totfield.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o totfield.o
LoadBalInner_r_Perm_totfield.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o totfield.o shunt.o
grad_r_reduce_totfield.o: mpole.o mpidat.o deriv.o energi.o virial.o totfield.o
emrecip1_1and2_3b_Polar_totfield.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
pmestuff_3bpolz_totfield_pt2.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
emrecip1_1_3b_Polar_totfield.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
emrecip1_1_3b_Polar_totfield_parl_513.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o chunks.o
emrecip1_2_3b_Polar_totfield.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o
pmestuff_3bpolz_totfield.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
erecip1_1mod_3b_totfieldnpole3b.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
erecip1_2mod_3b_totfieldnpole3b.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
udirect1_reg.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
erecip1_3b_totfieldnpole3b.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o 
erecip1_1_3b_totfieldnpole3b.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
erecip1_2_3b_totfieldnpole3b.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
erecip1_1_3b_totfieldnpole.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
erecip1_2_3b_totfieldnpole.o: sizes.o atoms.o chgpot.o mpole.o polar.o boxes.o ewald.o ewreg3bpolz.o math.o
stotfield_dyn_rMPIOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o limits.o
totfield_dyn_rMPIOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o limits.o aprx.o
totfieldloadbalsmooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
stotfieldloadbalsmooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
s1btotfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
s1btotfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
s1btotfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
sAll2bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
s2btotfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
s2btotfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
s2btotfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
print3Npoletotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
Printtotfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
tot1bNpolefield_dyn_rMPIOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
printloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
totfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
dfield0c_totfield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o totfield.o pme.o fft.o shunt.o
ereal1c1_3b_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o
ereal1c2_3b_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o
ereal1c1_3b_totfieldnpole3b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o
ereal1c2_3b_totfieldnpole3b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o 
empole1c_3bPolar_totfield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o totfield.o chunks.o
All2bAll3bNew2smooth_innerloop_empole1c_3b_totfield.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
All2bAll3bNew2smooth_innerloop_empole1c_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
pmestuff_3bpolz_pt2.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o 
sAll2bAll3bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
stotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o totfield.o
stotfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
stotfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
stotfield3NpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
tot1bfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
tot1bNpolefieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
All2bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
tot1b3Npolefieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
totfield1body_w1body_rtrndeptemp2.o:  sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfieldNpole1body_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfield3Npole1body_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfield1b_smooth_gradpolar.o: mpole.o mpidat.o deriv3b.o
totfield3Npole2body_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
empole1a_3bPolar_totfieldnpole3.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o 
empole1a1_3b_Polar_totfieldnpole_mod3.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole_mod3.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
BigMscratch2.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
tot2b3Npolefieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
tot2bNpolefieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
tot2bfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
totfield2b_smooth_gradpolar.o: mpole.o mpidat.o deriv3b.o
totfieldNpole2body_w1body.o:  sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfield2body_w1body_rtrndeptemp2.o:  sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
printtotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o 
printNpoletotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
PrinttotfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
printtotfield_smooth_gradpolar.o: mpole.o mpidat.o deriv3b.o 
empole1a2_3b_Polar_totfieldnpole_mod.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a1_3b_Polar_totfieldnpole_mod.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
totfieldloadbalsmoothInner_1a_3bPolar_rtrndeptemp2.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
empole1a2_3b_Polar_totfieldnpole3bmod_dol1_2.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a_3bPolar_totfieldnpole3b_rtrndeptemp2.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
totfieldloadbalsmoothInner_1a_3bPolar_w1body_rtrndeptemp2.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
empole1a2_3b_Polar_totfieldnpole3bmod_doall.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
totfieldloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfieldNpoleloadbalsmoothInner_1a_3bPolar_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o 
totfieldloadbalsmoothInner_1a_2bPolar_w1body.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
totfieldNpole4body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
induce0a_3b_scf_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o 
empole1a_mod.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o inter.o mplpot.o
totfieldNpole2body_w1body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfieldNpole2body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfield4body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
Scratch4body.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
totfieldloadbalsmoothInner_1a_2bPolar.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
dfield0a_totfield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
induce0a_3b_PolelecOnly_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a1_3b_Polar_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a_3bPolar_totfieldnpole.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
totfieldNpoleloadbalsmoothInner_1a_3bPolar.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole3bmod2_test.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole3bmod_nohalf.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole3bmod_dol1.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
udirtotfieldloadbalsmoothInner_1a_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
empole1a_3bPolar_subtrdir.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
induce0a_3b_PolelecOnly_subtrdir.o: sizes.o atoms.o polar.o
empole1b_dir.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o totfield.o
udirtotfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
torque3_Polar.o: atoms.o deriv.o mpole.o sizes.o
torque_3b_new_npole.o: atoms.o deriv.o mpole.o sizes.o
empole1a1_3b_Polar_totfieldnpole3bmod.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole3bmod2.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
dfield0a_totfieldnpole3bmod.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole3bmod.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
dfield0a_totfieldnpole3b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a1_3b_Polar_totfieldnpole3b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a2_3b_Polar_totfieldnpole3b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o
empole1a_3bPolar_totfieldnpole3b.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o totfield.o
totfield.o:
dfield0b_totfield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o totfield.o
empole1a_3bPolar_totfield.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o totfield.o
totfieldloadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o cell.o limits.o
totfield_smooth_gradpolar.o: mpole.o mpidat.o deriv3b.o totfield.o aprx.o
totfieldloadbalsmoothInner_1a_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
pme3bsmalldyn_rOMP_kBcast.o:  sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
nosmoothpme3bsmalldyn_rOMP_kBcast.o:  sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
fft3d_3b.o: fft.o openmp.o pme.o sizes.o
kewald3b.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
pmestuff_3bpolz.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
emrecip1_3b_Polar.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o fft.o
ew3bsmalldyn_rOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
loadbalsmoothInner_1c_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o 
nosmoothew3bsmalldyn_rOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
loadbalNOsmoothInner_1c_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
ewsmooth_gradpolar_no3lst.o: mpole.o mpidat.o deriv3b.o
testvirpolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cell.o
nocut_ewsmooth_gradpolar.o: mpole.o mpidat.o deriv3b.o
allpolz3lst_innerloop_empole1c_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
nocutew3b3lstsmalldyn_rOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
ewreg3bpolz.o: sizes.o
induce0c_field_amat.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o totfield.o
empole1c_3bPolar.o: sizes.o bound.o boxes.o ewald.o math.o mpole.o ewreg3bpolz.o atoms.o couple.o group.o polgrp.o polpot.o polar.o ewald.o cell.o limits.o chgpot.o fft.o
kewaldreg3b.o: limits.o sizes.o math.o boxes.o ewald.o ewreg3bpolz.o 
ewsmooth_gradpolar.o: mpole.o mpidat.o deriv3b.o 
3lst_innerloop_empole1c_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
ewpol_loadbal.o: mpidat.o neigh3b.o molcul.o
ew3b3lstsmalldyn_rOMP_kBcast.o:  sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
pol_loadbal3b.o: mpidat.o neigh3b.o molcul.o
3List.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o cobar.o
3lstsmalldyn_rOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
3lst_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
uind_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
cobar.o:
mdstat_propertiesDaveDeTo.o: sizes.o atoms.o bath.o bound.o boxes.o inform.o inter.o iounit.o limits.o math.o mdstuf.o molcul.o moment.o units.o usage.o warp.o
UindAll2bAll3bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
All2bAll3bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
Scratch_All2bAll3bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
PrintAll2bAll3bNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
UindloadbalsmoothInner_1a_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
Uindempole1a_3bPolar.o: sizes.o atoms.o cell.o chgpot.o couple.o group.o mpole.o polar.o polgrp.o polpot.o shunt.o usage.o
anintempole1a_3bPolar.o: sizes.o atoms.o cell.o chgpot.o couple.o group.o mpole.o polar.o polgrp.o polpot.o shunt.o usage.o
Uindsmooth_gradpolar.o: mpole.o mpidat.o deriv3b.o
uindpolzloadbalsmoothdynamic_3b_parallel_simult_master_too2_2_emrealOMPOnly_emrecipMPIbcst.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
printsmalldyn3b_emrealOMPOnly_emrecipBcast_allpolz.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smalldyn3b_emrealOMPOnly_emrecipBcast_allpolz.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o limits.o
alloc3b_uind.o: sizes.o atoms.o deriv3b.o
moments_onlyDipole.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o sizes.o solute.o units.o usage.o deriv3b.o
get_ntask_r.o: keys.o mpidat.o
New_r_loadbal.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o 
r_loadbal.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
smalldyn_rOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o tuckerman06.o
polz2bloadbalsmoothdynamic_3b_parallel_simult_master_too2_2_emrealOMPOnly_emrecipMPIbcst.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
sz2elstloadbalsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPIbcst64.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
noifdyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
dyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smalldyn_rMPIOMP_kBcast_rloadbal.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
emreal_load_balance.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
pol_loadbal.o: mpidat.o neigh3b.o molcul.o
paremneigh.o:
sendmlist.o: sizes.o mpidat.o neigh.o paremneigh.o 
sendmlist_load_balanced.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
sendmlist_load_balanced2.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
bcastmlist_load_balanced.o: sizes.o mpidat.o neigh.o paremneigh.o mpole.o
LoadBalInner_r_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o
LoadBalInnerloop_ereal1d_3b_Perm_fulllist2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o mpidat.o
gradient_covalent_vdwTest_novdw.o: bound.o group.o sizes.o atoms.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o vdwpot.o virial.o deriv3b.o mpole.o mpidat.o couple.o molcul.o usage.o boxes.o polar.o polgrp.o neigh3b.o
bcast_vdw.o: mpidat.o sizes.o atoms.o vdw.o vdwpot.o limits.o potent.o mutant.o
mech_par.o: atoms.o inform.o iounit.o limits.o potent.o vdwpot.o molcul.o
mechanic_pme.o: inform.o iounit.o limits.o potent.o vdwpot.o sizes.o atoms.o
mech_pme3.o: inform.o iounit.o limits.o potent.o vdwpot.o sizes.o atoms.o
mech_pme4.o: inform.o iounit.o limits.o potent.o vdwpot.o sizes.o atoms.o
mechanic_ewaldreg.o: inform.o iounit.o limits.o potent.o vdwpot.o sizes.o atoms.o
mechanic_ewaldreglist.o: inform.o iounit.o limits.o potent.o vdwpot.o sizes.o atoms.o
PermElec.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
PermElec2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o neigh2clust.o
emrecip1_3b_Perm2NoQPole.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
emrecip1_3b_Perm2NoQPoleNoDPole.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
emrecip1PermReduce.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o ewreg.o math.o mpole.o polar.o polpot.o units.o virial.o
parallel_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o ewreg.o math.o mpole.o polar.o polpot.o units.o virial.o
parallel2_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o ewreg.o math.o mpole.o polar.o polpot.o units.o virial.o
parallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o ewreg.o math.o mpole.o polar.o polpot.o units.o virial.o
Newmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o
Newmpiparallel3_erecip1_3b_Perm_master.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o
CondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o
NoQList2CondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o ewaldneigh.o
List2CondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o ewaldneigh.o
SendList2CondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o ewaldneigh.o
ListCondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o ewaldneigh.o
RmndrListCondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o ewaldneigh.o
CondensedLoopNewmpiparallel3_erecip1_3b_Perm_master.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o
NoQPCondensedLoopNewmpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o
NoQPCondensedLoopNewmpiparallel3_erecip1_3b_Perm_master.o: sizes.o atoms.o boxes.o chgpot.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o mpidat.o
mpiparallel3_erecip1_3b_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o virial.o mpidat.o
mpiparallel3_erecip1_3b_Perm_master1.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o virial.o mpidat.o
mpiparallel3_erecip1_3b_Perm_master2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o parewreg.o math.o mpole.o polar.o polpot.o units.o virial.o mpidat.o
Innerloop_ereal1d_3b_Perm_fulllist2.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
Inner_r_Perm.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o paremneigh.o
innerloop_pairwise_additive_fulllist_ereal1d_only.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
innerloop_pairwise_additive_fulllist_single_ereal1d_only.o: sizes.o atoms.o boxes.o chgpot.o deriv.o energi.o ewald.o inter.o math.o mpole.o polar.o polpot.o virial.o bound.o pme.o potent.o couple.o neigh.o limits.o mplpot.o potent.o
nblist3b.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o
par_nblist_2b.o: sizes.o atoms.o bound.o boxes.o iounit.o neigh3b.o molcul.o atomid.o neigh.o light.o cell.o
innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
New2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
anintloadbalsmoothInner_1a_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o
loadbalsmoothInner_1a_3bPolar.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o cobar.o
b2test.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o
LoadBalcom_smooth2b_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o mpidat.o
AllatomNew2smooth_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
innerloop_empole1a_NoCut2b3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
innerloop_empole1a_NoCut2bonly.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
innerloop_empole1a_2bonly.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
innerloop_empole1a_3bonly.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
printR123_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
err_nosqrt_printL3_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
nosqrt_printL3_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
serial_printL3_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
printL3_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
L3ongTaperR123smooth3b2b_2_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
NewL3ongTaperR123smooth3b2b_2_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
NoOMPL3ongTaperR123smooth3b2b_2_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
com_smooth2b_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
nolistcom_smooth2b_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
com_smooth2b3pt5Angst_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
com_NoSmooth2b_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
com_NOSmooth2b_innerloop_empole1a_3b.o: sizes.o atoms.o atomid.o molcul.o mpole.o neigh3b.o
empole1a_3bPolar.o: sizes.o atoms.o chgpot.o couple.o group.o mpole.o polar.o polgrp.o polpot.o usage.o
gradient_3b.o: bound.o group.o sizes.o atoms.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o vdwpot.o virial.o deriv3b.o mpole.o mpidat.o couple.o molcul.o usage.o boxes.o polar.o polgrp.o neigh3b.o cobar.o limits.o totfield.o aprx.o cho.o pcg.o
smooth_gradpolar.o: mpole.o mpidat.o deriv3b.o aprx.o mpidat3.o
printL_gradient_polar.o: mpole.o mpidat.o deriv3b.o
smooth2b_gradient_polar.o: mpole.o mpidat.o deriv3b.o
grad_r_reduce.o: mpole.o mpidat.o deriv.o energi.o virial.o
grad_k_reduce.o: mpole.o mpidat.o deriv.o energi.o virial.o
Modgrad_3b.o: bound.o group.o sizes.o atoms.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o vdwpot.o virial.o deriv3b.o mpole.o mpidat.o couple.o molcul.o usage.o boxes.o polar.o polgrp.o neigh3b.o parewreg.o
test_grad.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o
test_grad2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o
dynamic_3b_parallel.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o 
just_polz_wlist_wbcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
just_serial2_wlist_wbcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
Newjust_pairwiseadd.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o freeze.o
Newjust_vdw_cov.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o freeze.o
just_serial2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
just_polz.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
just_list.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
just_covalent.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
just_bcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
just_setup.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
dynamic_3b_parallel_simult.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
dynamic_3b_parallel_simult_master_too.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o
dynamic_3b_parallel_simult_master_too2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
dynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
Hopperdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o
loadbalsmooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o usolve.o
loadbalsmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o
nosmooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o
noewald_nosmooth_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o 
noewald_nosmooth_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o
smooth2b_virial_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o boxes.o
printL_findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smooth_findiff_test_nompipolz.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
findiff_test.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smooth2bdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smooth2b3pt5Angstdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPI.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
loadbalsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPI.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
loadbalsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPbcast_emrecipMPI.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smallpolz2bloadbalsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPIbcst64.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
dyn_rMPIOMP_kReduce.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smalldyn_rMPIOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
dyn_rMPIOMP_kBcast.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_2b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPI.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPIbcst.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPIsend.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend_emrecipMPI_notbetterpar.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
ewaldregsmoothdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
parewaldregsmoothdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o
par2ewaldregsmoothdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o
listpar2ewaldregsmoothdynamic_3b_parallel_simult_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
listpar2ewaldregsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMP.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
sendlistpar2ewaldregsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMP.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
sendlistpar2ewaldregsmoothdynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMPsend.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
just_ew_cov_vdw.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
just_polz_cov_vdw.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
just_ereal_cov_vdw.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o parewreg.o ewaldneigh.o
parallel_gradient_erecip.o: sizes.o atoms.o mpole.o mpidat.o deriv.o energi.o virial.o
listparallel_gradient_erecip.o: sizes.o atoms.o mpole.o mpidat.o deriv.o energi.o virial.o ewaldneigh.o parewreg.o
ewaldneigh.o: 
smoothdynamic_testcov.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
just_covvdwpme_wlistwbcast_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
serial_smoothdynamic_testcov.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_testcov_bcastlattice.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_testcov_bcastlattice_emrealmaster1.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_testcov_bcastlattice_emrealmaster1_covvdwonly.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_testcov_bcastlattice_emrealmaster1_covvdwonlymaster1.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_testcov_bcastlattice_emrealmaster1_covvdwonlymaster1_emrealmaster.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
smoothdynamic_testcov_bcastlattice_emrealmaster1_covvdwonlymaster1_emrealmaster1.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o vdw.o
smoothdynamic_testcov_bcastlattice_mlistonlymaster1.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
dynamic_3b_parallel_simult_master_too2_2_emreal_MPI_OMP.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
just_setup_master_too2_2.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
just_setup_master_too2_3.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o
dynamic_3b_parallel_simult_master_too3.o: sizes.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o solute.o stodyn.o usage.o mpidat.o moldyn.o virial.o deriv3b.o molcul.o neigh3b.o mpole.o energi.o polar.o polgrp.o couple.o group.o rgddyn.o files.o atomid.o units.o deriv.o neigh.o pme.o chunks.o
nose_par.o: sizes.o atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o units.o usage.o virial.o mpidat.o boxes.o cell.o
verlet_par.o: sizes.o atomid.o atoms.o freeze.o moldyn.o units.o usage.o
beeman_pt1_pt2_master.o: sizes.o atomid.o atoms.o freeze.o mdstuf.o moldyn.o units.o usage.o
mdinit_par.o: atomid.o atoms.o bath.o bound.o files.o freeze.o group.o inform.o iounit.o keys.o mdstuf.o molcul.o moldyn.o mpole.o rgddyn.o rigid.o sizes.o stodyn.o units.o uprior.o usage.o
gradslow_serial.o: limits.o potent.o
action.o:
active.o: atoms.o inform.o iounit.o keys.o sizes.o usage.o
alchemy.o: analyz.o atoms.o energi.o files.o inform.o iounit.o katoms.o mutant.o potent.o sizes.o units.o usage.o
align.o:
analysis.o: analyz.o atoms.o bound.o energi.o group.o inter.o iounit.o limits.o potent.o sizes.o vdwpot.o
analyz.o:
analyze.o: action.o analyz.o angang.o angbnd.o angpot.o angtor.o atomid.o atoms.o bitor.o bndstr.o bound.o boxes.o charge.o chgpot.o couple.o dipole.o energi.o ewald.o fields.o files.o improp.o imptor.o inform.o inter.o iounit.o korbs.o ktrtor.o kvdws.o limits.o math.o molcul.o moment.o mpole.o opbend.o opdist.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o potent.o sizes.o solute.o strbnd.o strtor.o tors.o tortor.o units.o urey.o vdw.o vdwpot.o virial.o
angang.o:
angbnd.o:
angles.o: angbnd.o atmlst.o atoms.o couple.o iounit.o sizes.o
angpot.o:
angtor.o:
anneal.o: atomid.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o mdstuf.o potent.o sizes.o solute.o usage.o warp.o
archive.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o usage.o
argue.o:
ascii.o:
atmlst.o:
atomid.o: sizes.o
atoms.o: sizes.o
attach.o: atoms.o couple.o iounit.o sizes.o
bar.o: atoms.o energi.o files.o inform.o iounit.o keys.o sizes.o units.o
basefile.o: ascii.o files.o
bath.o:
beeman.o: atomid.o atoms.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o
bicubic.o:
bitor.o:
bitors.o: angbnd.o atoms.o bitor.o couple.o iounit.o sizes.o
bndpot.o:
bndstr.o:
bonds.o: atmlst.o atoms.o bndstr.o couple.o iounit.o sizes.o
born.o: atomid.o atoms.o bath.o chgpot.o couple.o deriv.o inform.o iounit.o math.o pbstuf.o sizes.o solute.o virial.o
bound.o:
bounds.o: atomid.o atoms.o boxes.o molcul.o sizes.o
boxes.o:
bussi.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o
calendar.o:
cell.o: sizes.o
center.o: align.o sizes.o
charge.o:
chgpot.o:
chkpole.o: atoms.o mpole.o sizes.o
chkring.o: couple.o sizes.o
chkxyz.o: atoms.o iounit.o sizes.o
cholesky.o:
chrono.o:
chunks.o:
clock.o: chrono.o
cluster.o: atomid.o atoms.o bound.o group.o inform.o iounit.o keys.o limits.o molcul.o sizes.o
column.o: sizes.o
command.o: argue.o
connect.o: atoms.o couple.o sizes.o zclose.o zcoord.o
connolly.o: atoms.o faces.o inform.o iounit.o math.o sizes.o
control.o: argue.o inform.o keys.o output.o sizes.o
correlate.o: ascii.o atomid.o atoms.o files.o inform.o iounit.o sizes.o
couple.o: sizes.o
crystal.o: atomid.o atoms.o bound.o boxes.o couple.o files.o iounit.o math.o molcul.o sizes.o
cspline.o: iounit.o
cutoffs.o: atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o sizes.o tarray.o neigh3b.o cobar.o totfield.o aprx.o neigh2clust.o cho.o pcg.o
cutoffs_vdw_nolist.o: atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o sizes.o tarray.o neigh3b.o
deflate.o: iounit.o
delete.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
deriv.o:
diagq.o:
diffeq.o: atoms.o iounit.o math.o sizes.o warp.o
diffuse.o: atomid.o atoms.o bound.o inform.o iounit.o molcul.o sizes.o usage.o
dipole.o:
disgeo.o:
distgeom.o: angbnd.o atomid.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o kvdws.o math.o refer.o restrn.o sizes.o tors.o
dma.o:
document.o: iounit.o
domega.o:
dynamic.o: atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o sizes.o solute.o stodyn.o usage.o
eangang.o: angang.o angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangang1.o: angang.o angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangang2.o: angang.o angbnd.o angpot.o atoms.o bound.o deriv.o group.o hessn.o math.o sizes.o
eangang3.o: action.o analyz.o angang.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangle.o: angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangle1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangle2.o: angbnd.o angpot.o atoms.o bound.o deriv.o group.o hessn.o math.o sizes.o
eangle3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangtor.o: angbnd.o angtor.o atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o
eangtor1.o: angbnd.o angtor.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o
eangtor2.o: angbnd.o angtor.o atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o
eangtor3.o: action.o analyz.o angbnd.o angtor.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o
ebond.o: atoms.o bndpot.o bndstr.o bound.o energi.o group.o sizes.o usage.o
ebond1.o: atoms.o bndpot.o bndstr.o bound.o deriv.o energi.o group.o sizes.o usage.o virial.o
ebond2.o: atmlst.o atoms.o bndpot.o bndstr.o bound.o couple.o group.o hessn.o sizes.o
ebond3.o: action.o analyz.o atomid.o atoms.o bndpot.o bndstr.o bound.o energi.o group.o inform.o iounit.o sizes.o usage.o
ebuck.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o iounit.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
ebuck1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
ebuck2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o iounit.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
ebuck3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
echarge.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o iounit.o light.o limits.o math.o neigh.o pme.o shunt.o sizes.o usage.o warp.o
echarge1.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o light.o limits.o math.o molcul.o neigh.o pme.o shunt.o sizes.o usage.o virial.o warp.o
echarge2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o ewald.o group.o hessn.o limits.o math.o shunt.o sizes.o warp.o
echarge3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o warp.o
echgdpl.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
echgdpl1.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
echgdpl2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o group.o hessn.o shunt.o sizes.o units.o
echgdpl3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
edipole.o: atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
edipole1.o: atoms.o bound.o cell.o chgpot.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
edipole2.o: atoms.o bound.o cell.o chgpot.o dipole.o group.o hessn.o shunt.o sizes.o units.o
edipole3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
egauss.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egauss1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
egauss2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o warp.o
egauss3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egeom.o: atomid.o atoms.o bound.o energi.o group.o math.o molcul.o restrn.o sizes.o usage.o
egeom1.o: atomid.o atoms.o bound.o deriv.o energi.o group.o inter.o math.o molcul.o restrn.o sizes.o usage.o virial.o
egeom2.o: atomid.o atoms.o bound.o deriv.o group.o hessn.o math.o molcul.o restrn.o sizes.o
egeom3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o inter.o iounit.o math.o molcul.o restrn.o sizes.o usage.o
ehal.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
ehal1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
ehal1c_noswitch_bcast.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
ehal2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
ehal3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
eimprop.o: atoms.o bound.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o
eimprop1.o: atoms.o bound.o deriv.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o virial.o
eimprop2.o: atoms.o bound.o group.o hessn.o improp.o math.o sizes.o torpot.o
eimprop3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o improp.o inform.o iounit.o math.o sizes.o torpot.o usage.o
eimptor.o: atoms.o bound.o energi.o group.o imptor.o sizes.o torpot.o usage.o
eimptor1.o: atoms.o bound.o deriv.o energi.o group.o imptor.o sizes.o torpot.o usage.o virial.o
eimptor2.o: atoms.o bound.o group.o hessn.o imptor.o sizes.o torpot.o
eimptor3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o imptor.o inform.o iounit.o math.o sizes.o torpot.o usage.o
elj.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
elj1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
elj2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
elj3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
embed.o: angbnd.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o keys.o light.o math.o minima.o output.o refer.o restrn.o sizes.o tors.o units.o
emetal.o: atomid.o atoms.o couple.o energi.o kchrge.o sizes.o
emetal1.o: atomid.o atoms.o couple.o deriv.o energi.o kchrge.o sizes.o
emetal2.o:
emetal3.o: action.o analyz.o atomid.o atoms.o energi.o kchrge.o sizes.o
emm3hb.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o light.o limits.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
emm3hb1.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
emm3hb2.o: atmlst.o atomid.o atoms.o bndstr.o bound.o cell.o chgpot.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
emm3hb3.o: action.o analyz.o atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
empole.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o group.o limits.o math.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
empole1.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole1a_Perm.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
torque_Perm.o: atoms.o deriv.o mpole.o sizes.o 
empole2.o: atoms.o bound.o boxes.o chgpot.o couple.o deriv.o group.o hessn.o limits.o molcul.o mplpot.o mpole.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
empole3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
energi.o:
energy.o: bound.o energi.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o
eopbend.o: angbnd.o angpot.o atoms.o bound.o energi.o fields.o group.o math.o opbend.o sizes.o usage.o
eopbend1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o opbend.o sizes.o usage.o virial.o
eopbend2.o: angbnd.o angpot.o atoms.o bound.o deriv.o group.o hessn.o math.o opbend.o sizes.o
eopbend3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o opbend.o sizes.o usage.o
eopdist.o: angpot.o atoms.o bound.o energi.o group.o opdist.o sizes.o usage.o
eopdist1.o: angpot.o atoms.o bound.o deriv.o energi.o group.o opdist.o sizes.o usage.o virial.o
eopdist2.o: angpot.o atoms.o bound.o group.o hessn.o opdist.o sizes.o usage.o
eopdist3.o: action.o analyz.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o opdist.o sizes.o usage.o
epitors.o: atoms.o bound.o energi.o group.o pitors.o sizes.o torpot.o usage.o
epitors1.o: atoms.o bound.o deriv.o energi.o group.o pitors.o sizes.o torpot.o usage.o virial.o
epitors2.o: angbnd.o atoms.o bound.o deriv.o group.o hessn.o pitors.o sizes.o torpot.o usage.o
epitors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o pitors.o sizes.o torpot.o usage.o
erf.o: iounit.o math.o
erxnfld.o: atoms.o chgpot.o energi.o mpole.o rxnfld.o rxnpot.o shunt.o sizes.o usage.o
erxnfld1.o: atoms.o deriv.o energi.o sizes.o
erxnfld2.o:
erxnfld3.o: action.o analyz.o atomid.o atoms.o chgpot.o energi.o inform.o iounit.o mpole.o shunt.o sizes.o usage.o
esolv.o: atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o kvdws.o limits.o math.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
esolv1.o: atomid.o atoms.o bound.o boxes.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o virial.o warp.o
esolv2.o: atoms.o charge.o chgpot.o hessn.o math.o potent.o shunt.o sizes.o solute.o warp.o
esolv3.o: action.o analyz.o atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inform.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
estrbnd.o: angbnd.o angpot.o atoms.o bndstr.o bound.o energi.o group.o math.o sizes.o strbnd.o usage.o
estrbnd1.o: angbnd.o angpot.o atoms.o bndstr.o bound.o deriv.o energi.o group.o math.o sizes.o strbnd.o usage.o virial.o
estrbnd2.o: angbnd.o angpot.o atoms.o bndstr.o bound.o group.o hessn.o math.o sizes.o strbnd.o
estrbnd3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strbnd.o usage.o
estrtor.o: atoms.o bndstr.o bound.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o
estrtor1.o: atoms.o bndstr.o bound.o deriv.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o virial.o
estrtor2.o: atoms.o bndstr.o bound.o group.o hessn.o sizes.o strtor.o torpot.o tors.o
estrtor3.o: action.o analyz.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strtor.o torpot.o tors.o usage.o
etors.o: atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o warp.o
etors1.o: atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o warp.o
etors2.o: atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o warp.o
etors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o warp.o
etortor.o: atomid.o atoms.o bitor.o bound.o couple.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
etortor1.o: atoms.o bitor.o bound.o deriv.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o virial.o
etortor2.o: atoms.o bitor.o bound.o group.o hessn.o ktrtor.o math.o sizes.o torpot.o tortor.o units.o
etortor3.o: action.o analyz.o atoms.o bitor.o bound.o energi.o group.o inform.o iounit.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
eurey.o: atoms.o bound.o energi.o group.o sizes.o urey.o urypot.o usage.o
eurey1.o: atoms.o bound.o deriv.o energi.o group.o sizes.o urey.o urypot.o usage.o virial.o
eurey2.o: atoms.o bound.o couple.o group.o hessn.o sizes.o urey.o urypot.o
eurey3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o sizes.o urey.o urypot.o usage.o
evcorr.o: atoms.o bound.o boxes.o limits.o math.o mutant.o shunt.o sizes.o vdw.o vdwpot.o
ewald.o:
ewreg.o: sizes.o
parewreg.o: sizes.o
extra.o: energi.o
extra1.o: atoms.o deriv.o energi.o sizes.o
extra2.o: atoms.o hessn.o sizes.o
extra3.o: action.o analyz.o atoms.o energi.o sizes.o
faces.o:
fatal.o: iounit.o
fft.o: sizes.o
fft3d.o: fft.o openmp.o pme.o sizes.o
fftpack.o: math.o
field.o: keys.o potent.o sizes.o
fields.o: sizes.o
files.o:
final.o: align.o analyz.o angang.o angbnd.o angtor.o atmlst.o bitor.o bndstr.o charge.o chunks.o couple.o deriv.o dipole.o disgeo.o domega.o faces.o fracs.o freeze.o group.o hessn.o hpmf.o improp.o imptor.o inform.o iounit.o light.o merck.o molcul.o moldyn.o mpole.o mutant.o neigh.o nonpol.o opbend.o opdist.o orbits.o paths.o pdb.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o qmstuf.o refer.o restrn.o rgddyn.o rigid.o ring.o sizes.o socket.o solute.o stodyn.o strbnd.o strtor.o syntrn.o tarray.o tors.o tortor.o uprior.o urey.o usage.o usolve.o vdw.o vibs.o warp.o ewaldneigh.o totfield.o dEtensor.o chunklight.o deriv1bmat.o neigh2clust.o deriv2bmat.o uind1bmatrix.o uind2bmat.o
flatten.o: atoms.o fields.o inform.o iounit.o keys.o sizes.o warp.o
fracs.o:
freeunit.o: iounit.o
freeze.o:
gda.o: atoms.o files.o iounit.o minima.o potent.o sizes.o vdwpot.o warp.o
geometry.o: atoms.o math.o sizes.o
getint.o: atoms.o inform.o iounit.o output.o sizes.o
getkey.o: argue.o files.o iounit.o keys.o openmp.o
getmol.o: files.o iounit.o
getmol2.o: files.o iounit.o
getnumb.o: ascii.o
getpdb.o: iounit.o
getprm.o: files.o iounit.o keys.o params.o
getref.o: atomid.o atoms.o couple.o files.o refer.o sizes.o titles.o
getstring.o: ascii.o
gettext.o: ascii.o
getword.o: ascii.o
getxyz.o: inform.o iounit.o output.o
ghmcstep.o: atomid.o atoms.o bath.o freeze.o iounit.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
gkstuf.o: sizes.o
gradient.o: atoms.o bound.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o virial.o
gradrgd.o: atoms.o group.o rigid.o sizes.o
gradrot.o: atoms.o deriv.o domega.o omega.o potent.o rotbnd.o sizes.o
group.o:
groups.o: group.o sizes.o
grpline.o: atomid.o atoms.o group.o rgddyn.o sizes.o
gyrate.o: atoms.o sizes.o usage.o
hescut.o:
hessian.o: atoms.o bound.o couple.o hescut.o hessn.o inform.o iounit.o limits.o mpole.o potent.o rigid.o sizes.o usage.o vdw.o vdwpot.o
hessn.o:
hessrgd.o: atoms.o group.o rigid.o sizes.o
hessrot.o: math.o omega.o sizes.o zcoord.o
hpmf.o:
hybrid.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o charge.o couple.o dipole.o imptor.o inform.o iounit.o kangs.o katoms.o kbonds.o kchrge.o kdipol.o kitors.o kstbnd.o ksttor.o ktorsn.o kvdws.o math.o mutant.o sizes.o strbnd.o strtor.o tors.o vdw.o vdwpot.o
image.o: boxes.o cell.o sizes.o
impose.o: align.o inform.o iounit.o sizes.o
improp.o:
imptor.o:
induce.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o chunks.o fft.o totfield.o
inertia.o: atomid.o atoms.o iounit.o math.o sizes.o
inform.o:
initatom.o: ptable.o
initial.o: align.o atoms.o bath.o bound.o boxes.o cell.o files.o group.o inform.o iounit.o keys.o linmin.o minima.o molcul.o mutant.o neigh.o openmp.o output.o params.o pdb.o precis.o rigid.o scales.o sequen.o sizes.o socket.o warp.o zclose.o
initprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o math.o merck.o mplpot.o polpot.o rxnpot.o sizes.o solute.o torpot.o units.o urypot.o vdwpot.o
initres.o: resdue.o
initrot.o: atoms.o couple.o group.o inform.o iounit.o math.o omega.o potent.o restrn.o rotbnd.o sizes.o usage.o zcoord.o
insert.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
intedit.o: atomid.o atoms.o files.o iounit.o katoms.o sizes.o zcoord.o
inter.o:
intxyz.o: files.o iounit.o titles.o
invbeta.o:
invert.o: iounit.o
iounit.o:
jacobi.o: iounit.o
kanang.o: sizes.o
kangang.o: angang.o angbnd.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kanang.o keys.o potent.o sizes.o tors.o
kangle.o: angbnd.o angpot.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kangs.o keys.o merck.o potent.o ring.o sizes.o usage.o
kangs.o:
kangtor.o: angtor.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kantor.o keys.o potent.o sizes.o tors.o
kantor.o:
katom.o: atomid.o atoms.o couple.o inform.o iounit.o katoms.o keys.o sizes.o
katoms.o: sizes.o
kbond.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kbonds.o keys.o merck.o potent.o sizes.o tors.o usage.o
kbonds.o:
kcharge.o: atomid.o atoms.o charge.o chgpot.o couple.o fields.o inform.o iounit.o kchrge.o keys.o merck.o potent.o sizes.o
kchrge.o: sizes.o
kdipol.o:
kdipole.o: atmlst.o atoms.o bndstr.o couple.o dipole.o inform.o iounit.o kdipol.o keys.o potent.o sizes.o
kewald.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
kewaldreg.o: limits.o math.o boxes.o ewald.o ewaldneigh.o mpidat.o parewreg.o 
kewald_nopmealloc.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
send_ewaldreglist.o: limits.o math.o boxes.o ewald.o ewaldneigh.o mpidat.o parewreg.o
sendnomaster_ewaldreglist.o: limits.o math.o boxes.o ewald.o ewaldneigh.o mpidat.o parewreg.o
mpidat.o: 
kextra.o:
keys.o:
kgeom.o: atomid.o atoms.o bound.o couple.o group.o iounit.o keys.o molcul.o potent.o restrn.o sizes.o
khbond.o:
kimprop.o: atomid.o atoms.o couple.o improp.o inform.o iounit.o keys.o kiprop.o potent.o sizes.o tors.o
kimptor.o: atomid.o atoms.o couple.o imptor.o inform.o iounit.o keys.o kitors.o math.o potent.o sizes.o tors.o
kinetic.o: atomid.o atoms.o bath.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o usage.o
kiprop.o:
kitors.o:
kmetal.o:
kmpole.o: atoms.o couple.o inform.o iounit.o keys.o kmulti.o mpole.o polar.o polgrp.o potent.o sizes.o units.o
kmulti.o:
kopbend.o: angbnd.o angpot.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kopbnd.o merck.o opbend.o potent.o sizes.o usage.o
kopbnd.o:
kopdist.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o kopdst.o opdist.o potent.o sizes.o
kopdst.o:
korbit.o: atomid.o atoms.o bndstr.o inform.o iounit.o keys.o korbs.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
korbs.o: sizes.o
kpitor.o:
kpitors.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o keys.o kpitor.o pitors.o potent.o sizes.o tors.o
kpolar.o: atoms.o couple.o inform.o iounit.o keys.o kpolr.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o usolve.o
kpolr.o: sizes.o
ksolv.o: angbnd.o atmlst.o atomid.o atoms.o bath.o bndstr.o chgpot.o couple.o gkstuf.o hpmf.o inform.o iounit.o keys.o kvdws.o math.o nonpol.o pbstuf.o potent.o sizes.o solute.o
kstbnd.o:
kstrbnd.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kstbnd.o merck.o potent.o ring.o sizes.o strbnd.o
kstrtor.o: atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o ksttor.o potent.o sizes.o strtor.o tors.o
ksttor.o:
ktors.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o ktorsn.o math.o merck.o potent.o ring.o sizes.o tors.o usage.o
ktorsn.o:
ktortor.o: atomid.o atoms.o bitor.o inform.o iounit.o keys.o ktrtor.o potent.o sizes.o tortor.o
ktrtor.o:
kurey.o: angbnd.o atomid.o atoms.o inform.o iounit.o keys.o kurybr.o potent.o sizes.o urey.o
kurybr.o:
kvdw.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o khbond.o kvdwpr.o kvdws.o math.o merck.o potent.o sizes.o vdw.o vdwpot.o
kvdwpr.o:
kvdws.o: sizes.o
lattice.o: boxes.o cell.o inform.o iounit.o math.o sizes.o
lbfgs.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o scales.o
light.o:
lights.o: bound.o boxes.o cell.o iounit.o light.o sizes.o
limits.o:
linmin.o:
makeint.o: atoms.o couple.o inform.o iounit.o math.o sizes.o zclose.o zcoord.o
makeref.o: atomid.o atoms.o couple.o files.o refer.o sizes.o titles.o
makexyz.o: atoms.o sizes.o zcoord.o
math.o:
maxwell.o: units.o
mdinit.o: atomid.o atoms.o bath.o bound.o files.o freeze.o group.o inform.o iounit.o keys.o mdstuf.o molcul.o moldyn.o mpole.o rgddyn.o rigid.o sizes.o stodyn.o units.o uprior.o usage.o
mdrest.o: atomid.o atoms.o bound.o group.o inform.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o
mdsave.o: atomid.o atoms.o bound.o boxes.o files.o group.o inform.o iounit.o mdstuf.o moldyn.o mpole.o output.o polar.o potent.o rgddyn.o sizes.o socket.o titles.o units.o
mdstat.o: atoms.o bath.o bound.o boxes.o inform.o inter.o iounit.o limits.o mdstuf.o molcul.o sizes.o units.o usage.o warp.o
mdstuf.o:
mechanic.o: inform.o iounit.o limits.o potent.o vdwpot.o
merck.o: sizes.o
merge.o: atomid.o atoms.o couple.o iounit.o refer.o sizes.o
minima.o:
minimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
minirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
minrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o sizes.o
molcul.o:
moldyn.o:
molecule.o: atomid.o atoms.o couple.o molcul.o sizes.o
molxyz.o: files.o iounit.o titles.o
moment.o:
moments.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o sizes.o solute.o units.o usage.o
monte.o: atoms.o files.o inform.o iounit.o omega.o output.o sizes.o units.o usage.o zcoord.o
mplpot.o:
mpole.o:
mutant.o:
mutate.o: atomid.o atoms.o charge.o inform.o iounit.o katoms.o keys.o mpole.o mutant.o polar.o potent.o sizes.o
nblist.o: atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
neigh.o:
newton.o: atoms.o files.o inform.o iounit.o keys.o sizes.o usage.o
newtrot.o: files.o hescut.o inform.o iounit.o keys.o math.o omega.o sizes.o zcoord.o
nextarg.o: argue.o
nexttext.o:
nonpol.o:
nose.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o virial.o
nspline.o:
nucleic.o: atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o nucleo.o output.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
nucleo.o: sizes.o
number.o: inform.o iounit.o
numeral.o:
numgrad.o: atoms.o sizes.o
ocvm.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o potent.o scales.o sizes.o
omega.o: sizes.o
opbend.o:
opdist.o:
openend.o:
openmp.o:
optimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
optirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
optrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o sizes.o
optsave.o: atoms.o files.o iounit.o math.o omega.o output.o scales.o sizes.o socket.o usage.o zcoord.o
orbital.o: atomid.o atoms.o bndstr.o couple.o iounit.o keys.o piorbs.o potent.o sizes.o tors.o
orbits.o:
orient.o: atomid.o atoms.o group.o math.o rigid.o sizes.o
orthog.o:
output.o:
overlap.o: units.o
params.o:
path.o: align.o atomid.o atoms.o files.o inform.o iounit.o linmin.o minima.o output.o paths.o sizes.o
paths.o:
pbstuf.o: sizes.o
pdb.o:
pdbxyz.o: atomid.o atoms.o couple.o fields.o files.o inform.o iounit.o katoms.o pdb.o resdue.o sequen.o sizes.o titles.o
phipsi.o: sizes.o
picalc.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
piorbs.o:
pistuf.o:
pitors.o:
pme.o: sizes.o
pmestuff.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmestuffNoQPole.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmestuffNoQPoleNoDPole.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmpb.o: iounit.o
polar.o:
polarize.o: atoms.o inform.o iounit.o molcul.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
poledit.o: atomid.o atoms.o couple.o dma.o files.o iounit.o keys.o kpolr.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
polgrp.o:
polpot.o:
polymer.o: atoms.o bndstr.o bound.o boxes.o iounit.o keys.o sizes.o
potent.o:
potential.o: atomid.o atoms.o bndstr.o charge.o chgpot.o couple.o dipole.o files.o inform.o iounit.o katoms.o kchrge.o kdipol.o keys.o kmulti.o kpolr.o math.o minima.o moment.o mpole.o neigh.o output.o polar.o potent.o potfit.o refer.o sizes.o titles.o units.o
potfit.o: sizes.o
precis.o:
precise.o:
pressure.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o sizes.o units.o usage.o virial.o
ptest_parallel.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o sizes.o units.o usage.o virial.o
prmedit.o: angpot.o bndpot.o iounit.o math.o params.o sizes.o urypot.o vdwpot.o
prmkey.o: angpot.o bndpot.o chgpot.o fields.o mplpot.o polpot.o potent.o rxnpot.o sizes.o torpot.o urypot.o vdwpot.o
promo.o: iounit.o
protein.o: atomid.o atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o output.o phipsi.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
prtdyn.o: atoms.o boxes.o files.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o titles.o
prterr.o: files.o output.o
prtint.o: atomid.o atoms.o files.o inform.o sizes.o titles.o zclose.o zcoord.o
prtmol2.o: atomid.o atoms.o bndstr.o couple.o files.o iounit.o sizes.o titles.o
prtpdb.o: files.o pdb.o sequen.o sizes.o titles.o
prtprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o mplpot.o polpot.o sizes.o urypot.o vdwpot.o
prtseq.o: files.o sequen.o sizes.o
prtxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o sizes.o titles.o
pss.o: atoms.o files.o hescut.o inform.o iounit.o math.o omega.o refer.o sizes.o tree.o warp.o zcoord.o
pssrigid.o: atoms.o files.o group.o inform.o iounit.o math.o minima.o molcul.o refer.o rigid.o sizes.o warp.o
pssrot.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o refer.o sizes.o warp.o zcoord.o
ptable.o:
qmstuf.o:
quatfit.o: align.o sizes.o
radial.o: argue.o atomid.o atoms.o bound.o boxes.o files.o inform.o iounit.o limits.o math.o molcul.o potent.o sizes.o
random.o: inform.o iounit.o keys.o math.o sizes.o
rattle.o: atomid.o atoms.o freeze.o group.o inform.o iounit.o moldyn.o sizes.o units.o usage.o virial.o
readdyn.o: atoms.o boxes.o files.o group.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o
readgau.o: ascii.o iounit.o qmstuf.o sizes.o units.o
readint.o: atomid.o atoms.o files.o inform.o iounit.o sizes.o titles.o zclose.o zcoord.o
readmol.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o sizes.o titles.o
readmol2.o: atomid.o atoms.o couple.o files.o iounit.o sizes.o titles.o
readpdb.o: files.o inform.o iounit.o pdb.o resdue.o sequen.o titles.o
readprm.o: fields.o iounit.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o merck.o params.o sizes.o
readseq.o: files.o iounit.o resdue.o sequen.o
readxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o
refer.o: sizes.o
replica.o: bound.o boxes.o cell.o inform.o iounit.o sizes.o
resdue.o:
respa.o: atomid.o atoms.o freeze.o limits.o moldyn.o potent.o sizes.o units.o usage.o virial.o
restrn.o:
rgddyn.o:
rgdstep.o: atomid.o atoms.o bound.o group.o iounit.o rgddyn.o sizes.o units.o virial.o
rigid.o:
ring.o:
rings.o: angbnd.o atoms.o bitor.o bndstr.o couple.o inform.o iounit.o ring.o tors.o
rmsfit.o: align.o sizes.o
rotbnd.o: sizes.o
rotlist.o: atoms.o couple.o iounit.o molcul.o rotbnd.o sizes.o zclose.o
rotpole.o: atoms.o mpole.o sizes.o
rxnfld.o:
rxnpot.o:
saddle.o: atoms.o inform.o iounit.o keys.o linmin.o minima.o sizes.o syntrn.o titles.o zcoord.o
scales.o: sizes.o
scan.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o output.o sizes.o zcoord.o
sdstep.o: atomid.o atoms.o bath.o couple.o freeze.o kvdws.o math.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
search.o: linmin.o math.o sizes.o
sequen.o: sizes.o
server.o:
shakeup.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o bound.o couple.o freeze.o keys.o math.o molcul.o ring.o sizes.o usage.o
shunt.o:
sigmoid.o:
sizes.o:
sizes.o:
sktstuff.o: atomid.o atoms.o charge.o couple.o deriv.o fields.o files.o inform.o iounit.o keys.o moldyn.o mpole.o polar.o potent.o sizes.o socket.o
sniffer.o: atoms.o files.o inform.o iounit.o linmin.o math.o minima.o output.o scales.o sizes.o usage.o
socket.o: sizes.o
solute.o:
sort.o:
spacefill.o: atomid.o atoms.o files.o inform.o iounit.o kvdws.o math.o sizes.o usage.o
spectrum.o: files.o iounit.o math.o units.o
square.o: inform.o iounit.o keys.o minima.o sizes.o
stodyn.o:
strbnd.o:
strtor.o:
suffix.o: ascii.o
superpose.o: align.o atomid.o atoms.o bound.o files.o inform.o iounit.o sizes.o titles.o
surface.o: atoms.o inform.o iounit.o math.o sizes.o usage.o
surfatom.o: atoms.o iounit.o math.o sizes.o
switch.o: limits.o nonpol.o shunt.o sizes.o
sybylxyz.o: files.o iounit.o titles.o
syntrn.o:
tarray.o:
temper.o: atomid.o atoms.o bath.o group.o mdstuf.o molcul.o moldyn.o rgddyn.o sizes.o units.o usage.o
testgrad.o: atoms.o deriv.o energi.o files.o inform.o inter.o iounit.o sizes.o solute.o usage.o
testhess.o: atoms.o files.o hescut.o inform.o iounit.o sizes.o usage.o
testpair.o: atoms.o deriv.o energi.o inform.o iounit.o light.o limits.o neigh.o polpot.o potent.o sizes.o tarray.o vdwpot.o
testpol.o: atoms.o bound.o inform.o iounit.o limits.o polar.o polpot.o potent.o rigid.o sizes.o units.o
testrot.o: domega.o energi.o inform.o iounit.o math.o omega.o sizes.o zcoord.o
timer.o: atoms.o hescut.o inform.o iounit.o limits.o sizes.o
timerot.o: iounit.o limits.o omega.o sizes.o
titles.o:
tncg.o: atoms.o hescut.o inform.o iounit.o keys.o linmin.o math.o minima.o output.o piorbs.o potent.o sizes.o
torphase.o:
torpot.o:
torque.o: atoms.o deriv.o mpole.o sizes.o
tors.o:
torsfit.o: atomid.o atoms.o files.o inform.o iounit.o keys.o ktorsn.o math.o output.o potent.o qmstuf.o restrn.o scales.o sizes.o tors.o usage.o
torsions.o: atoms.o bndstr.o couple.o iounit.o sizes.o tors.o
tortor.o:
tree.o:
trimtext.o:
unitcell.o: bound.o boxes.o iounit.o keys.o sizes.o
units.o:
uprior.o:
urey.o:
urypot.o:
usage.o:
usolve.o:
valence.o: angbnd.o angpot.o atomid.o atoms.o bndpot.o bndstr.o couple.o files.o hescut.o inform.o iounit.o kangs.o kbonds.o keys.o kopbnd.o kstbnd.o ktorsn.o kurybr.o kvdws.o linmin.o math.o minima.o opbend.o output.o potent.o qmstuf.o scales.o sizes.o strbnd.o torpot.o tors.o units.o urey.o urypot.o usage.o valfit.o vdwpot.o
valfit.o:
vdw.o:
vdwpot.o:
verlet.o: atomid.o atoms.o freeze.o moldyn.o sizes.o units.o usage.o
version.o: iounit.o output.o
vibbig.o: atomid.o atoms.o bound.o couple.o files.o hescut.o hessn.o inform.o iounit.o keys.o limits.o mpole.o potent.o rigid.o sizes.o units.o usage.o vdw.o vdwpot.o vibs.o
vibrate.o: atomid.o atoms.o files.o hescut.o iounit.o math.o sizes.o units.o usage.o
vibrot.o: iounit.o omega.o sizes.o
vibs.o:
virial.o:
volume.o: atoms.o iounit.o math.o sizes.o
warp.o:
xtalfit.o: atomid.o atoms.o bound.o boxes.o charge.o dipole.o files.o fracs.o iounit.o kvdws.o math.o molcul.o potent.o sizes.o vdw.o xtals.o
xtalmin.o: atoms.o boxes.o files.o inform.o iounit.o keys.o math.o scales.o sizes.o
xtals.o:
xyzatm.o: atoms.o inform.o iounit.o math.o sizes.o
xyzedit.o: atomid.o atoms.o bound.o boxes.o couple.o files.o iounit.o limits.o math.o molcul.o refer.o sizes.o titles.o units.o usage.o
xyzint.o: files.o iounit.o titles.o
deriv3b.o:
xyzpdb.o: atomid.o atoms.o couple.o fields.o files.o inform.o molcul.o pdb.o resdue.o sequen.o sizes.o
xyzsybyl.o: files.o iounit.o sizes.o titles.o
zatom.o: angbnd.o atomid.o atoms.o bndstr.o fields.o iounit.o kangs.o katoms.o kbonds.o sizes.o zclose.o zcoord.o
neigh3b.o:
zclose.o: sizes.o
zcoord.o: sizes.o
