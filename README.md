# 3-AMOEBA 

## Table of Contents

* [Literature References] (#literature-references)
* [Code Layout](#code-layout)
* [Theory and implementation](#theory-and-implementation)
* [Compiling from source](#compiling-from-source)
* [Running the code](#running-the-code)

## Literature References

Account of the convergence (accuracy) of the many-body expansion of classical polarization for
the energy, force, and virial.  Speedups of the parallel implementation for the 
various versions of the 3-body approximation under different levels of electrostatic
embedding (Referred to in Reference 2 as 3m- and 3M-AMOEBA.):

1. Demerdash O, Head-Gordon T. Convergence of the many-body expansion for energy and forces
for classical polarizable models in the condensed phase. Journal of Chemical Theory and
Computation, 2016, 12(8): 3884-3893.

Account of the parallel implementations of the 3m- and 3M-AMOEBA models, where 3m- and 3M-
are, respectively, the non-embedded and the 'quasi'-embedded versions:

2. Demerdash O, Head-Gordon T. Parallel implementation of approximate atomistic models of the
AMOEBA polarizable model. Chemical Physics Letters, 2016, 664: 191-198.

## Code Layout

This is the `master` branch of the 3-AMOEBA implementation in Tinker
7. The directories in this repository are:

* `Src` - where the main source distribution is located.
* `Example` - example for the code.
* `thirdparty` - contains a local version of `fftw` used by the main code.


## Theory and implementation

The 3-body approximation of the polarization energy/gradient/virial
has now been implemented in the molecular dynamics routine (for the
Verlet and Nose-Hoover methods) using a hybrid MPI/OpenMP scheme for
both the 3m- and 3M-AMOEBA models. The motivation for this work is the
desire for a computationally cheaper, yet accurate, approximation for 
calculating the most expensive component of the AMOEBA energy/gradient/virial
calculation, which is the contribution due to many-body polarization,
the sole many-body, non-pairwise additive component of the AMOEBA
force field.  The approach towards a computationally efficient calculation of 
polarization pursued here rests on the many-body expansion, which states
that any many-body potential U (here only polarization) may be expanded as
         U = U(1) + Delta U(2) + Delta U(3) + ...

   where:
                  __
           U(1) = \
                  /  Ui
                  --
                   i
                       __ __
           deltaU(2) = \  \
                       /  /  Uij - Ui - Uj
                       -- --
                       i  j>i

                       __ __ __
           deltaU(3) = \  \  \
                       /  /  /   Uijk - Delta U(i,j) - Delta U(i,k) - Delta U(j,k) - Ui -Uj - Uk
                       -- -- --
                       i  j>i k>j
                       __ __ __
                     = \  \  \
                       /  /  /   Uijk - (Uij - Ui - Uj) - (Uik - Ui - Uk) - (Ujk - Uj - Uk) - Ui -Uj - Uk
                       -- -- --
                       i  j>i k>j

                       __ __ __
                     = \  \  \
                       /  /  /   Uijk - Uij - Uik - Ujk + Ui + Uj + Uk
                       -- -- --
                       i  j>i k>j

and where:

* `Ui` is the energy of the ith body or fragment, e.g., a single water molecule,collection of molecules, molecule fragment from a larger polymeric structure.
* `Uij` is the energy associated with the interaction of body/fragment `i` with body/fragment `j`, and
* `Uijk` is the energy associated with the interaction among bodies/fragments `i`, `j`, and `k`.
 
The utility of the truncated many-body expansion relies on the fact that 
successive terms are progressively smaller and may be neglected without sacrificing
accuracy.  Here, models truncated at the 2- and 3-body level have been implemented.  The
computational efficiency rests on the fact that solving for the induced dipoles 
of the smaller subsystems is computationally cheaper than the full N-atom
system, and that these subsystems may be calculated in parallel.  Moreover,
the decay of 2- and 3-body polarization with distance enables the use of cutoffs
and neighbor lists, reducing a potentially O(N^3) computation
to one that is O(N), where N here is the number of 3-body subsystems to be
computed.  

Models have been implemented using a hybrid MPI/OpenMP code using the 
technique of atom-decomposition, or replicated data, as spatial decomposition cannot
be realized due to the spatially diffuse nature of the 3-body interactions, whereas
with pairwise-additive potentials, a spatial decomposition is more straightforward
due to the localized nature of the interactions.  Despite the use of potentially
expensive collective MPI communication calls in this approach, we have seen during
development that these costs are always negligible compared to computation.  Moreover,
the cost of communication can be hidden in many cases through the use of non-blocking
collectives.  Moreover, efficient load balancing has been implemented, both in terms of a specific
interaction, such as distributing the number of polarization calculations evenly among tasks,
but also in terms of task parallelism, allowing different classes of interactions, e.g. vdW 
and permanent electrostatics, to be executed on separate tasks simultaneously. Due to Amdahls
Law, the next most expensive non-covalent interactions in AMOEBA, the multipole permanent 
electrostatics and the van der Waals, were distributed among multiple tasks in a hybrid MPI/OpenMP
implementation, as is the polarization. Specifics of the parallel implementation are 
described in Reference 2.  

Key to defining the accuracy of physical properties of interest, 
such as the radial distribution function of water, are the accuracy 
of the gradient of the polarization energy and the polarization 
contribution of the virial.  It is explained in Reference 1 that 
a seemingly rapid convergence of the many-body expansion for 
polarization energy does not correspond to rapid convergence in 
the accuracy of the force, whose RMS error is an order of magnitude
greater than that of the energy.  Enhancing the accuracy entails a 
better description of the permanent electrostatic environment
giving rise to many-body polarization in a procedure called 
electrostatic embedding, wherein higher-order effects are
captured implicitly.  A more computationally efficient form of 
embedding is 3M-AMOEBA where the "body" is a k-means cluster of 
molecules/fragments instead of a single molecule/fragment, which
is the simple 3m-AMOEBA model.  The 2M-AMOEBA model, using the same 
clustering approach as 3M-AMOEBA but with the 2-body approximation
instead, has been found to be accurate also. In practice, we have used
 k-means clustering, but in principle, other forms, such as hierarchical
clustering, are likely suitable.  


There are hybrid MPI/OpenMP versions of:

### I. Simple 3-body approximation of polarization without electrostatic embedding (3m-AMOEBA):
       The MPI/OpenMP version of Tinkers 'dynamic.x' that runs using the 3-body approximation of
       polarization is entitled 'dynamic_3mAMOEBA.x'. There
       are a few notes on the execution that are distinct from the standard OpenMP-only 
       'dynamic.x' (which runs only the standard, full N-body mutually polarizable AMOEBA model).
        
       Parallel execution: 1.Due to task parallelism,wherein all pairwise-additive interactions
                             are executed simultaneously on separate tasks, at least 4 MPI tasks 
                             (processes) must be requested.

                           2.An efficient parallelization of the pairwise-additive, real-space,
                             non-covalent interactions presumes an O(N) implementation.  Therefore,
                             this implementation necessitates that the user specify the use of
                             neighbor lists for the van der Waals and the real-space permanent
                             electrostatic interactions.  In the .key file, the following flags
                             must be used: 'mpole-list' and 'vdw-list'.

       
       Specific flags:     1. Cutoff distances for the 2- and 3-body interactions must be specified.
                              For 2-body interactions, these are distances between oxygen atoms for 
                              pairs of water molecules.  For the 3-body interaction, the inter-body
                              distance is defined as the sum of the 2 shortest distances out of the 3
                              sets of distances corresponding to the 3 pairs of water molecules in 
                              the trimer subsystem (the distance of the pair is defined identically
                              to that of the 2-body interaction).  In addition to the cutoff distances,
                              to conserve energy, the interactions must be smoothly turned off to zero
                              starting at a distance about 1-2 Angstroms less than the cutoff distance, the 
                              tapering distance. The flags for these distances should be in the .key file.
                              A typical example follows with values for cutoff and tapering distances that
                              yield the errors reported in Reference 1 for water interactions:

                                   # 3-body cutoff distance
                                   COBARCUT3B  9.0
                                   # 3-body tapering distance
                                   TAPER3B     7.0
                                   # 2-body cutoff distance
                                   CUT2B  6.0
                                   # 2-body tapering distance
                                   TAPER2B  5.0

                           2. Truncation at the 3-body level is the default mode for running
                              'dynamic_3mAMOEBA.x', and is the
                              level of truncation of the many-body expansion which ensures errors
                              in the range of a few percent.  It is possible to run at the level of
                              the 2-body approximation (2m-AMOEBA), but this is not recommended due
                              to the large error with respect to the full, N-body mutually 
                              polarizable model.  2m-AMOEBA can be invoked with the following flag:
                                   APRXMODE 2BODYMODE 
                              Without this flag, 3m-AMOEBA is run.                                             


### II. 2- and 3-body approximations under modified electrostatic embedding (2M- and 3M-AMOEBA):
        The 2- or 3-body approximations of polarization may be calculated under definitions of a body
        wherein the body is not a single molecule or fragment.  It has been found (See Ref. 1) that
        a model wherein the body is redefined as a cluster of water molecules yields better agreement
        of the polarization contribution to the energy, force, and virial, recovering much of the 
        computational accuracy afforded by full electrostatic embedding (See Reference 1.). 
        'dynamic_2M_AMOEBA.x' runs the 2-body
        approximation where the body is now a cluster of molecules/fragments. 
        'dynamic_3M_AMOEBA.x' runs the 3-body
        approximation.  In general, for bulk water systems of typical simulation size
        (~4000-1,000,000 atoms), a decimation of the system into 10-20 clusters achieves 
        reasonable accuracy in the polarization contribution to the force.


       Parallel execution: 1.Due to task parallelism,wherein all pairwise-additive interactions
                             are executed simultaneously on separate tasks, at least 4 MPI tasks
                             (processes) must be requested.  An additional requirement is that a 
                             number of tasks equal to the number of clusters must be specified should
                             this number be greater than 4, which is typically the case if a fast
                             model is desired.  
      

                           2.An efficient parallelization of the pairwise-additive, real-space,
                             non-covalent interactions presumes an O(N) implementation.  Therefore,
                             this implementation necessitates that the user specify the use of
                             neighbor lists for the van der Waals and the real-space permanent
                             electrostatic interactions.  In the .key file, the following flags
                             must be used: 'mpole-list' and 'vdw-list'.

       Specific flags:     1. Cutoff distances for the 2- and 3-body interactions must be specified.
                              For 2-body interactions, these are distances between the geometric 
                              centers of the clusters.  For the 3-body interaction, the inter-body
                              distance is defined as the sum of the 2 shortest distances out of the 3
                              sets of distances corresponding to the 3 pairs of clusters in
                              the trimer subsystem (the distance of the pair is defined identically
                              to that of the 2-body interaction).  Smoothening has not yet been
                              implemented for 2M- or 3M-AMOEBA.  For a decimation of a typical
                              simulation box, cutoffs in the following range are adequate:

                                   # 3-body cutoff distance
                                   COBARCUT3B  40.0 # A more aggressive cutoff for 3-body may be used for 3M-AMOEBA
                                   # 2-body cutoff distance
                                   CUT2B  40.0

                           2. 2M- and 3M-AMOEBA require an additional command line argument
                              which is a text file prescribing the clustering.  On the first
                              line of the file is the number of clusters, followed by the
                              size of the largest cluster.  Each subsequent line then is
                              simply the cluster to which each water molecule in the system
                              belongs:
                              20 1921 #20 clusters in the system, 1921 molecules in largest clust
                              20      #1st molecule in cluster index 20
                              13      #2nd moleculee in cluster index 13
                              ...
                              ...,etc.

                           3. Additional flags for running 2M- and 3M-AMOEBA are as follows:
                              CLUSTLIST          # Use neighbor lists for the clusters
                              APRXMODE 2BODYMODE # 2BODYMODE for 2M-AMOEBA and 3BODYMODE for 3M-AMOEBA
                              PCG                # Use Preconditioned conjugate gradient SCF solver
                              KMEANSCLUST        # Use k-means clusters to define bodies
                              NOOMPOUTER3B       # Dont call routines where OpenMP is at 'outer' loops
                              MLISTKMEANSCLUST   # Use neighborlist for molecules withing the subsystems
                              EWALDCLUST         # Use Ewald in the clustering method
                              POLAR-PREDICT     # Use values of induced dipoles of prior timesteps to form initial guess
 
### III. Tensor method for direct polarization with the iAMOEBA model (iAMOEBA-tensor).  In
         contrast to mutual polarization, one may define a polarization model wherein induced dipoles
         result only from the field generated by the permanent multipoles, avoiding the computationally
         expensive SCF.  Such a model must be reparameterized due to the missing contribution from 
         mutual induction, and has been done in the iAMOEBA water model.  
         Here, direct polarization has been implemented in a hybrid MPI/OpenMP manner and benefits from
         enhanced efficiency of real-space computations.  Specifically,the permanent field 
         and the gradient of the field are calculated and saved as a vector and tensor, respectively, 
         together with the permanent electrostatic contribution to the energy, gradient, and virial.
         Subsequently, the induced dipoles are calculated from the product of the polarizabilities and
         the saved field, the polarization energy is calculated from the inner product of the 
         induced dipoles and permanent field, and finally the polarization
         gradient and virial are calculated from the induced dipoles and saved field gradient.
         Use the executable entitled 'dynamic_iAMOEBA_tensor.x'

       Parallel execution: 1.Due to task parallelism,wherein all pairwise-additive interactions
                             are executed simultaneously on separate tasks, at least 4 MPI tasks
                             (processes) must be requested.


                           2.An efficient parallelization of the pairwise-additive, real-space,
                             non-covalent interactions presumes an O(N) implementation.  Therefore,
                             this implementation necessitates that the user specify the use of
                             neighbor lists for the van der Waals and the real-space permanent
                             electrostatic interactions.  In the .key file, the following flags
                             must be used: 'mpole-list' and 'vdw-list'.

       Specific flags:     1. In the .key file, 'polarization DIRECT' and 'UZE', to indicate that
                              an Particle Mesh Ewald contribution to direct polarization must be
                              calculated separately from the real-space component, calculated w/ tensor.                             
                              Additionally, since running a direct polarization-only model requires
                              reparameterization, a different parameter file must be specified. This is
                              not yet provided in the Example directory, but it will be obtained from
                              the authors and provided soon.

## Compiling from source 

The `Makefile` is currently modified for compiling with MPI on a
Linux platform (`F77 = mpif77`, where `module load gcc48` is
invoked prior to compiling using gfortran on our local cluster).  
The relevant make commands are as follows:

1. If the `fftw` library in `thirdparty/fftw` jhas not been compiled
   and installed you need to do this first. For instructions on doing
   this refer to the document in `thirdparty/fftw/0README`.

1. To remove all previously made object files (`.o`) and executables
(`.x`) files:

            make clean

1. To build all executables and any object files whose corresponding
source has been modified:

           make all

## Running the code

          Please refer to the Example directory 

Each of the `Example` directory contains a bash shell script that can
be customised to run the code for a specific platform.

### Inputs and Outputs

An explanation of the relevant input:

* `<KEYFILE>`: Specifies the relevant parameters. For this
               calculation,these are the periodic box size (keyword:
               `<a,b,c>-axis`), the energetic parameters (keyword:
               `parameters`), the specification to run mutual
               polarization (keyword: `polarization`), specification of
               Ewald electrostatics for the permanent, fixed-charge
               electrostatics (keyword: `ewald`), and the real-space
               cutoff for Ewald electrostatics (keyword:
               `ewald-cutoff`).

* `<COORDINATE_FILE>`: Cartesian atomic coordinates in Tinker `.xyz`
                       format.


