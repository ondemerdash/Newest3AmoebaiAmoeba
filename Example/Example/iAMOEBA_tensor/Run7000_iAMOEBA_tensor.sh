#!/bin/bash -l
#SBATCH -p debug
#SBATCH -D ./
#SBATCH -N 150 
#SBATCH -t 00:05:00
#SBATCH -J MPI7000_3600
#SBATCH -o MPI7000_3600.o%j
#SBATCH -e MPI7000_3600.e%j

module unload darshan
module unload xt-shmem
export CRAY_ROOTFS=DSL

# Run iAMOEBA using tensor method for 1000 time steps, 1.0 fs time step, update restart every 5.0 ps, Use NVT ensemble (2), Use Temp. of 300 K
srun -n 300 -c 12 ../../Src/dynamic_iAMOEBA_tensor.x -k ../inputfiles_allmodels/7000_iAMOEBA_tensor.key ../inputfiles_allmodels/7000_cluster.xyz 1000 1.0 5.0 2 300.0  > log7000_direct
