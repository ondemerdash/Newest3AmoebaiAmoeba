#!/bin/bash -l
#SBATCH -p debug
#SBATCH -D ./
#SBATCH -N 200 
#SBATCH -t 00:05:00
#SBATCH -J 20clustMPI7000_4800
#SBATCH -o 20clustMPI7000_4800.o%j
#SBATCH -e 20clustMPI7000_4800.e%j

module unload darshan
module unload xt-shmem
export CRAY_ROOTFS=DSL

# Run 3M-AMOEBA for 1000 time steps, 1.0 fs time step, update restart every 5.0 ps, Use NVT ensemble (2), Use Temp. of 300 K
srun -n 400 -c 12 ../../Src/dynamic_3M_AMOEBA.x -k ../inputfiles_allmodels/7000_3MAMOEBA.key ../inputfiles_allmodels/7000_cluster.xyz ../inputfiles_allmodels/7000_20clust.txt2 1000 1.0 5.0 2 300.0  > log7000_20clust 
