#!/bin/bash -x
#SBATCH --job-name=s_alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --time=96:00:00
#SBATCH --partition=batch

module purge
module load gcc/10.3.0
module load fftw/gcc/single/sse/3.3.9
module load gromacs/nompi/cpu/gcc/single/2021.1


gfortran alpha.f90
./a.out
