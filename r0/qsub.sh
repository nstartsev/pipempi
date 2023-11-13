#!/bin/sh

#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -p hpc4-3d
#SBATCH -D /s/ls4/users/startsevnikolay/Pipe_modified/r0
#SBATCH -t 14:15:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node 16

module load intel-parallel-studio/2017

hostname

mpirun -np 32 ./Pipe
