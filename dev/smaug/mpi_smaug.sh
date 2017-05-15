#!/bin/bash
#$ -l rmem=12G
#$ -l mem=12G
#$ -l gpu=1

module load libs/CUDA/8.0.44/binary                         
module load mpi/openmpi/1.10.4/gcc-4.9.4

mpirun -np 1 bin/smaug
