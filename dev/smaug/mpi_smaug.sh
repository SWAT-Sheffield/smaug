#!/bin/bash
#$ -l rmem=4G
#$ -l mem=4G
#$ -l gpu=4

module load libs/CUDA/8.0.44/binary                         
module load mpi/openmpi/1.10.4/gcc-4.9.4-TESTING

mpirun -np 4 bin/smaug
