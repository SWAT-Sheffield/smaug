#!/bin/bash
#$ -l h_rt=1:00:00
# Change 4 to the number of slots you want
#$ -pe mpi 4
# 8Gb per slot
#$ -l rmem=8G
#$ -l mem=8G
#$ -l gpu=1

module load libs/CUDA/7.5.18/binary                         
module load mpi/openmpi/1.10.4/gcc-4.9.4-TESTING

mpirun -np 4 bin/smaug
