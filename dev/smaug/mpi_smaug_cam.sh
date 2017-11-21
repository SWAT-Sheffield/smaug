#!/bin/bash
#$ -l rmem=12G
#$ -l mem=12G
#$ -l gpu=1

module load cuda/8.0                       
module load openmpi/gnu/1.10.3

mpirun -np $1 bin/smaug
