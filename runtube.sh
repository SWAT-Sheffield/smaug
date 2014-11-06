#!/bin/bash
#$ -j y
#$ -l arch=intel*
#$ -N sactubet1part3
#$ -l mem=6G
#$ -l h_rt=168:00:00


module add compilers/pgi/10.2
./vac < vac_tube_128_128_128_t6.par
