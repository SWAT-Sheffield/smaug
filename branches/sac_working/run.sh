#!/bin/bash
#$ -l h_rt=99:00:00
#$ -cwd
#$ -l arch=intel*
#$ -l mem=6G

#$ -N S2D



module add compilers/pgi/14.4


echo Starting vac now.
./vac < vac.par


