#!/bin/bash
#$ -j y
#$ -l arch=intel*
#$ -l gpu=1
#$ -l gpu_arch=nvidia-m2070
##$ -l gpu=1,gpu_arch=nvidia-k40m
##$ -P mhd
#$ -N tube1_128
#$ -l mem=12G
#$ -l rmem=12G
#$ -l h_rt=168:00:00
module load libs/cuda/6.5.14



cd include
#cp iosmaugparams_ot_1020.h iosmaugparams.h
#cp iosmaugparams_ot_256.h iosmaugparams.h
cp iosmaugparams_tube1_128.h iosmaugparams.h

cd ..

cd src
cp usersource_tube1_128.cu usersource.cu
cp boundary_3d.cu boundary.cu
make clean
#make -f Makefile_k40 smaug
#make smaug
make -f Makefile_3d smaug

cd ..



export TIMECOUNTER=0
source timeused
bin/smaug a
source timeused

