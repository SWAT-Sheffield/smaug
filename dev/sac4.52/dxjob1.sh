#!/bin/bash
#$ -cwd
#$ -l h_cpu=96:00:00
#$ -j y



export DXPROCESSORS=1
dx -script vnn_prepdat6.net
