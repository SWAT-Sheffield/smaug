#!/bin/bash
#$ -j y
#$ -V
export TMPDIR="/scratch"

./vac < vac_nm_bach3D_t1.par
