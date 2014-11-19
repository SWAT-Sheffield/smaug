#!/bin/bash
#$ -l arch=intel-e5-2650v2
#$ -l mem=12G
#$ -l rmem=12G
#$ -P cs-test



source /data/cs1mkg/tools/vapor/vapor-2.2.2/bin/vapor-setup.sh
export IDL_STARTUP=runconvert_spic
idl
