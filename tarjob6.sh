#!/bin/bash

#$ -l arch=intel*
#$ -l h_rt=168:00:00

#$ -N smaug_datapack
##test

tar -zcvf spic4b0_3d.tgz spic4b0_3d/*zerospic1*
tar -zcvf spic4b0_b1_3d.tgz spic4b0_b1_3d/*zerospic1*
tar -zcvf spic4b0_1_3d.tgz spic4b0_1_3d/*zerospic1*
tar -zcvf spic4b0_2_3d.tgz spic4b0_2_3d/*zerospic1*
tar -zcvf spic4b0_3_3d.tgz spic4b0_3_3d/*zerospic1*
tar -zcvf spic4b1_1_3d.tgz spic4b1_1_3d/*zerospic1*
tar -zcvf spic4b1_2_3d.tgz spic4b1_2_3d/*zerospic1*


