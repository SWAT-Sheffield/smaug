#!/bin/bash
#$ -cwd
#$ -l h_cpu=96:00:00
#$ -j y
nsteps=415
a=401
date 
while [ $a -le $nsteps ]
do
rep="s/%step%/"$a"/g"
echo $rep > subst.txt
sed -f subst.txt convertdata2.txt > tmp1.txt
sed -f subst.txt convertdata2.txt > convtemp1.txt
./convertdata < convtemp1.txt
a=$(($a+1))

done

date
