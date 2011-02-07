#!/bin/bash
nsteps=210
a=200
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
