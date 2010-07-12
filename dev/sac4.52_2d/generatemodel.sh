#!/bin/bash

#Set the parameters for compiling vac and vacini
#user provides a string containing the switches for vac
#switches are as follows

#./setvac -on=mpi
#./setvac -g=1976,44

#Setting for demo is
#./setvac -d=22 -phi=0 -z=0 -g=1976,44 -p=mhd -u=sim1 -on=cd,rk,mpi -off=mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist

dpar="22"
phipar="0"
zpar="0"
gpar="1976,44"
ppar="mhd"
upar="sim1"
onpar="cd,rk,mpi"
offpar="mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist"

cd src

./setvac -d=$dpar -phi=$phipar -z=$zpar -p=$ppar -u=$upar -on=$onpar -f=$offpar -s

#search for mpi in $onpar
if 
  #on


else
  #off
  onpar=0

fi

make vac
make vacini

cp vac ../vac
cp vacini ../vacini

cd ..
./vacini < vacini/vacini.par

#if we are parallel 
#distribute the model 
# execute /data/distribution to get help
if $onpar -eq 1
	
data/distribution -D -s=0 /data/vacinifile.ini /data/newdistributedinfile

fi




