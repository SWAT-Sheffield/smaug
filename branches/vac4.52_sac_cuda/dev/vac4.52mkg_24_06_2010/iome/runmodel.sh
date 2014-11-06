#!/bin/bash

#Set the parameters for compiling vac and vacini
#user provides a string containing the switches for vac
#switches are as follows

#./setvac -on=mpi
#./setvac -g=1976,44

#Setting for demo is
#./setvac -d=22 -phi=0 -z=0 -g=1976,44 -p=mhd -u=sim1 -on=cd,rk,mpi -off=mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist

echo $1

#load the default simulation params
if [ $# -gt 0 ]
then
	source $1
else
	source simulationparams.sh
fi


echo $2
#load the default model params
if [ $# -gt 1 ]
then
	source $2
else
	source modelparams.sh
fi

#search for mpi in $onpar
#using bash string replacement
testvar=${onpar/mpi/xyzxxx}


if [ $testvar == $onpar ]
then
  #off  -  mpi not found in $testvar
  mpion="0"


else
  #on  -  mpi not found in $testvar
  mpion="1"

fi

if [ $runmake == "1" ]
then




	cd src

	#show the current settings
	echo "current settings"
	./setvac -s


	./setvac -d=$dpar -g=$gpar -phi=$phipar -z=$zpar -p=$ppar -u=$upar -on=$onpar -f=$offpar -s

	#show the current settings
	echo "new settings for making vac"
	./setvac -s

	make clean

	make vac -f $vacmakefile

	./setvac -d=$dpar -g=$gparvacini -phi=$phipar -z=$zpar -p=$ppar -u=$upar -on=$onparvacini -f=$offparvacini -s

	#show the current settings
	echo "new settings for making vacini"
	./setvac -s
	
	make vacini -f $vacinimakefile

#	cp vac ../vac
#	cp vacini ../vacini

	cd ..

fi
#end of the make vac loop

if [ $runvacini == "1" ]
then
	#./vacini < vacini/vacini.par
	./vacini < $vaciniparfile
fi

#if we are parallel 
#distribute the model 
# execute /data/distribution to get help
if [ $mpion == "1" ]
then
	#the distributed infile must be of the form
        #filename_npFFLL.ini
        #FF and LL are integeres
        #FF is the first processor
        #LL is the last processor
        #example
        #filename_np0110.ini
   
	data/distribution -D -s=0 $inifile $distribinifile

fi

if [ $runvac == "1" ]
then
	#./vacini < vacini/vacini.par
	./vac < $vacparfile
fi




