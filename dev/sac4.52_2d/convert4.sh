#!/bin/bash
export VACMATFILENAME="../testcont2.out"


nsteps=413
a=412
date 
while [ $a -le $nsteps ]
do
	export VACMATNUMPICT=$a
	export VACMATOUTFILENAME="../testascmat_"$VACMATNUMPICT"_rub.out"
	matlab -nosplash -nodisplay -r mygetpict

	a=$(($a+1))

done

date

