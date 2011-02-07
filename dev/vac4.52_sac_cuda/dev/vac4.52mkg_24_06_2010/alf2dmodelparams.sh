#!/bin/sh 


simulationparams="alf2dsimparams.sh"
simulationfile=$simulationname".xml"
rundate=`date`

modelname="alfven2d_0"
modelnotes="test alfven2d eta 0"

vaciniparfile="alf2dvacini.par"
vacparfile="alf2dvac.par"

#distributed ini file configured for 10 processor iceberg
	#the distributed infile must be of the form
        #filename_npFFLL.ini
        #FF and LL are integeres
        #FF is the first processor
        #LL is the last processor
        #example
        #filename_np0110.ini
inifile="data/shearalfven2d.ini"
distribinifile="data/vacinifile_np0110.ini"

#vacpar parameters
#file must contain the line starting #vacparlist enclosed by %%
#     this line contains a space separated list of the variables
#     the first character of this variable name is its type 
#     d=double, i=integer, v=vector, s=string
#%vacparlist% suplogfile supoutfile
par[0]="shearalfven2d_0.log"
par[1]="shearalfven2d_0.out"


#vacinipar parameters
#file must contain the line starting #vaciniparlist enclosed by %%
#     this line contains a space separated list of the variables
#     the first character of this variable name is its type 
#     d=double, i=integer, v=vector, s=string
#%vaciniparlist% supgamma supeta supg1 supg2 suprho supv1 supv2 supp supb1 supb2
ipar[0]="2"
ipar[1]="0.0"
ipar[2]="0.2"
ipar[3]="0.3"
ipar[4]="0.1"
ipar[5]="0.1"
ipar[6]="0.1"
ipar[7]="0.1"
ipar[8]="0.1"
ipar[9]="0.1"
#ipar1=""
#ipar2=""
#ipar3=""
