#!/bin/sh 


projectname="alfven2d_2"
simulationname="defaultvacsim"
simulationnotes="alfven2d regenerate ini file using shock tube to enable pulse propagation"

author="mikeg"
vacversion="vac4.52_mkg_23_06_2010"
vacdownload="http://www.scispace.com/action/file/download?file_guid=66566"
creationdate=`date`
computeresource=`echo $HOSTNAME`

dpar="22"
phipar="0"
zpar="0"
gpar="260,260"
gparvacini="260,260"
ppar="mhd"
upar="gravity"
onpar="cd,tvdlf,tvd,poisson,resist"
onparvacini="cd,tvdlf,tvd,poisson,resist"
offpar="mc,fct,impl,ct,gencoord,rk,mpi"
offparvacini="mc,fct,impl,ct,gencoord,rk,mpi"
vacinimakefile="Makefile"
vacmakefile="Makefile"
mpion="0"
runmake="0"
runvacini="1"
runvac="1"
runprescript1="1"
runprescript2="1"
runpostscript1="1"
prescript1="tmp/testscript1.sh"
prescript2="scripts/runvacinipro.sh"
postscript1="scripts/runconvalf2d_2pro.sh"
genreplacementfiles="testscript1.sh testscript2.sh"
templatevaciniparfile="templates/tempalf2d_1vacini.par"
templatevacparfile="templates/tempalf2d_1vac.par"

modelparams="../data/vacinputs/scripts/alf2d_2modelparams.sh"


inputdirectory="../data/vacinputs/"
outputdirectory="../data/vacoutputs/"

modelname="alf2d_1"
modelnotes="test alfven2d eta 0"
simulationfile=$inputdirectory"xml/"$modelname".xml"
rundate=`date`



vaciniparfile=$inputdirectory"tmp/alf2dvacini.par"
vacparfile=$inputdirectory"tmp/alf2dvac.par"

#distributed ini file configured for 10 processor iceberg
	#the distributed infile must be of the form
        #filename_npFFLL.ini
        #FF and LL are integeres
        #FF is the first processor
        #LL is the last processor
        #example
        #filename_np0110.ini
inifile=$inputdirectory"modshearalfven2d.ini"
distribinifile=$inputdirectory"vacinifile_np0110.ini"

#vacpar parameters
#file must contain the line starting #vacparlist enclosed by %%
#     this line contains a space separated list of the variables
#     the first character of this variable name is its type 
#     d=double, i=integer, v=vector, s=string
#%vacparlist% suplogfile supoutfile
par[0]="shearalfven2d_3.log"
par[1]="shearalfven2d_3.out"



#vacinipar parameters
#file must contain the line starting #vaciniparlist enclosed by %%
#     this line contains a space separated list of the variables
#     the first character of this variable name is its type 
#     d=double, i=integer, v=vector, s=string
#%vaciniparlist% supgamma supeta supg1 supg2 suprho supv1 supv2 supp supb1 supb2
ipar[0]="2"
ipar[1]="1.0"
ipar[2]="0.0"
ipar[3]="0.0"
ipar[4]="0.1"
ipar[5]="0.1"
ipar[6]="0.1"
ipar[7]="0.1"
ipar[8]="0.1"
ipar[9]="0.1"

#ipar1=""
#ipar2=""
#ipar3=""

#general file replacement parameters
#file must contain the line starting #vaciniparlist enclosed by %%
#     this line contains a space separated list of the variables
#     the first character of this variable name is its type 
#     d=double, i=integer, v=vector, s=string
#%genparlist% sup1 sup2 sup3
gpar[0]="2"
gpar[1]="0.0"
gpar[2]="0.2"


