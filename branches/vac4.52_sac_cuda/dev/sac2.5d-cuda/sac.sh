#!/bin/bash

if [ $# -gt 0 ]; then
  IOME_SIMNAME=$1
else
  IOME_SIMNAME="mysim"
fi

if [ $# -gt 1 ]; then
  SIMFILE=$2
  echo "simfile is" $SIMFILE
fi



iogs initiome null $IOME_SIMNAME null >& iogs.err &
sleep 3
INPUT=`cat ${IOME_SIMNAME}0_port.txt`
IOME_WSPORT=$(echo $INPUT | cut -d' ' -f1 )
echo port is $IOME_WSPORT


if [ $# -gt 1 ]; then
#  iogs readsimulation simfile.xml 0 $IOME_WSPORT localhost
  ./iosac $IOME_SIMNAME $SIMFILE
else
  ./iosac $IOME_SIMNAME
fi

iogs exitiome 0 $IOME_WSPORT
#./killio.sh $IOME_SIMNAME > /dev/null
