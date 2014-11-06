#!/bin/bash

cp ../data/vacinputs/scripts/runconvalf2d_1.pro Idl/runconvalf2d_1.pro
cp ../data/vacinputs/tmp/convalf2d_1.pro Idl/convalf2d_1.pro
cd Idl
rm -rf vtk
mkdir vtk
mkdir %supmodelname%
gdl<runconvalf2d_1.pro
mv -f vtk/*.* %supmodelname%

cd ..
