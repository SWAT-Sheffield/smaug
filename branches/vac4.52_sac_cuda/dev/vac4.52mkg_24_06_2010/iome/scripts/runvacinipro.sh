#!/bin/bash

cp ../data/vacinputs/scripts/runvacinipro.pro Idl/runvacinipro.pro
cp ../data/vacinputs/tmp/genalfven2d_2vacini.pro Idl/genalfven2d_2vacini.pro
cd Idl
gdl<runvacinipro.pro

cd ..
