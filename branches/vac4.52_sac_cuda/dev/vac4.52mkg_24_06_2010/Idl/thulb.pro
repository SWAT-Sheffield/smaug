;^CFG COPYRIGHT VAC_UM
; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='../data/htmhd22.out';
npict=40;
.r getpict

pm=w(*,*,4)*w(*,*,4)+w(*,*,5)*w(*,*,5)
pu=w(*,*,1)*w(*,*,1)+w(*,*,2)*w(*,*,2)




