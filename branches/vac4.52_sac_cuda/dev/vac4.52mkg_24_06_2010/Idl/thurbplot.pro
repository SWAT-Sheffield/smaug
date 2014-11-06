;^CFG COPYRIGHT VAC_UM
; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='../data/shearalfven2d.ini';
npict=1;
.r getpict

sx=15
window,3,xsize=320,ysize=320
window,4,xsize=320,ysize=320
window,5,xsize=320,ysize=320
window,6,xsize=320,ysize=320

wset,3 & plot,w(sx,*,0)
wset,4 & plot,w(sx,*,3)
wset,5 & plot,pm(sx,*)
wset,6 & plot,pu(sx,*)




