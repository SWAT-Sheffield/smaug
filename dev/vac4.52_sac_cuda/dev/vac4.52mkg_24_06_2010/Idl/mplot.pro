;^CFG COPYRIGHT VAC_UM
; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE

getpict

window,1,xsize=320,ysize=320
window,2,xsize=320,ysize=320
wset,1 & surface,w(*,*,0)
wset,2 & vector,w(*,*,1),w(*,*,2)



