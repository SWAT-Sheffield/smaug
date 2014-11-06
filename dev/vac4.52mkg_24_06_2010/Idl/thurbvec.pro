;^CFG COPYRIGHT VAC_UM
; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='../data/shearalfven2d.out';
npict=10;
.r getpict

;window,1,xsize=320,ysize=320
;window,2,xsize=320,ysize=320


;wset,1 & surf,w(*,*,0)
;wset,1 & vector,w(*,*,1),w(*,*,2)
;wset,2 & vector,w(*,*,4),w(*,*,5)






