; This file plots the 10-th picture from data/example22.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='../../data/test10_390_asc.out'
npict=10
getpict
headline='Example PostScript Plot'
func='h v1^v2'
plotmode='contfill vel'
multiplot=1
plottitle='h V'
set_plot,'PS'
device,filename='example.ps',/color,bits=8,xsize=18,ysize=18,$
	xoffset=1,yoffset=2
loadct,3
plotfunc
device,/close
set_plot,'X'
