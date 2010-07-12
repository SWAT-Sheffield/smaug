; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='../data/shearalfven2d.out';
npict=15;
.r getpict
vtkfile='mom';
vac2vtk,npict,w(*,*,*),x(*,*),1,2,vtkfile
