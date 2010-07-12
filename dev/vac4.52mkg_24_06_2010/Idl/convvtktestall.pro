; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='../data/shearalfven2d_1.out';

for i=1,2 do begin;
npict=i;
.r getpict
vtkfile='mom';
vac2vtk,npict,w(*,*,*),x(*,*),1,2,vtkfile

endfor
end

