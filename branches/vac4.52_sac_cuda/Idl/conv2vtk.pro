;^CFG COPYRIGHT VAC_UM
;===========================================================================
;    Read the npict-th picture from an ascii or binary ini or out file 
;
;    Usage: 
;
; .r getpict
;
;    "getpict" will prompt you for "filename(s)" and "npict"
;    unless they are already set. Previous settings can be erased by 
;
; .r defaults
;
;    or modified explicitly, e.g.:
;
; filename='data/example.ini'
; npict=1
;
;    The "x" and "w" arrays and the header info will be read from the file. 
;
;    If a file is read with generalized coordinates, "gencoord=1" is set,
;    and the original data is transformed according to the "transform"
;    string variable into "xreg" and "wreg".
;
;    The same npict-th snapshot can be read from 2 or 3 files by e.g. setting
;
; filename='data/file1.ini data/file2.out'
;
;    In this case the data is read into x0,w0 and x1,w1 for the two files,
;    and possibly transformeed into wreg0,wreg1.
;
;    To plot a variable, type e.g.:
;
; surface,w(*,*,2)
;
;    or 
;
; .r plotfunc
;
;===========================================================================
;filename='../data/example23.out'
filename='../data/otmhd22.out'
directory='../data/'
nfile=0
npictinfile=1
headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)
nn=0
close,1
openr,1,'/home/cs1mkg/data/vac3d/3D_tube_196_100_100_new.ini',/f77_unf

for npict=1,1 do begin
;*****************************************************

readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

n1=nx(0)
n2=nx(1)
n3=nx(2)

x=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)

wi=dblarr(n1,n2,n3)

readu,1,x
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor
print, n1,n2,n3

;*************************************************





vtkfile='ne';
vacscalar2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),4,1,vtkfile

vtkfile='neb';
vacscalar2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),8,1,vtkfile

vtkfile='nrhob';
vacscalar2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),9,1,vtkfile

vtkfile='nrho';
vacscalar2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),0,1,vtkfile

vtkfile='nmom';
vac2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),1,3,vtkfile

vtkfile='nb';
vac2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),5,3,vtkfile

vtkfile='nbb';
vac2vtk3d,npict,w(*,*,*,*),x(*,*,*,*),10,3,vtkfile




;****************************************************
;write the fields here
outfile=directory+'ascdat'+strtrim(string(npict),2)+'.out'
openw,3,outfile



printf,3,npict
;for i=0,nx(0)-1 do begin
;for j=0,nx(1)-1 do begin
; ix=x(i,j,0)
; iy=x(i,j,1)

;  printf,3,j,i,w(i,j,0),w(i,j,1),w(i,j,2),w(i,j,3);,format='(i),(X),(i),(X)'
;  printf,3,i,j,w(i,j,0),w(i,j,1),w(i,j,2),w(i,j,3),w(i,j,4),w(i,j,5),w(i,j,6),w(i,j,7);,format='(i),(X),(i),(X)'
;  printf,3,i,j,w(i,j,0),w(i,j,1),w(i,j,2),w(i,j,3),w(i,j,4),w(i,j,5),w(i,j,6),w(i,j,7),format='(i,x,i,x,f,x,f,x,f,x,f,x,f,x,f,f,x,f)'
 
;writeu,1,wi
;endfor
;endfor


 
close,3



endfor

end
