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
close,1
;openr,1,'./grav2_test0.out',/f77_unf
openr,1,'/home/mikeg/data/vac3d/3D_tube_196_100_100_multidriver_test.out',/f77_unf
for npict=1,50 do begin

readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2,n3)
readu,1,x
for iw=0,nw-1 do begin
  
 readu,1,wi
 for i=0,n2-1 do w(*,i,*,iw)=reform(wi(*,i,*))
endfor


vtkfile='e';
vacscalar2vtk3d,npict,w(*,*,*,*),x(*,*,*),4,1,vtkfile

vtkfile='rho';
vacscalar2vtk3d,npict,w(*,*,*,*),x(*,*,*),0,1,vtkfile

vtkfile='mom';
vac2vtk3d,npict,w(*,*,*,*),x(*,*,*),1,3,vtkfile

vtkfile='b';
vac2vtk3d,npict,w(*,*,*,*),x(*,*,*),5,3,vtkfile

endfor

end
