;^CFG COPYRIGHT VAC_UM
; routine for converting output from sac2.5d to
; vtk format
; modify inputfilename and the number of files
;===========================================================================
filename='../out/test1.out'

;***************
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
openr,1,filename
;readf,1,headline
;readu,1,it,time,ndim,neqpar,nw
ndim=2
nw=8
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
nx(0)=150
nx(1)=150

if (ndim gt 0)then begin
  n1=nx(0)
endif
if (ndim gt 1) then begin
  n2=nx(1)
endif
if (ndim gt 2) then begin
  n3=nx(2)
endif

if (ndim eq 3) then begin
  x=dblarr(n1,n2,n3,ndim)
  w=dblarr(n1,n2,n3,nw)
endif else if (ndim eq 2) then begin
  x=dblarr(n1,n2,ndim)
  w=dblarr(n1,n2,nw)
endif else if (ndim eq 1) then begin
  x=dblarr(n1,ndim)
  w=dblarr(n1,nw)
endif


nfile=0
npictinfile=1

;for npict=1,23 do begin
for npict=1,99 do begin

readf,1,step
print,'step is ',step
for i=0,nx(0)-1 do begin
for j=0,nx(1)-1 do begin
 ;print, iw
 readf,1,iy,ix,rho,m1,m2,m3,energ,b1,b2,b3
 ; print,ix,iy
 x(i,j,0)=ix
 x(i,j,1)=iy
w(i,j,0)=rho
w(i,j,1)=m1
w(i,j,2)=m2
w(i,j,3)=m3
w(i,j,4)=energ
w(i,j,5)=b1
w(i,j,6)=b2
w(i,j,7)=b3



; print,i,j,rho,m1,m2,m3,energ,b1,b2,b3

endfor
endfor

print,'step is ',step



vtkfile='e';
vacscalar2vtk,npict,w(*,*,*),x(*,*),4,1,vtkfile

vtkfile='rho';
vacscalar2vtk,npict,w(*,*,*),x(*,*),0,1,vtkfile

vtkfile='mom';
vac2vtk,npict,w(*,*,*),x(*,*),1,3,vtkfile

vtkfile='b';
vac2vtk,npict,w(*,*,*),x(*,*),5,3,vtkfile

endfor

close,1

end
