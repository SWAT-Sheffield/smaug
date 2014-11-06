; Reading initialisation file used by versatile advection code
; read the vacini file which is all zeros, set the 
; fields as required and write to the new ini file
; 
; user provide valid vac initialisation file with correct header information
; modify file and directory as required
;
; generates shock tube configuration
directory='/home/mike/proj/sparc/data/vacinputs/ini/'
file='modshearalfven2d.ini'
;file='shearalfven.ini'

infile=directory+file
outfile=directory+'new'+file
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


openr,1,infile,/f77_unf

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
  x_code=dblarr(n1,n2,n3,ndim)
  w=dblarr(n1,n2,n3,nw)
  wi=dblarr(n1,n2,n3)
endif else if (ndim eq 2) then begin
  x_code=dblarr(n1,n2,ndim)
  w=dblarr(n1,n2,nw)
  wi=dblarr(n1,n2)
endif else if (ndim eq 1) then begin
  x_code=dblarr(n1,ndim)
  w=dblarr(n1,nw)
  wi=dblarr(n1)
endif

readu,1,x_code
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
 if (ndim eq 3) then begin
    w(*,*,*,iw)=wi
 endif else if (ndim eq 2) then begin
    w(*,*,iw)=wi
 endif else if (ndim eq 1) then begin
    w(*,iw)=wi
 endif
endfor
print, n1

close,1

;modify the fields here
;**************************************************
; reproduce 1d keppens tutorial example 2d version

rho=1
m1=0
m2=0
e=0.5
b1=1.0
b2=0.0

w(*,*,0)=rho
w(*,*,1)=0
w(*,*,2)=0
w(*,*,3)=0.5
w(*,*,4)=1.0
w(*,*,5)=0.0

;determine m2 below
width=10
starty=120

seg1=62
seg2=125
seg3=187

m2max=0.001
for i2=0,n2-1 do begin
  for i1=starty,starty+width do begin

   if (i1 gt seg1) then begin
    if (i1 lt seg2) then begin
      w(i1,i2,2)=m2max
    endif
   endif

   if (i1 gt seg2) then begin
    if (i1 lt seg3) then begin
      w(i1,i2,2)=m2max*(i1-seg2)/(seg3-seg2)
    endif
   endif

   if (i1 ge seg3) then begin
      w(i1,i2,2)=m2max*(n2-i1)/(n2-seg3)
   endif

  endfor
endfor



;****************************************************
;write the fields here

openw,1,outfile,/f77_unf
writeu,1,headline
writeu,1,it,time,ndim,neqpar,nw
writeu,1,nx
writeu,1,eqpar
writeu,1,varname
writeu,1,x_code
for iw=0,nw-1 do begin
 if (ndim eq 3) then begin
    wi=w(*,*,*,iw)
 endif else if (ndim eq 2) then begin
    wi=w(*,*,iw)
 endif else if (ndim eq 1) then begin
    wi=w(*,iw)
 endif
writeu,1,wi
endfor


 
close,1

;wwww :
print, 'complete'



end


