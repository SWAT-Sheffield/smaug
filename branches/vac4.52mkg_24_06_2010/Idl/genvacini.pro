; Reading initialisation file used by versatile advection code
; read the vacini file which is all zeros, set the 
; fields as required and write to the new ini file
; 
; user provide valid vac initialisation file with correct header information
; modify file and directory as required
;
directory='/home/mike/proj/sparc/vac4.52work/data/'
file='newshearalfven2d.ini'
;file='newshearalfven.ini'

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


