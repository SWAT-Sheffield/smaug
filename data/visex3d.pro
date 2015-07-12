;if (n_elements(size(w)) ne 7) then begin

;print,'read file?'
;read,as

loadct,4
gamma=1.666667

tmas=fltarr(1)
tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
avga=fltarr(1)
tena=fltarr(1)


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


;close,1
;openr,1,'./grav2_test1.out',/f77_unf

;readu,1,headline
;readu,1,it,time,ndim,neqpar,nw
;gencoord=(ndim lt 0)
;tarr=[tarr,time]
;ndim=abs(ndim)
;nx=lonarr(ndim)
;readu,1,nx
;eqpar=dblarr(neqpar)
;readu,1,eqpar
;readu,1,varname

;n1=nx(0)
;n2=nx(1)
;n3=nx(2)
;x=dblarr(n1,n2,n3,ndim)
;w1=dblarr(n1,n2,n3,nw)   ;was n1,n2,nw
;wi=dblarr(n1,n2,n3)
;readu,1,x
;for iw=0,nw-1 do begin
; readu,1,wi
; for i=0,n2-1 do w1(*,i,*,iw)=reform(wi(*,i,*))
;endfor



close,1
openr,1,'./grav2_test0.out',/f77_unf

!P.MULTI=[0,3,4]
;device,font_size=10

while not(eof(1)) do begin


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

;w=w-w1


zmax=max(x(*,*,*,2))/1E8   ;0
zmin=min(x(*,*,*,2))/1E8   ;0
xmax=max(x(*,*,*,1))/1E8   ;0
xmin=min(x(*,*,*,1))/1E8   ;0
ymax=max(x(*,*,*,0))/1E8   ;1
ymin=min(x(*,*,*,0))/1E8   ;1


rho=w(*,*,*,0)+w(*,*,*,9)

ttt=(gamma-1.0)*(w(*,*,*,4)+w(*,*,*,8)-0.5*(w(*,*,*,5)^2.0+w(*,*,*,6)^2.0+w(*,*,*,7)^2.0)-0.5*(w(*,*,*,1)^2+w(*,*,*,2)^2+w(*,*,*,3)^2)/w(*,*,*,0))*1.25/(8.31E7*w(*,*,*,0))


tvframe,reform(w(n1/2,*,*,0)),btitle='rho',charsize=2.0,/bar,/sample,xtitle='y',ytitle='z'

tvframe,reform(w(*,n2/2,*,0)),btitle='rho',charsize=2.0,/bar,/sample,xtitle='x',ytitle='z'

tvframe,reform(w(*,*,n3/2,0)),btitle='rho',charsize=2.0,/bar,/sample,xtitle='x',ytitle='y'

tvframe,reform(w(n1/2,*,*,1)/rho(n1/2,*,*)),btitle='vz',charsize=2.0,/bar;,/sample,xtitle='y',ytitle='z'

tvframe,reform(w(*,n2/2,*,1)/rho(*,n2/2,*)),btitle='vz',charsize=2.0,/bar,/sample,xtitle='x',ytitle='z'

tvframe,reform(w(*,*,n3/2,1)/rho(*,*,n3/2)),btitle='vz',charsize=2.0,/bar,/sample,xtitle='x',ytitle='y'

tvframe,reform(w(n1/2,*,*,2)/rho(n1/2,*,*)),btitle='vx',charsize=2.0,/bar,/sample,xtitle='y',ytitle='z'

tvframe,reform(w(*,n2/2,*,2)/rho(*,n2/2,*)),btitle='vx',charsize=2.0,/bar,/sample,xtitle='x',ytitle='z'

tvframe,reform(w(*,*,n3/2,2)/rho(*,*,n3/2)),btitle='vx',charsize=2.0,/bar,/sample,xtitle='x',ytitle='y'

tvframe,reform(w(n1/2,*,*,3)/rho(n1/2,*,*)),btitle='vy',charsize=2.0,/bar,/sample,xtitle='y',ytitle='z'

tvframe,reform(w(*,n2/2,*,3)/rho(*,n2/2,*)),btitle='vy',charsize=2.0,/bar,/sample,xtitle='x',ytitle='z'

tvframe,reform(w(*,*,n3/2,3)/rho(*,*,n3/2)),btitle='vy',charsize=2.0,/bar,/sample,xtitle='x',ytitle='z'


print,time

endwhile



tarr=reform(tarr[1:n_elements(tarr)-2])
maxa=reform(maxa[1:n_elements(maxa)-2])
mina=reform(mina[1:n_elements(mina)-2])
avga=reform(avga[1:n_elements(avga)-2])
tena=reform(tena[1:n_elements(tena)-2])

end
