fleft=30
fright=220

nj=512/2

dsl=dindgen(nj)-fleft
dsr=dindgen(nj)-fright

ww0=dindgen(2*nj)


dsl=atan(dsl)
dsr=-atan(dsr)

maxdsl=max(dsl)
mindsl=min(dsl)

maxdsr=max(dsr)
mindsr=min(dsr)

   dsl=(dsl-mindsl)/(maxdsl-mindsl)
   dsr=(dsr-mindsr)/(maxdsr-mindsr)

for i=0,2*nj-1 do begin
 if i lt nj then ww0(i)=dsl(i) else ww0(i)=dsr(i-nj) 
endfor  

www=dblarr(144,2*nj)

for i=0, 144-1 do begin
   www(i,*)=1.d0-ww0(*)
endfor

tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

tmas=fltarr(1)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, 0, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
WINDOW, 1, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

loadct, 4

mass=dblarr(1)
Vzmax=dblarr(1)

ia=1.0


dd=40

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

ii=1
close,1




aa=0

if aa eq 0 then begin

openr,1,'/data/ap1vf/3MHDCD4mpi16512_test.out',/f77_unf


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
w=dblarr(n1,n2,n3,nw)
wi=dblarr(n1,n2,n3)

readu,1,x


xout=dblarr(2)
yout=dblarr(2)
zout=dblarr(2)

for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,*,iw)=wi
endfor

endif

!p.multi =  0

xzstart=0
xzend=511

xystart=0
xyend=143
kk=72


w0=dblarr(n2,n3)
; for i=xzstart,xzend do begin
;   w0(*,i)=1.d0-exp(-((i-(xzend-xzstart)/2.d0)/ 1.0/(xzend-xzstart) )^2.d0)
; endfor

;maxw0=max(w0)
;minw0=min(w0)

; for i=xzstart,xzend do begin
;   w0(*,i)=(w0(*,i)-minw0)/(maxw0-minw0)
; endfor
 
w0=www

window, 0,xsize=500,ysize=500

tvframe, w0, /bar

stop

wF=dblarr(n2,n3)
wSin=dblarr(n2,n3)


xyi=findgen(xyend+1)
xzi=findgen(xzend+1)


for i=xystart,xyend do begin
 for j=xzstart,xzend do begin
   wSin(i,j)=sin(xyi(i)*!Pi/(xyi(xyend)-xyi(xystart)))* $
                sin(xzi(j)*!Pi/(xzi(xzend)-xzi(xzstart)))
 endfor
endfor

window, 1,xsize=500,ysize=500
  
tvframe,w(kk,*,*,3)/w(kk,*,*,0),/bar, /sample

stop

wF(*,*)=(w(kk,*,*,3)/w(kk,*,*,0)) ;*wSin(*,*)

wF=FFT(wF(*,*),-1)


stop

wF1=ABS(wF)^2.d0
wF1=reform(wF1)

ff=total(wF1,1)

window, 2,xsize=500,ysize=500

tvframe, alog10(wF1),/bar

;plot,ff
stop


wF3=reform(wF)
wF4=w0*wF3


wF4=double(FFT(wF4,1)) ;/(wSin+0.00001))

window, 3, xsize=500,ysize=500

tvframe, wF4,/bar,/sample
stop


wF3=reform(wF)

wF4=ABS(wF3*w0)^2.d0

tvframe, alog10(wF4), /bar

stop





end
