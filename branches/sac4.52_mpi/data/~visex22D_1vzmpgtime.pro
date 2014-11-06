tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


ii=1

if (ii eq 1) then begin
;loadct,4
mixct
endif else begin
loadct,0
tek_color
endelse




mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)

ia=1.0

headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)

; Open an MPEG sequence: 
;mpegID = MPEG_OPEN([700,1200],FILENAME='myMovie.mpg') 

window, 0,xsize=600,ysize=600,XPOS = 350, YPOS = 300 


nn=0
nn_i=0

close,1
close,2

openr,1,'/data/ap1vf/zero3B.out',/f77_unf


j=0

xend=1800

Vzt=dblarr(509,xend)


while (j lt xend) do begin
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




xout=dblarr(2)
yout=dblarr(2)


n1=nx(0)
n2=nx(1)
x=dblarr(n1,n2,ndim)
if (nn eq 0) then w=dblarr(n2,n1,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2)
;e2=dblarr(n1,n2)
readu,1,x
;readu,1,e2
;e2=rotate(e2,1)
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=rotate(wi,1)
endfor


Vt=dblarr(n1,n2)
B=dblarr(n1,n2)
B_bg=dblarr(n1,n2)

p=dblarr(n1,n2,1)
;e2=dblarr(n1,n2)

print,time


;wset,0
;!p.multi = 0
if (ii eq 1) then begin

Vzt(*,j)=w(55,*,1)/(w(55,*,0)+w(55,*,7))
;tvframe,Vzt,/sample, /bar,title='vz', xtitle='x', ytitle='y',charsize=1.0  
print,j
j=j+1

endif else begin

endelse

endwhile

tvframe,Vzt,/sample, /bar,title='vz', xtitle='x', ytitle='y',charsize=1.0  
 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device

device, /close
set_plot, 'x'


end
