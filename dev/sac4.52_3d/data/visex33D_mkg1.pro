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

;if (ii eq 1) then begin
;loadct,4
;mixct
;endif 
;else begin
;loadct,0
;tek_color
;endelse




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

window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
;window, 1,xsize=1500,ysize=600,XPOS = 500, YPOS = 80



nn=0
np=0
kkk=4

nn_i=0

close,1
close,2




;openr,1,'/data/ap1vf/3D_509_36_36_300s.out',/f77_unf
;openr,1,'/data/ap1vf/3D_tube_196_100_100.ini',/f77_unf

openr,1,'/home/mikeg/data/vac3d/3D_tube_196_100_100_multidriver_test.out',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60t.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100t.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_120s_full.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_200s_puls.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf


readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


print, 'tuta1'
xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

wi=dblarr(n1,n2,n3)

w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

xx=dblarr(n2)
yy=dblarr(n3)
zz=dblarr(n1)


xx(*)=x(1,*,1,1)
yy(*)=x(1,1,*,2)
zz(*)=x(*,1,1,0)

Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)


mu=4.0*!PI/1.0e7

print,'******************* time = ', time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

scale=1.d6

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=99
ystart=0
yend=99

pp=50 ;x
kk=5  ;y

wset,0
!p.multi = [0,4,4,0,1]


zstart=0
zend=185

wt=dblarr(zend-zstart+1,xend-xstart+1,iw)
wt=reform(w(zstart:zend,xstart:xend,pp,*))

wy=dblarr(n1,n3,iw)
wy=reform(w(zstart:zend,pp,*,*))

wt(*,*,3)=reform(w(zstart:zend,xstart:xend,pp+2,3))

wt(*,*,12)=reform(w(zstart:zend,pp,ystart:yend,12))

saeb=dblarr(zend-zstart+1,xend-xstart+1)
sarho_t=dblarr(zend-zstart+1,xend-xstart+1)
sabz_t=dblarr(zend-zstart+1,xend-xstart+1)
sabx_t=dblarr(zend-zstart+1,xend-xstart+1)
saby_t=dblarr(zend-zstart+1,xend-xstart+1)

saeb(*,*)=wt(*,*,8)
sarho_t(*,*)=wt(*,*,0)+wt(*,*,9)
sabz_t(*,*)=wt(*,*,10)
sabx_t(*,*)=wt(*,*,11)
saby_t(*,*)=wt(*,*,12)

vt=dblarr(n1,n2,n3)
vvt=dblarr(n2,n3)
;vt(*,*,*)=sqrt(w(*,*,*,1)^2.d0+w(*,*,*,2)^2.d0+w(*,*,*,3)^2.d0)/(w(*,*,*,0)+w(*,*,*,9))

;****************** Pressure background begin ********************
TP=saeb
TP=TP-(sabx_t^2.0+saby_t^2.0+sabz_t^2.0)/2.0
TP=(gamma-1.d0)*TP
;****************** Pressure background end ********************

;tvframe,rotate(wt(*,*,0),1), /bar,title='rho',$ 
;        /sample, xtitle='x', ytitle='y',charsize=2.0, $
;	xrange=[xx[xstart]/scale, xx[xend]/scale], $
;	yrange=[zz[zstart]/scale, zz[zend]/scale]

end
