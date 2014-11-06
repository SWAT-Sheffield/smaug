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

istart=1

if (ii eq 1) then begin
;loadct,4
;mixct
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

window, 0,xsize=1100,ysize=450,XPOS = 350, YPOS = 300 


nn=799
kkk=4

nn_i=0

close,1
close,2

print, 'tuta'
openr,1,'/data/ap1vf/2_6Mnzx1976400B41cont.out',/f77_unf

print, 'tuta 1'

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




xout=dblarr(2)
yout=dblarr(2)


n1=nx(0)
n2=nx(1)
x=dblarr(n1,n2,ndim)
if (nn eq 0) then w=dblarr(n2,n1,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2)
readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=rotate(wi,1)
endfor

zmin=min(x(*,1,0))
zmax=max(x(*,1,0))

xmin=min(x(1,*,1))
xmax=max(x(1,*,1))

scale=1.0d6

Vt=dblarr(n1,n2)
B=dblarr(n1,n2)
B_bg=dblarr(n1,n2)

p=dblarr(n1,n2,1)
;e2=dblarr(n1,n2)


mu=4.0*!PI/1.0e7

print,time

kk=200

label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

arho=w(*,*,0)
arho_t=w(*,*,7)+w(*,*,0)
avz=w(*,*,1)/arho_t
avx=w(*,*,2)/arho_t
amz=w(*,*,1)
amx=w(*,*,2)
ae=w(*,*,3)
aeb=w(*,*,6)
arho=w(*,*,0)
abz=w(*,*,4)
abx=w(*,*,5)
abz_t=w(*,*,8)+w(*,*,4)
abx_t=w(*,*,9)+w(*,*,5)
abz0=w(*,*,8)
abx0=w(*,*,9)

R=8.3e+003
mu_gas=0.6
gamma=1.66667

zstart=0
zend=1975

xstart=0
xend=399



;zstart=0
;zend=300

;xstart=0
;xend=799


n11=zend-zstart+1
n22=xend-xstart+1

Bxx=dblarr(n22,n11)
Bzz=dblarr(n22,n11)

xx=dblarr(n22)
zz=dblarr(n11)

xx(*)=x(1,xstart:xend,1)
zz(*)=x(zstart:zend,1,0)


sarho_t=dblarr(n22,n11)
savz=dblarr(n22,n11)
savx=dblarr(n22,n11)
samz=dblarr(n22,n11)
samx=dblarr(n22,n11)
sae=dblarr(n22,n11)
saeb=dblarr(n22,n11)
sarho=dblarr(n22,n11)
sabz=dblarr(n22,n11)
sabx=dblarr(n22,n11)
sabz_t=dblarr(n22,n11)
sabx_t=dblarr(n22,n11)
sabz0=dblarr(n22,n11)
sabx0=dblarr(n22,n11)

sarho_t=arho_t(xstart:xend,zstart:zend)
savz=avz(xstart:xend,zstart:zend)
savx=avx(xstart:xend,zstart:zend)
samz=amz(xstart:xend,zstart:zend)
samx=amx(xstart:xend,zstart:zend)
sae=ae(xstart:xend,zstart:zend)
saeb=aeb(xstart:xend,zstart:zend)
sarho=arho(xstart:xend,zstart:zend)
sabz=abz(xstart:xend,zstart:zend)
sabx=abx(xstart:xend,zstart:zend)
sabz_t=abz_t(xstart:xend,zstart:zend)
sabx_t=abx_t(xstart:xend,zstart:zend)
sabz0=abz0(xstart:xend,zstart:zend)
sabx0=abx0(xstart:xend,zstart:zend)



btot=sqrt(sabx_t^2.0+sabz_t^2.0)


Cosph=savx*sabx0+savz*sabz0 

AbsV=sqrt(savx^2.d0+savz^2.d0)
AbsB=sqrt(sabx0^2.d0+sabz0^2.d0)

Vpl=Cosph/AbsB

Vpl(where(abs(sabz0) le 1.0d-260))=0.d0

Cosph=savx*sabz0-savz*sabx0 

Vpd=Cosph/AbsB

Vpd(where(abs(sabz0) le 1.0d-260))=0.d0


;****************** Pressure begin ********************
T=sae+saeb
T=T-(savx^2.0+savz^2.0)*sarho_t/2.0
T=T-(sabx_t^2.0+sabz_t^2.0)/2.0
T=(gamma-1.d0)*T
;****************** Pressure end ********************


wset,0
!p.multi = [0,2,0,0,1]
tvframe,Vpl,/sample, /bar,title='!3V!I| |!N!3 [m/s]', xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]


xx(*)=x(1,xstart:xend,1)
zz(*)=x(zstart:zend,1,0)


nxz=30
nxnz = [nxz,nxz]

xxi=interpol(xx,nxz)
zzi=interpol(zz,nxz)

avxi=congrid(savx,nxz,nxz)
avzi=congrid(savz,nxz,nxz)

;avxi=interp2d(savx,xx,zz,xxi,zzi,nxnz, /grid)
;avzi=interp2d(savz,xx,zz,xxi,zzi,nxnz, /grid)
;stop

tek_color
VELOVECT, avxi, avzi,xxi/scale,zzi/scale, /OVERPLOT ;, color=6

line2_b_classic, sabx_t,sabz_t,xx, zz, 1, scale ; 0 - plot, 1 - oplot

beta=(sabx_t^2.0+sabz_t^2.0)/2.d0/T
xx(*)=x(1,xstart:xend,1)
zz(*)=x(zstart:zend,1,0)

tek_color
contour, beta,xx/scale,zz/scale, LEVELS = [0.001,0.01,0.1,1.0, 10.0], $
         C_Annotation = ['1000.0','100.0','10.0','1.0','0.1'], /overplot, /FOLLOW	  

;************************** V pd ***********************************
 
tvframe,Vpd,/sample, /bar,title='V!I!Mx!N!3 [m/s]',xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale] 

	  
tek_color
VELOVECT, avxi, avzi,xxi/scale,zzi/scale, /OVERPLOT ;, color=6
  
	  
xx(*)=x(1,xstart:xend,1)
zz(*)=x(zstart:zend,1,0)

line2_b_classic, sabx_t,sabz_t,xx, zz, 1, scale ; 0 - plot, 1 - oplot
xx(*)=x(1,xstart:xend,1)
zz(*)=x(zstart:zend,1,0)

tek_color
contour, beta,xx/scale,zz/scale, LEVELS = [0.001,0.01,0.1,1.0, 10.0], $
         C_Annotation = ['1000.0','100.0','10.0','1.0','0.1'], /overplot, /FOLLOW

 ss='time ='+strTrim(string(time),1)
 
 xyouts,20,2, ss, /device, color=200
 

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/REAL/VpVd_P30_R100_A500_B45_Dx4_Dz4_6_1976_400/'+indexss+'.png',image_p, red,green, blue


nn=nn+1



endwhile


end


