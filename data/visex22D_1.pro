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
loadct,4
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

window, 0,xsize=1200,ysize=900,XPOS = 1050, YPOS = 500 

window, 1,xsize=1200,ysize=200,XPOS = 20, YPOS = 80 

window, 2,xsize=600,ysize=585,XPOS = 20, YPOS = 400 

;window, 3,xsize=600,ysize=300,XPOS = 20, YPOS = 720 

nn=0
kkk=4

nn_i=0

close,1
close,2

;openr,1,'/data/ap1vf/2_6MmBBz252.out',/f77_unf
openr,1,'/data/cs1mkg/VAC_NN/testcont2.out',/f77_unf
;openr,1,'/data/ap1vf/2_6Mnzx1023400B41visk03.out',/f77_unf
;openr,1,'/data/ap1vf/2_6Mnzx1976400B41.ini',/f77_unf

;openr,1,'/data/ap1vf/2_6Mnzx1023400B41cont.out',/f77_unf

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
;e2=dblarr(n1,n2)
readu,1,x
;readu,1,e2
;e2=rotate(e2,1)
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=rotate(wi,1)
endfor

if istart eq 1 then TT=dblarr(n2,n1)

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

R=8.3e+003

;mu_gas=0.6d0
mu_gas=1.2d0

gamma=1.66667

xstart=0
xend=511


wset,0
!p.multi = [0,4,4,0,1]
if (ii eq 1) then begin



;xstart=330
;xend=470

zstart=0
;zend=1022

zstart=0
zend=1975

xstart=0
xend=399

;stop
n11=zend-zstart+1
n22=xend-xstart+1

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


tvframe,alog10(sarho_t), /bar,title='log rhob+rho',/sample, xtitle='x [Mm]', ytitle='z [Mm]',charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]

tvframe,savz,/sample, /bar,title='Vz [m/s*rho]', xtitle='x [Mm]', ytitle='z [Mm]',charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]


tvframe,savx,/sample, /bar,title='Vx [m/s*rho]',xtitle='x [Mm]', ytitle='z [Mm]',charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]

tvframe,sae,/bar,/sample, title='e', xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]
	
tvframe,saeb,/bar,/sample, title='eb',  xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]
	
tvframe,sarho,/bar,/sample, title='rho kg/m!E3!N ', xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]

tvframe,sabz*sqrt(mu)*1.0e4,/bar,/sample, title='b_z [Gauss]',  xtitle='x [mm]', ytitle='z [Mm]', charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]
	
tvframe,sabx*sqrt(mu)*1.0e4,/bar,/sample, title='b_x [Gauss]',  xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]

tvframe,sabz_t*sqrt(mu)*1.0e4,/bar,/sample, title='bT_z [Gauss]', xtitle='x [Mm]', ytitle='z [Mm]', $ 
        charsize=2.0,         xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]
	
tvframe,sabx_t*sqrt(mu)*1.0e4,/bar,/sample, title='bT_x [Gauss]', $
         xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
        xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]

endif

Va=sqrt((abz_t^2.d0+abx_t^2.d0)/arho_t)



tvframe, Va, title='Va [m/s]',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample


T=ae+aeb
T=T-(avx^2.0+avz^2.0)*arho_t/2.0

T=T-(abx_t^2.0+abz_t^2.0)/2.0

T=(gamma-1.d0)*T


Cs=sqrt(gamma*T/arho_t)

tvframe, Cs, title='Cs [m/s]',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample

tvframe, sqrt(Cs^2.d0+Va^2.d0), title='sqrt(Cs^2.d0+Va^2.d0) [m/s]',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample


print, '########## min sqrt(Cs^2.d0+Va^2.d0)', min(sqrt(Cs^2.d0+Va^2.d0))


;tvframe, alog10(T),title='log p',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample


beta=(abx_t^2.0+abz_t^2.0)/2.d0/T

contour, beta(xstart:xend, zstart:zend), LEVELS = [0.01,0.1,1.0, 10.0], C_Annotation = ['100.0','10.0','1.0','0.1'] 

;tvframe,beta,/bar,/sample, title='1/beta',  xtitle='x', ytitle='z', charsize=2.0

;plot, zz/scale,alog10(T(kk,*)),title='log(p)', xtitle='x', ytitle='z',charsize=2.0 

T=mu_gas*T/R/arho_t


TTT=T

if (istart eq 1) then begin
   TT=T
   istart=2
;      tvframe, alog10(T(*,*)),title='T [K]',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample,$
;      xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]

   endif else begin
;      tvframe, T(*,*)-TT(*,*),title='delta T [K]',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample,$
;      xrange=[min(xx)/scale, max(xx)/scale],yrange=[min(zz)/scale, max(zz)/scale]
 ;     tvframe, TT(*,*),title='TT',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample  

   endelse



 tvframe, alog10(T(*,*)),title='T [K]',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample
 
tek_color

;plot,w(kk,*,8)+w(kk,*,4),title='bT_z',  xtitle='x', ytitle='z', charsize=2.0, /ys, /xs
 

 
 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200
 



indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
;write_png,'/data/ap1vf/png/REAL/All_P30_R100_A500_B45_Dx4_Dz4_1976_400/'+indexss+'.png',image_p, red,green, blue

;write_png,'/data/cs1mkg/VAC_NN/plots/plot1/'+indexss+'.png',image_p, red,green, blue
print, '******************************', nn

;nn=nn+1

if (ia eq 1) then begin
 tm(0)=time
 mass(0)=total(w(*,*,7)+w(*,*,0)) 
 egas(0)=total(w(*,*,3)+w(*,*,6))
 ia=2.0
endif else begin
 tm=[tm,time]
 mass=[mass,total(w(*,*,7)+w(*,*,0))] 
 egas=[egas,total(w(*,*,3)+w(*,*,6))]
endelse


;goto, jump1 

sst=0
sen=1022

wset,1
!p.multi = [0,2,1,0,1]
;plot,zz/scale,savz(100,150:221),title='mx', xtitle='x', ytitle='y',charsize=1.0  
plot,savx(200,sst:sen),title='vz', xtitle='x', ytitle='y',charsize=1.0  
;oplot,zz/scale, savz(100,150:221),psym=4, color=100
oplot, savx(200,sst:sen),psym=4, color=100

tek_color
plot, zz/scale, alog10(T(kk,*)),title='log(T) tube', xtitle='x',ytitle='z',charsize=1.0, yrange=[3.4,5.6d0]
;oplot, zz/scale,alog10(T(kk,*)), color=200, psym=4


;oplot, alog10(T(kk,*)), psym=3, color=200


;if (ii eq 1) then loadct,4

;tvframe,alog10(T(*,*)), /bar,title='T',xtitle='x',/sample, ytitle='z',charsize=2.0 



Bxx=dblarr(n22,n11)
Bzz=dblarr(n22,n11)


Temp=dblarr(n22,n11)

Temp=TTT(xstart:xend,zstart:zend)

Bxx(*,*)=reform((w(xstart:xend,zstart:zend,5)+w(xstart:xend,zstart:zend,9)))*SQRT(mu)*1.0e4
Bzz(*,*)=reform((w(xstart:xend,zstart:zend,4)+w(xstart:xend,zstart:zend,8)))*SQRT(mu)*1.0e4

btot=sqrt(Bxx^2.0+Bzz^2.0)

wset,2
!P.multi=0

label_bt='B_T'
bb=sqrt(Bzz^2.d0+Bxx^2.d0)
;tvframe,bb*sqrt(mu)*1.0e4,/bar, yrange=[min(zz)/scale,max(zz)/scale],$
;        xrange=[min(xx)/scale,max(xx)/scale], title=label_bt,charsize=ch_size,$  
;	xtitle='Horizontal distance [Mm]', ytitle='Height [Mm]'

tvframe,avz(xstart:xend,zstart:zend), yrange=[min(zz)/scale,max(zz)/scale],$
        xrange=[min(xx)/scale,max(xx)/scale], title='V_z [m/s]',charsize=ch_size,$  
	xtitle='Horizontal distance [Mm]', ytitle='Height [Mm]', /bar, /sample


tek_color
contour, beta,xx/scale,zz/scale, LEVELS = [0.001,0.01,0.1,1.0, 10.0], $
         C_Annotation = ['1000.0','100.0','10.0','1.0','0.1'], /overplot, /FOLLOW

nxz=30
nxnz = [nxz,nxz]

xxi=interpol(xx,nxz)
zzi=interpol(zz,nxz)

avxi=congrid(savx,nxz,nxz)
avzi=congrid(savz,nxz,nxz)

VELOVECT, avxi, avzi,xxi/scale,zzi/scale, /OVERPLOT ;

line2_b_classic, Bxx,Bzz,xx, zz, 1, scale ; 0 - plot, 1 - oplot

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200


indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
;write_png,'/data/ap1vf/png/REAL/Vz_P30_R100_A500_B45_Dx4_Dz2_6_1023_400/'+indexss+'.png',image_p, red,green, blue

write_png,'/data/cs1mkg/VAC_NN/plots/Vz/'+indexss+'.png',image_p, red,green, blue


nn=nn+1

if (ia eq 1) then begin
 tm(0)=time
 mass(0)=total(w(*,*,7)+w(*,*,0)) 
 egas(0)=total(w(*,*,3)+w(*,*,6))
 ia=2.0
endif else begin
 tm=[tm,time]
 mass=[mass,total(w(*,*,7)+w(*,*,0))] 
 egas=[egas,total(w(*,*,3)+w(*,*,6))]
endelse

;plot,mass, charsize=2.0, ystyle=1, title='mass'

;plot,egas, charsize=2.0, ystyle=1, title='egas'

;maxa=[maxa,max(w(*,*,2)/(w(*,*,0)+w(*,*,7)))]
maxa=[maxa,max(Vt)]


jump1 :

endwhile


end
