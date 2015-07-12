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

window, 0,xsize=900,ysize=725,XPOS = 350, YPOS = 300 

window, 1,xsize=1200,ysize=200,XPOS = 20, YPOS = 80 

nn=0
kkk=4

nn_i=0

close,1
close,2

;openr,1,'/data/ap1vf/goriz1024100Bsm.out',/f77_unf
;openr,1,'/data/ap1vf/gorz16300_1024.out',/f77_unf
;openr,1,'/data/ap1vf/8192.ini',/f77_unf 
;openr,1,'/data/ap1vf/4096.out',/f77_unf
;openr,1,'/data/ap1vf/zero_np0104_001.out',/f77_unf
openr,1,'/data/ap1vf/zero_bach.out',/f77_unf
;openr,1,'/data/ap1vf/100.out',/f77_unf


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


Vt=dblarr(n1,n2)
B=dblarr(n1,n2)
B_bg=dblarr(n1,n2)

p=dblarr(n1,n2,1)
;e2=dblarr(n1,n2)






mu=4.0*!PI/1.0e7

print,time

kk=50

label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'


R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=511

Vt(*,*)=sqrt((w(*,*,1)/(w(*,*,0)+w(*,*,7)))^2.0+(w(*,*,2)/(w(*,*,0)+w(*,*,7)))^2.0)

B(*,*)=sqrt((w(*,*,4)*SQRT(mu))^2.0+(w(*,*,5)*SQRT(mu))^2.0)
;B_bg(*,*,*)=sqrt((w1(*,*,*,5)*SQRT(mu))^2.0+(w1(*,*,*,6)*SQRT(mu))^2.0+(w1(*,*,*,7)*SQRT(mu))^2.0)

wset,0
!p.multi = [0,4,4,0,1]
if (ii eq 1) then begin




tvframe,alog10(w(*,*,7)+w(*,*,0)), /bar,title='log rhob+rho',/sample, xtitle='x', ytitle='y',charsize=2.0  
;tvframe,w(*,0:40,1)/(w(*,0:40,0)+w(*,0:40,7)), /bar,title='v1', xtitle='x', ytitle='y',charsize=2.0  

;tvframe,w(*,0:40,1) ;/(w(*,0:40,0)+w(*,0:40,7)), /bar,title='v1', xtitle='x', ytitle='y',charsize=2.0  
tvframe,w(*,*,1)/(w(*,*,0)+w(*,*,7)),/sample, /bar,title='v1', xtitle='x', ytitle='y',charsize=2.0  

;plot,w(50,110:119,1)/(w(50,110:119,0)+w(50,110:119,7)), psym=10, /ys ;, /bar,title='v1', xtitle='x', ytitle='y',charsize=2.0  
;tvframe,w(*,0:40,2)/(w(*,0:40,0)+w(*,0:40,7)),/sample, /bar,title='v2',xtitle='x', ytitle='z',charsize=2.0 

tvframe,w(*,*,2)/(w(*,*,0)+w(*,*,7)),/sample, /bar,title='v2',xtitle='x', ytitle='z',charsize=2.0 

tvframe,w(*,*,3),/sample, title='e', xtitle='x', ytitle='z', charsize=2.0                                                                                                   

tvframe,w(*,*,6),/sample, title='eb',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,0),/bar,/sample, title='rho', xtitle='x', ytitle='z', charsize=2.0

tvframe,w(*,*,4),/bar,/sample, title='b_z',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,5),/bar,/sample, title='b_x',  xtitle='x', ytitle='z', charsize=2.0

tvframe,w(*,*,8)+w(*,*,4),/bar,/sample, title='bT_z',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,9)+w(*,*,5),/bar,/sample, title='bT_x',  xtitle='x', ytitle='z', charsize=2.0

;bbx = reform(w(*,*,2)/(w(*,*,0)+w(*,*,7))) 
;bbz = reform(w(*,*,1)/(w(*,*,0)+w(*,*,7))) 

;bbx = reform(w(*,*,9)+w(*,*,5))
;bbz = reform(w(*,*,8)+w(*,*,4))
bbx = reform(w(*,*,5))
bbz = reform(w(*,*,4))

;line2_b, n2,n1, bbx,bbz,max(x(1,*,1)),max(x(*,1,0))
;bbx=congrid(bbx,100,50)
;bbz=congrid(bbz,100,50)


;velovect, bbx,bbz, length=10

endif else begin

;plot, (w(kk1:kk2,*,1)/(w(kk1:kk2,*,0)+w(kk1:kk2,*,7))),title='Vx', xtitle='x', ytitle='z',charsize=2.0 ;, psym=3
;plot, (w(kk1:kk2,*,2)/(w(kk1:kk2,*,0)+w(kk1:kk2,*,7))),title='Vy', xtitle='x', ytitle='z',charsize=2.0   


;surface,w(*,*,0)

x1=0
x2=255

kk1=97
kk2=97

plot, (w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7)),title='rho+rhoB', xtitle='x', ytitle='y',charsize=2.0, /ys  ,/xs,/ylog
oplot,(w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7)),psym=4,color=3

plot, (w(kk1:kk2,x1:x2,0)),title='rho', xtitle='x', ytitle='y',charsize=2.0, /ys  ,/xs
oplot,(w(kk1:kk2,x1:x2,0)),psym=4,color=7

;(w(kk1:kk2,x1:x2,1)/(w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7)))
plot, w(kk1:kk2,x1:x2,1)/(w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7)),title='v1', xtitle='x', ytitle='y',charsize=2.0 ,/ys ,/xs
oplot,w(kk1:kk2,x1:x2,1)/(w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7)),psym=4,color=2
plot, (w(kk1:kk2,x1:x2,2)/(w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7))), title='v2',xtitle='x', ytitle='z',charsize=2.0 ,/ys,/xs
oplot,(w(kk1:kk2,x1:x2,2)/(w(kk1:kk2,x1:x2,0)+w(kk1:kk2,x1:x2,7))),psym=4,color=6


plot, (w(kk1:kk2,x1:x2,3)),title='e', xtitle='x', ytitle='y',charsize=2.0,/xs ;,psym=3
oplot,(w(kk1:kk2,x1:x2,3)),psym=4,color=6

plot, w(kk1:kk2,x1:x2,6)+w(kk1:kk2,x1:x2,3), title='eb+e', xtitle='x', ytitle='z', charsize=2.0   ,/xs 
oplot, w(kk1:kk2,x1:x2,6)+w(kk1:kk2,x1:x2,3),psym=5,color=5


plot, (w(kk1:kk2,x1:x2,7)), title='rhob', xtitle='x', ytitle='z', charsize=2.0, /ys    ,/xs

plot, w(kk1:kk2,x1:x2,6)+w(kk1:kk2,x1:x2,3), title='e_T', xtitle='x', ytitle='z', charsize=2.0   ,/xs 

endelse




T=dblarr(n2,n1)


T[*,*]=w[*,*,3]+w[*,*,6]


T[*,*]=T[*,*]-(w[*,*,1]^2.0+w[*,*,2]^2.0)/(w[*,*,0]+w[*,*,7])/2.0

T[*,*]=T[*,*]-((w[*,*,4]+w[*,*,8])^2.0+(w[*,*,5]+w[*,*,9])^2.0)/2.0


beta=dblarr(n2,n1)

beta[*,*]=(((w[*,*,4]+w[*,*,8])*sqrt(mu)*1.0e4)^2.0+((w[*,*,5]+w[*,*,9])*sqrt(mu)*1.0e4)^2.0)/2.0/((gamma-1.d0)*T[*,*])

tvframe,beta,/bar,/sample, title='1/beta',  xtitle='x', ytitle='z', charsize=2.0



;VEL, reform(w(*,*,4)), reform(w(*,*,5))

;tek_color
;CONTOUR, beta, LEVELS = [0.1,0.5,1.5], C_ANNOTATION = ['0.1', '0,5','1,5'], /Overplot

;plot, beta[kk,*],title='1/beta',xtitle='x', ytitle='z',charsize=2.0 



T[*,*]=(gamma-1.d0)*T[*,*]
tvframe, (T(*,*)),title='p', xtitle='x', ytitle='z',charsize=2.0,/sample 


T[*,*]=mu_gas*T[*,*]/R/(w[*,*,0]+w[*,*,7])



plot,w(kk,*,8)+w(kk,*,4),title='bT_z',  xtitle='x', ytitle='z', charsize=2.0, /ys, /xs

wset,1
!p.multi = [0,2,1,0,1]
plot,w(kk,*,1)/(w(kk,*,0)+w(kk,*,7)),title='v1', xtitle='x', ytitle='y',charsize=1.0  
oplot,w(kk,*,1)/(w(kk,*,0)+w(kk,*,7)),psym=4, color=100

tek_color
;plot, alog10(T(kk,*)),title='log(T)', xtitle='x',yrange=[3.d0, 6.d0], ytitle='z',charsize=1.0 
plot, alog10(T(kk,*)), psym=10, color=200,yrange=[3.d0, 6.d0]


;if (ii eq 1) then loadct,4

;tvframe,alog10(T(*,*)), /bar,title='T',xtitle='x',/sample, ytitle='z',charsize=2.0 



a=''
;read,a



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

 time=time/31536000.d0
 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200
 
print,  kkk
if (kkk eq 4) then begin

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/test2/'+indexss+'.png',image_p, red,green, blue


nn=nn+1

kkk=1

endif

kkk=kkk+1
;if (ii eq 2) then read,a

;maxa=[maxa,max(w(*,*,2)/(w(*,*,0)+w(*,*,7)))]
maxa=[maxa,max(Vt)]




endwhile


end
