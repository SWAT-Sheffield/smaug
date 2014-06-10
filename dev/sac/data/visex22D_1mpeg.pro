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
;mpegID = MPEG_OPEN([1100,825],BITRATE=104857200, $
;         MOTION_VEC_LENGTH=1,IFRAME_GAP=10000, QUALITY=100, FILENAME='bx_T300.mpg')
  
nn=0
nni=4


window, 0,xsize=650,ysize=460
;window, 0,xsize=800,ysize=300
nn=0
kkk=10

nn_i=1

close,1
close,2


;openr,1,'/data/ap1vf/gorz1604n.out',/f77_unf
;openr,1,'/data/ap1vf/opozmf.out',/f77_unf
;openr,1,'/data/ap1vf/test.out',/f77_unf
openr,1,'/data/ap1vf/goriz1024100sm.out',/f77_unf
;openr,1,'/data/ap1vf/tube256.out',/f77_unf


zstart=250
zend=350

timeend=2400
timeout=dblarr(timeend,zend-zstart+1)



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


Vt=dblarr(n2,n1)

B=dblarr(n2,80)
B_bg=dblarr(n1,n2)

p=dblarr(n1,n2,1)
;e2=dblarr(n1,n2)




!p.multi = [0,2,3,0,1]
;!p.multi = [0,2,0,0,1]

mu=4.0*!PI/1.0e7

print,time

kk=128

label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'


R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

; gauss Tesla
factor=100000.0

xstart=0
xend=511


label_Vz='Velocity '+'!19V!X!DZ'+'!N'+' (m/s)'
label_Vx='Velocity '+'!19V!X!DX'+'!N'+' (m/s)'
label_Vt='Total velocity '+'!19V!X!DT'+'!N'+' (m/s)'

label_bz='Magnetic field '+'!19b!X!DZ'+'!N'+' (Gauss)'
label_bx='Magnetic field '+'!19b!X!DX'+'!N'+' (Gauss)'
label_bt='Total magnetic field '+'!19B!X!DT'+'!N'+' (Gauss)'

ch_size=1.6

Vt(*,*)=sqrt((w(*,*,1)/(w(*,*,0)+w(*,*,7)))^2.0+(w(*,*,2)/(w(*,*,0)+w(*,*,7)))^2.0)

B(*,*)=sqrt((w(*,80:159,8)+w(*,80:159,4))^2.0+(w(*,80:159,9)*w(*,80:159,5))^2.0)



if (ii eq 1) then begin

tvframe,w(*,0:890,1)/(w(*,0:890,0)+w(*,0:890,7)), /bar,title=label_Vz, /sample, charsize=ch_size,$ 
        xtitle='Horizontal distance (Mm)', xrange=[0,16],ytitle='Height (Mm)', yrange=[0,7.0],$
	YTICKFORMAT='(I6)'
tvframe,w(*,0:890,2)/(w(*,0:890,0)+w(*,0:890,7)), /bar,title=label_Vx,/sample,charsize=ch_size, $
        xtitle='Horizontal distance (Mm)', xrange=[0,16], yrange=[0,7.0], ytitle='Height (Mm)' 

tvframe,Vt(*,0:890),/bar,/sample, title=label_Vt, charsize=ch_size,$
                xtitle='Horizontal distance (Mm)', xrange=[0,16],  ytitle='Height (Mm)', yrange=[0,7.0]

tvframe,w(*,0:890,4)*SQRT(mu)*factor,/bar,/sample, title=label_bz,charsize=ch_size, $
        xtitle='Horizontal distance (Mm)', xrange=[0,16], ytitle='Height (Mm)', yrange=[0,7.0]
tvframe,w(*,0:890,5)*SQRT(mu)*factor,/bar,/sample, title=label_bx,   charsize=ch_size, $
        xtitle='Horizontal distance (Mm)', xrange=[0,16], ytitle='Height (Mm)', yrange=[0, 7.0]




;zstart=0
;zend=255

n11=zend-zstart+1
n22=n2

Bxx=dblarr(n22,n11)
Bzz=dblarr(n22,n11)

xx=max(x(1,*,1))
zz=max(x(zend,1,0))



Bxx(*,*)=reform(w(*,zstart:zend,5)+w(*,zstart:zend,9))
Bzz(*,*)=reform(w(*,zstart:zend,4)+w(*,zstart:zend,8))

;Bxx(*,*)=reform(w(*,zstart:zend,5))
;Bzz(*,*)=reform(w(*,zstart:zend,4))

;Bxx(*,*)=reform(w(*,zstart:zend,2)/(w(*,zstart:zend,0)+w(*,zstart:zend,7)))
;Bzz(*,*)=reform(w(*,zstart:zend,1)/(w(*,zstart:zend,0)+w(*,zstart:zend,7)))


btot=sqrt(Bxx^2.0+Bzz^2.0)*SQRT(mu)*factor


tvframe,btot,$
        /bar,/sample, title=label_bt,charsize=ch_size,$  
	xtitle='Horizontal distance (Mm)', xrange=[0,16], ytitle='Height (Mm)', yrange=[1.96,2.75]
tek_color	
;oplot, [3.57,3.57],[1.9,3.0], linestyle=6, thick=2	

;bbx = reform(w(*,*,5)+w(*,*,9))
;bbz = reform(w(*,*,4)+w(*,*,8))
;bbx = reform(w(*,*,1)/(w(*,*,0)+w(*,*,7)))
;bbz = reform(w(*,*,2)/(w(*,*,0)+w(*,*,7)))
print, n1,n2

tek_color

line2_b, n22,n11, Bxx,Bzz,xx,zz,0.d0,x(zstart,1,0)



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
beta[*,*]=(((w[*,*,4]+w[*,*,8])^2.0+(w[*,*,5]+w[*,*,9])^2.0))/2.0/((gamma-1.d0)*T[*,*])

;tvframe,beta,/bar,/sample, title='beta',  xtitle='x', ytitle='z', charsize=2.0


;plot, beta[kk,*],title='1/beta',xtitle='x', ytitle='z',charsize=2.0 


T[*,*]=mu_gas*(gamma-1.d0)*T[*,*]/R/(w[*,*,0]+w[*,*,7])


;plot, alog10(T(kk,*)),title='T', xtitle='x',yrange=[3.d0, 6.d0], ytitle='z',charsize=2.0 

tek_color

if (ii eq 1) then loadct,4

;tvframe,alog10(T(*,*)), /bar,title='T',xtitle='x',/sample, ytitle='z',charsize=2.0 

 ss0=string(time,format="(F10.1)")

 ss='time ='+strTrim(ss0,1)  ;+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
; s0 = string(format="(1x,a,F.2)"
 ;FORMAT = '("Value: ", I0)', 23
 xyouts,2,2, ss, /device, color=200

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()

;WRITE_PPM, '/data/ap1vf/png/opoz16/'+indexss+'.png',image_p

if (nni eq 4) then begin

;  WRITE_PPM, '/data/ap1vf/png/finalppm/'+indexss+'.ppm',image_p
  write_png,'/data/ap1vf/png/final/'+indexss+'.png',image_p, red,green, blue
  
  nn=nn+1
  nni=0
endif 

nni=nni+1
nn_i=nn_i+1

;endif
endwhile
; Save the MPEG sequence in the file myMovie.mpg: 
;MPEG_SAVE, mpegID 
 
; Close the MPEG sequence: 
;MPEG_CLOSE, mpegID 



end
