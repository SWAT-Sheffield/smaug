tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


ii=2

if (ii eq 1) then begin
loadct,4
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
;mpegID = MPEG_OPEN([900,900],FILENAME='myMovie.mpg') 

nn=0


nn_i=0

close,1
close,2


openr,1,'/data/ap1vf/zeroBW.out',/f77_unf


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



;!p.multi = [0,1,0,0,1]
!p.multi = [0,4,4,0,1]
;if nn eq 1 then stop

;read, aa




mu=4.0*!PI/1.0e7

print,time

kk=64

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



if (ii eq 1) then begin

;tvframe,reform(w(*,*,1)/(w(*,*,0)+w(*,*,7))),/bar,/sample, title='Vx', xtitle='x', ytitle='z',charsize=2.0    
;tvframe,reform(w(*,*,2)/(w(*,*,0)+w(*,*,7))),/bar,/sample, title='Vy', xtitle='x', ytitle='z',charsize=2.0   



tvframe,w(*,*,7)+w(*,*,0), /bar,title='log rho_b',/sample, xtitle='x', ytitle='y',charsize=2.0  
tvframe,w(*,*,1)/(w(*,*,7)+w(*,*,0)), /bar,title='v1',/sample, xtitle='x', ytitle='y',charsize=2.0  
tvframe,w(*,*,2)/(w(*,*,7)+w(*,*,0)), /bar,title='mv2',xtitle='x',/sample, ytitle='z',charsize=2.0 

tvframe,w(*,*,3),/bar,/sample, title='e', xtitle='x', ytitle='z', charsize=2.0                                                                                                   

tvframe,w(*,*,6),/bar,/sample, title='eb',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,0),/bar,/sample, title='rho',  xtitle='x', ytitle='z', charsize=2.0

tvframe,w(*,*,4),/bar,/sample, title='b_z',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,5),/bar,/sample, title='b_x',  xtitle='x', ytitle='z', charsize=2.0

tvframe,w(*,*,8),/bar,/sample, title='bg_z',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,9),/bar,/sample, title='bg_x',  xtitle='x', ytitle='z', charsize=2.0



endif else begin

;plot, (w(kk1:kk2,*,1)/(w(kk1:kk2,*,0)+w(kk1:kk2,*,7))),title='Vx', xtitle='x', ytitle='z',charsize=2.0 ;, psym=3
;plot, (w(kk1:kk2,*,2)/(w(kk1:kk2,*,0)+w(kk1:kk2,*,7))),title='Vy', xtitle='x', ytitle='z',charsize=2.0   


;surface,w(*,*,0)

x1=0
x2=799

kk1=4
kk2=4

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


plot, w(kk1:kk2,x1:x2,5),title='By', xtitle='x', ytitle='y',charsize=2.0, /ys  ,/xs
oplot,w(kk1:kk2,x1:x2,5),psym=4,color=3


endelse




T=dblarr(n2,n1)


T[*,*]=w[*,*,3]+w[*,*,6]


T[*,*]=T[*,*]-(w[*,*,1]^2.0+w[*,*,2]^2.0)/(w[*,*,0]+w[*,*,7])/2.0

T[*,*]=T[*,*]-(w[*,*,4]^2.0+w[*,*,5]^2.0)/2.0


B=dblarr(n2,n1)

B[*,*]=sqrt(((w[*,*,8]+w[*,*,4])^2.0+(w[*,*,9]+w[*,*,5])^2.0))

tvframe,B,/bar,/sample, title='B',  xtitle='x', ytitle='z', charsize=2.0


;plot, beta[kk,*],title='1/beta',xtitle='x', ytitle='z',charsize=2.0 


;T[*,*]=mu_gas*(gamma-1.d0)*T[*,*]/R/(w[*,*,0]+w[*,*,7])


;plot, alog10(T(kk,*)),title='T', xtitle='x', ytitle='z',charsize=2.0 

tek_color

if (ii eq 1) then loadct,4

;tvframe,reform(T(*,*)), /bar,title='T',xtitle='x',/sample, ytitle='z',charsize=2.0 



a=''

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

plot,mass, charsize=2.0, ystyle=1, title='mass'
plot,egas, charsize=2.0, ystyle=1, title='egas'

indexs=strtrim(nn,2)

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200
;if (ii eq 2) then read,a

;maxa=[maxa,max(w(*,*,2)/(w(*,*,0)+w(*,*,7)))]
maxa=[maxa,max(Vt)]


endwhile


end
