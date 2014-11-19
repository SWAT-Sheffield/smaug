function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
;for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2)
res(1)=res(2)
res(nel-1)=res(nel-3)
res(nel-2)=res(nel-3)
return,res
end

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

ii=1

if (ii eq 1) then begin
loadct,4
mixct
endif else begin
loadct,0
tek_color
endelse


window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
window, 1,xsize=800,ysize=800,XPOS = 500, YPOS = 80

wset,0
!p.multi = [0,4,4,0,1]

print,'******************* time = ', time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

;scale=1.d6
scale=1.0;

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=511
ystart=0
yend=511

pp=255 ;x
kk=5  ;y

wset,0
!p.multi = [0,4,4,0,1]


zstart=0
zend=255

wt=dblarr(zend-zstart+1,xend-xstart+1,iw)
wt=reform(w(zstart:zend,xstart:xend,pp,*))

wy=dblarr(n1,n3,iw)
wy=reform(w(zstart:zend,pp,*,*))

wt(*,*,3)=reform(w(zstart:zend,xstart:xend,pp,3))

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
vt(*,*,*)=sqrt(w(*,*,*,1)^2.d0+w(*,*,*,2)^2.d0+w(*,*,*,3)^2.d0)/(w(*,*,*,0)+w(*,*,*,9))


;****************** Pressure background begin ********************
TP=saeb
TP=TP-(sabx_t^2.0+saby_t^2.0+sabz_t^2.0)/2.0
TP=(gamma-1.d0)*TP
;****************** Pressure background end ********************

if (ii eq 1) then begin


tvframe,rotate(wt(*,*,0),1), /bar,title='rho',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]		
	
st=strTrim(it,1)	

ww=dblarr(xend-xstart+1,zend-zstart+1)
xxs=xx[xstart]
xxe=xx[xend]
zzs=zz[zstart]
zze=zz[zend]
	
;close,10
;openw,10,'/data/ap1vf/data_line_prof/rho/rho_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate((wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10
	

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz',$
        xtitle='x', ytitle='y',charsize=2.0; title='dPdT'


;close,10
;openw,10,'/data/ap1vf/data_line_prof/vz/vz_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10

tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='x', ytitle='z',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]; CT='dPdT'

;close,10
;openw,10,'/data/ap1vf/data_line_prof/vx/vx_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10	

tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'
	
;close,10
;openw,10,'/data/ap1vf/data_line_prof/vy/vy_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10	
	

tvframe,rotate(wy(*,*,2)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'
	
tvframe,rotate(wt(*,*,4),1),/bar, /sample, title='e', xtitle='x', ytitle='z', $
        charsize=2.0

tvframe,rotate(wt(*,*,5),1)*sqrt(mu)*1.0e4,/bar,/sample, title='bz', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,11),1)*sqrt(mu)*1.0e4,/bar,/sample, title='Bx_b', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,12),1)*sqrt(mu)*1.0e4,/bar,/sample, title='By_b', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,8),1),/bar,/sample, title='eb', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,9),1),/bar,/sample, title='rho_b', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,10),1)*sqrt(mu)*1.0e4,/bar,/sample, title='Bz_b', $
        xtitle='x', ytitle='z', charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]	

;stop
T=mu_gas*TP/R/sarho_t
;goto, jump20

tvframe,rotate(T(*,*),1),/bar,/sample, title='T', $
        xtitle='x', ytitle='z', charsize=2.0
	
save, filename='Temp_VALIIIC_old.sav', T
save, filename='z_VALIIIC_old.sav', zz
	
;close,10
;openw,10,'/data/ap1vf/data_line_prof/T/T_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(T(*,*),1)
;writeu,10,ww
;close, 10	

;tvframe,rotate(wt(*,*,12),1)*sqrt(mu)*1.0e4,/bar,/sample, title='By_b', $
;        xtitle='y', ytitle='z', charsize=2.0

Va=dblarr(zend-zstart+1,xend-xstart+1)
Vap=dblarr(zend-zstart+1,xend-xstart+1)

Va(*,*)=sqrt((wt(*,*,10)^2.d0+wt(*,*,11)^2.d0+wt(*,*,12)^2.d0)/wt(*,*,9))*sqrt(mu)*1.0e4

Vap(*,*)=sqrt((wt(*,*,5)^2.d0+wt(*,*,6)^2.d0+wt(*,*,7)^2.d0)/(wt(*,*,9)+wt(*,*,0)))*sqrt(mu)*1.0e4

tvframe, rotate(Va(*,*),1)/1000.d0, title='V!DA!N!3 [km/s]',/bar,/sample, $
        xtitle='y', ytitle='z', charsize=2.0


plot, rotate(Va(10,*),1)/1000.d0  
oplot, rotate(Va(10,*),1)/1000.d0 , psym=4 
oplot, rotate(Vap(10,*),1)/1000.d0
oplot, rotate(Vap(10,*),1)/1000.d0, psym=5


Cs=dblarr(zend-zstart+1,xend-xstart+1)

Cs(*,*)=sqrt(gamma*TP(*,*)/wt(*,*,9))

tvframe, rotate(Cs(*,*)/1000.d0,1),  title='C!DS!N!3 [km/s]',/bar, charsize=2.0,/sample

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200	


 
indexs=strtrim(np,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
;write_png,'/data/ap1vf/png/3D/tube/test_200_puls/all/'+indexss+'.png',image_p, red,green, blue


np=np+1

;goto, jump10
wset,1
!p.multi = [0,3,2,0,1]

tvframe,rotate(wy(*,*,1)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'

tvframe,rotate(wy(*,*,2)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'


tvframe,rotate(wy(*,*,3)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vz', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'



;for hh=0,n1-1 do begin

hh=150

vvt(*,*)=vt(hh,*,*)
cs=1.2

hxmin=20
hymin=20

hxmax=80
hymax=80

wv=reform(w(hh,*,*,*))

savx=dblarr(hxmax-hxmin+1,hymax-hymin+1)
savy=savx
sabx=savx
saby=savx



savx(*,*)=wv(hxmin:hxmax,hymin:hymax,2)/(wv(hxmin:hxmax,hymin:hymax,0)+wv(hxmin:hxmax,hymin:hymax,9))
savy(*,*)=wv(hxmin:hxmax,hymin:hymax,3)/(wv(hxmin:hxmax,hymin:hymax,0)+wv(hxmin:hxmax,hymin:hymax,9))

sabx(*,*)=wv(hxmin:hxmax,hymin:hymax,6)
saby(*,*)=wv(hxmin:hxmax,hymin:hymax,7)

nxy=50
nxny = [nxy,nxy]

hxx=dblarr(hxmax-hxmin+1)
hyy=dblarr(hymax-hymin+1)

hxx=xx(hxmin:hxmax)
hyy=yy(hymin:hymax)

xxi=interpol(hxx,nxy)
yyi=interpol(hyy,nxy)

avxi=congrid(savx,nxy,nxy)
avyi=congrid(savy,nxy,nxy)

abxi=congrid(sabx,nxy,nxy)
abyi=congrid(saby,nxy,nxy)


hight=strTrim(string(hh),1)



jump22:



endif else begin

endelse


;jump1 :



indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

;image_p = TVRD_24()
;write_png,'/data/ap1vf/png/3D/tube/test_200_puls/all/'+indexss+'.png',image_p, red,green, blue


nn=nn+1

;goto, jump

;wset,1
;!p.multi = [0,2,1,0,1]


jump10 :

jump :

goto, jump20

jump20:
stop
;endwhile
;endfor


end
