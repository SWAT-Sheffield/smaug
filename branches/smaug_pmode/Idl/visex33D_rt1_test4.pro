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
mixct
endif else begin
loadct,0
tek_color
endelse

n1=128
wdt1v=dblarr(373,n1,4)
window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
window, 1,xsize=800,ysize=800,XPOS = 500, YPOS = 80

n1=128
n2=128
wdt1v=dblarr(281,n1,4)
wdt1h=dblarr(281,n2,4,2)
window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
window, 1,xsize=800,ysize=800,XPOS = 500, YPOS = 80

;for ipic=pic,489000 do begin
for ipic=280,489 do begin
;for ipic=1,371 do begin
;for ipic=1,281 do begin
;for ipic=1,3 do begin

;while not(eof(1)) do begin

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




nn=0
np=0
kkk=4

nn_i=0

close,1
close,2




;openr,1,'/data/ap1vf/3D_509_36_36_300s.out',/f77_unf
;openr,1,'/data/ap1vf/3D_tube_modif_200_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_vert_driver/3D_tube_196_100_100_phi.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_vertical/3D_tube_200_100_100_vertical.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_multidriver_lower.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100t.out',/f77_unf



;openr,1,'/data/ap1vf/vxx.040000',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60t.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_196_100_100_test.out',/f77_unf

;openr,1,'/data/ap1vf/3D_vert_driver/3D_tube_196_100_100_vert.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_modif_200_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_200s_puls.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf

;openr,1,'/home/mikeg/proj/sac2.5d-cuda/test_OT.out'
;directory='/home/mikeg/proj/sac2.5d-cuda/out_OT_withhyper/'
;directory='../out/'
;directory='/nobackup/shemkg/modes3/sp2rc/3demt1/'
;directory='/nobackup/shemkg/modes3/f1/'
;directory='/fastdata/cs1mkg/smaug/spicule4b0_3d/'
;directory='/data/cs1mkg/smaug_spicule1/spicule5b0_3d/'
;directory='/fastdata/cs1mkg/smaug/spicule5b0_3d/'
directory='/data/cs1mkg/smaug_spicule1/spicule5b0_3d/'
;directory='/data/cs1mkg/smaug_spicule1/spicule4b0_b1_3d/'

;directory='/nobackup/shemkg/modes/'
;pic=999
;name='zeroOT_'
;name='3D_em_t1_bin_np010808_00'
;name='3D_em_f1_bin_'
name='zerospic1__'
;ndim=2
;n1=800
;n2=6


;picid=ipic*5+4
picid=ipic*1000L
;picid=ipic;
outfile=directory+name+strtrim(string(picid),2)+'.out'
print,'ipic=',ipic
;openr,1,outfile,/f77_unf
openr,1,outfile

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

;ndim=3
;n1=256
;n2=512
;n3=512
;nw=13

x=dblarr(n1,n2,n3,ndim)

print,'about to allocate'
print,n1,n2,n3,nw

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

dt=0.001

myrange=[x(0,0,0)/1.0d6,x(n1-1,0,0)/1.0d6]
mxrange=[x(0,0,1)/1.0d6,x(0,n2-1,1)/1.0d6]
mtrange=[0,picid*dt]
mv1range=[-0.008,0.008]
mv2range=[-0.008,0.008]




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

;scale=1.d6
scale=1.0;

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=127
ystart=0
yend=127

pp=63 ;x
kk=5  ;y



zstart=0
zend=127

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

nsect1=63
wdt1v(ipic,*,0)=wt(nsect1,*,0)
wdt1v(ipic,*,1)=wt(nsect1,*,1)/(wt(nsect1,*,0)+wt(nsect1,*,9))
wdt1v(ipic,*,2)=wt(nsect1,*,2)/(wt(nsect1,*,0)+wt(nsect1,*,9))
wdt1v(ipic,*,3)=wt(nsect1,*,3)


wdt1h(ipic,*,0,0)=wt( *,nsect1,0)
wdt1h(ipic,*,1,0)=wt( *,nsect1,1)/(wt( *,nsect1,0)+wt( *,nsect1,9))
wdt1h(ipic,*,2,0)=wt( *,nsect1,2)/(wt( *,nsect1,0)+wt( *,nsect1,9))
wdt1h(ipic,*,3,0)=wt( *,nsect1,3)


nsect1=20
wdt1h(ipic,*,0,1)=wt( *,nsect1,0)
wdt1h(ipic,*,1,1)=wt( *,nsect1,1)/(wt( *,nsect1,0)+wt( *,nsect1,9))
wdt1h(ipic,*,2,1)=wt( *,nsect1,2)/(wt( *,nsect1,0)+wt( *,nsect1,9))
wdt1h(ipic,*,3,1)=wt( *,nsect1,3)

;****************** Pressure background begin ********************
TP=saeb
TP=TP-(sabx_t^2.0+saby_t^2.0+sabz_t^2.0)/2.0
TP=(gamma-1.d0)*TP
;****************** Pressure background end ********************


goto, jump30

wset,0
!p.multi = [0,4,4,0,1]

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

tvframe,rotate(wt(*,*,6),1)*sqrt(mu)*1.0e4,/bar,/sample, title='bx', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,7),1)*sqrt(mu)*1.0e4,/bar,/sample, title='by', $
        xtitle='x', ytitle='z', charsize=2.0




;tvframe,rotate(wt(*,*,11),1)*sqrt(mu)*1.0e4,/bar,/sample, title='Bx_b', $
;        xtitle='x', ytitle='z', charsize=2.0

;tvframe,rotate(wt(*,*,12),1)*sqrt(mu)*1.0e4,/bar,/sample, title='By_b', $
 ;       xtitle='x', ytitle='z', charsize=2.0

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
jump30:
wset,1
!p.multi = [0,3,2,0,1]

;;;tvframe,rotate(wy(*,*,1)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vx', $
 ;;;       xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

;;;tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
 ;;;       xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'

;;;tvframe,rotate(wy(*,*,2)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vy', $
 ;;;       xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

;;;tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
 ;;;       xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'


;;;tvframe,rotate(wy(*,*,3)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vz', $
 ;;;       xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

;;;tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz', $
 ;;;       xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'





tvframe,rotate(wy(*,*,1)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='y', ytitle='z',charsize=2.0, xrange=mxrange, yrange=myrange, brange=mv1range  
   ;, title='dPdT'

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='x', ytitle='z',charsize=2.0, xrange=mxrange, yrange=myrange, brange=mv2range  
; CT='dPdT'

tvframe,rotate(wy(*,*,2)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='y', ytitle='z',charsize=2.0, xrange=mxrange, yrange=myrange, brange=mv1range  
;, title='dPdT'

tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='x', ytitle='z',charsize=2.0, xrange=mxrange, yrange=myrange, brange=mv2range  
; CT='dPdT'


tvframe,rotate(wy(*,*,3)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vz', $
        xtitle='y', ytitle='z',charsize=2.0, xrange=mxrange, yrange=myrange, brange=mv1range  
;, title='dPdT'

tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz', $
        xtitle='x', ytitle='z',charsize=2.0, xrange=mxrange, yrange=myrange, brange=mv2range  
; CT='dPdT'




 
;stop

jump22:

	
;wait, 0.5
;endfor

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200


endif else begin

endelse


T=dblarr(n1,n2,n3)


T[*,*,*]=(w[*,*,*,4]+w[*,*,*,8])


T[*,*,*]=T[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/(w[*,*,*,0]+w[*,*,*,9])/2.0

T[*,*,*]=T[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0+(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.d0


beta=dblarr(n1,n2,n3)

beta[*,*,*]=(((w[*,*,*,5]+w[*,*,*,10])*sqrt(mu)*1.0e4)^2.0+((w[*,*,*,6]+w[*,*,*,11])*sqrt(mu)*1.0e4)^2.0+$
              ((w[*,*,*,7]+w[*,*,*,12])*sqrt(mu)*1.0e4)^2.0)/2.0/((gamma-1.d0)*T[*,*,*])

;tvframe,rotate(beta(*,*,pp),1),/bar,/sample, title='1/beta',  xtitle='x', ytitle='z', charsize=2.0



;plot, alog10(1.d0/beta[*,kk,pp]),title='1/beta',xtitle='x', ytitle='z',charsize=2.0 

T[*,*,*]=(gamma-1.d0)*T[*,*,*]


C=dblarr(n1,n2,n3)

C[*,*,*]=sqrt(gamma*T[*,*,*]/(reform(w[*,*,*,0]+w[*,*,*,9])))



indexs=strtrim(picid,2)
;indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='00000000'+indexs                                             
 2:indexss='0000000'+indexs                                              
 3:indexss='000000'+indexs                                               
 4:indexss='00000'+indexs
 5:indexss='0000'+indexs                                             
 6:indexss='000'+indexs                                              
 7:indexss='00'+indexs                                               
 8:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/cs1mkg/smaug_spicule1/Idl/images/spicule5b0_3d/s5b0_3d'+indexss+'.png',image_p, red,green, blue


nn=nn+1



jump10 :

jump :


;stop
;endwhile
endfor


tvframe,wdt1v(*,*,1), /bar,title='v1',/sample, xtitle='dt (sec)', ytitle='Height (Mm)',charsize=1.1, xrange=mtrange, yrange=myrange 
tvframe,wdt1v(*,*,2), /bar,title='v2',/sample, xtitle='dt (sec)', ytitle='Height (Mm)',charsize=1.1, xrange=mtrange, yrange=myrange 




image_p = TVRD_24()
write_png,'/data/cs1mkg/smaug_spicule1/Idl/images/spicule5b0_3d/pm300s_sindrv_n0_5b03d_vdt.png',image_p, red,green, blue


save,wdt1har,filename='wdt1var_5b03d.dat'


tvframe,wdt1h(*,*,1,0), /bar,title='v1',/sample, xtitle='dt (sec)', ytitle='Height (Mm)',charsize=1.1, xrange=mtrange, yrange=myrange 
tvframe,wdt1h(*,*,2,0), /bar,title='v2',/sample, xtitle='dt (sec)', ytitle='Height (Mm)',charsize=1.1, xrange=mtrange, yrange=myrange 
write_png,'/data/cs1mkg/smaug_spicule1/Idl/images/spicule5b0_3d/pm300s_sindrv_n0_5b03d_h0dt.png',image_p, red,green, blue




image_p = TVRD_24()





tvframe,wdt1h(*,*,1,1), /bar,title='v1',/sample, xtitle='dt (sec)', ytitle='Height (Mm)',charsize=1.1, xrange=mtrange, yrange=myrange 
tvframe,wdt1h(*,*,2,1), /bar,title='v2',/sample, xtitle='dt (sec)', ytitle='Height (Mm)',charsize=1.1, xrange=mtrange, yrange=myrange 
write_png,'/data/cs1mkg/smaug_spicule1/Idl/images/spicule5b0_3d/pm300s_sindrv_n0_5b03d_h1dt.png',image_p, red,green, blue




image_p = TVRD_24()


save,wdt1har,filename='wdt1har_5b03d.dat'






end
