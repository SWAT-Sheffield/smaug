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

window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
;window, 1,xsize=800,ysize=800,XPOS = 500, YPOS = 80
window, 1,xsize=1025,ysize=1025,XPOS = 950, YPOS = 200 


nn=0
np=0
kkk=4

nn_i=0

close,1
close,2




;openr,1,'/data/ap1vf/3D_509_36_36_300s.out',/f77_unf
;openr,1,'/data/ap1vf/3D_tube_modif_200_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_vert_driver/3D_tube_196_100_100_phi.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_vertical/3D_tube_modif_200_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_multidriver_lower.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_150_100_100_vigeesh_sss.ini',/f77_unf



;openr,1,'/data/ap1vf/vxx.040000',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60t.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_196_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_vert_driver/3D_tube_196_100_100_vert.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100.ini',/f77_unf



;openr,1,'/data/ap1vf/3D/torsional_driver_puls_long/3D_tube_196_100_100.out',/f77_unf

;openr,1,'/fastdata/cs1mkg/VAC_NN_tests/3D_tube_128_128_128.ini',/f77_unf
;openr,1,'/fastdata/cs1mkg/VAC_NN_tests/3D_tube_128_128_128_nodrv_nobge.out',/f77_unf
openr,1,'/fastdata/cs1mkg/VAC_NN_tests/3D_tubet1_128_128_128.out',/f77_unf
;openr,1,'/data/ap1vf/3D/torsional_driver/3D_tube_196_100_100_puls_one.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_200s_puls.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf

for pic=0,14 do begin
;while not(eof(1)) do begin
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

;stop
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

;xstart=0
;xend=99
;ystart=0
;yend=99

xstart=0
xend=127
ystart=0
yend=127

pp=49 ;x
pp=64 ;x

kk=5  ;y

wset,0
!p.multi = [0,4,4,0,1]


zstart=0
zend=127

;zend=150

wt=dblarr(zend-zstart+1,xend-xstart+1,iw)
;wt=reform(w(zstart:zend,xstart:xend,pp,*))

wy=dblarr(n1,n3,iw)
wy=reform(w(zstart:zend,pp,*,*))

;wt(*,*,3)=reform(w(zstart:zend,xstart:xend,pp,3))

;wt(*,*,12)=reform(w(zstart:zend,pp,ystart:yend,12))
;wt(*,*,11)=reform(w(zstart:zend,pp,ystart:yend,11))


for iw=0,nw-1 do begin 
 wt(*,*,iw)=reform(w(zstart:zend,pp,ystart:yend,iw))
endfor

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
	

	

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz',$
        xtitle='x', ytitle='y',charsize=2.0;, CT='dPdT'





tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='x', ytitle='z',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]




tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='x', ytitle='z',charsize=2.0;, CT='dPdT'
	

	

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

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200	

image_p = TVRD_24()

;write_png,vz_slice+indexss+'.png',image_p, red,green, blue
;write_png,'tubeactslice'+indexss+'.png',image_p, red,green, blue
np=np+1

;goto, jump10
wset,1
!p.multi = [0,4,4,0,1]


;for hh=0,n1-1 do begin

hh=10

vvt(*,*)=vt(hh,*,*)
cs=1.2

hxmin=20
hymin=20

hxmax=80
hymax=80

wv=reform(w(hh,*,*,*))

savx=dblarr(hxmax-hxmin+1,hymax-hymin+1)
savy=savx




jump22:

 
indexs=strtrim(np,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

;image_p = TVRD_24()

;write_png,vz_slice+indexss+'.png',image_p, red,green, blue
;write_png,'tubeslice'+indexss+'.png',image_p, red,green, blue


endif else begin

endelse


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



jump10 :

jump :


;endwhile
endfor

close,1
close,2


;****************************************************************************************************
;display GPU generated

;openr,1,'/data/ap1vf/3D_tube_196_100_100_200s_puls.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf
;openr,1,'/data/cs1mkg/sac_cuda/out_withdrv_withbge/3D_tube_128_128_128_5.out'
;openr,1,'/data/cs1mkg/sac_cuda/out_nodrv_withbge/3D_tube_128_128_128_10.out'
;openr,1,'/data/cs1mkg/sac_cuda/out_nodrv_nobge/3D_tube_128_128_128_80.out'
;openr,1,'/data/cs1mkg/sac_cuda/out/3D_tube_128_128_128_10.out'
;openr,1,'/data/cs1mkg/sac_cuda/out_nodrv_nobge_nograv_nohyp/3D_tube_128_128_128_10.out'
openr,1,'/fastdata/cs1mkg/sac_cuda/out_driver_hyp_tube/3D_atubet1_128_128_128_150.out'
;;;;while not(eof(1)) do begin
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
time=time+6.795570135

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

;stop
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

;xstart=0
;xend=99
;ystart=0
;yend=99

xstart=0
xend=127
ystart=0
yend=127

pp=49 ;x
pp=64 ;x

kk=5  ;y


wset,1
!p.multi = [0,4,4,0,1]

zstart=0
zend=127

;zend=150

gwt=dblarr(zend-zstart+1,xend-xstart+1,iw)
gwt=reform(w(zstart:zend,xstart:xend,pp,*))

dwt=dblarr(zend-zstart+1,xend-xstart+1,iw)
dwt=reform(w(zstart:zend,xstart:xend,pp,*))

wy=dblarr(n1,n3,iw)
wy=reform(w(zstart:zend,pp,*,*))

;gwt(*,*,3)=reform(w(zstart:zend,xstart:xend,pp,3))

;gwt(*,*,12)=reform(w(zstart:zend,pp,ystart:yend,12))

;dwt(*,*,i)=(gwt(*,*,*)-wt(*,*,*))*(gwt(*,*,*)-wt(*,*,*))

for iw=0,nw-1 do begin 
 gwt(*,*,iw)=reform(w(zstart:zend,pp,ystart:yend,iw))
 dwt(*,*,iw)=sqrt((gwt(*,*,iw)-wt(*,*,iw))*(gwt(*,*,iw)-wt(*,*,iw)))
endfor

saeb=dblarr(zend-zstart+1,xend-xstart+1)
sarho_t=dblarr(zend-zstart+1,xend-xstart+1)
sabz_t=dblarr(zend-zstart+1,xend-xstart+1)
sabx_t=dblarr(zend-zstart+1,xend-xstart+1)
saby_t=dblarr(zend-zstart+1,xend-xstart+1)

saeb(*,*)=gwt(*,*,8)
sarho_t(*,*)=gwt(*,*,0)+gwt(*,*,9)
sabz_t(*,*)=gwt(*,*,10)
sabx_t(*,*)=gwt(*,*,11)
saby_t(*,*)=gwt(*,*,12)

vt=dblarr(n1,n2,n3)
vvt=dblarr(n2,n3)
vt(*,*,*)=sqrt(w(*,*,*,1)^2.d0+w(*,*,*,2)^2.d0+w(*,*,*,3)^2.d0)/(w(*,*,*,0)+w(*,*,*,9))


tvframe,rotate(dwt(*,*,0),1), /bar,title='rho',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]	

tvframe,rotate(dwt(*,*,1),1), /bar,title='Vz',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(dwt(*,*,4),1), /bar,title='e',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]


tvframe,rotate(dwt(*,*,2),1), /bar,title='Vx',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(dwt(*,*,3),1), /bar,title='Vy',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]


tvframe,rotate(dwt(*,*,5),1), /bar,title='bz',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]


tvframe,rotate(dwt(*,*,7),1), /bar,title='bx',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(dwt(*,*,6),1), /bar,title='by',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

image_p = TVRD_24()

;write_png,vz_slice+indexss+'.png',image_p, red,green, blue
write_png,'diftubeactslice'+indexss+'.png',image_p, red,green, blue
	

end
