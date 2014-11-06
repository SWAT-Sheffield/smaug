tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)



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

window, 0,xsize=800,ysize=500,XPOS = 950, YPOS = 300 
window, 1,xsize=800,ysize=450,XPOS = 500, YPOS = 80

window, 2,xsize=800,ysize=500,XPOS = 1050, YPOS = 300 
window, 3,xsize=800,ysize=450,XPOS = 800, YPOS = 80



nn=0
kkk=4

nn_i=0

close,1
close,2


;openr,1,'/data/ap1vf/3D_509_36_36_300s.out',/f77_unf
openr,1,'/data/ap1vf/3D_tube_196_100_100_120s_full.out',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60t.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_120s_full.out',/f77_unf
;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf

while not(eof(1)) do begin
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
pm=50 ;x
kk=5  ;y

wset,0
!p.multi = [0,3,2,0,1]


zstart=0
zend=195

wt=dblarr(zend-zstart+1,xend-xstart+1,iw)
wtm=dblarr(zend-zstart+1,xend-xstart+1,iw)


wt=reform(w(zstart:zend,xstart:xend,pp,*))
wtm=reform(w(zstart:zend,xstart:xend,pm,*))

wy=dblarr(n1,n3,iw)
wy=reform(w(zstart:zend,pp,*,*))

wt(*,*,3)=reform(w(zstart:zend,xstart:xend,pp+2,3))

wt(*,*,12)=reform(w(zstart:zend,pp,ystart:yend,12))

saeb=dblarr(zend-zstart+1,xend-xstart+1)
sabz_t=dblarr(zend-zstart+1,xend-xstart+1)
sabx_t=dblarr(zend-zstart+1,xend-xstart+1)
saby_t=dblarr(zend-zstart+1,xend-xstart+1)

saeb(*,*)=wt(*,*,8)
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


; *************** VxVyVzbxbybz

tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx [m/s]', $
        xtitle='x [Mm]', ytitle='z [Mm]',charsize=2.0, CT='dPdT', $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(wtm(*,*,6),1)*sqrt(mu)*1.0e4,/bar,/sample, title='bx [Gauss]', $
        xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy [m/s]', $
        xtitle='x [Mm]', ytitle='z [Mm]',charsize=2.0, CT='dPdT', $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(wtm(*,*,7),1)*sqrt(mu)*1.0e4,/bar,/sample, title='by [Gauss]', $
        xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz [m/s]',$
        xtitle='x [Mm]', ytitle='z [Mm]',charsize=2.0, CT='dPdT', $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]		

tvframe,rotate(wtm(*,*,5),1)*sqrt(mu)*1.0e4,/bar,/sample, title='bz [Gauss]', $
        xtitle='x [Mm]', ytitle='z [Mm]', charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]

	

ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)
 xyouts,20,5, ss, /device, color=200	


 
indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/slice/VxVyVzbxbybz/'+indexss+'.png',image_p, red,green, blue


;*********************************** VrVpVzbrbpbz

wset,1
!p.multi = [0,3,2,0,1]

;for hh=0,n1-1 do begin

hh=12

vvt(*,*)=vt(hh,*,*)
cs=2

hxmin=20
hymin=20

hxmax=80
hymax=80

wv=reform(w(hh,*,*,*))


hight=strTrim(string(FORMAT='(6F10.2)',x(hh,1,1,0)/1.0d6),1)+' [Mm]'




r=dblarr(n2,n3)
vr=dblarr(n2,n3)
vphi=dblarr(n2,n3)

br=dblarr(n2,n3)
bphi=dblarr(n2,n3)


for i=0,n2-1 do begin
 for j=0,n3-1 do begin
 r[i,j]=sqrt((xx[i]-xx[n2/2])^2.d0+(yy[j]-yy[n3/2])^2.d0)

 vr[i,j]= wv(i,j,2)/(wv(i,j,0)+wv(i,j,9))*(xx[i]-xx[n2/2])/r[i,j]+$
          wv(i,j,3)/(wv(i,j,0)+wv(i,j,9))*(yy[j]-yy[n3/2])/r[i,j]
  
 vphi[i,j]=-wv(i,j,2)/(wv(i,j,0)+wv(i,j,9))*(yy[j]-yy[n3/2])/r[i,j]+$
          wv(i,j,3)/(wv(i,j,0)+wv(i,j,9))*(xx[i]-xx[n2/2])/r[i,j]

 br[i,j]= wv(i,j,6)*(xx[i]-xx[n2/2])/r[i,j]+$
          wv(i,j,7)*(yy[j]-yy[n3/2])/r[i,j]
  
 bphi[i,j]=-wv(i,j,6)*(yy[j]-yy[n3/2])/r[i,j]+$
          wv(i,j,7)*(xx[i]-xx[n2/2])/r[i,j]
 br[n2/2,n3/2]=0.d0	  
 bphi[n2/2,n3/2]=0.d0
 endfor
 
endfor


	

tvframe,vr(hxmin:hxmax,hymin:hymax),/bar,/sample, title='Vr [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,br(hxmin:hxmax,hymin:hymax)*sqrt(mu)*1.0e4,/bar,/sample, title='br [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	
	

tvframe,vphi(hxmin:hxmax,hymin:hymax),/bar,/sample, title='Vphi [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	


tvframe,bphi(hxmin:hxmax,hymin:hymax)*sqrt(mu)*1.0e4,/bar,/sample, title='bphi [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	
	
	
tvframe,wv(hxmin:hxmax,hymin:hymax,1)/(wv(hxmin:hxmax,hymin:hymax,0)+wv(hxmin:hxmax,hymin:hymax,9)),$
        /bar,/sample, title='Vz [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,wv(hxmin:hxmax,hymin:hymax,5)*sqrt(mu)*1.0e4,/bar,/sample, $
        title='bz [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	



ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)
 xyouts,20,5, ss, /device, color=200	



indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/slice/VrVpVzbrbpbz/'+indexss+'.png',image_p, red,green, blue


;*********************************** VxVyVzVrVpbz

wset,2
!p.multi = [0,3,2,0,1]

;for hh=0,n1-1 do begin


vvt(*,*)=vt(hh,*,*)
cs=2


wv=reform(w(hh,*,*,*))


hight=strTrim(string(FORMAT='(6F10.2)',x(hh,1,1,0)/1.0d6),1)+' [Mm]'




r=dblarr(n2,n3)
vr=dblarr(n2,n3)
vphi=dblarr(n2,n3)

br=dblarr(n2,n3)
bphi=dblarr(n2,n3)


for i=0,n2-1 do begin
 for j=0,n3-1 do begin
 r[i,j]=sqrt((xx[i]-xx[n2/2])^2.d0+(yy[j]-yy[n3/2])^2.d0)

 vr[i,j]= wv(i,j,2)/(wv(i,j,0)+wv(i,j,9))*(xx[i]-xx[n2/2])/r[i,j]+$
          wv(i,j,3)/(wv(i,j,0)+wv(i,j,9))*(yy[j]-yy[n3/2])/r[i,j]
  
 vphi[i,j]=-wv(i,j,2)/(wv(i,j,0)+wv(i,j,9))*(yy[j]-yy[n3/2])/r[i,j]+$
          wv(i,j,3)/(wv(i,j,0)+wv(i,j,9))*(xx[i]-xx[n2/2])/r[i,j]

 endfor
 
endfor


	
tvframe,wv(hxmin:hxmax,hymin:hymax,2)/(wv(hxmin:hxmax,hymin:hymax,0)+wv(hxmin:hxmax,hymin:hymax,9)),$
        /bar,/sample, title='Vx [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,vr(hxmin:hxmax,hymin:hymax),/bar,/sample, title='Vr [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,wv(hxmin:hxmax,hymin:hymax,3)/(wv(hxmin:hxmax,hymin:hymax,0)+wv(hxmin:hxmax,hymin:hymax,9)),$
        /bar,/sample, title='Vy [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	
	
tvframe,vphi(hxmin:hxmax,hymin:hymax),/bar,/sample, title='Vphi [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

	
tvframe,wv(hxmin:hxmax,hymin:hymax,1)/(wv(hxmin:hxmax,hymin:hymax,0)+wv(hxmin:hxmax,hymin:hymax,9)),$
        /bar,/sample, title='Vz [m/s], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,wv(hxmin:hxmax,hymin:hymax,5)*sqrt(mu)*1.0e4,/bar,/sample, $
        title='bz [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	
	
	



ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)
 xyouts,20,5, ss, /device, color=200	



indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/slice/VxVyVzVrVpbz/'+indexss+'.png',image_p, red,green, blue


;*********************************** bxbybzbrbpbz

wset,3
!p.multi = [0,3,2,0,1]

;for hh=0,n1-1 do begin

vvt(*,*)=vt(hh,*,*)
cs=2


wv=reform(w(hh,*,*,*))


hight=strTrim(string(FORMAT='(6F10.2)',x(hh,1,1,0)/1.0d6),1)+' [Mm]'




r=dblarr(n2,n3)
vr=dblarr(n2,n3)
vphi=dblarr(n2,n3)

br=dblarr(n2,n3)
bphi=dblarr(n2,n3)


for i=0,n2-1 do begin
 for j=0,n3-1 do begin
 r[i,j]=sqrt((xx[i]-xx[n2/2])^2.d0+(yy[j]-yy[n3/2])^2.d0)

 br[i,j]= wv(i,j,6)*(xx[i]-xx[n2/2])/r[i,j]+$
          wv(i,j,7)*(yy[j]-yy[n3/2])/r[i,j]
  
 bphi[i,j]=-wv(i,j,6)*(yy[j]-yy[n3/2])/r[i,j]+$
          wv(i,j,7)*(xx[i]-xx[n2/2])/r[i,j]
 br[n2/2,n3/2]=0.d0	  
 bphi[n2/2,n3/2]=0.d0
 endfor
 
endfor


	

tvframe,wv(hxmin:hxmax,hymin:hymax,6)*sqrt(mu)*1.0e4,/bar,/sample, title='bx [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,br(hxmin:hxmax,hymin:hymax)*sqrt(mu)*1.0e4,/bar,/sample, title='br [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	
	

tvframe,wv(hxmin:hxmax,hymin:hymax,7)*sqrt(mu)*1.0e4,/bar,/sample, title='by [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	


tvframe,bphi(hxmin:hxmax,hymin:hymax)*sqrt(mu)*1.0e4,/bar,/sample, title='bphi [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	
	
	
tvframe,wv(hxmin:hxmax,hymin:hymax,5)*sqrt(mu)*1.0e4,/bar,/sample, $
        title='bz [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	

tvframe,wv(hxmin:hxmax,hymin:hymax,5)*sqrt(mu)*1.0e4,/bar,/sample, $
        title='bz [Gauss], h='+hight, $
        xtitle='x [Mm]', ytitle='y [Mm]', charsize=cs, $
	xrange=[xx(hxmin)/scale, xx(hxmax)/scale], yrange=[yy(hymin)/scale, yy(hymax)/scale]	


ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)
 xyouts,20,5, ss, /device, color=200	



indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/slice/bxbybzbrbpbz/'+indexss+'.png',image_p, red,green, blue





nn=nn+1


endwhile



end
