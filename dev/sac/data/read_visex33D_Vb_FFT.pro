
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1060,ysize=760,XPOS = 950, YPOS = 300 

!p.multi = [0,2,3,0,1]

;filename1='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver/slice/VrVpVzbrbpbz_1Mm/Vphi_h1Mm.dat'
;close,1

;filename2='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver/slice/VrVpVzbrbpbz_05Mm/Vphi_h05Mm.dat'
;close,2

filename1='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver_lower/normalised/slice/VrVpVzbrbpbz_1Mm/Vphi_h1Mm.dat'
close,1

filename2='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver_lower/normalised/slice/VrVpVzbrbpbz_05Mm/Vphi_h05Mm.dat'
close,2

;nn=258 ; for 30, 60,90, 120
nn=1060 ; for 120,180, 240,300,350

ni=findgen(nn)


ww1=dblarr(100,100,nn)
ww05=dblarr(100,100,nn)

openr,1,filename1,/f77_unf
openr,2,filename2,/f77_unf

readu,1,ww1
readu,2,ww05
close, 1
close, 2


;for i=0, nn do begin
;wset,0
;tvframe, ww1(*,*,i), /bar,$
;         xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0,$
;	 xrange=[xstart/scale, xend/scale], $
;	 yrange=[ystart/scale, yend/scale]
;
;wset,1

;tvframe, ww05(*,*,i), /bar,$
;         xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0,$
;	 xrange=[xstart/scale, xend/scale], $
;	 yrange=[ystart/scale, yend/scale]
;endfor


;stop


out=dblarr(nn,10000)

maxarr_1=dblarr(100,100)
maxom_1=dblarr(100,100)

maxarr_05=dblarr(100,100)
maxom_05=dblarr(100,100)


wfft_re=dblarr(100,100)
wfft_Im=dblarr(100,100)

abs_ww=dblarr(100,100)
omega_ww=dblarr(100,100)

k=0


for i=0,99 do begin
 for j=0,99 do begin
 a=dblarr(nn)
  a(*)=ww1(i,j,*)
  
  a = FFT(a, -1,/OVERWRITE)

 b=dblarr(nn)
  b(*)=ww05(i,j,*)
  b = FFT(b, -1,/OVERWRITE)

dtime=1.4d0 ; for 30, 60,90,120  
dtime=0.7d0 ; for lower !!!!!


maxarr_05[i,j]=max(abs(b(*)), p)
maxom_05[i,j]=p/dtime/nn*1000.d0

maxarr_1[i,j]=max(abs(a(*)),p)


maxom_1[i,j]=p/dtime/nn*1000.d0

  ;wait, 0.02
 
 print, i,j
 endfor

;tvframe, out(0:20,*), /bar
 
endfor

omega=8.3
eps=0.1

maxpower_omega=dblarr(100,100)

for i=0, 99 do begin
for j=0, 99 do begin

omratio=maxom_05[i,j]/omega-1.d0
if (abs(omratio) le eps) then maxpower_omega[i,j]=maxarr_05[i,j] else maxpower_omega[i,j]=0.d0   

endfor
endfor



wset,0

tvframe, maxom_05, /bar, CT='dPdT', title='0.5 Mm' 

tvframe, maxom_1, /bar, CT='dPdT', title='1 Mm', charsize=1.5

sm_maxom_05=smooth(maxom_05,3)
tvframe, sm_maxom_05, /bar, CT='dPdT', title='smooth 0.5 Mm'


xstart=10000.d0
xend=1990000.d0
ystart=10000.d0
yend=1990000.d0
scale=1000000.d0

	 
tvframe, maxpower_omega, /bar, title='0.5 Mm', $
         xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0,$
	 xrange=[xstart/scale, xend/scale], $
	 yrange=[ystart/scale, yend/scale]
	 
	 
stop	 
	 
	 
;*************************************** begin B background ************************
tarr=dblarr(1)
mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)



headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)


close,1

openr,1,'/data/ap1vf/3D_tube_196_100_100_multidriver_lower.out',/f77_unf


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

mu=4.0*!PI/1.0e7

print,'******************* time = ', time



scale=1.d6

mu=1.257E-6

xstart=0
xend=99
ystart=0
yend=99

pp=50 ;x
kk=5  ;y


zstart=0
zend=185

wt=dblarr(zend-zstart+1,xend-xstart+1,iw)
wt=reform(w(zstart:zend,xstart:xend,pp,*))


;tvframe,rotate(wt(*,*,10),1)*sqrt(mu)*1.0e4,/bar,/sample, title='Bz_b', $
;        xtitle='x', ytitle='z', charsize=2.0, $
;	xrange=[xx[xstart]/scale, xx[xend]/scale], $
;	yrange=[zz[zstart]/scale, zz[zend]/scale]	


;endwhile

;xstart=40
;xend=58
;ystart=40 
;yend=58

xstart=32
xend=66  
ystart=32 
yend=66  

aa=dblarr(xend-xstart+1)

maxom_05(50,50)=maxom_05(50,51)


aa(*)=maxom_05(50,xstart:xend)

ab05=smooth(aa,4)
plot,ab05/max(ab05), charsize=1.5, yrange=[0,1.1], title='0.5 Mm' 
oplot, ab05/max(ab05),psym=4
oplot, rotate(wt(60,xstart:xend,10),1)/max(wt(60,*,10)), color=200


aa=dblarr(xend-xstart+1)

maxom_1(50,50)=maxom_1(50,51)

aa(*)=maxom_1(50,xstart:xend)

ab1=smooth(aa,4)
;ab1=smooth(aa,4)
;ab1=smooth(aa,4)
;ab1=smooth(aa,4)

plot,ab1/max(ab1), charsize=1.5, yrange=[0,1.1], title='1.0 Mm'
oplot, ab1/max(ab1),psym=4
oplot, rotate(wt(121,xstart:xend,10),1)/max(wt(121,*,10)), color=200


;*************************************** end B background ************************


;******************************* save PS files ************************* 



window, 0,xsize=300,ysize=300,XPOS = 950, YPOS = 300 
!p.multi = [0,1,1,0,1]

set_plot, 'ps'
device, filename='Omega_Mag_field_1Mm_350_300_240_180_120_2D.ps', ysize=12, xsize=14, $
        xoffset=1, yoffset=4, /color
!p.thick = 5
!x.thick = 1
!y.thick = 1
!z.thick = 1
!p.font = 1.0


scale=1000000.d0

tvframe, maxom_1(xstart:xend,ystart:yend), /bar, CT='dPdT', title='h=1.0 Mm', $
         xtitle='x [Mm]', ytitle='y [Mm]', $
         xrange=[xx[xstart]/scale, xx[xend]/scale], $
	 yrange=[zz[zstart]/scale, zz[zend]/scale], BTITLE='omega [mHz]', charsize=1.2	
 
device, /close
set_plot, 'x'


window, 0,xsize=300,ysize=300,XPOS = 950, YPOS = 300 
!p.multi = [0,1,1,0,1]

set_plot, 'ps'
device, filename='Omega_Mag_field_1Mm_350_300_240_180_120.ps', ysize=12, xsize=14,$
        xoffset=1, yoffset=4, /color
!p.thick = 5
!x.thick = 1
!y.thick = 1
!z.thick = 1
!p.font = 1.0


scale=1000000.d0


;plot,xx(xstart:xend)/scale,ab1/max(ab1), charsize=1.5, yrange=[0,1.1], $
;     title='1.0 [Mm]', xtitle='x [Mm]',  ytitle='normalised frequency, Bz' 
;oplot, xx(xstart:xend)/scale, ab1/max(ab1),psym=4
;oplot, xx(xstart:xend)/scale, rotate(wt(121,xstart:xend,10),1)/max(wt(121,*,10)), color=200

plot,xx(xstart:xend)/scale,ab05/max(ab05), charsize=1.5, yrange=[0,1.1], $
     title='0.5 [Mm]', xtitle='x [Mm]',  ytitle='normalised frequency, Bz' 
oplot, xx(xstart:xend)/scale, ab05/max(ab05),psym=4
oplot, xx(xstart:xend)/scale, rotate(wt(60,xstart:xend,10),1)/max(wt(60,*,10)), color=200


device, /close
set_plot, 'x'


end
