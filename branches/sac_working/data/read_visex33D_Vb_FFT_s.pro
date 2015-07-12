
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

;filename1='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver_lower/normalised/slice/VrVpVzbrbpbz_05Mm/Vphi_h05Mm.dat'
;close,1


;********************** begin main information *****************************
title='1.0'
;zb=60 ; for 0.5 Mm
zb=121 ; for 1.0 Mm

; for 30, 60,90, 120
;nn=258
;dtime=1.4d0 
;filename_out='Omega_Mag_field_'+title+'Mm_120_90_60_30'


; for 120,180,240,300,350
nn=1060 
dtime=0.7d0
filename_out='Omega_Mag_field_'+title+'Mm_350_300_240_180_120'

;********************* end main information ********************************

ni=findgen(nn)


ww=dblarr(100,100,nn)
openr,1,filename1,/f77_unf
readu,1,ww

close, 1

out=dblarr(nn,10000)

maxarr=dblarr(100,100)
maxom=dblarr(100,100)

wfft_re=dblarr(100,100)
wfft_Im=dblarr(100,100)

abs_ww=dblarr(100,100)
omega_ww=dblarr(100,100)

k=0


for i=0,99 do begin
 for j=0,99 do begin
 a=dblarr(nn)
  a(*)=ww(i,j,*)
  
  a = FFT(a, -1,/OVERWRITE)

maxarr[i,j]=max(abs(a(*)),p)
maxom[i,j]=p/dtime/nn*1000.d0
 
 print, i,j
 endfor
 
endfor


wset,0

maxom_s=smooth(maxom,4)

tvframe, maxom_s, /bar, CT='dPdT', title='h='+title+ '[Mm]', charsize=1.5
	 
	 
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
;*************************************** end B background ************************
xstart=40
xend=58
ystart=40 
yend=58

;xstart=40
;xend=68
;ystart=40 
;yend=68

xstart=40 ; 32
xend=62   ;66
ystart=40 ;32
yend=62   ;66


xstart=15 ; 32
xend=84   ;66
ystart=15 ;32
yend=84   ;66



aa=dblarr(xend-xstart+1)
ab=dblarr(xend-xstart+1)

maxom(50,50)=maxom(50,51)


;aa(*)=maxom(xstart:xend,50)
aa(*)=maxom(50,ystart:yend)
;ab=aa
ab=smooth(aa,4)
ab=smooth(ab,4)
ab=smooth(ab,4)
ab=smooth(ab,4)

plot,ab/max(ab), charsize=1.5, yrange=[0,1.1], title='h='+title+ '[Mm]'
oplot, ab/max(ab),psym=4
oplot, rotate(wt(zb,xstart:xend,10),1)/max(wt(zb,*,10)), color=200



;******************************* save PS files ************************* 

window, 0,xsize=300,ysize=300,XPOS = 950, YPOS = 300 
!p.multi = [0,1,1,0,1]

set_plot, 'ps'
device, filename=filename_out+'2D_large_domain.ps', ysize=12, xsize=14, $
        xoffset=1, yoffset=4, /color
!p.thick = 5
!x.thick = 1
!y.thick = 1
!z.thick = 1
!p.font = 1.0


scale=1000000.d0
;
tvframe, maxom_s(xstart:xend,ystart:yend), /bar, CT='dPdT', title='h='+title+ '[Mm]', $
         xtitle='x [Mm]', ytitle='y [Mm]', $
         xrange=[xx[xstart]/scale, xx[xend]/scale], $
	 yrange=[zz[zstart]/scale, zz[zend]/scale], BTITLE='omega [mHz]', charsize=1.2	
 
device, /close
set_plot, 'x'


window, 0,xsize=300,ysize=300,XPOS = 950, YPOS = 300 
!p.multi = [0,1,1,0,1]

set_plot, 'ps'
device, filename=filename_out+'no_smoothing.ps', ysize=12, xsize=14,$
        xoffset=1, yoffset=4, /color


!p.thick = 5
!x.thick = 1
!y.thick = 1
!z.thick = 1
!p.font = 1.0

plot,xx(xstart:xend)/scale,ab/max(ab), charsize=1.5, yrange=[0,1.1], $
     title='h='+title+ '[Mm]', xtitle='y [Mm]',  ytitle='normalised frequency, Bz' 
oplot, xx(xstart:xend)/scale, ab/max(ab),psym=4
oplot, xx(xstart:xend)/scale, rotate(wt(zb,xstart:xend,10),1)/max(wt(zb,*,10)), color=200


device, /close
set_plot, 'x'

print, '*** DONE ***'
end
