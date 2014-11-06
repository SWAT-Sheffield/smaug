
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1000,ysize=1000,XPOS = 950, YPOS = 300 
window, 1,xsize=1000,ysize=1000,XPOS = 950, YPOS = 300 

!p.multi = [0,1,0,0,1]

filename1='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver/slice/VrVpVzbrbpbz_1Mm/Vphi_h1Mm.dat'
close,1

filename2='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196_multidriver/slice/VrVpVzbrbpbz_05Mm/Vphi_h05Mm.dat'
close,2

nn=258

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
  a(*)=ww1(i,j,*)
  a = FFT(a, -1,/OVERWRITE)

 b=dblarr(nn)
  b(*)=ww05(i,j,*)
  b = FFT(b, -1,/OVERWRITE)
  
;  wfft_re[i,j]=Real_part(max(a))
;  wfft_Im[i,j]=Imaginary(max(a))

;  abs_ww[i,j]=max(abs(a), p)-max(abs(b), p)
;   abs_ww[i,j]=max(abs(a), p)
  
 ; print, abs_ww[i,j], p
  
  
  
;  plot, a, xrange=[0,20], psym=4
 ;  omega_ww[i,j]=ni[p]/1.4d0/(nn*1.d0)*1000.d0 
;print, max(a)


;wset,0
;  plot, ni/1.4d0/(nn*1.d0)*1000.d0, abs(a), xrange=[0.d0,50.d0]
;  oplot, ni/1.4d0/(nn*1.d0)*1000.d0, abs(a), psym=4
  
;wset,1
;  plot, ni/1.4d0/(nn*1.d0)*1000.d0, abs(b), xrange=[0.d0,50.d0]
;  oplot, ni/1.4d0/(nn*1.d0)*1000.d0, abs(b), psym=4
  


 ; plot, ni/1.4d0/nn*1000.d0, (abs(a)+1.d0)/(abs(b)+1.d0), xrange=[0.d0,40.d0]
 ; oplot, ni/1.4d0/nn*1000.d0, (abs(a)+1.d0)/(abs(b)+1.d0), psym=4



;out(*,k) = (abs(a(*))+1.d0)/(abs(b(*))+1.d0)

;maxarr[i,j]=max((abs(a(*))+1.d0)/(abs(b(*))+1.d0), p)
;maxom[i,j]=p/1.4d0/nn*1000.d0

;maxarr[i,j]=max((abs(b(*))+1.d0)/(abs(a(*))+1.d0), p)
;maxom[i,j]=p/1.4d0/nn*1000.d0

maxarr[i,j]=max(abs(b(*)), p)
maxom[i,j]=p/1.4d0/nn*1000.d0

; k=k+1 


 
  ;wait, 0.02
 
 print, i,j
 endfor

;tvframe, out(0:20,*), /bar
 
endfor

wset,0
;tvframe,maxarr, /bar, CT='dPdT'
tvframe,maxom, /bar, CT='dPdT' ;, xrange=[20,80], yrange=[20,80]

stop
xstart=10000.d0
xend=1990000.d0
ystart=10000.d0
yend=1990000.d0
scale=1000000.d0

tvframe, wfft_re, /bar,$
         xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0,$
	 xrange=[xstart/scale, xend/scale], $
	 yrange=[ystart/scale, yend/scale]
	 
tvframe, omega_ww, /bar,$
         xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0,$
	 xrange=[xstart/scale, xend/scale], $
	 yrange=[ystart/scale, yend/scale]
	 

;filename='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/slice/VrVpVzbrbpbz_1Mm/FFT_h1Mm.png'

;image_p = TVRD_24()
;write_png,filename,image_p, red,green, blue
stop

for i=0,nn-1 do begin
tvframe, ww(*,*,i)
endfor

end
