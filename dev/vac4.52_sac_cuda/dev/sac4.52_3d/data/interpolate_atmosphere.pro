
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


!P.multi=0
close,11
nvaldata=13
harr=dblarr(nvaldata)
Tarr=dblarr(nvaldata)
parr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)

openr, 11, '../McWriter.dat'
 for i=0,nvaldata-1 do begin
		readf, 11, h, T, p
		
            print,h, T, p,i		
		harr(i)=h*1000.d0	 ; [m]
		Tarr(i)=T                ; [K]
		parr(i)=p*0.1d0          ; [H/m^2]
 endfor
close,11



mu=0.6d0
R=8.31e3

rhoarr=parr*mu/Tarr/R

for i=nvaldata-1,0,-1 do begin
  print,harr(i)/1.e3+300.d0,rhoarr(i)
endfor

plot, harr/1.0e6,alog10(rhoarr), charsize=2.0

stop

rhoarr_N=dblarr(nvaldata/2)
harr_N=dblarr(nvaldata/2)

for i=0,127 do begin
rhoarr_N[i]=rhoarr[2*i]
harr_N[i]=harr[2*i]
endfor

;plot, alog10(rhoarr_N),xrange=[0,255], psym=4, /ys
;oplot, alog10(rhoarr_N)

;stop

print,'----->',n_elements(rhoarr)



nvaldata=100

rhoarr_N = CONGRID(rhoarr_N, nvaldata,/interp) 
harr_N = CONGRID(harr_N, nvaldata, /interp) 

;plot, rhoarr_N ,xrange=[300,600], /ys

plot, alog10(rhoarr_N) ,xrange=[0,119], /ys
oplot,alog10(rhoarr_N),psym=4




print,'------------------------>',n_elements(rhoarr)



openw, 11, '../valrho_n100.dat'
 for i=0,nvaldata-1 do begin
		printf, 11, harr_N(i), rhoarr_N(i)
 endfor
close,11


end
