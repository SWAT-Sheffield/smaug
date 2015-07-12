
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

;set_plot, 'ps'


;device, filename='sol.ps', /color, BITS=8


!P.multi=0
close,11
nvaldata=21
harr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)


;openr, 11, '../valrhoVALMcnew_n256.dat'
;openr, 11, '../valrhoPhot.dat'
;openr, 11, '../valrhoPot12MM.dat'
openr, 11, '../valrhoPhot_short.dat'

 for i=0,nvaldata-1 do begin
		readf, 11, valh, valrho
		
            print,valh, valrho,i		
		harr(i)=valh
		rhoarr(i)=valrho
 endfor
close,11

nvaldata=256
harr1=dblarr(nvaldata)
rhoarr1=dblarr(nvaldata)


;openr, 11, '../valrho_n256.dat'
; for i=0,nvaldata-1 do begin
;		readf, 11, valh, valrho
		
;            print,valh, valrho,i		
;		harr1(i)=valh
;		rhoarr1(i)=valrho
; endfor
;close,11;


;plot, harr,alog10(rhoarr), charsize=2.0
;oplot, harr1,alog10(rhoarr1), psym=4
;device, /close
;set_plot, 'x'
;stop


;rhoarr_N=dblarr(nvaldata/2)
;harr_N=dblarr(nvaldata/2)

;for i=0,127 do begin
;rhoarr_N[i]=rhoarr[2*i]
;harr_N[i]=harr[2*i]
;endfor

;plot, alog10(rhoarr_N),xrange=[0,255], psym=4, /ys
;oplot, alog10(rhoarr_N)

;stop

;print,'----->',n_elements(rhoarr)



nvaldata=256

rhoarr = CONGRID(rhoarr, nvaldata,/interp) 
harr = CONGRID(harr, nvaldata, /interp) 

;plot, harr,rhoarr   /ys

plot, harr, alog10(rhoarr), /ys
oplot,harr,alog10(rhoarr),psym=4




;print,'------------------------>',n_elements(rhoarr)



openw, 11, '../valrhoPhotshort_n256.dat'
 for i=0,nvaldata-1 do begin
		printf, 11, harr(i), rhoarr(i)
 endfor
close,11


g=-273.0d0

p=dblarr(nvaldata-5)

for i=nvaldata-2,5, -1 do begin
 
 comi=abs(harr(i+1)-harr(i))
 p[i]=p[i+1]+rhoarr[i]*comi*g
 print, harr[i],p[i]
endfor


end
