
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1200,ysize=325,XPOS = 50, YPOS = 700 

!p.multi = [0,3,0,0,1]

close,11
nvaldata=252
harr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)
muarr=dblarr(nvaldata)

;x=9.09d6-findgen(nvaldata)*9.09d6/(nvaldata-1)

x=findgen(nvaldata)*9.09d6/(nvaldata-1)

;stop
;openr, 11, '../valrhoVALMc_n256.dat'
;openr, 11, '../valrho_n256.dat'
openr, 11, '../val.dat'
 for i=nvaldata-1,0,-1 do begin
		readf, 11, valh, valrho
		
            print,valh, valrho,i		
		harr(i)=valh*1000.d0
		rhoarr(i)=valrho
 endfor
close,11

openr, 11, '../mu_252.dat'
 for i=0,nvaldata-1 do begin
		readf, 11, valh, valmu
		
            print,valh, valmu,i		
		muarr(i)=valmu
 endfor
close,11


for i=0,nvaldata-1 do begin
print, harr[i],muarr[i],x[i],i
endfor
;stop
;rhoarr = INTERPOL(rhoarr, harr, x)



;openw, 11, '../valrhoVALMcnew_n256.dat'
; for i=0,nvaldata-1 do begin
;		printf, 11, x(i)/1000, rhoarr(i)
; endfor
;close,11

;stop
plot, x,alog10(rhoarr), charsize=2.0, psym=4
g=-274.0d0

parr=dblarr(nvaldata)

;parr[nvaldata-1]=789000*8.31e3*(9.0412855e-13)/0.6d0

parr[nvaldata-1]=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0

;parr[nvaldata-1]=1000000.0d0*8.31e3*(1.1790001e-11)/0.6d0

print,'************',parr[nvaldata-1]

for i=nvaldata-2,0, -1 do begin
 
 comi=-abs(x[i+1]-x[i])
 parr[i]=parr[i+1]+(rhoarr[i]+rhoarr[i+1])*comi*g/2.d0
 print, x[i],parr[i], rhoarr[i],i
endfor


for i=2,nvaldata-3 do begin
 
 rhoarr[i]=-1.d0/g*(1.d0/(12.D0*(x[i+1]-x[i])))*(parr[i+2]-8.d0*parr[i+1]+8.d0*parr[i-1]-parr[i-2])

endfor


;oplot, x,alog10(rhoarr), psym=4
stop

plot, x,alog10(parr), charsize=2.0

mu=0.6
R=8.31e3


Tarr=dblarr(nvaldata)

Tarr=parr*mu/rhoarr/R


plot, x,alog10(Tarr), charsize=2.0
;oplot, x,alog10(Tarr), psym=4

print,'*********** Temp *****************'
for i=0,nvaldata-1 do begin
print, Tarr[i],parr[i],mu,x[i]/1000.d0
endfor



;openw, 11, '../rho_p_252.dat'
; for i=0,nvaldata-1 do begin
;		printf, 11, x(i), rhoarr(i), parr(i)
; endfor
;close,11

print,'DONE'

end
