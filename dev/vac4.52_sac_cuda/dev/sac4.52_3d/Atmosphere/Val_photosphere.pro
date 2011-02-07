
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1200,ysize=625,XPOS = 50, YPOS = 700 

!p.multi = [0,3,2,0,1]

close,11
add=0

nvaldata=49+add
harr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)
parr=dblarr(nvaldata)
Tarr=dblarr(nvaldata)
PgParr=dblarr(nvaldata)
mu=dblarr(nvaldata)


openr, 11, 'VALIIIC.dat'
 for i=nvaldata-1-add,0,-1 do begin
		readf, 11, valh, valrho, valp, valT, valPgP
		
            print,valh, valrho, valp, valT, i		
		harr(i)=valh*1000.d0
		rhoarr(i)=valrho
		parr(i)=valp*valPgP
		Tarr(i)=valT
 endfor
close,11




plot, harr, alog10(rhoarr), charsize=1.5, title='log rho'
plot, harr, rhoarr, charsize=1.5, title='rho'

plot, harr, alog10(Tarr), charsize=1.5, title='log T'
plot, harr, Tarr, charsize=1.5, title='T'

plot, harr, alog10(parr), charsize=1.5, title='log p'
plot, harr, parr, charsize=1.5, title='p'




stop


R=8.31e7
;mu=Tarr*R*rhoarr/parr

;ii=1.d0
;for i=nvaldata-add, nvaldata-1 do begin
;  harr[i]=(980.0+ii*100.d0) *1000.d0
;  print, harr[i],i
;  mu[i]=0.6d0
;  ii=ii+1
;endfor

;plot, harr,mu, charsize=2.0
;oplot, harr,mu, psym=4

print,'********** mu ****************'
for i=0,nvaldata-1 do begin
print, harr[i],mu[i],i
endfor

stop


nvaldata=252

x=findgen(nvaldata)*2.543d6/(nvaldata-1)

mu = INTERPOL(mu, harr, x, /spline)

for i=0,nvaldata-1 do begin
print, x[i],mu[i],i
endfor

plot, x,mu, charsize=2.0
oplot, x,mu, psym=4

for i=0,nvaldata-1 do begin
print, x[i]/1000.0,mu[i]
endfor



openw, 11, '../mu_252.dat'
 for i=0,nvaldata-1 do begin
		printf, 11, x(i), mu(i)
 endfor
close,11

print, 'DONE'

end
