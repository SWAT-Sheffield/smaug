!P.multi=0
close,11
nvaldata=505
harr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)

openr, 11, '../Atmosphere/VALMc_rho_505.dat'
 for i=0,nvaldata-1 do begin
		readf, 11, valh, valrho
		
            print,valh, valrho,i		
		harr(i)=valh
		rhoarr(i)=valrho
 endfor
close,11

;rhoarr_N=dblarr(nvaldata/2)
;harr_N=dblarr(nvaldata/2)

;for i=0,127 do begin
;rhoarr_N[i]=rhoarr[2*i]
;harr_N[i]=harr[2*i]
;endfor

;plot, alog10(rhoarr_N),xrange=[0,255], psym=4, /ys
;oplot, alog10(rhoarr_N)

;stop

print,'----->',n_elements(rhoarr)



nvaldata=1024

rhoarr = CONGRID(rhoarr, nvaldata,/interp) 
harr = CONGRID(harr, nvaldata, /interp) 

;plot, rhoarr_N ,xrange=[300,600], /ys

plot, alog10(rhoarr) ,xrange=[0,1023], /ys
oplot,alog10(rhoarr),psym=4




print,'------------------------>',n_elements(rhoarr)



openw, 11, '../Atmosphere/VALMc_rho_1023.dat'
 for i=0,nvaldata-2 do begin
		printf, 11, harr(i), rhoarr(i)
 endfor
close,11


end
