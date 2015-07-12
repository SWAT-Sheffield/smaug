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

PRO gfunct, X, A, F, pder  
  bx = EXP(A[1] * X)  
  F = A[0] * bx + A[2]  
  
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF N_PARAMS() GE 4 THEN $  
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]  
END  

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1400,ysize=1025,XPOS = 700, YPOS = 700 

!p.multi = [0,4,3,0,1]


;****************** VALIIIF atmosphere begin ***************
close,11
add=0

nvaldata=49-add
harr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)
parr=dblarr(nvaldata)
Tarr=dblarr(nvaldata)
PgParr=dblarr(nvaldata)



openr, 11, 'VAL_IIIF.dat'
 for i=nvaldata-1-add,0,-1 do begin
		readf, 11, valh, valrho, valp, valT, valPgP
		
            print,valh, valrho, valp, valT, i		
		harr(i)=valh*1000.d0
		rhoarr(i)=valrho
		parr(i)=valp*valPgP
		Tarr(i)=valT
 endfor
close,11
;****************** VALIIIF atmosphere end ***************

!P.multi=0


nvaldata=200
hmax=1.6e3*1000.d0


x=findgen(nvaldata)*hmax/(nvaldata-1)

rhoarr = INTERPOL(rhoarr, harr, x, /spline)
Tarr = INTERPOL(Tarr, harr, x, /spline)



rhoarr1=dblarr(2*nvaldata)

x1=findgen(2*nvaldata)*2.d0*hmax/(2*nvaldata-1)

for i=0, 2*nvaldata-1 do begin
 if (i le nvaldata-1) then begin
  rhoarr1[i]=rhoarr[i]
 endif else begin
  rhoarr1[i]=rhoarr[nvaldata-1]
 endelse
endfor 

rhoarr1 = INTERPOL(rhoarr1, x1, x1, /spline)

stop


plot, x/1.0e6, alog10(rhoarr), charsize=1.5, title='log rho'
;plot, x,deriv(x,rhoarr)

rho=rhoarr1
;T=Tarr

stop
rho = SMOOTH( rhoarr1, 2) ;, /EDGE_TRUNCATE) 
for i=0,100 do rho = SMOOTH( rho, 2)

!P.multi=0

plot, x1/1.0d6, alog10(rho), psym=4
oplot, x1/1.0d6, alog10(rhoarr1)

stop

plot, x1,deriv1(rho,x1), psym=4
oplot, x1,deriv1(rhoarr1,x1)
stop


plot, x1,deriv1(deriv1(rhoarr1,x1),x1)
oplot, x1,deriv1(deriv1(rho,x1),x1), psym=4

stop

openw, 11, 'VALIIIF_rho_200_1_6Mm_test.dat'
 for i=2*nvaldata-2,3,-1 do begin
		printf, 11, x1[i], rho[i]*1.0e3		;
 endfor
close,11

 


print, '*DONE'

end
