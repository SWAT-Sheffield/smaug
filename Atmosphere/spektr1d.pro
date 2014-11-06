;****************** read file begin *****************

nn=249

harr=dblarr(nn)
rhoarr=harr

openr, 11, 'VALMc_rho_249.dat'
 for i=nn-1,0,-1 do begin
		readf, 11, h,rho
		harr(i)=h
		rhoarr(i)=rho
 endfor
close,11

;****************** read file end *****************



;****************** filter begin *******************
fleft=30
fright=110

nj=nn/2

dsl=dindgen(nj)-fleft
dsr=dindgen(nj)-fright

ww0=dindgen(2*nj)


dsl=atan(dsl)
dsr=-atan(dsr)

maxdsl=max(dsl)
mindsl=min(dsl)

maxdsr=max(dsr)
mindsr=min(dsr)

   dsl=(dsl-mindsl)/(maxdsl-mindsl)
   dsr=(dsr-mindsr)/(maxdsr-mindsr)

for i=0,2*nj-1 do begin
 if i lt nj then ww0(i)=dsl(i) else ww0(i)=dsr(i-nj) 
endfor  

;****************** filter end *******************
stop

wF=FFT(rhoarr,-1)
stop

wF1=ABS(wF)^2.d0
stop
wF4=ww0*wF1


wF4=double(FFT(wF4,1)) ;/(wSin+0.00001))



wF4=ABS(wF4*ww0)^2.d0

plot, alog10(wF4)





end
