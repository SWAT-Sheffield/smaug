function omega, tt

;delta=260.d0
;f=2.d0*!Pi/120.0+2.d0*!Pi/exp(-t*t/delta/delta)/460.0

f120=2.d0*!Pi/120.d0
f30=2.d0*!Pi/30.d0

t120=120.d0
t90=90.d0
t60=60.d0
t30=30.d0

om120=1.d0/t120*1000.d0
om90=1.d0/t90*1000.d0
om60=1.d0/t60*1000.d0
om30=1.d0/t30*1000.d0

print, '************************************'
print, om120,' ',om90,' ',om60,' ',om30
print, '************************************'
;ff=(f120-f30)/(t120-t30)*(tt-t30)+f30

ff=sin(2.d0*!Pi/t30*tt)+sin(2.d0*!Pi/t60*tt)+sin(2.d0*!Pi/t90*tt)+sin(2.d0*!Pi/t120*tt)

;print, ff,tt


res=ff

return, res
end

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=800,ysize=800,XPOS = 1100, YPOS = 400 
nn=10000

!p.multi = [0,1,2,0,1]

delta=0.1d0
t=findgen(nn)*delta
Om=dblarr(nn)
omega_ww=dblarr(nn)


Om_120=findgen(nn)
Om_30=findgen(nn)

for i=0, nn-1 do Om_120[i]=2.d0*!Pi/120.d0
for i=0, nn-1 do Om_30[i]=2.d0*!Pi/30.d0

Om=omega(t)

print, max(Om)
stop

plot, t, Om, charsize=3.0 ;, psym=4
;oplot, t, Om_120, psym=4
;oplot, t, Om_30, psym=3


  Om = FFT(Om, -1,/OVERWRITE)


for i=0, nn-1 do omega_ww[i]=i*1.d0/delta/nn*1000.d0 

plot, omega_ww, abs(Om), xrange=[00,100], charsize=3.0




end
