DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
;--------------------
; Make a vector of 16 points, A[i] = 2pi/16: 
A = FINDGEN(17) *  (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
USERSYM, 0.6*COS(A), 0.6*SIN(A), /fill
;--------------------

Window, 1, xsize=1400,ysize=300, $
           XPOS = 850, YPOS = 500

dtime=1.04d0
!P.multi=[0,3,0,0,1]

RESTORE, 'T_Vphi_s.saw' ; T_Vphi_s, Vr, phi_s, Vphi_test


nn=n_elements(Vphi_test)
ni=findgen(nn)

plot, Vphi_test, charsize=1.2
oplot, Vphi_test, psym=8
a = FFT(Vphi_test, -1)

plot, ni(1:50)/dtime/nn*1000.d0, a(1:50), charsize=1.2
oplot,ni(1:50)/dtime/nn*1000.d0, a(1:50), psym=8

plot, ni/dtime/nn*1000.d0, charsize=2
oplot, ni/dtime/nn*1000.d0, psym=4
end
