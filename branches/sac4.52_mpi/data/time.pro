
set_plot, 'ps'
zstart=250
zend=350

;timeend=2400
;timeout1=dblarr(timeend,zend-zstart+1)

;openr, 11, 'data.dat'
;    readf, 11, timeout1,btot
;close,11
;print, timeout1
;stop


loadct,4

device, filename='sol_01.ps', /color, BITS=8
DEVICE, XSIZE=10, YSIZE=3, /INCHES
!p.thick = 4
!x.thick = 4
!y.thick = 4
!z.thick = 4
!p.font = 1.0



tvframe,btot,$
        /bar, title=label_bt, $ ;charsize=1.0,$  
	xtitle='Horizontal distance (Mm)', xrange=[0,16], ytitle='Height (Mm)', yrange=[1.96,2.75]
	
tek_color

line2_b, n22,n11, Bxx,Bzz,xx,zz,0.d0,x(zstart,1,0)
	
oplot, [5,5],[1.9,3.0], linestyle=2, thick=2, color=230	

tvframe, timeout, xrange=[0,907], yrange=[1.96,2.75], ytitle='Height (Mm)',xtitle='Time (s)',$
         charsize=1.0, title='Time-distance'


device, /close
set_plot, 'x'


end
