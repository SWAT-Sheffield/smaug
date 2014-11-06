
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


zstart=250
zend=350

;timeend=2400
;timeout1=dblarr(timeend,zend-zstart+1)

;openr, 11, 'data.dat'
;    readf, 11, timeout1,btot
;close,11
;print, timeout1
;stop

window, 0,xsize=1000,ysize=600

!p.multi = [0,2,0,0,1]



loadct,4


tvframe,btot,/bar, /sample, title=label_bt, charsize=1.2,$  
	xtitle='Horizontal distance (Mm)', xrange=[0,16], ytitle='Height (Mm)', yrange=[1.96,2.75]
	
tek_color

line2_b, n22,n11, Bxx,Bzz,xx,zz,0.d0,x(zstart,1,0)

	
oplot, [5,5],[1.9,3.0], linestyle=2, thick=2, color=230	

tvframe, timeout, /sample, xrange=[0,907], yrange=[1.96,2.75], ytitle='Height (Mm)',xtitle='Time (s)',$
         charsize=1.2, title='Time-distance'


image_p = TVRD_24()

  write_png,'/data/ap1vf/png/test/out.png',image_p, red,green, blue



end
