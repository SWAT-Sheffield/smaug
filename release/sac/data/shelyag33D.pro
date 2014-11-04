
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


; Open an MPEG sequence: 
;mpegID = MPEG_OPEN([700,1200],FILENAME='myMovie.mpg') 

window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 

close,1 

openr,1,'/data/ap1vf/vxx.040000'
wx=fltarr(480,100,480)
readu, 1, wx
close,1 
openr,1,'/data/ap1vf/vxx.040000'
wy=fltarr(480,100,480)
readu, 1, wy
close,1 
openr,1,'/data/ap1vf/vxx.040000'
wz=fltarr(480,100,480)
readu, 1, wz
close,1 

wset,0
!p.multi = [0,1,3,0,1]


tvframe,wx(*,50,*)
tvframe,wy(*,50,*)	
tvframe,wz(*,50,*)


end
