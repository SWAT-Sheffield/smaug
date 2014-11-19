cdtarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)




DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

;window, 0,xsize=1200,ysize=700,XPOS = 1000, YPOS = 300 ; ZOOM
;window, 0,ysize=1600,xsize=800,XPOS = 100, YPOS = 300
;window, 0,ysize=500,xsize=800,XPOS = 100, YPOS = 300
window, 0,ysize=1000,xsize=1600,XPOS = 100, YPOS = 300

ii=1



if (ii eq 1) then begin
loadct,4
endif else begin
loadct,0
tek_color
endelse

ipic=long(1)

n1=2048
n2=1024
nsect1=512
nsect2=512

pic=521
;pic=pic/2

;!p.multi = [0,1,0,0,1]
;!p.multi = [0,3,2,0,1]
;!p.multi = [0,2,2,0,1]
!p.multi = [0,4,2,0,1]



wdtv1=dblarr(pic+1,2)
wdt=dblarr(pic+1,n1,4)
wdt1=dblarr(pic+1,n1,4)

wdt1v=dblarr(pic+1,n1,4)
wdt2v=dblarr(pic+1,n1,4)



wfft1=dblarr(pic+1,n1,4)
wfft2=dblarr(pic+1,n1,4)


restore,'fft/wdt1h_7b0_521s.dat'
restore,'fft/wdt2h_7b0_521s.dat'

dt=0.001

myrange=[0,6]
mxrange=[0,15]

mv1range=[-1,1]
mv2range=[-1,1]




for i=1,n1-1 do begin

   wfft1(*,i,0)=fft(wdt1v(*,i,0))
   wfft1(*,i,1)=fft(wdt1v(*,i,1))
   wfft1(*,i,2)=fft(wdt1v(*,i,2))
   wfft1(*,i,3)=fft(wdt1v(*,i,3))

   wfft2(*,i,0)=fft(wdt2v(*,i,0))
   wfft2(*,i,1)=fft(wdt2v(*,i,1))
   wfft2(*,i,2)=fft(wdt2v(*,i,2))
   wfft2(*,i,3)=fft(wdt2v(*,i,3))


endfor

tvframe,wfft1(0:50,*,1),/bar ,/sample , yrange=myrange, title='fft centre v1', charsize=2.0
tvframe,wfft1(0:50,*,2),/bar ,/sample, yrange=myrange, title='fft centre v2', charsize=2.0

tvframe,wfft2(0:50,*,1),/bar ,/sample , yrange=myrange, title='fft off centre v1', charsize=2.0
tvframe,wfft2(0:50,*,2),/bar ,/sample, yrange=myrange, title='fft off centre v2', charsize=2.0


tvframe,wdt1v(*,*,1),/bar ,/sample , yrange=myrange, title='centre v1', charsize=2.0
tvframe,wdt1v(*,*,2),/bar ,/sample ,  yrange=myrange, title='centre v2', charsize=2.0

tvframe,wdt2v(*,*,1),/bar ,/sample ,  yrange=myrange, title='off centre v1', charsize=2.0
tvframe,wdt2v(*,*,2),/bar ,/sample ,  yrange=myrange, title='off centre v1', charsize=2.0



end




