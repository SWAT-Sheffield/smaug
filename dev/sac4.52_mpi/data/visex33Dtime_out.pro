
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
;PRINT, 'Date:      ', systime(0)
;PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
;PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


;window, 0,xsize=800,ysize=425,XPOS = 950, YPOS = 300 
;!p.multi = [0,2,1,0,1]

ns1=intarr(1)
tst1=double(1)
ted1=double(1)

nt1=intarr(1)
nz1=intarr(1)

close,3
;openr,3,'/data/ap1vf/3D_509_36_36_300s_e_1616_t.out',/f77_unf
openr,3,'/data/ap1vf/3D_396_60_60_t.out',/f77_unf
readu,3,tst1,ted1,ns1
readu,3,nz1
print, ns1

timearr1=dblarr(ns1[0],nz1[0])
beta1=dblarr(ns1[0],nz1[0])
zz1=dblarr(nz1[0])

readu,3,timearr1
;readu,3,beta1
readu,3,zz1

close,3

xx1=dblarr(ns1[0])

for i=0,ns1[0]-1 do xx1[i]=tst1+(ted1-tst1)/(ns1[0]*1.d0)*(i*1.d0)

scale=1.0d6

; ****************** ps, eps begin ****************************************
;   xs=19.d0
;   k=15.d0/19.d0
   
;   SET_PLOT,'ps'  

;   device, filename='/data/ap1vf/ps/time.eps', $
;   BITS=8, /color, xsize=xs, ysize=k*xs, /encap
   
;!p.thick = 2
;!x.thick = 2
;!y.thick = 2
;!z.thick = 2


;timearr(where(abs(timearr) le 0.5d0))=0.d0

;title='(!7q!N-!7q!I0!N)/!7q!I0!N!3'
tit_pd='V!I!Mx!N!3 [m/s]'
tit_pl='!3V!I| |!N!3 [km/s]'
tit_T='T-T!I0!N [K]'
tit_e='e'
tit_emag='e_mag'

zs=80

tvframe, timearr1[*,0:196], /sample, /bar, title=tit_e, $
        xtitle='Time [s]', ytitle='z [Mm]',charsize=1.0, CT='mixct16', $
        xrange=[xx1[0], xx1[ns1[0]-1]], yrange=[zz1[zs]/scale, max(zz1)/scale]


mixct10
;loadct,3
;contour, timearr,xx,zz/scale, LEVELS = [0.01], /overplot, /FOLLOW, thick=2.0, color=200

;contour, beta,xx,zz/scale, LEVELS = [0.01,0.1,1.0, 10.0], $
;	 C_Thick = [2.0,2.0,4.0,2.0],  $
;	 /overplot, /FOLLOW, thick=4.0, color=128, charsize=4, $
;         C_LABELS = [100, 10, 0, 1], $
;         C_Annotation = ['100.0','10.0','1.0','0.1']


;device, /close
;set_plot, 'x'

end


