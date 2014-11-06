
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
;PRINT, 'Date:      ', systime(0)
;PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
;PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


;window, 0,xsize=800,ysize=425,XPOS = 950, YPOS = 300 
!p.multi = [0,2,2,0,1]

ns1=intarr(1)
tst1=double(1)
ted1=double(1)

nt1=intarr(1)
nz1=intarr(1)

close,3
openr,3,'/data/ap1vf/3D_509_36_36_300s_te_t.out',/f77_unf

readu,3,tst1,ted1,ns1
readu,3,nz1
print, ns1

timearr1=dblarr(ns1[0],nz1[0])
beta1=dblarr(ns1[0],nz1[0])
zz1=dblarr(nz1[0])

readu,3,timearr1
readu,3,beta1
readu,3,zz1

close,3

xx1=dblarr(ns1[0])

for i=0,ns1[0]-1 do xx1[i]=tst1+(ted1-tst1)/(ns1[0]*1.d0)*(i*1.d0)

; **********************************************************

ns2=intarr(1)
tst2=double(1)
ted2=double(1)

nt2=intarr(1)
nz2=intarr(1)

close,3
openr,3,'/data/ap1vf/3D_509_36_36_300s_mf_te_t.out',/f77_unf

readu,3,tst2,ted2,ns2
readu,3,nz2
print, ns2

timearr2=dblarr(ns2[0],nz2[0])
beta2=dblarr(ns2[0],nz2[0])
zz2=dblarr(nz2[0])

readu,3,timearr2
readu,3,beta2
readu,3,zz2

close,3

xx2=dblarr(ns2[0])

for i=0,ns2[0]-1 do xx2[i]=tst2+(ted2-tst2)/(ns2[0]*1.d0)*(i*1.d0)

; **********************************************************
; **********************************************************

ns3=intarr(1)
tst3=double(1)
ted3=double(1)

nt3=intarr(1)
nz3=intarr(1)

close,3
openr,3,'/data/ap1vf/3D_509_36_36_300s_te_t.out',/f77_unf

readu,3,tst3,ted3,ns3
readu,3,nz3
print, ns3

timearr3=dblarr(ns3[0],nz3[0])
beta3=dblarr(ns3[0],nz3[0])
zz3=dblarr(nz3[0])

readu,3,timearr3
readu,3,beta3
readu,3,zz3

close,3

xx3=dblarr(ns3[0])

for i=0,ns3[0]-1 do xx3[i]=tst3+(ted3-tst3)/(ns3[0]*1.d0)*(i*1.d0)

; **********************************************************
; **********************************************************

ns4=intarr(1)
tst4=double(1)
ted4=double(1)

nt4=intarr(1)
nz4=intarr(1)

close,3
openr,3,'/data/ap1vf/3D_509_36_36_300s_mf_te_t.out',/f77_unf

readu,3,tst4,ted4,ns4
readu,3,nz4
print, ns4

timearr4=dblarr(ns4[0],nz4[0])
beta4=dblarr(ns4[0],nz4[0])
zz4=dblarr(nz4[0])

readu,3,timearr4
readu,3,beta4
readu,3,zz4

close,3

xx4=dblarr(ns4[0])

for i=0,ns4[0]-1 do xx4[i]=tst4+(ted4-tst4)/(ns4[0]*1.d0)*(i*1.d0)

; **********************************************************

scale=1.0d6

; ****************** ps, eps begin ****************************************
   xs=19.d0
   k=12.d0/19.d0
   
   

;timearr(where(abs(timearr) le 0.5d0))=0.d0

;title='(!7q!N-!7q!I0!N)/!7q!I0!N!3'
tit_pd='V!I!Mx!N!3 [m/s]'
tit_pl='!3V!I| |!N!3 [km/s]'
tit_T='T-T!I0!N [K]'
tit_e='Vx'
tit_emag='e_mag'

zs=100
chsz=0.6
aa=zz3[zs]/scale

;stop

nup=300

;tvframe, timearr3[*,zs:nup], /sample, /bar, title=tit_e, $
;        xtitle='Time [s]', ytitle='z [Mm]',charsize=chsz, CT='mixct16', $
;        xrange=[xx3[0], xx3[ns3[0]-1]], yrange=[aa, max(zz3)/scale]


;tvframe, timearr4[*,zs:nup], /sample, /bar, title=tit_e, $
;        xtitle='Time [s]', ytitle='z [Mm]',charsize=chsz, CT='mixct16', $
;        xrange=[xx4[0], xx4[ns4[0]-1]], yrange=[zz4[zs]/scale, max(zz4)/scale]

tvframe, timearr1[*,zs:nup], /sample, /bar, title=tit_e, $
        xtitle='Time [s]', ytitle='z [Mm]',charsize=chsz, CT='mixct16', $
        xrange=[xx2[0], xx2[ns2[0]-1]], yrange=[zz1[zs]/scale, zz1[nup]/scale]

tvframe, timearr2[*,zs:nup], /sample, /bar, title=tit_e, $
        xtitle='Time [s]', ytitle='z [Mm]',charsize=chsz, CT='mixct16', $
        xrange=[xx2[0], xx2[ns2[0]-1]], yrange=[zz2[zs]/scale, zz2[nup]/scale]




end


