tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


ii=1

if (ii eq 1) then begin
;loadct,4
;mixct
endif else begin
loadct,0
tek_color
endelse




mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)

ia=1.0

headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)

; Open an MPEG sequence: 
;mpegID = MPEG_OPEN([700,1200],FILENAME='myMovie.mpg') 

window, 0,xsize=500,ysize=500,XPOS = 950, YPOS = 300 



nn=0
np=0
kkk=4

nn_i=0

close,1
close,2


;openr,1,'/data/ap1vf/3D_509_36_36_300s.out',/f77_unf
;openr,1,'/data/ap1vf/3D_tube_196_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60t.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

openr,1,'/data/ap1vf/3D_tube_196_100_100_120s_full.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_200s_puls.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf

while not(eof(1)) do begin
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


print, 'tuta1'
xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

wi=dblarr(n1,n2,n3)

w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

xx=dblarr(n2)
yy=dblarr(n3)
zz=dblarr(n1)


xx(*)=x(1,*,1,1)
yy(*)=x(1,1,*,2)
zz(*)=x(*,1,1,0)

Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)


mu=4.0*!PI/1.0e7

print,'******************* time = ', time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

scale=1.d6

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=99
ystart=0
yend=99

;xstart=3
;xend=65
;ystart=39
;yend=65


hh=111 ;z

wset,0
!p.multi = [0,1,1,0,1]


zstart=0
zend=180

wt=dblarr(xend-xstart+1,yend-ystart+1,iw)
wt=reform(w(hh,xstart:xend,ystart:yend,*))

bx=dblarr(xend-xstart+1,yend-ystart+1)
by=dblarr(xend-xstart+1,yend-ystart+1)
bz=dblarr(xend-xstart+1,yend-ystart+1)

surf=dblarr(xend-xstart+1,yend-ystart+1)

bx(*,*)=(rotate(wt(*,*,5),1)+rotate(wt(*,*,10),1))*sqrt(mu)*1.0e4
by(*,*)=(rotate(wt(*,*,6),1)+rotate(wt(*,*,11),1))*sqrt(mu)*1.0e4
bz(*,*)=(rotate(wt(*,*,7),1)+rotate(wt(*,*,12),1))*sqrt(mu)*1.0e4



surf=sqrt(bx^2.d0+by^2.d0+bz^2.d0)

;tvframe,bx,/bar,/sample, title='bz', $
;        xtitle='x', ytitle='z', charsize=2.0

;tvframe,by,/bar,/sample, title='bx', $
;        xtitle='x', ytitle='z', charsize=2.0

;tvframe,bz,/bar,/sample, title='by', $
;        xtitle='x', ytitle='z', charsize=2.0

contour, surf, /follow ;,levels=[15.1,15.4,15.8,16.2,16.6,17.0,17.5]

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200	


 
indexs=strtrim(np,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,'/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/area/h111_B/'+indexss+'.png',image_p, red,green, blue


np=np+1



endwhile



end
