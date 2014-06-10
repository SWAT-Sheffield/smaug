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

window, 0,xsize=900,ysize=725,XPOS = 350, YPOS = 300 

window, 1,xsize=600,ysize=200,XPOS = 20, YPOS = 80 

nn=0
kkk=4

nn_i=0

close,1
close,2

;openr,1,'/data/ap1vf/goriz1024100Bsm.out',/f77_unf
;openr,1,'/data/ap1vf/gorz16300_1024.out',/f77_unf
;openr,1,'/data/ap1vf/8192.ini',/f77_unf 
;openr,1,'/data/ap1vf/4096.out',/f77_unf
;openr,1,'/data/ap1vf/zero_np0104_001.out',/f77_unf
openr,1,'/data/ap1vf/zero_bach16.out',/f77_unf
;openr,1,'/data/ap1vf/100.out',/f77_unf


while not(eof(1)) do begin
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname




xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)
if (nn eq 0) then w=dblarr(n3,n2,n1,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2,n3)
;e2=dblarr(n1,n2,n3)
readu,1,x
;readu,1,e2
;e2=rotate(e2,1)
for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
; w(*,*,*,iw)=rotate(wi,1)
endfor


Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)
;e2=dblarr(n1,n2)






mu=4.0*!PI/1.0e7

print,time

kk=64

label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'


R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=511

pp=63

wset,0
!p.multi = [0,4,4,0,1]
if (ii eq 1) then begin




tvframe,alog10(w(*,pp,*,9)+w(*,pp,*,0)), /bar,title='log rhob+rho',/sample, xtitle='x', ytitle='y',charsize=2.0  

tvframe,w(*,pp,*,1)/(w(*,pp,*,0)+w(*,pp,*,9)),/sample, /bar,title='v1', xtitle='x', ytitle='y',charsize=2.0  


tvframe,w(pp,*,*,2)/(w(pp,*,*,0)+w(pp,*,*,9)),/sample, /bar,title='v2',xtitle='x', ytitle='z',charsize=2.0 

tvframe,w(*,pp,*,3)/(w(*,pp,*,0)+w(*,pp,*,9)),/sample, /bar,title='v3',xtitle='x', ytitle='z',charsize=2.0 

tvframe,w(pp,*,*,4),/sample, title='e', xtitle='x', ytitle='z', charsize=2.0                                                                                                   

tvframe,w(*,*,pp,8),/sample, title='eb',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,pp,0),/bar,/sample, title='rho', xtitle='x', ytitle='z', charsize=2.0

endif else begin

endelse



wset,1
!p.multi = [0,1,1,0,1]
plot,w(kk,*,pp,1)/(w(kk,*,pp,0)+w(kk,*,pp,9)),title='v1', xtitle='x', ytitle='y',charsize=1.0  
oplot,w(kk,*,pp,1)/(w(kk,*,pp,0)+w(kk,*,pp,9)),psym=4, color=100

 time=time/31536000.d0
 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200
 
print,  kkk
if (kkk eq 4) then begin

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
;write_png,'/data/ap1vf/png/test2/'+indexss+'.png',image_p, red,green, blue


nn=nn+1

kkk=1

endif

kkk=kkk+1
;if (ii eq 2) then read,a

;maxa=[maxa,max(w(*,*,2)/(w(*,*,0)+w(*,*,7)))]
maxa=[maxa,max(Vt)]




endwhile





;data=reform(alog10(w(*,*,*,0)+w(*,*,*,9))) 

data=reform(w(*,*,*,4)) 

;data=reform(w(*,*,*,3)/(w(*,*,*,0)+w(*,*,*,9)))

datax=reform(x(*,*,*,0)) 
datay=reform(x(*,*,*,1)) 
dataz=reform(x(*,*,*,2)) 

;data=(data-min(data))/(max(data)-min(data))

;ivolume, data

hData = PTR_NEW(data, /NO_COPY) 
 

SLICER3, hdata, DATA_NAMES='Dave' 

end
