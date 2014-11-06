tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


loadct,0
tek_color


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

window, 1,xsize=1200,ysize=200,XPOS = 20, YPOS = 80 

nn=0
kkk=4

nn_i=0

close,1
close,2

openr,1,'/data/ap1vf/hasan/mod1.20000.out',/f77_unf

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

wi=dblarr(n1,n2,n3)

w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
endfor


Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)


mu=4.0*!PI/1.0e7

print,time


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

pp=50  ;x
kk=50  ;y

wset,0
!p.multi = [0,4,4,0,1]

wt=dblarr(n1,n2,iw)

wt=reform(w(*,*,pp,*))



tvframe,alog10(w(*,*,0)), /bar,title='log rho',/sample, xtitle='x', ytitle='y',charsize=2.0  

tvframe,wt(*,*,1)/(wt(*,*,0)),/sample, /bar,title='v1', xtitle='x', ytitle='y',charsize=2.0  


tvframe,wt(*,*,2)/(wt(*,*,0)),/sample, /bar,title='v2',xtitle='x', ytitle='z',charsize=2.0 

tvframe,wt(*,*,3)/(wt(*,*,0)),/sample, /bar,title='v3',xtitle='x', ytitle='z',charsize=2.0 

tvframe,wt(*,*,4),/sample, title='e', xtitle='x', ytitle='z', charsize=2.0                                                                                                   

tvframe,wt(*,*,5),/bar,/sample, title='Bx', xtitle='x', ytitle='z', charsize=2.0

tvframe,wt(*,*,6),/bar,/sample, title='By', xtitle='x', ytitle='z', charsize=2.0

tvframe,wt(*,*,7),/bar,/sample, title='Bz', xtitle='x', ytitle='z', charsize=2.0

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200

endwhile


end
