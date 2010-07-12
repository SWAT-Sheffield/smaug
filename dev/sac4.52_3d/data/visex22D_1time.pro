tarr=dblarr(1)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


loadct,4
;mixct
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

nn=0
nni=4

window, 0,xsize=1000,ysize=300
nn=0
kkk=10

nn_i=1

close,1
close,2


;openr,1,'/data/ap1vf/gorz1604n.out',/f77_unf
;openr,1,'/data/ap1vf/opozmf.out',/f77_unf
;openr,1,'/data/ap1vf/test.out',/f77_unf
openr,1,'/data/ap1vf/goriz1024100sm.out',/f77_unf
;openr,1,'/data/ap1vf/tube256.out',/f77_unf


zstart=250
zend=350

timeend=2400
timeout=dblarr(timeend,zend-zstart+1)



while not(eof(1)) do begin

;set_plot, 'ps'
;loadct, 3

;device, filename='sol_.ps', /color, BITS=8
;DEVICE, XSIZE=12, YSIZE=3, /INCHES

;!p.thick = 4
;!x.thick = 4
;!y.thick = 4
;!z.thick = 4
;!p.font = 1.0

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




xout=dblarr(2)
yout=dblarr(2)


n1=nx(0)
n2=nx(1)
x=dblarr(n1,n2,ndim)
if (nn eq 0) then w=dblarr(n2,n1,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2)
;e2=dblarr(n1,n2)
readu,1,x
;readu,1,e2
;e2=rotate(e2,1)
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=rotate(wi,1)
endfor


Vt=dblarr(n2,n1)

B=dblarr(n2,80)
B_bg=dblarr(n1,n2)

p=dblarr(n1,n2,1)
;e2=dblarr(n1,n2)




;!p.multi = [0,2,3,0,1]
!p.multi = [0,2,0,0,1]

mu=4.0*!PI/1.0e7

print,time

kk=128

label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'


R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

; gauss Tesla
factor=100000.0

xstart=0
xend=511


label_Vz='Velocity '+'!19V!X!DZ'+'!N'+' (m/s)'
label_Vx='Velocity '+'!19V!X!DX'+'!N'+' (m/s)'
label_Vt='Velocity Total '+'!19V!X!DT'+'!N'+' (m/s)'

label_bz='Magnetic Field '+'!19b!X!DZ'+'!N'+' (Gauss)'
label_bx='Magnetic Field '+'!19b!X!DX'+'!N'+' (Gauss)'
label_bt='Total magnetic field '+'!19B!X!DT'+'!N'+' (Gauss)'



n11=zend-zstart+1
n22=n2

Bxx=dblarr(n22,n11)
Bzz=dblarr(n22,n11)

xx=max(x(1,*,1))
zz=max(x(zend,1,0))



Bxx(*,*)=reform(w(*,zstart:zend,5)+w(*,zstart:zend,9))
Bzz(*,*)=reform(w(*,zstart:zend,4)+w(*,zstart:zend,8))



btot=sqrt(Bxx^2.0+Bzz^2.0)*SQRT(mu)*factor

for j=1,zend-zstart do begin
timeout(nn_i,j)=btot(35,j)
endfor



tvframe,btot,$
        /bar,/sample, title=label_bt, $ ;charsize=1.0,$  
	xtitle='Horizontal distance (Mm)', xrange=[0,16], ytitle='Height (Mm)', yrange=[1.96,2.75]
tek_color	

oplot, [5,5],[1.9,3.0], linestyle=2, color=255 ;, thick=2	

print, n1,n2

tek_color

line2_b, n22,n11, Bxx,Bzz,xx,zz,0.d0,x(zstart,1,0)

tvframe, timeout, xrange=[0,907], yrange=[1.96,2.75], ytitle='Height (Mm)',xtitle='Time (s)',$
         charsize=1.0, title='Time-distance'


 ss='time ='+strTrim(string(time),1)  ;+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   
stop

if (nni eq 4) then begin

image_p = TVRD_24()
 ; write_png,'/data/ap1vf/png/time2/'+indexss+'.png',image_p, red,green, blue
  nn=nn+1
  nni=0

endif 

nni=nni+1
nn_i=nn_i+1
print,'nni=',nni

;device, /close
;set_plot, 'x'

endwhile




end
