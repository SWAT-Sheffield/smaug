DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW

window, 0,xsize=600,ysize=900,XPOS = 1300, YPOS = 300 

loadct, 4

!p.multi = [0,2,3,0,1]

n1=100
n2=100
x=dblarr(n1,n2)
y=dblarr(n1,n2)

for i=0, n1-1 do begin
 for j=0, n2-1 do begin
  x[i,j]=2.d0*i/100.d0-1.d0
  y[i,j]=2.d0*j/100.d0-1.d0
endfor
endfor

Vr=dblarr(n1,n2)
Vphi=dblarr(n1,n2)

vvv=dblarr(n1,n2)

Vx=dblarr(n1,n2)
Vy=dblarr(n1,n2)

r=sqrt(x^2.d0+y^2.d0)

deltar=0.1d0
r0=0.1

deltax=0.1d0
x0=0.1
deltay=0.1d0
y0=0.1

s_period=30.d0
delta_qt=0.5
A=10.0
qt=0.d0




while (qt le 100) do begin

for i=0, n1-1 do begin
 for j=0, n2-1 do begin
   Vr[i,j]=0.d0
   Vphi[i,j]=A*exp(-(r[i,j]-r0)^2.d0/deltar^2.d0)*sin(qt*2.d0*!pi/s_period)
 endfor
endfor

for i=0, n1-1 do begin
 for j=0, n2-1 do begin

 Vx[i,j]=A*exp(-(r[i,j]-r0)^2.d0/deltar^2.d0)*sin(qt*2.d0*!pi/s_period) 
 Vy[i,j]=-A*exp(-(r[i,j]-r0)^2.d0/deltar^2.d0)*cos(qt*2.d0*!pi/s_period) 

;   Vr1=0.d0
;   Vphi1=A*exp(-(r[i,j]-r0)^2.d0/deltar^2.d0)*sin(qt*2.d0*!pi/s_period)
 
 
;  if (r[i,j] ne 0.0d0) then begin
;   Vx[i,j]=(Vr1*x[i,j]-Vphi1*y[i,j])/r[i,j]
;   Vy[i,j]=(Vr1*y[i,j]+Vphi1*x[i,j])/r[i,j]
;  endif else begin
;   Vx[i,j]=0.d0
;   Vy[i,j]=0.d0
;  endelse 
 endfor
endfor




tvframe, Vx, /bar, title='Vx', charsize=2.d0
tvframe, Vy, /bar, title='Vy', charsize=2.d0
tvframe, sqrt(Vx^2.d0+Vy^2.d0), /bar, title='V', charsize=2.d0
tvframe, Vr, /bar, title='Vr', charsize=2.d0
tvframe, Vphi, /bar, title='Vphi', charsize=2.d0
hxmin=0
hymin=0

hxmax=n1-1
hymax=n2-1


savx=dblarr(n1,n2)
savy=savx

savx=Vx
savy=Vy


nxy=10
nxny = [nxy,nxy]

hxx=dblarr(n1)
hyy=dblarr(n2)

hxx=x(*,0)
hyy=y(0,*)

xxi=interpol(hxx,nxy)
yyi=interpol(hyy,nxy)

avxi=congrid(savx,nxy,nxy)
avyi=congrid(savy,nxy,nxy)
tek_color
VELOVECT, avxi, avyi, xxi, yyi, charsize=2 ;, /OVERPLOT ;, color=6

 ss='time ='+strTrim(string(qt),1)
 xyouts,50,2, ss, /device, color=200
 
qt=qt+delta_qt
wait, 0.05
endwhile

end
