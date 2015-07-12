DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW

window, 0,xsize=600,ysize=900,XPOS = 1300, YPOS = 300 

loadct, 4

!p.multi = [0,2,3,0,1]

n1=100
n2=100
x=dblarr(n1,n2)
y=dblarr(n1,n2)

xi=dblarr(n1)
yi=dblarr(n2)

for i=0, n1-1 do begin
 for j=0, n2-1 do begin
  x[i,j]=2.d0*i/100.d0-1.d0
  y[i,j]=2.d0*j/100.d0-1.d0
endfor
endfor

xi(*)=x(*,1)
yi(*)=y(1,*)


Vr=dblarr(n1,n2)
Vphi=dblarr(n1,n2)

xx=dblarr(n1,n2)
yy=dblarr(n1,n2)

xx_p=dblarr(n1,n2)
yy_p=dblarr(n1,n2)


x0n=dblarr(2)
y0n=dblarr(2)


x0n[0]=0.05
y0n[0]=0.05

x0n[1]=-0.1
y0n[1]=-0.1

xm=dblarr(1)
ym=dblarr(1)

xm1=dblarr(1)
ym1=dblarr(1)

xm=x0n[0]
ym=y0n[0]

xm1=x0n[1]
ym1=y0n[1]

xx_p=x
yy_p=y

dt=0.01d0

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




while (qt le 1000) do begin

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
; tvframe, Vr, /bar, title='Vr', charsize=2.d0
; tvframe, Vphi, /bar, title='Vphi', charsize=2.d0
hxmin=0
hymin=0

hxmax=n1-1
hymax=n2-1


for i=0, n1-1 do begin
 for j=0, n2-1 do begin
   xx[i,j]=xx_p[i,j]+Vx[i,j]*dt
   yy[i,j]=yy_p[i,j]+Vy[i,j]*dt
 endfor
endfor



 tvframe, xx, /bar, title='x', charsize=2.d0

xx_p=xx
yy_p=yy


;***************** begin start NEW positions ***************************
 for ii=0, n_elements(x0n)-1 do begin 
     xc=interpol(dindgen(n1),xi, x0n[ii])
     yc=interpol(dindgen(n2),yi, y0n[ii])

     Vvx=interpolate(Vx,xc,yc)
     Vvy=interpolate(Vy,xc,yc)
     
     x0n[ii]=x0n[ii]+Vvx*dt
     y0n[ii]=y0n[ii]+Vvy*dt
     
     print, '*', x0n[ii], y0n[ii], Vvx, VVy, ii       

endfor

   xm=[xm,x0n[0]]
   ym=[ym,y0n[0]]

   xm1=[xm1,x0n[1]]
   ym1=[ym1,y0n[1]]

PLOT, xm,ym, psym=3, charsize=2, XRANGE=[-1.0, 1.0], YRANGE=[-1.0, 1.0]
oplot, xm1,ym1, psym=3, color=200 
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
