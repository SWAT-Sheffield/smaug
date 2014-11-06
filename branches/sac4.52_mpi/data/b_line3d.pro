function B_x,x_i,y_i,z_i
 RETURN,z_i
end  
function B_y,x_i,y_i,z_i
 RETURN,1
end  

function B_z,x_i,y_i,z_i
 RETURN,x_i
end  



print, 'data or function ? (data "1", function "2")'
read,aa


;**********  begin data from file
if aa eq 1 then begin
close,1
;openr,1,'../../../../../data/ap1vf/MHD33_b.ini',/f77_unf
openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

tarr=dblarr(1)
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


readu,1,headline
print,'head = ', headline
readu,1,it,time,ndim,neqpar,nw
print,'it, time', it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
print, '*',ndim 
nx=lonarr(ndim)
readu,1,nx
print, 'nx = ', nx
eqpar=dblarr(neqpar)
readu,1,eqpar
print,eqpar
readu,1,varname
print,varname

n1=nx(0) 
n2=nx(1)
n3=nx(2)

xx=dblarr(n1,n2,n3,ndim)

readu,1,xx

close,1

endif else begin

n1=100
n2=100
n3=100

endelse

bx=dblarr(n1,n2,n3)
by=dblarr(n1,n2,n3)
bz=dblarr(n1,n2,n3)

;if aa eq 1 then begin

;************************** B E G I N   F I L E *******************************

; openr,1,'B_t_3.dat',/f77_unf
;  readu,1,bx,by,bz
; close,1

;************************** E N D   F I L E *******************************

;endif

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW

window, xsize=600,ysize=600

loadct, 12

;mixct1

!p.multi=[0,1,0,0,1]


x=dblarr(n1)
y=dblarr(n2)
z=dblarr(n3)


;****************end data from file 


x=dblarr(n1)
y=dblarr(n2)
z=dblarr(n3)

for i=0,n1-1 do begin
  if aa eq 1 then begin
    x(i)=xx(i,0,0,0)
  endif else begin
   x(i)=(i-n1/2.0)
  endelse 
 endfor

for j=0,n2-1 do begin
  if aa eq 1 then begin
    y(j)=xx(0,j,0,1)
  endif else begin
    y(j)=(j-n2/2.0)
  endelse 
 endfor

for k=0,n3-1 do begin
  if aa eq 1 then begin
    z(k)=xx(0,0,k,2)
  endif else begin
    z(k)=(k-n3/2.0)
  endelse  
endfor

if aa ne 1 then begin
  for i=0,n1-1 do begin
   for j=0,n2-1 do begin
    for k=0,n3-1 do begin 
      bx(i,j,k)=B_x(x(i),y(j),z(k))
      by(i,j,k)=B_y(x(i),y(j),z(k))
      bz(i,j,k)=B_z(x(i),y(j),z(k))  
    endfor
   endfor
  endfor
endif 

xyz_color=dblarr(1)

mi_x=min(x)
ma_x=max(x)

mi_y=min(y)
ma_y=max(y)

mi_z=min(z)
ma_z=max(z)

b_total=sqrt(bx(*,*,0)^2.0+by(*,*,0)^2.0+bz(*,*,0)^2.0)

tvframe,b_total, /bar

max_b=max(b_total)
min_b=min(b_total)

db=max_b-min_b
print,db


; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL 


;******************************************************
nn=10.0

;******************************************************



xm=dblarr(1)
ym=dblarr(1)
zm=dblarr(1)


print,'mi_x, ma_x, mi_y,ma_y,mi_z,ma_z = ', mi_x, ma_x, mi_y,ma_y,mi_z,ma_z


x0=mi_x
y0=mi_y
z0=20.0 ;mi_z

x0_s=x0
y0_s=y0
z0_s=z0

xm(0)=x0
ym(0)=y0
zm(0)=z0

xyz_color(0)=1.0

color_min=100.0
color_max=256.0


dLx=1.1 ;500.0
dLy=1.0 ;501.0
dLz=0.9 ;3009.0

df=1.0
f=1.0

C=50.0


delta_x=(ma_x-mi_x)/nn
delta_y=(ma_y-mi_y)/nn
delta_z=(ma_z-mi_z)/4.0


SURFACE, DIST(5), /NODATA, /SAVE, XRANGE=[mi_x, ma_x], $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z], XSTYLE=1, $
   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, xtitle='x',ytitle='y',title='Field Lines', $
   ztitle='z', POSITION=[0.2, 0.1, 0.95, 0.95, 0.1, 0.95], az=30.0, ax=30.0  




;while z0_s le ma_z do begin

while y0_s le ma_y do begin

while x0_s le ma_x do begin


   
while (x0 ge mi_x and x0 le ma_x and y0 ge mi_y and y0 le ma_y and z0 ge mi_z and z0 le ma_z) do begin 

    i=0
    j=0
    k=0
     while (x0-x(i))*(x0-x(i+1)) gt 0.0 do begin     
      i=i+1
     endwhile    
    ii=i+(x0-x(i))/(x(i+1)-x(i))
  ;  print,'i',i, x0,x(i),x(i+1),ii, (x0-x(i))*(x0-x(i+1))
 
     while (y0-y(j))*(y0-y(j+1)) gt 0.0 do begin         
     j=j+1
     endwhile
    jj=j+(y0-y(j))/(y(j+1)-y(j))
;    print,'j',j , y0,y(j),y(j+1),jj,(y0-y(j))*(y0-y(j+1))     

     while (z0-z(k))*(z0-z(k+1)) gt 0.0 do begin         
     k=k+1
     endwhile
    kk=k+(z0-z(k))/(z(k+1)-z(k))
;    print,'k',k , z0,z(k),z(k+1),kk,(z0-z(k))*(z0-z(k+1))     

;read,aa

    Bx_inter=INTERPOLATE(bx,ii,jj,kk)
    By_inter=INTERPOLATE(by,ii,jj,kk)
    Bz_inter=INTERPOLATE(bz,ii,jj,kk)
;print,'Bx,By,Bz', Bx_inter,Bz_inter,Bz_inter
   sq=sqrt(Bx_inter^2.0+By_inter^2.0+Bz_inter^2.0)

cl=color_max-(max_b-sq)/(max_b-min_b)*(color_max-color_min)
   
   xyz_color=[xyz_color,cl]
   
   dx=-df*dLx*Bx_inter/sq
   dy=-df*dLy*By_inter/sq
   dz=-df*dLz*Bz_inter/sq   

;print, dx,dy,dz
   x0=x0+dx
   y0=y0+dy
   z0=z0+dz
   
;print,'f= ',f
 if f eq 1.0 then begin
   x0_pre=x0
   y0_pre=y0
   z0_pre=z0   
   f=-1.0
 endif 
   
;read,aa
   xm=[xm,x0]
   ym=[ym,y0]
   zm=[zm,z0]
  
   
; PLOTS, xm,ym,zm, psym=3, /T3D 
endwhile


 
 if df eq -1.0 then begin
   f=1.0
 endif
  
;   print,'######*****',x0_s,y0_s
     
 if f eq -1.0 then begin
  df=-1.0   
    x0=x0_pre
    y0=y0_pre
    z0=z0_pre  
    
 endif else begin 
  df=1.0 
;*********************************************
;i=0
;j=0
;k=0
;     while (x0_s-x(i))*(x0_s-x(i+1)) gt 0.0 do begin     
;      i=i+1
;     endwhile    
;    ii=i+(x0_s-x(i))/(x(i+1)-x(i))
 
;     while (y0_s-y(j))*(y0_s-y(j+1)) gt 0.0 do begin         
;     j=j+1
;     endwhile
;    jj=j+(y0_s-y(j))/(y(j+1)-y(j))


;     while (z0_s-z(k))*(z0_s-z(k+1)) gt 0.0 do begin         
;     k=k+1
;     endwhile
;    kk=k+(z0_s-z(k))/(z(k+1)-z(k))

;    Bx_inter=INTERPOLATE(bx,ii,jj,kk)
;    By_inter=INTERPOLATE(by,ii,jj,kk)
;    Bz_inter=INTERPOLATE(bz,ii,jj,kk)
;print,'Bx,By,Bz', Bx_inter,Bz_inter,Bz_inter
;   sq=Bx_inter^2.0+By_inter^2.0+Bz_inter^2.0
;    
;    a_delta_x=C*Bx_inter/sq 
;    a_delta_y=C*By_inter/sq
     
;    if abs(a_delta_x) ge 0.5 then begin
;      delta_x=abs(a_delta_x)
;    endif  
;    
;    if abs(a_delta_y) ge 0.5 then begin
;      delta_y=abs(a_delta_y)
;    endif  
    
;    print,' delta_x = ' , delta_x, '  delta_y=', delta_y, ' x0_s=', x0_s,' y0_s=', y0_s,' sq=',sq, Bx_inter,By_inter 
    
;*********************************************  

    x0_s=x0_s+delta_x
;    y0_s=y0_s+delta_y
    x0=x0_s
    y0=y0_s 
    z0=z0_s
 endelse
 ;PLOTS, xm,ym,zm, color=xyz_color, psym=3, /T3D 
 
;   print,'######',x0,y0,f, df, delta_x 
;read , aa 

endwhile
  x0_s=mi_x
  x0=x0_s
  y0_s=y0_s+delta_y
  y0=y0_s

endwhile
; z0_s=z0_s+delta_z
; z0=z0_s
; x0_s=mi_x
; x0=x0_s
; y0_s=mi_y
;y0=y0_s
; print,'z =', z0

PLOTS, xm,ym,zm, color=xyz_color, psym=3, /T3D 

;endwhile
 ; oplot,xm,ym, psym=3.0   
 
end
