function BZxy, x, x0, y,y0,delta,A

res=A*exp(-3.d0*(x^2.d0/x0^2.d0+y^2.d0/y0^2.0d0))

;res=res-0.4*A*exp(-3.d0*(x-delta)^2.d0/(x0/2.d0)^2.d0)
;res=res-0.4*A*exp(-3.d0*(x+delta)^2.d0/(x0/2.d0)^2.d0)

return,res
end

; Define the Fxy function. 
FUNCTION Fxy, x, y 
   RETURN, y*COS(x^5) 
END 
 
; Define the limits of integration for y as a function of x: 
FUNCTION PQ_Limits, x 
   RETURN, [0.0, x^2] 
END 

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

; Define limits of integration for x: 
AB_Limits = [0.0, 2.0] 

window,0
window,1

;***************
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
close,1
openr,1,'/data/cs1mkg/VAC_NN/2_6Mnzx1976400.ini',/f77_unf
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname
n1=nx(0)
n2=nx(1)
n3=nx(2)
x_code=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)   ;was n1,n2,n3,nw
wi=dblarr(n1,n2,n3)
readu,1,x_code
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,*,iw)=wi
endfor

print,n1,n2,n3

;************** n1----z,  n2----x, n3----y
;stop


bx=dblarr(n1,n2,n3)
by=dblarr(n1,n2,n3)
bz=dblarr(n1,n2,n3)

phi=dblarr(n1,n2,n3)

y=reform(x_code(0,0,*,2))
x=reform(x_code(0,*,0,1))
z=reform(x_code(*,0,0,0))
rho=reform(w(*,*,*,0)+w(*,*,*,7))
rho1=dblarr(n1,n2,n3)

b0z=dblarr(n2,n3)

x=x-max(x)/2.d0
y=y-max(y)/2.d0

x0=2.0d6
y0=1.0d6
delta=3.0d6
A=0.5d6

for i=0,n2-1 do begin ;x
 for j=0,n3-1 do begin   ;y
   b0z(i,j)=BZxy(x(i),x0,y(j),y0,delta,A)
 endfor
endfor 

tvframe, b0z
stop

n_shift=30
n1l=n1+n_shift


zl=congrid(z,n1l, /interp)

xzl=dblarr(n1l,n2)
for k=0,n1l-1 do begin; z
for i=0,n2-1 do begin ;x
 for j=0,n3-1 do begin   ;y
 
     xzl(k,i,j)=int_tabulated_2D(x,y,b0z/sqrt((x(j)-x)^2.d0+(y(j)-y)^2.d0+zl(i)^2.d0))
 
 endfor
endfor 
tvframe, xzl(k,*,*)
endfor
   
end
