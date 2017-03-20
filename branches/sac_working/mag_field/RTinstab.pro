;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Create a sac domain for the 2D RT instability
;Created by Ben Snow
;Last modified: 2012/03/17
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
;WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


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

;It isn't necessary to load a file but this was easier 
;than writing the headline etc from scratch
;*****************************************************
openr,1,'../configs/zero1_ot_bin_256.ini',/f77_unf

readu,1,headline
readu,1,it,time,ndim,neqpar,nw
print,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Change a few of the key parameters, 
eqpar(0)=1.4         ;gamma
eqpar(3)=-0.1        ;gravity
print,eqpar
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Specify the number of cells in each dimension
nx(0)=512
nx(1)=512

n1=nx(0)
n2=nx(1)

x_code=dblarr(n1,n2,ndim)
w=dblarr(n1,n2,nw)

wi=dblarr(n1,n2)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Create variable x_code to define grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for i=0,n1-1 do begin
	for j=0,n2-1 do begin
		x_code(i,j,0)=1.0*i/n1
		x_code(i,j,1)=1.0*j/n2
	endfor
endfor

print, n1,n2

;*************************************************

gamma=1.4d0
mu=4.d0*!PI/1.0d7


y=reform(x_code(0,*,1))
x=reform(x_code(*,0,0))

print,max(x),min(x)
print,max(y),min(y)


;*************** Start bg density ******************
;Specify a smooth density profile between 1 and 10
;Variable 'grad' determines the steepness of the density change
rho=dblarr(n1,n2)
rho(*,*)=1.0
grad=80.0

for i=0,n2-1 do begin
;	if (x(i) LT 0.5) then rho(i,*)=1.0
;	if (x(i) GE 0.5) then rho(i,*)=2.0
	rho(i,*)=rho(i,*)+9.0/(1.0+exp(-grad*(x(i)-0.8)))
endfor

;*************** End bg density ******************


;*************** Start bg pressure ******************

p=dblarr(n1,n2)

for i=0,n2-1 do begin
	p[i,*]=2.5-0.1*rho[i,*]*x(i)
endfor


;*************** End bg pressure ******************

;*************** Begin bg magnetic field ******************

bx=dblarr(n1,n2)
by=dblarr(n1,n2)

bx(*,*)=1.0e-4*0.0125
;by(*,*)=1.0e-4*0.0125

;bx(*,*)=0.0
by(*,*)=0.0

;convert to VAC magnetic field

bx(*,*)=bx(*,*)/sqrt(mu)
by(*,*)=by(*,*)/sqrt(mu)

;*************** End bg magnetic field ******************

;*************** Begin bg energy ******************

e=p/(gamma-1.d0)+0.5d0*(bx*bx+by*by)

;*************** End bg energy ******************

;*************** Begin momentum perturbation ******************
;Generate a 2D random profile
A=randomn(5,n1,n2)
;Specify the maximum velocity magnuitude
vmax=0.01

for i=0,n2-1 do begin 

;use the sin functions to centre the maximum velocity perturbations
;Not really necessary but can help prevent boundary effects
	w(*,i,1)=A(*,i)*vmax*rho(*,i)*sin(3.14*x(*))*sin(3.14*y(i))

endfor

loadct,70
contour,reform(w(*,*,1)),/fill,nlevels=101

;*************** End momentum perturbation ******************

;'h m1 m2 e b1 b2 eb rhob bg1 bg2'

;w(*,*,1)=0.0
w(*,*,2)=0.0

;Energy perturbation
w(*,*,3)=0.5*rho*w(*,*,1)^2.0

w(*,*,4)=0.0
w(*,*,5)=0.0

w(*,*,7)=rho
w(*,*,6)=e 

w(*,*,8)=bx
w(*,*,9)=by


close,1
openw,1,'../configs/RT512M.ini',/f77_unf
writeu,1,headline
writeu,1,it,time,ndim,neqpar,nw
writeu,1,nx
writeu,1,eqpar
writeu,1,varname
writeu,1,x_code
for iw=0,nw-1 do begin
wi=w(*,*,iw)
writeu,1,wi
endfor

print,max(w(*,*,1)),min(w(*,*,1))
 
close,1


print, 'complete'
end





