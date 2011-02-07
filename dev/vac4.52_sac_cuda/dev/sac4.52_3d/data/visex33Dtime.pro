tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)

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
kout=20

nn_i=0

close,1


nt=162
nxx=396

ns=intarr(1)
tst=double(1)
ted=double(1)


timearr=dblarr(nt,nxx)
beta=dblarr(nt,nxx)
Cs=dblarr(nt,nxx)

print, 'tuta'

;openr,1,'/data/ap1vf/3D_509_36_36_300s_mf.out',/f77_unf
openr,1,'/data/ap1vf/3D_396_60_60_tube1.out',/f77_unf

!p.multi = [0,1,0,0,1]

file_N=1


jump1 :

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

if nn eq 0 then tst=time

xout=dblarr(2)
yout=dblarr(2)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

w=dblarr(n1,n2,n3,nw)   ;was n1,n2,nw

wi=dblarr(n1,n2,n3)
readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,*,iw)=wi
endfor

zmin=min(x(*,1,1,0))
zmax=max(x(*,1,1,0))

xmin=min(x(1,*,1,1))
xmax=max(x(1,*,1,1))


xx=dblarr(zend)

xx(zstart:zend-1)=x(zstart:zend-1,1,1,0)

scale=1.0d6

Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)
;e2=dblarr(n1,n2,n3)


mu=4.0*!PI/1.0e7

print,time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

arho_t=w(*,*,*,9)+w(*,*,*,0)
arho0=w(*,*,*,9)

avz=w(*,*,*,1)/arho_t

avx=w(*,*,*,2)/arho_t
avy=w(*,*,*,3)/arho_t

amz=w(*,*,*,1)
amx=w(*,*,*,2)
amy=w(*,*,*,3)

ae=w(*,*,*,4)

aeb=w(*,*,*,8)

arho=w(*,*,*,0)

abz=w(*,*,*,5)
abx=w(*,*,*,6)
aby=w(*,*,*,7)

abz_t=w(*,*,*,5)+w(*,*,*,10)
abx_t=w(*,*,*,6)+w(*,*,*,11)
aby_t=w(*,*,*,7)+w(*,*,*,12)

abz0=w(*,*,*,10)
abx0=w(*,*,*,11)
aby0=w(*,*,*,12)

R=8.3e+003
mu_gas=0.6
gamma=1.66667

zstart=0
zend=196


tmin=0.d0
tmax=410

kx=34
ky=26
;aa=dblarr(n1,n3)

;for i=0,n2-1 do begin
; for j=0,n3-1 do begin
;   vv=(avx[*,i,j]^2.d0+avy[*,i,j]^2.d0+avz[*,i,j]^2.d0)*arho_t[*,i,j]/2.d0
;   timearr[nn,*]=timearr[nn,*]+vv*4.d6*4.d6/35.d0/35.d0 ;+aeb[*,i,j]  
; endfor
;endfor

;aa(*,*)=total(ae+aeb, 2)

;for i=0,n2-1 do begin
; for j=0,n3-1 do begin
;   timearr[nn,*]=timearr[nn,*]+avx[*,i,j]  
; endfor
;endfor


;timearr[nn,*]=(avx[*,kx,ky]^2.d0+avy[*,kx,ky]^2.d0+avz[*,kx,ky]^2.d0)*arho_t[*,kx,ky]/2.d0

;timearr[nn,*]=ae[*,kx,ky]
timearr[nn,zstart:zend]=avx[zstart:zend,kx,ky]



;****************** Pressure begin **********************
T=ae[*,kx,ky]+aeb[*,kx,ky]
T=T-(avx[*,kx,ky]^2.0+avy[*,kx,ky]^2.0+avz[*,kx,ky]^2.0)*arho_t[*,kx,ky]/2.0
T=T-(abx_t[*,kx,ky]^2.0+aby_t[*,kx,ky]^2.0+abz_t[*,kx,ky]^2.0)/2.0
T=(gamma-1.d0)*T

;TT=saeb
;TT=TT-(sabx0^2.0+sabz0^2.0)/2.0
;TT=(gamma-1.d0)*TT
;****************** Pressure end ************************

;timearr[nn,*]=T[*]

nn=nn+1

;beta[nn,*]= 0.d0 ;reform((abx_t[*,kx,ky]^2.0+aby_t[*,kx,ky]^2.0+abz_t[*,kx,ky]^2.0)/2.d0/T[*])


ted=time

ns[0]=nn
timearr_out=dblarr(ns[0], nxx)
beta_out=dblarr(ns[0],nxx)

for i=0,ns[0]-1 do timearr_out[i,*]=timearr[i,*]
;for i=0,ns[0]-1 do beta_out[i,*]=beta[i,*]

tvframe, timearr_out, /sample, /bar,  CT='mixct16'

endwhile

close,2
openw,2,'/data/ap1vf/3D_396_60_60_t.out',/f77_unf
writeu,2,tst,ted,ns
writeu,2,zend
writeu,2,timearr_out
;writeu,2,beta_out
writeu,2,xx
close,2




close, 1

;file_N=file_N+1

;if (file_N eq 2) then begin
; openr,1,'/data/ap1vf/3D_509_36_36_300s_mf_1.out',/f77_unf
; goto, jump1
;endif

;if (file_N eq 3) then begin
; openr,1,'/data/ap1vf/3D_509_36_36_300s_mf_2.out',/f77_unf
; goto, jump1
;endif


end


