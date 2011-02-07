; profiles for cut-off period omega and rho 8.01.2009



function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
;for i=0,nel-2 do res(i)=(1.d0/(x(i+1)-x(i)))*(f(i+1)-f(i))
res(0)=res(2)
res(1)=res(2)
res(nel-1)=res(nel-3)
res(nel-2)=res(nel-3)
return,res
end

function deriv16,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=3,nel-4 do res(i)=(1.d0/60.D0/(x(i+1)-x(i)))*(-f(i-3)+9.d0*f(i-2)-45.d0*f(i-1)+45.d0*f(i+1)-9.d0*f(i+2)+f(i+3))

res(0)=res(3)
res(1)=res(3)
res(2)=res(3)
res(nel-1)=res(nel-4)
res(nel-2)=res(nel-4)
res(nel-3)=res(nel-4)
return,res
end

tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)


DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW

!p.multi = [0,2,0,0,1]


!p.thick = 1
!x.thick = 1
!y.thick = 1
!z.thick = 1
!p.font = 1.0


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


;openr,2,'/data/ap1vf/3D_256_80_80.ini',/f77_unf
openr,2,'/data/ap1vf/3D_509_36_36_30s.out',/f77_unf

readu,2,headline
readu,2,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,2,nx
eqpar=dblarr(neqpar)
readu,2,eqpar
readu,2,varname


xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

wi=dblarr(n1,n2,n3)

w=dblarr(n1,n2,n3,nw)

readu,2,x
for iw=0,nw-1 do begin
 print, iw 
 readu,2,wi
  w(*,*,*,iw)=wi
endfor

close,2

mu=4.0*!PI/1.0e7

print,time


R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

pp=8  ;x
kk=8  ;y

T=dblarr(n1)

P=dblarr(n1)

rhob=dblarr(n1)
rhob=reform(w[*,pp,kk,9])
Bx=reform(w[*,pp,kk,10])
By=reform(w[*,pp,kk,11])
Bz=reform(w[*,pp,kk,12])

P=reform(w[*,pp,kk,8])


z=dblarr(n1)
z=reform(x(*,pp,kk,0))

H_non_mag=dblarr(n1)
; V initial=0

P=P-(Bx^2.0+By^2.0+Bz^2.0)/2.d0
P=(gamma-1.d0)*P

beta=dblarr(n1)

beta=((Bx*sqrt(mu)*1.0e4)^2.0+(By*sqrt(mu)*1.0e4)^2.0+(Bz*sqrt(mu)*1.0e4)^2.0)/2.0/P

;plot, beta,title='1/beta',xtitle='x', ytitle='z',charsize=2.0 

T=mu_gas*P/R/rhob

Cs=sqrt(gamma*P/rhob)

gg=274.d0

H_non_mag=R*T/mu_gas/gg



close,1



; 
;openr,1,'/data/ap1vf/3D_256_80_80_test.ini',/f77_unf

openr,1,'/data/ap1vf/3D_509_36_36_30s.out',/f77_unf

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

print, '1'

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

mu=4.0*!PI/1.0e7

print,time


R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

pp=8  ;x
kk=8  ;y

T=dblarr(n1)

P=dblarr(n1)

rhob=dblarr(n1)
rhob=reform(w[*,pp,kk,9])
Bx=reform(w[*,pp,kk,10])
By=reform(w[*,pp,kk,11])
Bz=reform(w[*,pp,kk,12])

P=reform(w[*,pp,kk,8])


z=dblarr(n1)
z=reform(x(*,pp,kk,0))

H=dblarr(n1)
; V initial=0

P=P-(Bx^2.0+By^2.0+Bz^2.0)/2.d0
P=(gamma-1.d0)*P

T=mu_gas*P/R/rhob

Cs=sqrt(gamma*P/rhob)

gg=274.d0

H=R*T/mu_gas/gg

Pa=1.d0/(gamma*gg/4.d0/!Pi/Cs*sqrt(1.d0+2.d0*deriv1(H,z)))

label_Pa='!3P!Da !N [s]!N'

label_H='log H / Mm'

label_omega='!3!4q!3!N [!N mHz]'

label_Z='Height [Mm]'

scale=1.e6
n_end=244


plot, z/scale, alog10(H_non_mag),  YSTYLE=8, ytitle=label_H, charsize=1.5, $
      xtitle=label_Z, YRANGE=[alog10(min(H)), alog10(max(H))]


AXIS, 8.01d0, Pa(4), $
      YRANGE=[min(Pa(4:n_end)), max(Pa(4:n_end))], $
      YAXIS=1, /SAVE, charsize=1.5, ytitle=label_Pa
oplot, z(4:n_end)/scale, Pa(4:n_end)

omega=1.d0/Pa*1.e3

oplot, z(4:n_end)/scale, omega(4:n_end)

AXIS, 10.01d0, max(omega[4:n_end]), $
      YRANGE=[min(omega[4:n_end]), max(omega[4:n_end])], $
      YAXIS=1, /SAVE, charsize=1.5, ytitle=label_omega


oplot, z(4:n_end)/scale, omega(4:n_end)
      
print,max(1.d0/Pa(4:n_end))*1.e3, min(1.d0/Pa(4:n_end))*1.e3       


set_plot, 'ps'

!p.multi = [0,1,1,0,1]

device, filename='H_Pa1.eps', /encap
!p.thick = 4
!x.thick = 4
!y.thick = 4
!z.thick = 4
!p.font = 1.0

label_Pa='P!Dac!N [s]!N'

label_H='log !9L!3 / Mm'

label_omega='!9n!3!Dac!3!N [!N mHz]'

label_Z='!3Height [Mm]'

pos = Aspect(1.0)

plot, z/scale, alog10(H_non_mag),  YSTYLE=8, ytitle=label_H, charsize=1.5, $
      xtitle=label_Z, YRANGE=[alog10(min(H)), alog10(max(H))],Position=pos


AXIS, 8.0d0, Pa(4), $
      YRANGE=[min(Pa(4:n_end)), max(Pa(4:n_end))], $
      YAXIS=1, /SAVE, charsize=1.5, ytitle=label_Pa

oplot, z(4:n_end)/scale, Pa(4:n_end), linestyle=2

omega=1.d0/Pa*1.e3

AXIS, 10.2d0, max(omega[4:n_end]), $
      YRANGE=[min(omega[4:n_end]), max(omega[4:n_end])], $
      YAXIS=1, /SAVE, charsize=1.5, ytitle=label_omega


oplot, z(4:n_end)/scale, omega(4:n_end), linestyle=3


xyouts, 9600,8500, 'P!Dac', /device, charsize=1.5
xyouts, 5200,3800, '!9L!3', /device, charsize=1.5
xyouts, 4900,6200, '!9n!3!Dac', /device, charsize=1.5

device, /close
set_plot, 'x'
stop


end
