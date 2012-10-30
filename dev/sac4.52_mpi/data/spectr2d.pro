fleft=30
fright=220

nj=512/2

dsl=dindgen(nj)-fleft
dsr=dindgen(nj)-fright

ww0=dindgen(2*nj)


dsl=atan(dsl)
dsr=-atan(dsr)

maxdsl=max(dsl)
mindsl=min(dsl)

maxdsr=max(dsr)
mindsr=min(dsr)

   dsl=(dsl-mindsl)/(maxdsl-mindsl)
   dsr=(dsr-mindsr)/(maxdsr-mindsr)

for i=0,2*nj-1 do begin
 if i lt nj then ww0(i)=dsl(i) else ww0(i)=dsr(i-nj) 
endfor  

www=dblarr(144,2*nj)

for i=0, 144-1 do begin
   www(i,*)=1.d0-ww0(*)
endfor

tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

tmas=fltarr(1)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, 0, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
WINDOW, 1, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

loadct, 4

mass=dblarr(1)
Vzmax=dblarr(1)

ia=1.0


dd=40

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

ii=1
close,1




aa=0

if aa eq 0 then begin

openr,1,'/data/ap1vf/3MHDCD4mpi16512_test.out',/f77_unf


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


n1=nx(0)
n2=nx(1)
n3=nx(2)


x=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)
wi=dblarr(n1,n2,n3)

readu,1,x


xout=dblarr(2)
yout=dblarr(2)
zout=dblarr(2)

for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,*,iw)=wi
endfor

endif

!p.multi =  0


;******************


MU=10.0



;set_plot, 'ps'


;device, filename='all_A02_gauss_collaps.ps', BITS=8, /color
;device, filename='all_A02_gauss_collaps.ps', BITS=8, /color

!p.thick = 1
!x.thick = 1
!y.thick = 1
!z.thick = 1
!p.font = 1

  loadct, 4

N_time=double(1)
Br=dblarr(2^MU)
Bi=dblarr(2^MU)
Ba=dblarr(2^MU)
FBr=dblarr(2^MU)
FBi=dblarr(2^MU)
FBa=dblarr(2^MU)
time=double(1)
NN=double(1)

N_t=139

N_end=128

A=dblarr(N_t+1,2^MU)
A_out=dblarr(N_t+1,N_end)
yy_out=dblarr(N_end)

a_time=dblarr(N_t+1)

yy=dindgen(2^MU)



ia=1
nn_i=0

print,'plot data'

!p.multi =  [0,1,1,0,1]

close,1
openr,1,'all_ln_5.dat',/f77_unf
readu,1, N_time
print, N_time
i=0


delta=4
ii=0

while (i lt N_t+1) do begin

readu,1, time
readu,1, Br
readu,1, Bi
readu,1, Ba
readu,1, FBr
readu,1, FBi
readu,1, FBa
readu,1, NN

ii=ii+1

if (ii eq delta) then begin
 print,'i=',i
 A(i,*)=Br
 a_time(i)=time 
 print,'time', time,i, n_elements(Br)
; read,aa
i=i+1
ii=0
endif

endwhile




print,'~~~~~',n_elements(a_time)
ij=i-1

delta_jj=8
i_start=0
jj=i_start
for i=0,N_end-1  do begin
  print,'i',i,jj
  a_out(*,i)=a(*,jj)  
  yy_out(i)=yy(jj)/2^MU
  jj=jj+delta_jj
  print,'~~~~~',n_elements(a_out(*,i))
  print,'~~~~~',n_elements(a_out(10,*))  
endfor 

  print,'~~~~~',n_elements(yy_out)







   xstart = (0.5/5.0)
   ystart = (0.4/3.0)
   xend = (3.5/5.0) + xstart
   yend = (2.2/3.0) + ystart   


;*****************










n2=n_elements(a_time)
n3=n_elements(yy_out)

xzstart=0
xzend=n_elements(yy_out)-1

xystart=0
xyend=n_elements(a_time)-1


w0=dblarr(n2,n3)
 
w0=www

window, 0,xsize=500,ysize=500

tvframe, w0, /bar

stop

wF=dblarr(n2,n3)
wSin=dblarr(n2,n3)


xyi=findgen(xyend+1)
xzi=findgen(xzend+1)


for i=xystart,xyend do begin
 for j=xzstart,xzend do begin
   wSin(i,j)=sin(xyi(i)*!Pi/(xyi(xyend)-xyi(xystart)))* $
                sin(xzi(j)*!Pi/(xzi(xzend)-xzi(xzstart)))
 endfor
endfor

window, 1,xsize=500,ysize=500
  
tvframe,A_out,/bar, /sample

stop

wF(*,*)=(w(kk,*,*,3)/w(kk,*,*,0)) ;*wSin(*,*)

wF=FFT(wF(*,*),-1)


stop

wF1=ABS(wF)^2.d0
wF1=reform(wF1)

ff=total(wF1,1)

window, 2,xsize=500,ysize=500

tvframe, alog10(wF1),/bar

;plot,ff
stop


wF3=reform(wF)
wF4=w0*wF3


wF4=double(FFT(wF4,1)) ;/(wSin+0.00001))

window, 3, xsize=500,ysize=500

tvframe, wF4,/bar,/sample
stop


wF3=reform(wF)

wF4=ABS(wF3*w0)^2.d0

tvframe, alog10(wF4), /bar

stop





end
