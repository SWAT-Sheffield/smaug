tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

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

kout=0

nn_i=0

close,1


nt=397

zstart=0
zend=140

xstart=0
xend=99

ystart=0
yend=99

nm=100.d0

dd=dblarr(nm,nm)

dd(*,*)=0.5

ns=intarr(1)
tst=double(1)
ted=double(1)

timearr_vr=dblarr(nt,xend)
timearr_vphi=dblarr(nt,xend)
timearr_Vvphi=dblarr(nm,nt)

asi_s=dblarr(nt)

beta=dblarr(nt,zend)
Cs=dblarr(nt,zend)

print, 'tuta'

openr,1,'/data/ap1vf/3D_tube_196_100_100_120s_full.out',/f77_unf

dr='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/Vphi_analiz_r02_kz50/'

!p.multi = [0,2,2,0,1]

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

tmin=0.d0
tmax=521.41

xx=dblarr(n2)
yy=dblarr(n3)
zz=dblarr(n1)


xx(*)=x(1,*,1,1)
yy(*)=x(1,1,*,2)
zz(*)=x(*,1,1,0)

scale=1.0d6

Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

Vap=dblarr(n1,n2,n3)
Cs=dblarr(n1,n2,n3)

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



vr=dblarr(n1,n2,n3)
vphi=dblarr(n1,n2,n3)

;br=dblarr(n2,n3)
;bphi=dblarr(n2,n3)


if kout eq 0 then begin

r=dblarr(n2,n3)
xxn=dblarr(n2)
yyn=dblarr(n3)

 for i=0,n2-1 do begin
  for j=0,n3-1 do begin
  
  xxn[i]=xx[i]-xx[n2/2]
  yyn[j]=yy[j]-yy[n3/2]
  r[i,j]=sqrt(xxn[i]^2.d0+yyn[j]^2.d0)
  endfor
 endfor  

kout=1
endif



;************************************************************************
kx=49
ky=49
kz=50
;************************************************************************



for k=0,n1-1 do begin
 for i=0,n2-1 do begin
  for j=0,n3-1 do begin

 vr[k,i,j]= avx[k,i,j]*xxn[i]/r[i,j]+avy[k,i,j]*yyn[j]/r[i,j]
 vphi[k,i,j]=-avx[k,i,j]*yyn[j]/r[i,j]+avy[k,i,j]*xxn[i]/r[i,j]

  endfor
 endfor
endfor ; k


rr=0.2*1.d6
phi=0.0d0

phi_i=dblarr(nm)
xc=dblarr(nm)
yc=dblarr(nm)

ii=0

deltaphi=2.d0*!Pi/(nm-1)
while phi le 2.d0*!Pi do begin

   xi=rr*cos(phi)
   yi=rr*sin(phi)
   
   
   xc[ii]=interpol(dindgen(xend-xstart+1),xx(xstart:xend)-xx(n2/2),xi)
   yc[ii]=interpol(dindgen(yend-ystart+1),yy(ystart:yend)-yy(n3/2),yi)
   
   timearr_Vvphi[ii,nn]=interpolate(Vphi, kz,xc[ii],yc[ii])
   
   phi_i[ii]=phi
   
   
   ii=ii+1
   phi=phi+deltaphi
   
endwhile
;--------------------
; Make a vector of 16 points, A[i] = 2pi/16: 
A = FINDGEN(17) *  (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
USERSYM, 0.6*COS(A), 0.6*SIN(A), /fill
;--------------------
print, nn


cs =1.0
tvframe,vphi(kz,*,*),/sample, /bar, $
       ; xrange=[min(xx)/scale, max(xx)/scale], yrange=[min(yy)/scale, max(yy)/scale], $    
        xtitle='x [px]', ytitle='y [px]',CT='dPdT', $
	BTITLE='V!D phi!N!3 [m/s]', charsize=cs

for i=0, ii-1 do begin ; kolichestvo tochek
   plots, xc[i],yc[i], color=200, psym=8
endfor

tvframe,timearr_Vvphi,/sample, /bar,  $ 
        xrange=[phi_i[0],phi_i[nm-2]],$  
        xtitle='Phi ',CT='dPdT', $
	TITLE='V!D!7 u!3!N [m/s]', charsize=cs
	
	
plot,timearr_Vvphi(*,nn), charsize=cs

asi=timearr_Vvphi(*,nn)

asi_s[nn]=total(asi)/nm
print,'ss=',asi_s[nn], '   time=', time
plot, asi_s, psym=4, charsize=cs




indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
write_png,dr+indexss+'.png',image_p, red,green, blue

nn=nn+1

endwhile
close, 1

end


