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

ns=intarr(1)
tst=double(1)
ted=double(1)

timearr_br=dblarr(nt,zend)
timearr_bphi=dblarr(nt,zend)

beta=dblarr(nt,zend)
Cs=dblarr(nt,zend)

print, 'tuta'

openr,1,'/data/ap1vf/3D_tube_196_100_100_120s_full.out',/f77_unf

!p.multi = [0,2,0,0,1]

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



;vr=dblarr(n1,n2,n3)
;vphi=dblarr(n1,n2,n3)

br=dblarr(n1,n2,n3)
bphi=dblarr(n1,n2,n3)


if kout eq 0 then begin

r=dblarr(n1,n2,n3)
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



for k=0,n1-1 do begin
 for i=0,n2-1 do begin
  for j=0,n3-1 do begin

; vr[k,i,j]= avx[k,i,j]*xxn[i]/r[i,j]+avy[k,i,j]*yyn[j]/r[i,j]
  
; vphi[k,i,j]=-avx[k,i,j]*yyn[j]/r[i,j]+avy[k,i,j]*xxn[i]/r[i,j]


 br[k,i,j]= abx_t[k,i,j]*xxn[i]/r[i,j]+aby_t[k,i,j]*yyn[j]/r[i,j]
  
 bphi[k,i,j]=-abx_t[k,i,j]*yyn[j]/r[i,j]+aby_t[k,i,j]*xxn[i]/r[i,j]


; br[i,j]= wv(i,j,6)*(xx[i]-xx[n2/2])/r[i,j]+$
;          wv(i,j,7)*(yy[j]-yy[n3/2])/r[i,j]
  
; bphi[i,j]=-wv(i,j,6)*(yy[j]-yy[n3/2])/r[i,j]+$
;          wv(i,j,7)*(xx[i]-xx[n2/2])/r[i,j]
; br[n2/2,n3/2]=0.d0	  
; bphi[n2/2,n3/2]=0.d0
  endfor
 
 endfor

endfor ; k

;***************************************************************************
;*************** sound speed************************************************
;****************** Pressure background begin ********************
Rr=8.3e+003
mu=1.257E-6
mu_gas=1.2d0
gamma=1.66667

TP=aeb
TP=TP-(abx_t^2.0+aby_t^2.0+abz_t^2.0)/2.0
TP=(gamma-1.d0)*TP
;****************** Pressure background end ********************

;***************************************************************************
;***************alfven speed************************************************
Vap=sqrt((abx_t^2.d0+aby_t^2.d0+abz_t^2.d0)/arho_t)*1.0e4/1.0e3


kx=49
ky=49

timearr_br[nn,zstart:zend-1]=br[zstart:zend-1,kx,ky]
timearr_bphi[nn,zstart:zend-1]=bphi[zstart:zend-1,kx,ky]

print, nn

; ****************** ps, eps begin ****************************************
   xs=26.d0
   k=10.d0
   
   SET_PLOT,'ps'  
dr='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/ps/time_distance/'

   device, filename=dr+'P120_Vx_kx49_ky49_br_bphi_time.eps', $
   BITS=8, /color, xsize=xs, ysize=k, /encap
   
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2

tvframe,timearr_br(*,*),/sample, /bar, $
        yrange=[zz[zstart]/scale, zz[zend]/scale], xrange=[0, 521], $  
        xtitle='Time [s]', ytitle='z [Mm]',charsize=1.0, CT='dPdT', $
	BTITLE='B!D r!N!3 [m/s]' 

tvframe,timearr_bphi(*,*),/sample, /bar, $
        yrange=[zz[zstart]/scale, zz[zend]/scale], xrange=[0, 521], $  
        xtitle='Time [s]', ytitle='z [Mm]',charsize=1.0, CT='dPdT', $
	BTITLE='B!D!7 u!3!N [m/s]'


;plot, zz(zstart:zend),Vap(zstart:zend,kx,ky) 
	
;tvframe,vr(12,*,*),/sample, /bar, $
;        xrange=[min(xx)/scale, max(xx)/scale], yrange=[min(yy)/scale, max(yy)/scale], $  
;        xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0, CT='dPdT', $
;	BTITLE='V!D r!N!3 [m/s]'

;tvframe,vphi(12,*,*),/sample, /bar, $
;        xrange=[min(xx)/scale, max(xx)/scale], yrange=[min(yy)/scale, max(yy)/scale], $    
;        xtitle='x [Mm]', ytitle='z [Mm]',charsize=1.0, CT='dPdT', $
;	BTITLE='V!D phi!N!3 [m/s]'

;zz(*)=x(zstart:zend,1,0)
;for i=0,nt-1 do xx[i]=tmin+(tmax-tmin)/(nt*1.d0)*(i*1.d0)
;mixct10

;contour, beta(0:nt-1,*),xx(0:673),zz/scale, LEVELS = [0.01,0.1,1.0, 10.0], $
;         C_Annotation = ['100.0','10.0','1.0','0.1'],  $
;	 C_Thick = [2.0,2.0,4.0,2.0],  $
;	 /overplot, /FOLLOW, thick=4.0, color=128	  


; ss='time ='+strTrim(string(time),1)
 
; xyouts,20,2, ss, /device, color=200

nn=nn+1

device, /close
set_plot, 'x'

endwhile
close, 1

end


