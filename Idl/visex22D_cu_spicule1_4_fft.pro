cdtarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)





ipic=long(1)

n1=2048
n2=1024
nsect1=512
nsect2=512

pic=521
;pic=pic/2

;wdtv1=dblarr(pic+1,2)
;wdt=dblarr(pic+1,n2,4)
;wdt1=dblarr(pic+1,n2,4)

;wdtv=dblarr(pic+1,n1,4)
wdt1v=dblarr(pic+1,n1,4)
wdt2v=dblarr(pic+1,n1,4)

dt=0.001

for ipic=0,pic do begin
;while not(eof(1)) do begin

mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)

VV=dblarr(1)
VVa=dblarr(1)

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

; Open an MPEG sequence: 
;mpegID = MPEG_OPEN([900,900],FILENAME='myMovie.mpg') 

nn=0


nn_i=0

close,1
close,2


;openr,1,'/data/ap1vf/4mm4mm200Nx255Nz.ini',/f77_unf
;openr,1,'/data/ap1vf/1_8Mnz4Mnx882400.out',/f77_unf

;openr,1,'/data/ap1vf/1_8Mnz4Mnx1983400.ini',/f77_unf

;openr,1,'/data/ap1vf/4Mnzx1976400.ini',/f77_unf

;openr,1,'/data/cs1mkg/VAC_NN_tests/zeroBW_p13.out',/f77_unf
;openr,1,'/home/mikeg/proj/sac2.5d-cuda/zero1_BW_bin.ini',/f77_unf
;openr,1,'/home/mikeg/proj/sac2.5d-cuda/zero1_BW.ini',/f77_unf
;openr,1,'/home/mikeg/proj/sac2.5d-cuda/test_OT.out'
;directory='/home/mikeg/proj/sac2.5d-cuda/out_OT_withhyper/'
;directory='../out/'
;directory='/fastdata/cs1mkg/smaug/spicule1/'
;directory='/fastdata/cs1mkg/smaug/spicule_nob1/'
;directory='/fastdata/cs1mkg/smaug/spicule_nosource_nohydros/'
;directory='/fastdata/cs1mkg/smaug/spicule_nohydros/'
;directory='/fastdata/cs1mkg/smaug/spicule7_nob/'
;directory='/fastdata/cs1mkg/smaug/spicule7a_nob/'
;directory='/fastdata/cs1mkg/smaug/spicule6a_nob/'

;directory='/data/cs1mkg/smaug_spicule1/out/spicule5b4/'
;directory='/fastdata/cs1mkg/smaug/em4b4_bhor120/'
;directory='/fastdata/cs1mkg/smaug/spicule8b0_nob/'
directory='/fastdata/cs1mkg/smaug/spicule7b0_nob/'
;pic=999
name='zerospic1__'
;name='zerospic1_'





;ndim=2
;n1=800
;n2=6


;picid=ipic*5+4
;picid=ipic
picid=ipic*1000L
;picid=ipic*2000L
;picid=ipic*500L


;if(picid lt 149500)then begin
;	name='zerospic1_
;endif 
;if( picid ge 149500 && picid le 297000) then begin
;	;name='zerospic1_restart_'
;        name='zerospic1_'       
;endif

;if (picid gt 297000) then begin
;	name='zerospic1_'
;endif


outfile=directory+name+strtrim(string(picid),2)+'.out'
openr,1,outfile
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
;tarr=[tarr,time]
;ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

print,'varname ',varname


xout=dblarr(2)
yout=dblarr(2)


n1=nx(0)
n2=nx(1)
x=dblarr(n1,n2,ndim)
if (nn eq 0) then w=dblarr(n2,n1,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2)
;e2=dblarr(n1,n2)
readu,1,x
;readu,1,e2
;e2=rotate(e2,1)
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=rotate(wi,1)
endfor


for zslice=0,n1-1 do begin

yp=nsect1
wdt1v(ipic,zslice,0)=w(yp,zslice,0)
wdt1v(ipic,zslice,1)=w(yp,zslice,1)/(w(yp,zslice,0)+w(yp,zslice,7))
wdt1v(ipic,zslice,2)=w(yp,zslice,2)/(w(yp,zslice,0)+w(yp,zslice,7))
wdt1v(ipic,zslice,3)=w(yp,zslice,3)

yp=nsect1+256
wdt2v(ipic,zslice,0)=w(yp,zslice,0)
wdt2v(ipic,zslice,1)=w(yp,zslice,1)/(w(yp,zslice,0)+w(yp,zslice,7))
wdt2v(ipic,zslice,2)=w(yp,zslice,2)/(w(yp,zslice,0)+w(yp,zslice,7))
wdt2v(ipic,zslice,3)=w(yp,zslice,3)


endfor


mu=4.0*!PI/1.0e7

print,ipic,time


close,1
endfor

save,wdt1v,filename='fft/wdt1h_7b0_521s.dat'
save,wdt2v,filename='fft/wdt2h_7b0_521s.dat'

end
