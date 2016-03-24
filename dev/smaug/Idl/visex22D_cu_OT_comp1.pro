tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)




DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1050,ysize=350,XPOS = 1000, YPOS = 300 ; ZOOM
ii=1
!p.multi = [0,3,1,0,1]



if (ii eq 1) then begin
loadct,4
endif else begin
loadct,0
tek_color
endelse


npic=1


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


nn=0


nn_i=0

close,1
close,2


directory='../out/'

;directory='/fastdata/cs1mkg/smaug/ot_512/'
;pic=999
name='zeroOT_'
;ndim=2
;n1=800
;n2=6

;while not(eof(1)) do begin

;picid=ipic*5+4
;picid=ipic*10
picid=ipic
outfile=directory+name+strtrim(string(picid),2)+'.out'

;***********************************************
;  modify these lines to set input file
;  sac index is 1 more than smaug
;**************************************************
openr,1,'/data/cs1mkg/smaugtest/out/zeroOT_2.out'  ; set the index in here for the smaug output
openr,2,'/data/cs1mkg/sac_working/results/zeroOT.out',/f77_unf
pic=3;  set picture number for sac output
;**************************************************

;openr,1,outfile
readu,1,headline
readu,1,it,time,ndim,neqpar,nw


gencoord=(ndim lt 0)
tarr=[tarr,time]
;ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

print,'varname ',varname


xout=dblarr(2)
yout=dblarr(2)
;mxt=dblarr(npic)

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





indexs=strtrim(ipic,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

;stop
;endwhile
close,1


;pic=5000L
for ipic=1,pic do begin
;for ipic=0L,1000L do begin



readu,2,headline
readu,2,it,time,ndim,neqpar,nw

readu,2,nx
readu,2,eqpar
readu,2,varname

xout=dblarr(2)
yout=dblarr(2)


n1=nx(0)
n2=nx(1)
sx=dblarr(n1,n2,ndim)
if (nn eq 0) then sw=dblarr(n2,n1,nw)   ;was n1,n2,nw
swi=dblarr(n1,n2)
;e2=dblarr(n1,n2)
readu,2,sx
;readu,2,e2
;e2=rotate(e2,1)
for iw=0,nw-1 do begin
 readu,2,swi
 sw(*,*,iw)=rotate(swi,1)
endfor










endfor  ; end of looping over sac file images

;without the ghost cells
tvframe,w(2:253,2:253,7)+w(2:253,2:253,0), /bar,title='rho',/sample, xtitle='x', ytitle='y',charsize=2.0
tvframe,sw(2:253,2:253,7)+sw(2:253,2:253,0), /bar,title='srho',/sample, xtitle='x', ytitle='y',charsize=2.0
tvframe,(w(2:253,2:253,7)+w(2:253,2:253,0))-(sw(2:253,2:253,7)+sw(2:253,2:253,0)), /bar,title='difrho',/sample, xtitle='x', ytitle='y',charsize=2.0

;with the  ghost cells
;tvframe,w(*,*,7)+w(*,*,0), /bar,title='rho',/sample, xtitle='x', ytitle='y',charsize=2.0
;tvframe,sw(*,*,7)+sw(*,*,0), /bar,title='srho',/sample, xtitle='x', ytitle='y',charsize=2.0
;tvframe,(w(*,*,7)+w(*,*,0))-(sw(*,*,7)+sw(*,*,0)), /bar,title='difrho',/sample, xtitle='x', ytitle='y',charsize=2.0

image_p = TVRD_24()
write_png,'sac_smaug_OT_denscomp.png',image_p, red,green, blue




end
