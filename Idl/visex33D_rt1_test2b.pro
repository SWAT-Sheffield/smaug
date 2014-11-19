function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
;for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2)
res(1)=res(2)
res(nel-1)=res(nel-3)
res(nel-2)=res(nel-3)
return,res
end


tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


ii=1

if (ii eq 1) then begin
loadct,4
mixct
endif else begin
loadct,0
tek_color
endelse

window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
window, 1,xsize=800,ysize=800,XPOS = 500, YPOS = 80

;for ipic=pic,489000 do begin
;for ipic=280,489 do begin
;while not(eof(1)) do begin
;for ipic=1,371 do begin
for ipic=pic,pic do begin
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

; Open an MPEG sequence: 
;mpegID = MPEG_OPEN([700,1200],FILENAME='myMovie.mpg') 




nn=0
np=0
kkk=4

nn_i=0

close,1
close,2




;openr,1,'/data/ap1vf/3D_509_36_36_300s.out',/f77_unf
;openr,1,'/data/ap1vf/3D_tube_modif_200_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_vert_driver/3D_tube_196_100_100_phi.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_vertical/3D_tube_200_100_100_vertical.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_multidriver_lower.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100t.out',/f77_unf



;openr,1,'/data/ap1vf/vxx.040000',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60t.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_196_100_100_test.out',/f77_unf

;openr,1,'/data/ap1vf/3D_vert_driver/3D_tube_196_100_100_vert.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_modif_200_100_100.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_200s_puls.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf

;openr,1,'/home/mikeg/proj/sac2.5d-cuda/test_OT.out'
;directory='/home/mikeg/proj/sac2.5d-cuda/out_OT_withhyper/'
;directory='../out/'
;directory='/nobackup/shemkg/modes3/sp2rc/3demt1/'
;directory='/nobackup/shemkg/modes3/f1/'
;directory='/fastdata/cs1mkg/smaug/spicule4b0_3d/'
directory='/data/cs1mkg/smaug_spicule1/spicule4b0_b1_3d/'
;directory='/data/cs1mkg/smaug_spicule1/spicule5b0_3d/'
;directory='/data/cs1mkg/smaug_spicule1/spicule4b0_b1_3d/'

;directory='/nobackup/shemkg/modes/'
;pic=999
;name='zeroOT_'
;name='3D_em_t1_bin_np010808_00'
;name='3D_em_f1_bin_'
name='zerospic1__'
;ndim=2
;n1=800
;n2=6


;picid=ipic*5+4
;picid=ipic*1000L
picid=ipic;
outfile=directory+name+strtrim(string(picid),2)+'.out'
print,'ipic=',ipic
;openr,1,outfile,/f77_unf
openr,1,outfile

readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


print, 'tuta1'
xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)

;ndim=3
;n1=256
;n2=512
;n3=512
;nw=13

x=dblarr(n1,n2,n3,ndim)

print,'about to allocate'
print,n1,n2,n3,nw

wi=dblarr(n1,n2,n3)

w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

xx=dblarr(n2)
yy=dblarr(n3)
zz=dblarr(n1)


xx(*)=x(1,*,1,1)
yy(*)=x(1,1,*,2)
zz(*)=x(*,1,1,0)

Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)


mu=4.0*!PI/1.0e7

print,'******************* time = ', time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

;scale=1.d6
scale=1.0;

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=127
ystart=0
yend=127

pp=63 ;x
kk=5  ;y

wset,0
!p.multi = [0,4,4,0,1]


zstart=0
zend=127

wt=dblarr(zend-zstart+1,xend-xstart+1,iw)
wt=reform(w(zstart:zend,xstart:xend,pp,*))

wy=dblarr(n1,n3,iw)
wy=reform(w(zstart:zend,pp,*,*))

wt(*,*,3)=reform(w(zstart:zend,xstart:xend,pp,3))

wt(*,*,12)=reform(w(zstart:zend,pp,ystart:yend,12))

;wv=dblarr(n1,n3);
;wv=wt(*,*,2);
;wv(where(abs(wv) le 1.0d-3))=0.0
;wt(*,*,2)=wv;
;wv=wt(*,*,3);
;wv(where(abs(wv) le 1.0d-3))=0.0
;wt(*,*,3)=wv;


saeb=dblarr(zend-zstart+1,xend-xstart+1)
sarho_t=dblarr(zend-zstart+1,xend-xstart+1)
sabz_t=dblarr(zend-zstart+1,xend-xstart+1)
sabx_t=dblarr(zend-zstart+1,xend-xstart+1)
saby_t=dblarr(zend-zstart+1,xend-xstart+1)

saeb(*,*)=wt(*,*,8)
sarho_t(*,*)=wt(*,*,0)+wt(*,*,9)
sabz_t(*,*)=wt(*,*,10)
sabx_t(*,*)=wt(*,*,11)
saby_t(*,*)=wt(*,*,12)

vt=dblarr(n1,n2,n3)
vvt=dblarr(n2,n3)
vt(*,*,*)=sqrt(w(*,*,*,1)^2.d0+w(*,*,*,2)^2.d0+w(*,*,*,3)^2.d0)/(w(*,*,*,0)+w(*,*,*,9))


;****************** Pressure background begin ********************
TP=saeb
TP=TP-(sabx_t^2.0+saby_t^2.0+sabz_t^2.0)/2.0
TP=(gamma-1.d0)*TP
;****************** Pressure background end ********************

if (ii eq 1) then begin


tvframe,rotate(wt(*,*,0),1), /bar,title='rho',$ 
        /sample, xtitle='x', ytitle='y',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]		
	
st=strTrim(it,1)	

ww=dblarr(xend-xstart+1,zend-zstart+1)
xxs=xx[xstart]
xxe=xx[xend]
zzs=zz[zstart]
zze=zz[zend]
	
;close,10
;openw,10,'/data/ap1vf/data_line_prof/rho/rho_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate((wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10
	

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz',$
        xtitle='x', ytitle='y',charsize=2.0; title='dPdT'


;close,10
;openw,10,'/data/ap1vf/data_line_prof/vz/vz_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10

tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='x', ytitle='z',charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]; CT='dPdT'

;close,10
;openw,10,'/data/ap1vf/data_line_prof/vx/vx_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10	

tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'
	
;close,10
;openw,10,'/data/ap1vf/data_line_prof/vy/vy_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10	
	

tvframe,rotate(wy(*,*,2)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'
	
tvframe,rotate(wt(*,*,4),1),/bar, /sample, title='e', xtitle='x', ytitle='z', $
        charsize=2.0

tvframe,rotate(wt(*,*,5),1)*sqrt(mu)*1.0e4,/bar,/sample, title='bz', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,6),1)*sqrt(mu)*1.0e4,/bar,/sample, title='bx', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,7),1)*sqrt(mu)*1.0e4,/bar,/sample, title='by', $
        xtitle='x', ytitle='z', charsize=2.0




;tvframe,rotate(wt(*,*,11),1)*sqrt(mu)*1.0e4,/bar,/sample, title='Bx_b', $
;        xtitle='x', ytitle='z', charsize=2.0

;tvframe,rotate(wt(*,*,12),1)*sqrt(mu)*1.0e4,/bar,/sample, title='By_b', $
 ;       xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,8),1),/bar,/sample, title='eb', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,9),1),/bar,/sample, title='rho_b', $
        xtitle='x', ytitle='z', charsize=2.0

tvframe,rotate(wt(*,*,10),1)*sqrt(mu)*1.0e4,/bar,/sample, title='Bz_b', $
        xtitle='x', ytitle='z', charsize=2.0, $
	xrange=[xx[xstart]/scale, xx[xend]/scale], $
	yrange=[zz[zstart]/scale, zz[zend]/scale]	

;stop
T=mu_gas*TP/R/sarho_t
;goto, jump20

tvframe,rotate(T(*,*),1),/bar,/sample, title='T', $
        xtitle='x', ytitle='z', charsize=2.0
	
save, filename='Temp_VALIIIC_old.sav', T
save, filename='z_VALIIIC_old.sav', zz
	
;close,10
;openw,10,'/data/ap1vf/data_line_prof/T/T_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(T(*,*),1)
;writeu,10,ww
;close, 10	

;tvframe,rotate(wt(*,*,12),1)*sqrt(mu)*1.0e4,/bar,/sample, title='By_b', $
;        xtitle='y', ytitle='z', charsize=2.0

Va=dblarr(zend-zstart+1,xend-xstart+1)
Vap=dblarr(zend-zstart+1,xend-xstart+1)

Va(*,*)=sqrt((wt(*,*,10)^2.d0+wt(*,*,11)^2.d0+wt(*,*,12)^2.d0)/wt(*,*,9))*sqrt(mu)*1.0e4

Vap(*,*)=sqrt((wt(*,*,5)^2.d0+wt(*,*,6)^2.d0+wt(*,*,7)^2.d0)/(wt(*,*,9)+wt(*,*,0)))*sqrt(mu)*1.0e4

tvframe, rotate(Va(*,*),1)/1000.d0, title='V!DA!N!3 [km/s]',/bar,/sample, $
        xtitle='y', ytitle='z', charsize=2.0


plot, rotate(Va(10,*),1)/1000.d0  
oplot, rotate(Va(10,*),1)/1000.d0 , psym=4 
oplot, rotate(Vap(10,*),1)/1000.d0
oplot, rotate(Vap(10,*),1)/1000.d0, psym=5


Cs=dblarr(zend-zstart+1,xend-xstart+1)

Cs(*,*)=sqrt(gamma*TP(*,*)/wt(*,*,9))

tvframe, rotate(Cs(*,*)/1000.d0,1),  title='C!DS!N!3 [km/s]',/bar, charsize=2.0,/sample

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200	


 
indexs=strtrim(np,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

image_p = TVRD_24()
;write_png,'/data/ap1vf/png/3D/tube/test_200_puls/all/'+indexss+'.png',image_p, red,green, blue


np=np+1

;goto, jump10
wset,1
!p.multi = [0,3,2,0,1]

tvframe,rotate(wy(*,*,1)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

tvframe,rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vx', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'

tvframe,rotate(wy(*,*,2)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

tvframe,rotate(wt(*,*,2)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vy', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'


tvframe,rotate(wy(*,*,3)/(wy(*,*,0)+wy(*,*,9)),1),/sample, /bar,title='Vz', $
        xtitle='y', ytitle='z',charsize=2.0;, title='dPdT'

tvframe,rotate(wt(*,*,3)/(wt(*,*,0)+wt(*,*,9)),1),/sample, /bar,title='Vz', $
        xtitle='x', ytitle='z',charsize=2.0; CT='dPdT'





 
;stop

jump22:

	
;wait, 0.5
;endfor

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200


endif else begin

endelse


T=dblarr(n1,n2,n3)


T[*,*,*]=(w[*,*,*,4]+w[*,*,*,8])


T[*,*,*]=T[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/(w[*,*,*,0]+w[*,*,*,9])/2.0

T[*,*,*]=T[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0+(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.d0


beta=dblarr(n1,n2,n3)

beta[*,*,*]=(((w[*,*,*,5]+w[*,*,*,10])*sqrt(mu)*1.0e4)^2.0+((w[*,*,*,6]+w[*,*,*,11])*sqrt(mu)*1.0e4)^2.0+$
              ((w[*,*,*,7]+w[*,*,*,12])*sqrt(mu)*1.0e4)^2.0)/2.0/((gamma-1.d0)*T[*,*,*])

;tvframe,rotate(beta(*,*,pp),1),/bar,/sample, title='1/beta',  xtitle='x', ytitle='z', charsize=2.0



;plot, alog10(1.d0/beta[*,kk,pp]),title='1/beta',xtitle='x', ytitle='z',charsize=2.0 

T[*,*,*]=(gamma-1.d0)*T[*,*,*]


C=dblarr(n1,n2,n3)

C[*,*,*]=sqrt(gamma*T[*,*,*]/(reform(w[*,*,*,0]+w[*,*,*,9])))

;tvframe, rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9))/C(*,*,pp),1),title='V1/C',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample 
;tvframe, rotate(T(*,*,pp),1),title='P',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample 

;tvframe, rotate(alog10(T(*,*,pp)),1),title='p',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample 

;plot, alog10(T(*,kk,pp)),title='log(p)', xtitle='x', ytitle='z',charsize=2.0 

;T[*,*,*]=mu_gas*T[*,*,*]/R/(w[*,*,*,0]+w[*,*,*,9])


;tvframe, rotate(T(*,*,pp),1),title='T',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample 



;ttt=rotate(T(*,*,pp),1)

;plot,alog10(ttt(60,*))

;tvframe, rotate(C(*,*,pp),1),title='C',/bar, xtitle='x', ytitle='z',charsize=2.0,/sample 

;jump1 :


indexs=strtrim(picid,2)
;indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='00000000'+indexs                                             
 2:indexss='0000000'+indexs                                              
 3:indexss='000000'+indexs                                               
 4:indexss='00000'+indexs
 5:indexss='0000'+indexs                                             
 6:indexss='000'+indexs                                              
 7:indexss='00'+indexs                                               
 8:indexss='0'+indexs                                               
endcase   

;image_p = TVRD_24()
;write_png,'/data/cs1mkg/smaug_spicule1/Idl/images/spicule4b0_3d/s4b0_3d'+indexss+'.png',image_p, red,green, blue


nn=nn+1

;goto, jump

;wset,1
;!p.multi = [0,2,1,0,1]

;plot,w(*,kk,pp,1)/(w(*,kk,pp,0)+w(*,kk,pp,9)),title='v1', xtitle='x', ytitle='y',charsize=1.0  , /ys
;oplot,w(*,kk,pp,1)/(w(*,kk,pp,0)+w(*,kk,pp,9)),psym=4, color=100


;tek_color
;plot, alog10(T(*,kk,pp)),title='log(T)', yrange=[3.0,6.0], xtitle='x', ytitle='z',charsize=1.0 
;plot,w(*,kk,pp,5)+w(*,kk,pp,10),title='Bz_T', xtitle='x', ytitle='z', charsize=1.0


;plot,wt(*,kk,1)/(wt(*,kk,0)+wt(*,kk,9))/C(*,kk,pp),title='V1/C', xtitle='x', ytitle='z', charsize=1.0

;tek_color
;oplot,wt(*,kk,1)/(wt(*,kk,0)+wt(*,kk,9))/C(*,kk,pp), psym=4
 
;print,  kkk
;if (kkk eq 4) then begin


;kkk=1

;endif

;kkk=kkk+1
;if (ii eq 2) then read,a

;maxa=[maxa,max(w(*,*,2)/(w(*,*,0)+w(*,*,7)))]
maxa=[maxa,max(Vt)]

;data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9)))
;hData = PTR_NEW(data, /NO_COPY) 
;Slicer_ps, hdata, time

omega=dblarr(n1)

cutoff=dblarr(n1)

gg=274.d0
omega(*)=sqrt(-gg/w(*,10,10,9)*deriv(x(*,10,10,0),w(*,10,10,9)))

cutoff(*)=gamma*gg/2.d0/C(*, 10,10) ;*sqrt(1+2.d0*lambda)

goto, jump

set_plot, 'ps'

!p.multi = [0,1,1,0,1]

device, filename='omega.ps' ;, /encap
!p.thick = 4
!x.thick = 4
!y.thick = 4
!z.thick = 4
!p.font = 1.0

ss='!9w!3'




plot, x(*,10,10,0)/1.d6, alog10(w(*,10,10,9)), xtitle='!3z (Mm)', ytitle=ss,$
     charsize=1.5,  YSTYLE=8 ;/ylog,
;plot, x(*,10,10,0)/1.d6, cutoff, xtitle='!3z (Mm)', ytitle=ss,$
;     charsize=1.5 ;, /ylog, YSTYLE=8
;plot, x(*,10,10,0)/1.d6, omega, xtitle='!3z (Mm)', ytitle=ss,$
;     charsize=1.5 ;, /ylog, YSTYLE=8
;oplot, x(*,10,10,0)/1.d6, cutoff     
     
AXIS, YAXIS=1, YRANGE=[min(T(*,10,10)), max(T(*,10,10))], $
      /SAVE, charsize=1.5, ytitle='log(T)'

oplot, x(*,10,10,0)/1.d6, alog10(T(*,10,10))

xyouts, 5000,8500, ss, /device, charsize=1.5
xyouts, 7500,7000, 'log(T)', /device, charsize=1.5

device, /close
set_plot, 'x'

jump10 :

jump :

goto, jump20
;****************** begin vorticity *****************
wset,1
!p.multi = [0,1,1,0,1]
hh=20
vort_z= dblarr(n2,n3)
vvx=dblarr(n2,n3)
vvy=dblarr(n2,n3)
vvz=dblarr(n2,n3)

vvx(*,*)=w(hh,*,*,1)/(w(hh,*,*,0)+w(hh,*,*,9))
vvy(*,*)=w(hh,*,*,2)/(w(hh,*,*,0)+w(hh,*,*,9))
vvz(*,*)=w(hh,*,*,3)/(w(hh,*,*,0)+w(hh,*,*,9))

dvvxdy=dblarr(n2,n3)
dvvydx=dblarr(n2,n3)


print,'dx'

for i=0,n3-1 do begin
 dvvydx(*,i)=deriv1(vvy(*,i),xx)
endfor

for i=0,n2-1 do begin
 dvvxdy(i,*)=deriv1(vvx(i,*),yy)
endfor


vort_z=dvvydx-dvvxdy

tvframe, rotate(vort_z,1), /bar, /sample

;****************** end  vorticity *****************

jump20:
;stop
;endwhile
endfor


end
