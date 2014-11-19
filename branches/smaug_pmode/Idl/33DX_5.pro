
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)



;openr,1,'/data/ap1vf/3D_509_36_36_300s_mf_1.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60.out',/f77_unf
;openr,1,'/data/ap1vf/3D_196_100_100.ini',/f77_unf

;openr,1,'/home/mikeg/proj/sac2.5d-cuda/test_OT.out'
;directory='/home/mikeg/proj/sac2.5d-cuda/out_OT_withhyper/'
;directory='../out/'
;directory='/nobackup/shemkg/modes3/f1/'
;directory='/fastdata/cs1mkg/smaug/spicule5b0_3d/'
directory='/fastdata/cs1mkg/smaug/spic4b1_1_3d/'
;directory='/data/cs1mkg/smaug_spicule1/spicule5b0_3d/'
;directory='/fastdata/cs1mkg/smaug/spicule6b0_3d/'
;directory='/nobackup/shemkg/modes/3demt1/'
;directory='/nobackup/shemkg/modes/'
;pic=999
;name='zeroOT_'
;name='3D_em_t1_bin_'
;name='3D_em_t1_bin_np010808_00'
;name='3D_em_f1_bin_'
name='zerospic1__'
;ndim=2
;n1=800
;n2=6
ipic=0L
pic=489L
pic=280L
pic=1043L



nn=0

;pic=394L
for ipic=1L,pic do begin
;for ipic=280L,699L do begin


headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)

close,1


;while not(eof(1)) do begin

picid=ipic*1000L
;picid=ipic*5+4
;picid=ipic;
outfile=directory+name+strtrim(string(picid),2)+'.out'
print,'ipic=',ipic
;openr,1,outfile,/f77_unf
openr,1,outfile




;while not(eof(1)) do begin
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
;neqpar=4
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


n1=nx(0)
n2=nx(1)
n3=nx(2)

;n1=128
;n2=128
;n3=128
;ndim=3
;nw=13


;ndim=3
;n1=256
;n2=512
;n3=512
;nw=13


x=dblarr(n1,n2,n3,ndim)
wi=dblarr(n1,n2,n3)
w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
;for iw=0,nw-3 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

indexs=strtrim(picid,2)

a = strlen(indexs)                                                  
case a of
 1:indexss='000000000'+indexs
 2:indexss='00000000'+indexs
 3:indexss='0000000'+indexs
 4:indexss='000000'+indexs
 5:indexss='00000'+indexs                                                           
 6:indexss='0000'+indexs                                             
 7:indexss='000'+indexs                                              
 8:indexss='00'+indexs                                               
 9:indexss='0'+indexs                                               
endcase 

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

T=dblarr(n1,n2,n3)
BT=dblarr(n1,n2,n3)
DT=dblarr(n1,n2,n3)
T[*,*,*]=(w[*,*,*,4]+w[*,*,*,8])
T[*,*,*]=T[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/(w[*,*,*,0]+w[*,*,*,9])/2.0
T[*,*,*]=T[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0+(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.d0

BT[*,*,*]=(w[*,*,*,8])
BT[*,*,*]=BT[*,*,*]-((w[*,*,*,10])^2.0+(w[*,*,*,11])^2.0+(w[*,*,*,12])^2.0)/2.d0


T[*,*,*]=(gamma-1.d0)*T[*,*,*]
BT[*,*,*]=(gamma-1.d0)*BT[*,*,*]
DT[*,*,*]=T[*,*,*]-BT[*,*,*]

C=dblarr(n1,n2,n3)
C[*,*,*]=sqrt(gamma*T[*,*,*]/(reform(w[*,*,*,0]+w[*,*,*,9])))
	 
;data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9))/C(*,*,*))	 
data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9)))
;data=reform(DT(*,*,*))


;tvframe, data(*,*,128)
;stop
time=ipic*1000*0.001
 hData = PTR_NEW(data, /NO_COPY) 
 a='Mx/'+indexss                                          
 Slicer_X, hdata, time


;for i=1,2 do begin
; data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9)))
; data=reform(w(*,*,*,i))
; hData = PTR_NEW(data, /NO_COPY) 

;st='/data/ap1vf/png/3D/30s/NoMagField/509_36_36/'

;case i of                                                           
; 1:a=st+'Vz/'+indexss+'.png'                                         
; 2:a=st+'Vx/'+indexss+'.png'                                             
; 3:a=st+'Vy/'+indexss+'.png'                                              ;

;endcase

; Slicer_ps, hdata, time, a
;Slicer_X, hdata, time, a
 
; Slicer3, hdata, DATA_NAMES='Dave'
 
; stop
;endfor
;window,2
image_p = TVRD_24()
write_png,'/data/cs1mkg/smaug_spicule1/Idl/images/spic4b1_1_3d/im'+indexss+'.png',image_p, red,green, blue


nn=nn+1
print, 'time - - - ', time

;end loop over each input file
;endwhile
endfor 
end
