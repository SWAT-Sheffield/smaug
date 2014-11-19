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


;DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
;WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


;window, 0,xsize=1200,ysize=800,XPOS = 700, YPOS = 900 
!p.multi = [0,3,2,0,1]

;*************** read old ini file ***************
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
;openr,1,'/fastdata/cs1mkg/smaug/spicule4_nob/zerospic1_290000.out'
;openr,1,'/fastdata/cs1mkg/smaug/spicule4a_nob/zerospic1_149000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/out/spicule5b4/zerospic1__555000.out'
;openr,1,'/fastdata/cs1mkg/smaug/spicule4_nob/zerospic1_100000.out'
;openr,1,'/fastdata/cs1mkg/smaug/em6b4_bhor120/zerospic1__84000.out'
;openr,1,'/fastdata/cs1mkg/smaug/em5b4_bhor120/zerospic1__86000.out'
;openr,1,'/fastdata/cs1mkg/smaug/spicule8b0_nob/zerospic1__584000.out'
;openr,1,'/fastdata/cs1mkg/smaug/spicule8b0_nob/zerospic1__643000.out'
;openr,1,'/fastdata/cs1mkg/smaug/spicule7b0_nob/zerospic1__491000.out'
;openr,1,'/fastdata/cs1mkg/smaug/scnstrh5b0b_nob/zerospic1__391000.out'

;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bin.ini',/f77_unf
;openr,1,'/fastdata/cs1mkg/smaug/spicule5b0_3d/zerospic1__279000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert_bin.ini',/f77_unf
;openr,1,'/data/cs1mkg/smaug_spicule1/spicule4b0_b1_3d/zerospic1__39000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert2_bin.ini',/f77_unf
;openr,1,'/data/cs1mkg/smaug_spicule1/spicule4b0_b2_3d/zerospic1__12000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/spicule5b0_3d/zerospic1__612000.out'
;openr,1,'/fastdata/cs1mkg/smaug/spic5b0_3d/zerospic1__79000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_btube1_bin.ini',/f77_unf;
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_btube2_bin.ini',/f77_unf;
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_btube3_bin.ini',/f77_unf;
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert10G_bin.ini',/f77_unf;
;openr,1,'/data/cs1mkg/smaug_spicule1/spicule4b0_bom500_3d/zerospic1__434000.out'
;openr,1,'/data/cs1mkg/vac3d/3D_196_100_100.ini',/f77_unf
;openr,1,'/fastdata/cs1mkg/smaug/spic4b0_3d/zerospic1__1077000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/spicule2p3a_0_3_3d/zerospic1__592000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/spicule4p3a_0_1_3d/zerospic1__1206000.out'
openr,1,'/data/cs1mkg/smaug_spicule1/spicule6p7a_0_0_3d/zerospic1__1400000.out'
;openr,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert10G_bin.ini',/f77_unf
;openr,1,'/data/cs1mkg/smaugtemp/smaug/out/sedov_118000.out'
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

n1=nx(0)
n2=nx(1)
n3=nx(2)

x_code=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)

wi=dblarr(n1,n2,n3)

readu,1,x_code
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor
print, n1,n2,n3

;*************** END read old ini file ***************



;tek_color

R=8.31e3
mu=1.2d0
gamma=1.66666667
T=dblarr(n1)

nz=n_elements(zax)



nn=200 ; the number of elements in new arrays 


;***********************************save new ini file
close,1
;openw,1,'/fastdata/cs1mkg/smaug/spicule4_nob/zerospic1_asc_290001.ini'
;openw,1,'/fastdata/cs1mkg/smaug/spicule5_nob/zerospic1_asc_167501.ini'
;openw,1,'/fastdata/cs1mkg/smaug/spicule4a_nob/zerospic1_asc_149001.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/out/spicule5b4/zerospic1_asc_555001.ini'
;openw,1,'/fastdata/cs1mkg/smaug/em6b4_bhor120/zerospic1_asc_84000.ini'
;openw,1,'/fastdata/cs1mkg/smaug/em5b4_bhor120/zerospic1_asc_86000.ini'
;openw,1,'/fastdata/cs1mkg/smaug/spicule8b0_nquotaob/zerospic1_asc_584001.ini'
;openw,1,'/fastdata/cs1mkg/smaug/spicule8b0_nob/zerospic1_asc_643001.ini'
;openw,1,'/fastdata/cs1mkg/smaug/scnstrh5b0b_nob/zerospic1_asc_391000.ini'
;openw,1,'/data/ap1vf/3D_modif_200_100_100.ini',/f77_unf
;openw,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_asc.ini'
;openw,1,'/fastdata/cs1mkg/smaug/spicule5b0_3d/zerospic1_asc_279001.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert2_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/spicule4b0_b1_3d/zerospic1_asc_39000.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/spicule4b0_b2_3d/zerospic1_asc_5000.ini'
;openw,1,'/fastdata/cs1mkg/smaug/spic5b0_3d/zerospic1_asc_612000.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert1_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/spicule4b0_bom500_3d/zerospic1_asc_434000.ini'
;openw,1,'/data/cs1mkg/vac3d/3D_196_100_100_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/spicule4b3_3_3d/zerospic1_asc_360000.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/configs/3D_128_spic_btube1_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/configs/3D_128_spic_btube2_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/configs/3D_128_spic_btube3_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/configs/3D_128_spic_btube4_asc.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/3D_128_spic_bvert10G_asc.ini'
openw,1,'/data/cs1mkg/smaug_spicule1/spicule6p7a_0_0_3d/zerospic1_asc_1400000.ini'
;openw,1,'/data/cs1mkg/smaug_spicule1/spicule4p3a_0_1_3d/zerospic1_asc_1206000.ini'

;openw,1,'//data/cs1mkg/smaugtemp/smaug/out/sedov_asc_118000.ini'
;openr,1,'/data/cs1mkg/smaugtemp/smaug/out/sedov_118000.out'

printf,1, FORMAT='(%"%s ")',headline
printf,1,FORMAT='(%"%d %g %d %d %d ")',it,time,ndim,neqpar,nw
printf,1,FORMAT='(%"%d %d %d ")',nx(0),nx(1),nx(2)
printf,1,FORMAT='(%"%g %g %g %g %g %g %g")',eqpar(0),eqpar(1),eqpar(2),eqpar(3),eqpar(4),eqpar(5),eqpar(6)
printf,1,FORMAT='(%"%s ")',varname
;printf,1,x_code

i1=long(1)
i2=long(1)
i3=long(1)

for i3=0, n3-1 do begin
for i2=0, n2-1 do begin
for i1=0, n1-1 do begin
       ; printf,1,x_code(i1,i2,i3,0),x_code(i1,i2,i3,1),w(i1,i2,i3,0),w(i1,i2,i3,1),w(i1,i2,i3,2),w(i1,i2,i3,3),w(i1,i2,i3,4),w(i1,i2,i3,5),w(i1,i2,i3,6),w(i1,i2,i3,7),w(i1,i2,i3,8),w(i1,i2,i3,9),FORMAT='(15F5)'
 printf,1, FORMAT='(%"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ")',x_code(i1,i2,i3,0),x_code(i1,i2,i3,1),x_code(i1,i2,i3,2),w(i1,i2,i3,0),w(i1,i2,i3,1),w(i1,i2,i3,2),w(i1,i2,i3,3),w(i1,i2,i3,4),w(i1,i2,i3,5),w(i1,i2,i3,6),w(i1,i2,i3,7),w(i1,i2,i3,8),w(i1,i2,i3,9),w(i1,i2,i3,10),w(i1,i2,i3,11),w(i1,i2,i3,12)
	;for iw=0,nw-1 do begin
	; wi=w(*,*,iw)
	; printf,1,w(i1,i2,i3,iw)
	;endfor
endfor
endfor
endfor


 
close,1

print, 'complete'
end





