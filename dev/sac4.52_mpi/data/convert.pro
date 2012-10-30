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
;loadct,4
;mixct
endif else begin
loadct,0
tek_color
endelse




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


nn=0


close,1

;openr,1,'/tslive/3D_396_60_60_tt.out',/f77_unf

openr,1,'/data/ap1vf/3D_tube_196_100_100t.out',/f77_unf

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



xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

wi=dblarr(n1,n2,n3)

w=fltarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
endfor


print,time, nn

zmin=0
zmax=196

;///////////////////////////////////////
;-- VAPOR
;//////////////////////////////////////
dim = [zmax-zmin+1,n2,n3]

num_levels = 0

mfd = vdf_create(dim,num_levels)
timesteps = 1

vdf_setnumtimesteps, mfd,timesteps

varnames = ['bx','by','bz','vx','vy','vz','rho','bt']

vdf_setvarnames, mfd, varnames

extents = [24120.603d0, 0.0, 0.0, 1591959.8d0, 2.0d6, 2.0d6]
vdf_setextents, mfd, extents

;
; Set a global comment
;

vdf_setcomment, mfd, 'This is my SAC data'

attribute_name = 'MyMetadata'
f = findgen(100)
Vdf_setdbl,mfd,attribute_name, f
vdffile = 'sac.vdf'
vdf_write, mfd, vdffile
vdf_destroy, mfd




vdffile = '/data/ap1vf/3D_tube_196_100_100t.vdf'
dfd = vdc_bufwritecreate(vdffile)


; Get the data volume that we wish to store.




sac_rho=reform(w(zmin:zmax,*,*,0)+w(zmin:zmax,*,*,9))

sac_vx=reform(w(zmin:zmax,*,*,3)/(w(zmin:zmax,*,*,0)+w(zmin:zmax,*,*,9)))
sac_vy=reform(w(zmin:zmax,*,*,1)/(w(zmin:zmax,*,*,0)+w(zmin:zmax,*,*,9)))
sac_vz=reform(w(zmin:zmax,*,*,2)/(w(zmin:zmax,*,*,0)+w(zmin:zmax,*,*,9)))

sac_bx=reform(w(zmin:zmax,*,*,6)+w(zmin:zmax,*,*,11))
sac_by=reform(w(zmin:zmax,*,*,7)+w(zmin:zmax,*,*,12))
sac_bz=reform(w(zmin:zmax,*,*,5)+w(zmin:zmax,*,*,10))

sac_bt=sqrt((w(zmin:zmax,*,*,5)+w(zmin:zmax,*,*,10))^2.0+  $
            (w(zmin:zmax,*,*,6)+w(zmin:zmax,*,*,11))^2.0+  $
	    (w(zmin:zmax,*,*,7)+w(zmin:zmax,*,*,12))^2.0)


bx = sac_bx
by = sac_by
bz = sac_bz
vx = sac_vx
vy = sac_vy
vz = sac_vz
bt = sac_bt
rho = sac_rho

; Prepare the data set for writing. We need to identify
; the time step and the name of the variable that
; we wish to store. In this case, the first time step,
; zero, and the variable named ÔvxÕ
;
print, time
dim= [zmax-zmin+1,n2,n3]

vdc_openvarwrite, dfd, nn, 'bx'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, bx[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'by'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, by[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'bz'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, bz[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vx'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vx[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vy'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vy[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vz'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vz[*,*,z]
endfor
vdc_closevar, dfd
                          ;
; Write (transform) the volume to the data set one
; slice at a time


;////////////////////////////////////////
;////////////////////////////////////////

nn=nn+1

endwhile


;-- CLOSE VAPOR

;An Overview of VAPOR Data Collections 12
; Close the currently opened variable/time-step. We're
; done writing to it
;
;
; Destroy the "buffered write" data transformation
; object. We're done with it.
;
vdc_bufwritedestroy, dfd

end
