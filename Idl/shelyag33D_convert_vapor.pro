
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 


dim_x=128
dim_y=128
dim_z=128

close,1 

;openr,1,'/data/ap1vf/vxx.040000'
;vx=fltarr(dim_x,dim_y,dim_z)
;readu, 1, vx
;close,1 
;openr,1,'/data/ap1vf/vyy.040000'
;vy=fltarr(dim_x,dim_y,dim_z)
;readu, 1, vy
;close,1 
;openr,1,'/data/ap1vf/vzz.040000'
;vz=fltarr(dim_x,dim_y,dim_z)
;readu, 1, vz
;close,1 


;///////////////////////////////////////
;-- VAPOR
;//////////////////////////////////////
dim = [dim_x,dim_y,dim_z]

num_levels = 0

mfd = vdf_create(dim,num_levels)
timesteps = 1

vdf_setnumtimesteps, mfd,timesteps

varnames = ['vx','vy','vz']

vdf_setvarnames, mfd, varnames

extents = [133333.33, 1953.1, 1953.1, 5955555.6e0, 4.0d6, 4.0d6]
vdf_setextents, mfd, extents

;
; Set a global comment
;

vdf_setcomment, mfd, 'This is my SAC data'

attribute_name = 'MyMetadata'
f = findgen(100)
Vdf_setdbl,mfd,attribute_name, f
vdffile = 'shelyag.vdf'
vdf_write, mfd, vdffile
vdf_destroy, mfd




vdffile = '/data/cs1mkg/smaug_spicule1/spicule4b0_b4_3d/vap/shelyag33D.vdf'
dfd = vdc_bufwritecreate(vdffile)


; Get the data volume that we wish to store.




; Prepare the data set for writing. We need to identify
; the time step and the name of the variable that
; we wish to store. In this case, the first time step,
; zero, and the variable named ÔvxÕ
;
dim= [dim_x,dim_y,dim_z]

vdc_openvarwrite, dfd, nn, 'bx'
for z = 0, dim_z-1 do begin
vdc_bufwriteslice, dfd, vx[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'by'
for z = 0, dim_z-1 do begin
vdc_bufwriteslice, dfd, vy[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'bz'
for z = 0, dim_z-1 do begin
vdc_bufwriteslice, dfd, vz[*,*,z]
endfor
vdc_closevar, dfd

                          ;
; Write (transform) the volume to the data set one
; slice at a time


;////////////////////////////////////////
;////////////////////////////////////////

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
