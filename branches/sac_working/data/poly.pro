; Create an empty, 3-D array:  

n=20
SPHERE = FLTARR(n,n,n)  


x=findgen(n)/(n*1.d0)*3.d0
y=findgen(n)/(n*1.d0)*3.d0
z=findgen(n)/(n*1.d0)*3.d0 
  
; Create the spherical dataset:  
FOR i=0,n-1 DO FOR j=0,n-1 DO FOR k=0,n-1 DO $  
   SPHERE(i, j, k) = SQRT((x[i]-1.5)^2 + (y[j]-1.5)^2 + (z[k]-1.5)^2)  
  
; Find the vertices and polygons for a density level of 8:  
SHADE_VOLUME, SPHERE, 1.0, V, P  

mi_x=-10
ma_x=10
mi_y=-10
ma_y=10
mi_z=-10
ma_z=10


SURFACE, DIST(4), /NODATA, /SAVE, XRANGE=[mi_x, ma_x], $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z], XSTYLE=1, $
   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, xtitle='x [Mm]',ytitle='y [Mm]',title='Field Lines', $
   ztitle='z [Mm]', POSITION=[0.2, 0.1, 0.95, 0.95, 0.1, 0.95], az=30.0, ax=30.0
   ;BACKGROUND=FSC_Color('ivory')   
stop   
; Set up an appropriate 3-D transformation so we can see the  
; sphere. This step is very important:  
SCALE3, XRANGE=[0,20.0], YRANGE=[0,20.0], ZRANGE=[0,20.0]  
  
; Render the image. Note that the T3D keyword has been set so that   
; the previously-established 3-D transformation is used:  
;image = POLYSHADE(V, P, /T3D)  
  
; Display the image:  
TV, POLYSHADE(V, P, /T3D)
end
