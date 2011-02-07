; docformat = 'rst'

;+
; Creates a line-integral convolution flow visualization of a 2D vector field.
;
; :Author: Michael Galloy
;
; :Examples: 
;    See the main level program associated with this file::
;
;       IDL> .run mg_lic
;-


;+
; Generate streamline from a given point.
;
; :Params:
;    i : in, required, type=long
;       column of streamline origin
;    j : in, required, type=long
;       row of streamline origin
;
; :Keywords:
;    origin : out, optional, type=fltarr(2)
;       origin of the streamline
;    fwd : out, optional, type="fltarr(2, nSegments)"
;       streamline in the forward direction
;    bwd : out, optional, type="fltarr(2, nSegments)"
;       streamline in the backward direction
;    n_segments : in, required, type=long
;       number of segments to compute
;    stepsize : in, required, type=float
;       stepsize 
;    grid : in, required, type="fltarr(2, xsize, ysize)"
;       vector field to visualize
;-
pro mg_lic_streamline, i, j, origin=origin, fwd=fwd, bwd=bwd, $
                       n_segments=nSegments, stepsize=stepsize, grid=grid
  compile_opt strictarr
  
  fwd = fltarr(2, nSegments)
  nFwd = 0L
  bwd = fltarr(2, nSegments)
  nBwd = 0L
  
  origin = [i + 0.5, j + 0.5] 
  f = origin
  b = origin

  fwdValid = 1B
  bwdValid = 1B
  
  for k = 0L, nSegments - 1L do begin
    if (fwdValid) then begin
      mg_lic_rk, f, stepsize=stepsize, grid=grid, valid=fwdValid
      fwd[*, k] = f
      nFwd++
    endif
    
    if (bwdValid) then begin
      mg_lic_rk, b, stepsize=-stepsize, grid=grid, valid=bwdValid
      bwd[*, k] = b
      nBwd++
    endif
  endfor
  
  fwd = fwd[*, 0:nFwd - 1L]
  bwd = bwd[*, 0:nBwd - 1L]
end


;+
; Look up a point in the vector data grid.
;
; :Params:
;    p : in, required, type=fltarr(2)
;       point to look up
;
; :Keywords:
;    grid : in, required, type="fltarr(2, xsize, ysize)"
;       vector data
;    valid : out, optional, type=boolean
;       true if the point is inside the vector data grid
;-
function mg_lic_grid, p, grid=grid, valid=valid
  compile_opt strictarr
  
  valid = 1B
  sz = size(grid, /dimensions)
  if (p[0] ge 0 && p[0] lt sz[1] && p[1] ge 0 && p[1] lt sz[2]) then begin
    return, reform(grid[*, p[0], p[1]])
  endif
  
  valid = 0B
  return, [!values.f_nan, !values.f_nan]
end


;+
; Calculates the next point in a streamline.
;
; :Params:
;    p : in, out, required, type=fltarr(2)
;       starting point
;
; :Keywords:
;    stepsize : in, required, type=float
;       stepsize
;    grid : in, required, type="fltarr(2, xsize, ysize)"
;       vector field to visualize
;-
pro mg_lic_rk, p, stepsize=stepsize, grid=grid, valid=valid
  compile_opt strictarr

  v = mg_lic_grid(p, grid=grid, valid=valid)
  if (~valid) then return
  if (~array_equal(v, [0.0, 0.0])) then v /= sqrt(total(v^2, /preserve_type))
  
  k1 = v * stepsize
  v = mg_lic_grid(p + k1 * 0.5, grid=grid, valid=valid)
  if (~valid) then return
  if (~array_equal(v, [0.0, 0.0])) then v /= sqrt(total(v^2, /preserve_type))
  
  k2 = v * stepsize
  v = mg_lic_grid(p + k2 * 0.5, grid=grid, valid=valid)
  if (~valid) then return
  if (~array_equal(v, [0.0, 0.0])) then v /= sqrt(total(v^2, /preserve_type))
  
  k3 = v * stepsize
  v = mg_lic_grid(p + k3, grid=grid, valid=valid)
  if (~valid) then return
  if (~array_equal(v, [0.0, 0.0])) then v /= sqrt(total(v^2, /preserve_type))
  
  k4 = v * stepsize
  p += k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
end


;+
; Converts a list of (x, y) points into a list of indices.
;
; :Returns: lonarr
;
; :Params:
;    pts : in, required, type="fltarr(2, n)"
;       list of points
;    array : in, required, type=2D array
;       array that the points refer to
;-
function mg_lic_getindices, pts, array
  compile_opt strictarr
  
  sz = size(array, /dimensions)
  lpts = long(pts)
  return, reform(lpts[0, *]) + reform(lpts[1, *]) * sz[1]
end


;+
; Pads the original vector field on each size and copies the edge rows out into 
; the padding.
; 
; :Returns: fltarr(2, xsize, ysize)
;
; :Params:
;    u : in, required, type=fltarr
;       original x-components of the vector field
;    v : in, required, type=fltarr
;       original y-components of the vector field
;
; :Keywords:
;    xsize : out, optional, type=long
;       new xsize of the vector field
;    ysize : out, optional, type=long
;       new ysize of the vector field
;-
function mg_lic_pad, u, v, padding=padding, xsize=xsize, ysize=ysize
  compile_opt strictarr
  
  ; dimensions of the problem
  sz = size(u, /dimensions)
  xsize = sz[0] + 2 * padding
  ysize = sz[1] + 2 * padding
  
  ; original u, v data
  new_u = fltarr(xsize, ysize)
  new_v = fltarr(xsize, ysize)
  new_u[padding, padding] = u
  new_v[padding, padding] = v
  
  if (padding gt 0) then begin
    ; left side
    new_u[0:padding - 1L, padding:ysize - padding - 1L] $
      = rebin(u[0, *], padding, ysize - 2 * padding)
    new_v[0:padding - 1L, padding:ysize - padding - 1L] $
      = rebin(v[0, *], padding, ysize - 2 * padding)
        
    ; right side
    new_u[xsize - padding:xsize-1L, padding:ysize - padding - 1L] $
      = rebin(u[xsize - 2 * padding - 1L, *], padding, ysize - 2 * padding)
    new_v[xsize-padding:xsize-1L, padding:ysize-padding-1L] $
      = rebin(v[xsize - 2 * padding - 1L, *], padding, ysize - 2 * padding)
            
    ; bottom
    new_u[padding:xsize - padding - 1L, 0L:padding - 1L] $
      = rebin(u[*, 0], xsize - 2 * padding, padding)
    new_v[padding:xsize - padding - 1L, 0L:padding - 1L] $
      = rebin(v[*, 0], xsize - 2 * padding, padding)
  
    ; top
    new_u[padding:xsize - padding - 1L, ysize - padding:ysize - 1L] $
      = rebin(u[*, ysize - 2 * padding - 1L], xsize - 2 * padding, padding)
    new_v[padding:xsize - padding - 1L, ysize - padding:ysize - 1L] $
      = rebin(v[*, ysize - 2 * padding - 1L], xsize - 2 * padding, padding)
  endif
  
  vecdata = fltarr(2, xsize, ysize)
  vecdata[0, *, *] = new_u
  vecdata[1, *, *] = new_v  
  
  return, vecdata
end


;+
; Removes padding and normalizes result image.
;
; :Returns: bytarr
;
; :Params:
;    idata : in, required, type=fltarr
;       line integral convolution data
;    hitdata : in, required, type=lonarr
;       number of hits per pixel
;
; :Keywords:
;    padding : in, required, type=long
;       size of padding on each side of the image
;-
function mg_lic_normalize, idata, hitdata, padding=padding
  compile_opt strictarr
  
  sz = size(idata, /dimensions)
  
  new_idata = idata[padding:sz[0] - padding - 1L, padding:sz[1] - padding - 1L]
  new_hitdata = hitdata[padding:sz[0] - padding - 1L, padding:sz[1] - padding - 1L]
  
  return, bytscl(new_idata / new_hitdata)
end


;+
; Creates a line-integral convolution flow visualization of the given vector
; field.
; 
; :Returns: bytarr(3, xsize, xsize)
;
; :Bugs:
;    This routine is quite slow; it should probably be re-written as a DLM.
;
; :Params:
;    u : in, required, type="fltarr(xsize, ysize)"
;       x-component of vector field
;    v : in, required, type="fltarr(xsize, ysize)"
;       y-component of vector field
;
; :Keywords:
;    padding : in, optional, type=long, default=10L
;-
function mg_lic, u, v, padding=padding
  compile_opt strictarr
  
  ; constants
  maxIterations = 10L
  maxStepSize = 0.5
  minNumHits = 12L
  
  myPadding = n_elements(padding) eq 0 ? 10L : padding
  
  vecdata = mg_lic_pad(u, v, padding=myPadding, xsize=xsize, ysize=ysize)
  
  maxVecMag = max(sqrt(u * u + v * v))
  
  hitData = lonarr(xsize, ysize) 
  idata = lonarr(xsize, ysize)
  texdata = 255.0 * randomu(seed, xsize, ysize)

  for i = 0L, xsize - 1L do begin
    for j = 0L, ysize - 1L do begin
      ; generate streamline
      mg_lic_streamline, i, j, origin=origin, fwd=fwd, bwd=bwd, $
                         n_segments=maxIterations + 10L, $
                         stepsize=maxStepSize, $
                         grid=vecdata
      
      ; compute I
      numValid = (n_elements(fwd) + n_elements(bwd)) / 2
      find = mg_lic_getindices(fwd, texdata)
      bind = mg_lic_getindices(bwd, texdata)
      idata[i, j] = total(texdata[find] + texdata[bind]) / numvalid
      hitdata[i, j]++
      
      ; compute Ifwd
      ; compute Ibwd
    endfor
  endfor
  
  ; normalize and remove padding
  return, mg_lic_normalize(idata, hitdata, padding=padding)
end

startTime = systime(/seconds)

scale = 4L

restore, filepath('galaxy.dat', subdir=['examples','data'])
; globalwinds.dat
;ivector, u, v, x, y

u = rebin(u, 128L * scale, 64L * scale)
v = rebin(v, 128L * scale, 64L * scale)
x = rebin(x, 128L * scale)
y = rebin(y, 64L * scale)

im = mg_lic(u, v, padding=10L)

window, xsize=128L * scale, ysize=64L * scale * 3, $
        /free, title='LIC for globalwinds.dat'
tv, im, 0

h = bytarr(128L * scale, 64L * scale)
mag = sqrt(u * u + v *v)
m = mag / max(mag)
v = im / 255.0

color_convert, h, m, v, r, g, b, /hsv_rgb
im2 = bytarr(3L, 128L * scale, 64L * scale)
im2[0, *, *] = r
im2[1, *, *] = g
im2[2, *, *] = b

tv, im2, 1, true=1

file = dialog_pickfile(title='Choose color table file')
if (file ne '') then begin
  xloadct, file=file, /block
endif else begin
  xloadct, /block
endelse
;loadct, 10

tvlct, ctr, ctg, ctb, /get

lic_image = bytarr(3L, 128L * scale, 64L * scale)
lic_image[0, *, *] = im
lic_image[1, *, *] = im
lic_image[2, *, *] = im

mscaled = bytscl(m)
m_image = fltarr(3L, 128L * scale, 64L * scale)
m_image[0, *, *] = ctr[mscaled] / 255.0
m_image[1, *, *] = ctg[mscaled] / 255.0
m_image[2, *, *] = ctb[mscaled] / 250.0
im3 = lic_image * m_image

tv, im3, 2, true=1

write_png, 'globalwinds-' + strtrim(scale, 2) + '.png', im
write_png, 'globalwinds-mag-' + strtrim(scale, 2) + '.png', im2
write_png, 'globalwinds-mag2-' + strtrim(scale, 2) + '.png', im3

print, format='(%"Total time to compute: %0.1f seconds")', $
       systime(/seconds) - startTime 
       
end
