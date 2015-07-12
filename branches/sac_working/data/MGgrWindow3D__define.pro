;+
; Get properties of the MGgrWindow3D.
; 
; @keyword eye_separation {out}{optional}{type=float} number of
;          degrees of the cone formed by drawing lines from each eye to the
;          origin of the view
; @keyword _ref_extra {out}{optional}{type=keywords} keyword of 
;          IDLgrWindow::getProperty
;-
pro mggrwindow3d::getProperty, eye_separation=eye_separation, _ref_extra=e
  compile_opt strictarr

  if (arg_present(eye_separation)) then begin
      self.converter->getProperty, eye_separation=eye_separation
  endif

  if (n_elements(e) gt 0) then begin
      self->IDLgrWindow::getProperty, _strict_extra=e
  endif
end


;+
; Set properties of the MGgrWindow3D. Must intercept DIMENSIONS property to set 
; the converter's buffer size correctly; otherwise, just pass along stuff to 
; IDLgrWindow's setProperty method.
;
; @keyword dimensions {in}{optional}{type=intarr(2)} dimensions of the window
; @keyword eye_separation {out}{optional}{type=float} number of
;          degrees of the cone formed by drawing lines from each eye to the
;          origin of the view
; @keyword _extra {in}{optional}{type=keywords} keywords of IDLgrWindow's
;          "setProperty" method are accepted
;-
pro mggrwindow3d::setProperty, dimensions=dimensions, $
                               eye_separation=eye_separation, _extra=e
    compile_opt strictarr

    self->IDLgrWindow::setProperty, _extra=e

    if (n_elements(dimensions) gt 0) then begin
        self->IDLgrWindow::setProperty, dimensions=dimensions
        self.converter->setProperty, dimensions=dimensions
    endif

    if (n_elements(eye_separation) gt 0) then begin
        self.converter->setProperty, eye_separation=eye_separation
    endif
end


;+
; Draw the picture in 3D.
;
; @param opicture {in}{optional}{type=obj ref} the view, viewgroup, or scene to
;        be drawn; if the GRAPHICS_TREE property is set to a valid picture,
;        then this argument must <em>not</em> be given
;-
pro mggrwindow3d::draw, opicture
    compile_opt strictarr
    on_error, 2

    self->getProperty, graphics_tree=graphics_tree
    ipicture = obj_valid(opicture) ? opicture : graphics_tree

    oview = self.converter->convert(ipicture)

    self->idlgrwindow::draw, oview
end


;+
; Free resources.
;-
pro mggrwindow3d::cleanup
    compile_opt strictarr

    self->idlgrwindow::cleanup
    obj_destroy, self.converter
end


;+
; Initialize Window3D.
;
; @returns 1 for success, o/w for failure
; @keyword eye_separation {in}{optional}{type=float}{default=4.0} number of
;          degrees of the cone formed by drawing lines from each eye to the
;          origin of the view
; @keyword _extra {in}{optional}{type=keywords} keywords to IDLgrWindows "init"
;          method are accepted
;-
 function mggrwindow3d::init, eye_separation=eye_separation, $
                              dimensions=dimensions, _extra=e
    compile_opt strictarr

    if (~self->IDLgrWindow::init(dimensions=dimensions, _strict_extra=e)) then return, 0

    if (n_elements(dimensions) eq 0) then begin
        case strlowcase(!version.os_family) of
            'unix' : begin
                dims = [pref_get('idl_gr_x_width'), pref_get('idl_gr_x_height')]
            end
            'windows' : begin
                dims = [pref_get('idl_gr_win_width'), pref_get('idl_gr_win_height')]
            end
        endcase
    endif else dims = dimensions

    self.converter = obj_new('MGgr3DConverter', $
                             eye_separation=eye_separation, $
                             dimensions=dims)

    return, 1
end


;+
; Destination for object graphics that automatically creates a 3d anaglyph
; appropriate to view with red-blue glasses.
;
; @field converter object which takes a view and converts to a 3D anaglyph
;-
pro mggrwindow3d__define
    compile_opt strictarr

    define = { mggrwindow3d, inherits idlgrwindow, $
        converter : obj_new() $
        }
end
