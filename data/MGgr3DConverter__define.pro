;+
; Creates a combined image from images from the left and right eyes where the
; left eye is "shaded" red and the right eye is "shaded" blue.
;
; @private
; @returns bytarr(3, xsize, ysize)
; @param left_image {in}{optional}{type=bytarr(3, xsize, ysize)} image from
;        left eye
; @param right_image {in}{optional}{type=bytarr(3, xsize, ysize)} image from
;        right eye
;-
function mggr3dconverter::combineImages, left_image, right_image
    compile_opt strictarr

    ; define combined_image to the correct size
    combined_image = left_image * 0B
    dims = size(left_image, /dimensions)

    temp_left = byte(total(fix(left_image), 1) / 3)
    temp_right = byte(total(fix(right_image), 1) / 3)

    combined_image[0, 0, 0] = Reform(temp_left, 1, dims[1], dims[2])
    combined_image[1, 0, 0] = Reform(temp_right, 1, dims[1], dims[2])
    combined_image[2, 0, 0] = Reform(temp_right, 1, dims[1], dims[2])

    return, combined_image
end


;+
; Rotates "top-level" models of the given picture by the given number of
; degrees about the y-axis.
;
; @private
; @param picture {in}{required}{type=obj ref} the view, viewgroup, or scene to
;        be drawn
; @param degrees {in}{required}{type=float} number of degrees to rotate
;        "top-level" models
;-
pro mggr3dconverter::rotateModels, picture, degrees
    compile_opt strictarr
    
    ; if picture is a model then rotate it, but don't rotate models inside it
    if (obj_isa(picture, 'IDLgrModel')) then begin
        picture->rotate, [0, 1, 0], degrees
        return
    endif

    if (obj_isa(picture, 'IDL_Container')) then begin
        items = picture->get(/all, count=count)
        for i = 0L, count - 1 do begin
            self->rotateModels, items[i], degrees        
        endfor
    endif
end


;+
; Converts a standard object graphics picture to a view containing a 3D image.
;
; @returns IDLgrView object reference
; @param opicture {in}{optional}{type=obj ref} the view, viewgroup, or scene to
;        be drawn; if the GRAPHICS_TREE property is set to a valid picture,
;        then this argument must <em>not</em> be given
;-
function mggr3dconverter::convert, opicture

    ; rotate "top-level" models for left eye
    self->rotateModels, opicture, self.eye_separation / 2.

    ; draw picture to left eye buffer
    self.obuffer->draw, opicture

    ; get data out of left eye buffer
    oleft_image = self.obuffer->read()
    oleft_image->getProperty, data=left_image
    obj_destroy, oleft_image

    ; rotate "top-level" models for right eye
    self->rotateModels, opicture, - self.eye_separation

    ; draw picture to right eye buffer
    self.obuffer->draw, opicture

    ; get data out of left eye buffer
    oright_image = self.obuffer->read()
    oright_image->getProperty, data=right_image
    obj_destroy, oright_image

    ; rotate "top-level" models back to center
    self->rotateModels, opicture, self.eye_separation / 2.

    combined_image = self->combineImages(left_image, right_image)

    self.oimage->setProperty, data=combined_image

    return, self.oview
end


;+
; Get properties of the converter.
; 
; @keyword eye_separation {out}{optional}{type=float} number of degrees of the
;          cone formed by drawing lines from each eye to the origin     
; @keyword dimensions {out}{optional}{type=intarr(2)} dimensions of the window
;-
pro mggr3dconverter::getProperty, eye_separation=eye_separation, $
                                  dimensions=dimensions
    compile_opt strictarr

    if (arg_present(eye_separation)) then begin
        eye_separation = self.eye_separation
    endif

    if (arg_present(dimensions)) then begin
        self.obuffer->getProperty, dimensions=dimensions
    endif

end


;+
; Set properties of the converter.
;
; @keyword eye_separation {in}{optional}{type=float} number of degrees of the
;          cone formed by drawing lines from each eye to the origin  
; @keyword dimensions {in}{optional}{type=intarr(2)} dimensions of the window
;-
pro mggr3dconverter::setProperty, eye_separation=eye_separation, $
                                  dimensions=dimensions
    compile_opt strictarr

    if (n_elements(eye_separation) gt 0) then begin
        self.eye_separation = eye_separation
    endif

    if (n_elements(dimensions) gt 0) then begin
        self.oview->setProperty, viewplace_rect=[0, 0, dimensions]
        self.obuffer->setProperty, dimensions=dimensions
    endif
end


;+
; Free resources.
;-
pro mggr3dconverter::cleanup
    compile_opt strictarr

    obj_destroy, [self.oview, self.obuffer]
end


;+
; Initialize Window3D.
;
; @returns 1 for success, o/w for failure
; @keyword eye_separation {in}{optional}{type=float}{default=4.0} number of
;          degrees of the cone formed by drawing lines from each eye to the
;          origin
; @keyword dimensions {in}{required}{type=intarr(2)} dimensions of the window
; @keyword picture {out}{optional}{type=IDLgrView} view which will contain a 3D
;          image; the same view is updated each time that "convert_3d_picture"
;          method is called
; @keyword _extra {in}{optional}{type=keywords} keywords to IDLgrWindows "init"
;          method are accepted
;-
function mggr3dconverter::init, eye_separation=eye_separation, $
                                dimensions=dimensions, picture=picture, _extra=e
    compile_opt strictarr

    self.eye_separation = n_elements(eye_separation) eq 0 ? 4.0 : eye_separation

    self.obuffer = obj_new('IDLgrBuffer', dimensions=dimensions)

    self.oview = obj_new('IDLgrView', viewplane_rect=[0, 0, dimensions])

    omodel = obj_new('IDLgrModel')
    self.oview->add, omodel

    self.oimage = obj_new('IDLgrImage')
    omodel->add, self.oimage

    return, 1
end


;+
; Helper object to transform a normal object graphics scene to a 3d picture.
;
; @field eye_separation number of degrees of the cone formed by drawing lines
;        from each eye to the origin
; @field obuffer IDLgrBuffer to send left and right eye images to and extract
; @field oview IDLgrView to contain the 3D image
; @field oimage IDLgrImage actually being displayed
; @author Michael Galloy
;-
pro mggr3dconverter__define
    compile_opt strictarr

    define = { mggr3dconverter, $
        eye_separation : 0.0, $
        obuffer : obj_new(), $
        oview : obj_new(), $
        oimage : obj_new() $
    }
end
