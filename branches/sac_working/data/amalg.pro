;+
; $Id: scc_euvi_pairs.pro,v 1.3 2007/10/11 16:39:41 reduce Exp $
;
; Project   : STEREO SECCHI
;                   
; Name      : SCC_EUVI_PAIRS
;               
; Purpose   : Create STEREO pairs for EUVI data
;               
; Explanation: Given a list of EUVI images from the A-spacecraft (ONLY!!!), will find the corresponding
;              B-spacecraft image and make the STEREO pairs. The user can either return two images (fully
;              scaled/repositioned) or an anaglyph image for viewing with red/blue glasses.
;              SEE ROUTINE HEADER FURTHER BELOW FOR MORE INFO AND INSTRUCTIONS
;               
; Use       : IDL> scc_euvi_pairs,filelist,outdir='mydir',outsize=1024,/color_on
;    
; Inputs    : (optional) filelist -- list of **SPACECRAFT A** EUVI images
;             (If no input is given, scclister is invoked)
;               
; Outputs   : STEREO pairs or anaglyph images as pngs
;
; Keywords  :   OUTDIR: output directory. Default is current working directory
;               OUTSIZE: output size of images. Default is full-res 2048x2048. Images can be ANY SQUARE SIZE!
;               COLOR_ON: if set, images will have appropriate color tables applied
;               NRLVIEWER: creates the filenames with _L and _R (specifically for the NRL stereo viewer setup
;                          but there's no reason why anyone can't use this option.)
;               ANAGLYPH: Create a single anaglyph png
;               LISTER: Force scclister to appear, regardless of whether a filelist is entered or not
;               SAVELIST: Create a text file containing the file names of the files to be processed (creates file
;                         in current working directory)
;               ROI:   Interactively define a sub-region of interest rather than full field. (Might not work too well yet)
;
; Common    : None.
;
; Calls     : SECCHI_PREP
;               
; Restrictions: !!!!!!!!! ONLY GIVE SPACECRAFT _A_ IMAGES AS INPUT !!!!!!!!
;               
; Side effects: Possible headache if you frequently use the "cross-eye" method to fuse the 3D pairs
;               
; Category    : SECCHI, Display
;               
; Prev. Hist. : None.
;
; Written     : Karl Battams, NRL/I2, May-Jul 2007
;               
; $Log: scc_euvi_pairs.pro,v $
; Revision 1.3  2007/10/11 16:39:41  reduce
; Changed something. Not sure what. I am sure it works well though. Karl
;
; Revision 1.2  2007/10/02 15:33:08  reduce
; Got rid of the ugly is_even function (thanks Bogdan!). Karl B.
;
; Revision 1.1  2007/07/26 15:34:48  reduce
; Initial release. Karl B.
;
;
;-

; This function makes the code a little nicer.  You should probably ignore this and
; jump down a few lines to the main section...


; ###################################################
function mkanaglyph,imr,iml

szr=size(imr,/dim)

img=bytarr(3,szr[0],szr[1])

img[0,*,*]=iml
img[1,*,*]=imr
img[2,*,*]=imr*0.1

return,img
end
; ###################################################

; */*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/ ;
; #################################################################################### ;
; ############################ Here is the main routine  ############################# ;
; #################################################################################### ;
; */*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/ ;


PRO scc_euvi_pairs,filelist,outdir=outdir,outsize=outsize,color_on=color_on, $
        nrlviewer=nrlviewer,anaglyph=anaglyph,lister=lister,savelist=savelist,roi=roi

; Given an optional input list of EUVI images from the "A" spacecraft, this routine will locate the 
; corresponding 'B' image, then process both into "pretty" pictures. If files are not specified,
; scclister is called. The user should be sure and pick only SC-A files from scclister.
; Note that the routine will resize and recenter the Sun in the B images so that it matches the
; Sun in te A images.
; Output is png files.
;
; KEYWORDS: 
;       OUTDIR: output directory. Default is current working directory
;       OUTSIZE: output size of images. Default is full-res 2048x2048. Images can be ANY SQUARE SIZE!
;       COLOR_ON: if set, images will have appropriate color tables applied
;       NRLVIEWER: creates the filenames with _L and _R (specifically for the NRL stereo viewer setup
;                  but there's no reason why anyone can't use this option.)
;       ANAGLYPH: Create a single anaglyph png
;       LISTER: Force scclister to appear, regardless of whether a filelist is entered or not
;       SAVELIST: Create a text file containing the file names of the files to be processed (creates file
;                 in current working directory)
;       ROI:   Interactively define a sub-region of interest rather than full field.

; EXAMPLES:
;    1. Create a set of 'pretty' stereo pair images from both spacecraft in 1024x1024 size
;    scc_euvi_pairs,myfiles,outdir='mydir',outsize=1024,/color_on
;
;    2. Ceate a set of red/blue anaglyph images in 2048x2048 size
;    scc_euvi_pairs,myfiles,outdir='mydir',/anaglyph
;
;    3. Create set of pretty pngs for NRL stereoviewer setup in 1024x1024 size
;    scc_euvi_pairs,myfiles,outdir='mydir',outsize=1024,/nrlviewer
;
;  NOTE: to do any of the above without specifying a list of files, just omit the 'myfiles' parameter.


; IMPORTANT REMINDERS
;    1) The "myfiles" parameter is a string array of the complete file path to the set of spacecraft-A images ONLY! 
;    ** DO NOT point to the spacecraft-B images -- it finds them itself!!! ** 
;    2) The roi function is perhaps a little rusty still. I make no promises as to it's usefulness.
;    3) After ~June 1st, 2007, the EUVI data starts to look bad when viewed in 3D. For June/early-July, you can probably get
;       away with using sub-regions, as long as they're not on the limb, but after that it gets ugly



n=n_elements(filelist)

IF ( n EQ 0 ) OR keyword_set(lister) then begin
    files=scclister()
    nf=n_elements(files.SC_A)
    if (nf LT 1) THEN BEGIN
        PRINT,'ERROR: No SC-A files returned from SCC lister'
        return
    endif
    filelist=strarr(nf)
    FOR k=0,nf-1 DO filelist[k]=files.SC_A[k]
    n=nf
ENDIF

IF keyword_set(savelist) THEN BEGIN
    openw,1,'euvi_pairs_filelist.txt',/append
    for i=0,n-1 DO printf,1,filelist[i]
    free_lun,1
ENDIF
;stop

if not keyword_set(outsize) then outsize=2048
resize_flag=0
if (outsize gt 2048) THEN BEGIN
    oversized=outsize
    outsize=2048
    resize_flag=1
    PRINT,''
    PRINT,'%%%%%%%%%% WARNING %%%%%%%%%%
    PRINT,'SECCHI_PREP ONLY OUTPUTS IMAGES UP TO 2048x2048!'
    PRINT,'IMAGES WILL GET RESCALED OUTSIDE OF SECCHI_PREP!'
    PRINT,'(There is nothing necessarily wrong with this) 
    PRINT,'%%%%%%%% END WARNING %%%%%%%%
    PRINT,''
endif

if not keyword_set(color_on) then loadct,0

if keyword_set(outdir) then outdir=outdir+'/' else outdir='./'

FOR i=0,n-1 DO BEGIN
    PRINT,'
    PRINT,'Processing file ',strcompress(string(i+1),/remove_all)
    
    ; find corresponding images to list-A
    alist=filelist
    blist=filelist
    tmp=blist[i]
    strput,tmp,'/b/',strpos(blist[i],'/a/')
    strput,tmp,'B.fts',strpos(blist[i],'A.fts')
    blist[i]=tmp
    IF not file_exist(blist[i]) THEN BEGIN
        PRINT,'Could not find file ',blist[i]
        goto,skip
    ENDIF

    ; Loop over all pairs and pump thru secchi prep
    
    if keyword_set(roi) then BEGIN
        secchi_prep,alist[i],ahdr,aim,outsize=2048,color_on=color_on,/smask_on,/calimg_off,/rotate_on,/precommcorrect
        secchi_prep,blist[i],bhdr,bim,outsize=2048,color_on=color_on,/smask_on,/calimg_off,/rotate_on,/precommcorrect
    ENDIF ELSE BEGIN
        secchi_prep,alist[i],ahdr,aim,outsize=outsize,color_on=color_on,/smask_on,/calimg_off,/rotate_on,/precommcorrect
        secchi_prep,blist[i],bhdr,bim,outsize=outsize,color_on=color_on,/smask_on,/calimg_off,/rotate_on,/precommcorrect
    ENDELSE
    
    ; resize the B data
    a_obs_dsun=ahdr.dsun_obs
    b_obs_dsun=bhdr.dsun_obs

    if keyword_set(roi) then BEGIN
        sizetmp=outsize  ; temporarily change the outsize if we're doing a ROI
        outsize=2048
    ENDIF
    
    ratio=b_obs_dsun/a_obs_dsun
    newsize=floor(ratio*outsize) ; get new array dimensions for B
    bim=congrid(bim,newsize, newsize)
    diff_in_size=newsize-outsize
    if ~(diff_in_size mod 2) THEN BEGIN
        half_diff=diff_in_size/2
        xmin=half_diff
        xmax=(newsize - 1)-half_diff
        ymin=half_diff
        ymax=(newsize - 1)-half_diff
        bim=bim[xmin:xmax,ymin:ymax]
    ENDIF ELSE BEGIN
        half_diff=(diff_in_size-1)/2
        xmin=half_diff+1
        xmax=(newsize - 1)-half_diff
        ymin=half_diff+1
        ymax=(newsize - 1)-half_diff
        bim=bim[xmin:xmax,ymin:ymax]
    ENDELSE
    
    if keyword_set(roi) then outsize=sizetmp  ; restore the original params
    
    ; now fix b to get the Sun in the center...
    axcntr=ahdr.crpix1
    aycntr=ahdr.crpix2
    bxcntr=bhdr.crpix1
    bycntr=bhdr.crpix2
    xdiff=(axcntr-bxcntr) ; we use '0-' because that prepares the values for "shift"
    ydiff=(aycntr-bycntr) ; we use '0-' because that prepares the values for "shift"
    bim=shift(bim,xdiff,ydiff)
    
    ; now bytescale.....
    aim=scc_bytscl(aim,ahdr)
    bim=scc_bytscl(bim,bhdr)
     
    ;stop
    
    IF keyword_set(roi) THEN BEGIN
        IF i EQ 0 THEN BEGIN ; only do this once
            window,retain=2,xsize=1024,ysize=1024
            tmpim=rebin(aim,1024,1024)
            tvscl,tmpim
            PRINT,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PRINT,'$$$$$$   -- PLEASE DEFINE YOUR REGION OF INTEREST. --  $$$$$$
            PRINT,'###### Note that you do not have to be exact -- the software
            PRINT,'###### will use the min/max of your ROI to form a rectangular
            PRINT,'###### ROI. If you want an irregular ROI... well,umm... 
            PRINT,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            junk=defroi(1024,1024,xverts,yverts)     
            roixmin=min(xverts)*2
            roixmax=max(xverts)*2
            roiymin=min(yverts)*2
            roiymax=max(yverts)*2
            if resize_flag then ywidth = oversized else ywidth = outsize
            xwidth=floor( (ywidth/(1.0*roiymax-roiymin)) * (roixmax-roixmin) )
        ENDIF
        aim=aim[roixmin:roixmax,roiymin:roiymax]
        bim=bim[floor(roixmin*1.15):floor(roixmax+(roixmin*0.15)),roiymin:roiymax]
    ENDIF  
    ; write out the pairs
    IF keyword_set(nrlviewer) THEN BEGIN
        fna=strmid(ahdr.filename,0,16)+strcompress(string(ahdr.wavelnth),/remove_all)+strmid(ahdr.filename,18,2)+'_R.png'
        fnb=strmid(bhdr.filename,0,16)+strcompress(string(bhdr.wavelnth),/remove_all)+strmid(bhdr.filename,18,2)+'_L.png'
    ENDIF ELSE BEGIN
        fna=strmid(ahdr.filename,0,16)+strcompress(string(ahdr.wavelnth),/remove_all)+strmid(ahdr.filename,18,4)+'png'
        fnb=strmid(bhdr.filename,0,16)+strcompress(string(bhdr.wavelnth),/remove_all)+strmid(bhdr.filename,18,4)+'png'
    ENDELSE
        
    IF (resize_flag)  OR (keyword_set(roi)) THEN BEGIN
        ; Resize data
        PRINT,''
        PRINT,'Resizing data arrays now...'
        PRINT,''
        IF keyword_set(roi) THEN BEGIN
            PRINT,'ROI will be resized so that the height is equal to the specified image size...'
            aim=congrid(aim,xwidth,ywidth)
            bim=congrid(bim,xwidth,ywidth)
        ENDIF ELSE BEGIN
            aim=congrid(aim,oversized,oversized)
            bim=congrid(bim,oversized,oversized)
        ENDELSE
    ENDIF
   ; stop
    PRINT,'Writing pngs...'
    
    aim=scc_add_datetime(aim,ahdr)
    aim=scc_add_logo(aim,ahdr)
    bim=scc_add_datetime(bim,ahdr)   ; NOTE: I am INTENTIONALLY passing the ahdr to the b-image
    bim=scc_add_logo(bim,ahdr)       ; NOTE: I am INTENTIONALLY passing the ahdr to the b-image
    
    tvlct,r,g,b,/get
   ; stop
    IF keyword_set(anaglyph) THEN BEGIN
        ;stop
        imanag=mkanaglyph(aim,bim)
        imanag[where(imanag LE 75)]=75 ; bit of a hack to remove compression noise from images
        tvlct,r,g,b,/get
        fn=strmid(ahdr.filename,0,16)+strcompress(string(ahdr.wavelnth),/remove_all)+strmid(ahdr.filename,18,2)+'.png'
        write_png,outdir+fn,imanag,r,g,b
    ENDIF ELSE BEGIN
        write_png,outdir+fna,aim,r,g,b  
        write_png,outdir+fnb,bim,r,g,b  
    ENDELSE
        
    skip:

ENDFOR

END
