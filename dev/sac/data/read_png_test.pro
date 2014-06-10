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

DEVICE, true_color=24, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW


dr='/data/ap1vf/png/3D/tube/P120_R100_A200_B1000_Dxy2_Dz1.6_Nxy100_Nz196/normalised/
traj=dr+'/ends_moves_three_full_vect/traj_three_shape/'

red = indgen(256)  
green = indgen(256)  
blue = indgen(256)  


names = [traj+'00100.png',traj+'00100.png']  
FOR i=0,N_ELEMENTS(names)-1 DO BEGIN  
         img = READ_PNG(names[i])  
         HELP,img 

window=i

    siz=size(img)
 stop   
    window,window,xsize=siz(1),ysize=siz(2),XPOS=i*siz(1), YPOS = 400

    tvframe, img
ENDFOR  


  
END




