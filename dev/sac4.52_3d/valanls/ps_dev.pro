;  =======================================================================
;
;  ps_dev.pro
;
;  unified SR to set the device befor writting the ps-file
;
;  =======================================================================


pro ps_dev, nam_ps_file, a_cpl, sty, phe,pwi,che,cwi, spn, ps_rat,$
            eps=eps

    ;    sty   p=portrait l=landscape
    ;    phe   % of image height
    ;    pwi   % if image width
    ;    che   % adjust height
    ;    cwi   % adjust widht 
    ;
    ;    spn   space (normal) on the left,right,lower and upper image boundary
    ;          i.e the lower left edge is (  -spn(0) ,  -spn(2) )
    ;                  upper right     is (1+ spn(1) , 1+spn(3) )
    ;    ps_rat   [im_wid,im_hei]
    ;
    ;    eps=eps  --> encapsulated ps-file


    hei=28.0
    wid=19.1

    nam=nam_ps_file

    n_enc = 0
    if keyword_set(eps) then begin
       ip1 = strpos(nam,'.')
       if ip1 eq -1 then begin
          nam = nam+'.eps'
       endif else begin
          nam = strmid(nam,0,ip1) +'.eps'
       endelse
       n_enc=1
    endif
    nam_ps_file = nam

    set_plot,'ps'
    
        if phe le 0. then phe=0.8
        if pwi le 0. then pwi=0.8
    if sty eq 'p' then begin
	im_hei  = phe*hei
        im_wid  = pwi*wid
        im_cen_hei = (hei - im_hei)*che
        im_cen_wid = (wid - im_wid)*cwi
    endif else begin
       	im_hei  = phe*wid
        im_wid  = pwi*hei
        im_cen_hei = (wid - im_hei)*che
        im_cen_wid = (hei - im_wid)*cwi
    endelse

    ps_rat=[im_wid,im_hei]

    if a_cpl eq 'y' then begin
       if sty eq 'l' then begin
           device,filename=nam,/color,bit=8,encapsulated=n_enc $	
            ,/landscape $
                 ,xoffset=     0.4 + im_cen_hei , xsize = im_wid $
	         ,yoffset= hei+0.8 - im_cen_wid , ysize = im_hei    

                  spn_le= (       im_cen_wid         ) / im_wid 
                  spn_ri= ( hei - im_cen_wid - im_wid) / im_wid
                  spn_do= (       im_cen_hei         ) / im_hei
                  spn_up= ( wid - im_cen_hei - im_hei) / im_hei

       endif else  if sty eq 'p' then begin
           device,filename=nam,/color,bit=8,encapsulated=n_enc $
                  ,/portrait $
                  ,xoffset = 0.8 + im_cen_wid , xsize = im_wid $ 
                  ,yoffset = 0.8 + im_cen_hei , ysize = im_hei

                  spn_le= (       im_cen_wid         ) / im_wid 
                  spn_ri= ( wid - im_cen_wid - im_wid) / im_wid
                  spn_do= (       im_cen_hei         ) / im_hei
                  spn_up= ( hei - im_cen_hei - im_hei) / im_hei

       endif else begin
           print,'Cannot choose device!'
           stop
       endelse
    endif else begin
       if sty eq 'l' then begin      
           device,filename=nam,bit=8,encapsulated=n_enc $	
  	         ,/landscape $
	         ,xoffset=     0.4 + im_cen_hei , xsize = im_wid $
	         ,yoffset= hei+0.8 - im_cen_wid , ysize = im_hei    

                  spn_le= (       im_cen_wid         ) / im_wid 
                  spn_ri= ( hei - im_cen_wid - im_wid) / im_wid
                  spn_do= (       im_cen_hei         ) / im_hei
                  spn_up= ( wid - im_cen_hei - im_hei) / im_hei

       endif else  if sty eq 'p' then begin 
           device,filename=nam, bit=8,encapsulated=n_enc  $	
                 ,/portrait $
                  ,xoffset = 0.8 + im_cen_wid , xsize = im_wid $ 
                  ,yoffset = 0.8 + im_cen_hei , ysize = im_hei

                  spn_le= (       im_cen_wid         ) / im_wid 
                  spn_ri= ( wid - im_cen_wid - im_wid) / im_wid
                  spn_do= (       im_cen_hei         ) / im_hei
                  spn_up= ( hei - im_cen_hei - im_hei) / im_hei
       endif else begin
           print,'Cannot choose device!'
           stop
       endelse
    endelse

    spn = [spn_le,spn_ri,spn_do,spn_up]

end





