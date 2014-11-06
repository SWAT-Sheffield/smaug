pro dist2d_pict, npict
default2D,npx,npy,np,xd,yd,xpmin,ypmin,zpmin,xpmax,$
        ypmax,zpmax,xmin,ymin,xpSmin,ypSmin,xpSmax,ypSmax

s_filename = strarr(np) 

fname='/data/ap1vf/goriz1024100sm_np0108_0'
endfname='.out'


filename='/data/ap1vf/goriz1024100sm'+endfname

for nn=0,np-1 do begin
indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0'+indexs                                             
 2:indexss=''+indexs                                              
endcase   

s_filename(nn)=fname+indexss+endfname
print, s_filename(nn) 

endfor



step_save=1
ii=step_save

i_begin=0
unit_all=100
openfilew,unit_all,filename,'binary'


unit_b=10
for unit=unit_b,unit_b+np-1 do begin
 openfile,unit,s_filename(unit-unit_b),'binary'
endfor
 
 
 

   for unit=unit_b,unit_b+np-1 do begin

    gotopict,unit,'binary',npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
    error

     gethead,unit,'binary',headline,physics,it,time,gencoord, $
             ndim,neqpar,nw,nx,eqpar,variables,varname,pictsize=pictsize


     getpict_bin,unit,ndim,nw,nx,x,w
 
     xpb=round(xd*(x(0,xpSmin,1)-xmin)+xpmin)
     xpe=round(xd*(x(0,xpSmax,1)-xmin)+xpmin)

     print,'###',xpSmin,xpSmax ,unit-10, ' ', it, ' ', time,' ', x(0,0,0),' ',x(0,0,1),xpb,xpe


if i_begin eq 0 then begin
    nxL=lonarr(ndim)
    nxL(1)=(xpSmax+1)*npx
    nxL(0)=nx(0)

    xL=dblarr(nxL(0),nxL(1),ndim)
    wL=dblarr(nxL(0),nxL(1),nw)
    wiL=dblarr(nxL(0),nxL(1))
   i_begin=1
endif	    
 
    
		  
i1=0

    for i=xpb,xpe do begin
      xL(*,i,*)=x(*,i1,*)
      wL(*,i,*)=w(*,i1,*)
      i1=i1+1
    endfor
    
   endfor
    
    
;*********** start save to file *************    

    while strlen(headline) lt 79 do begin
      headline=headline+' '
    endwhile

    while strlen(varname) lt 79 do begin
      varname=varname+' '
    endwhile
    
    puthead,unit_all,'binary',headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nxL,eqpar,varname
    writeu,unit_all,xL
    for iw=0,nw-1 do begin
            wiL(*,*)=wL(*,*,iw)      
            writeu,unit_all,wiL
    end

    ;*********** end save to file *************    
 
 
 
for unit=unit_b,unit_b+np-1 do begin 
 close,unit
endfor
 close, unit_all 

end
