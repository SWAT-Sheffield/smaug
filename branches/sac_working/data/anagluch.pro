;window,0,col=256
nx=100
ny=100

l_i=bytarr(nx,ny)
r_i=bytarr(nx,ny)

for i=0,nx-1 do begin
 for j=0, ny-1 do begin
   aa=abs(sin(i*!Pi/100.d0)*cos(j*!Pi/100.d0)*256.d0)
   l_i[i,j]=aa
   r_i[i,j]=aa
 endfor
endfor
 
;tv,l_i 


a_i_rgb=bytarr(nx,ny,3)
a_i_rgb(*,*,0)=r_i(*,*)
a_i_rgb(*,*,1)=l_i(*,*)
a_i_rgb(*,*,2)=l_i(*,*)
stop
a_i=color_quan(a_i_rgb,3,r,g,b,colors=256,get_translation=ctrans)
tvlct,r,g,b
tv,a_i ;, true=3
end
