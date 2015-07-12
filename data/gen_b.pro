function deriv1,x,f

deriv1=dblarr(n_elements(f))
for i=1,n_elements(f)-2 do deriv1[i]=(f(i+1)-f(i-1))/(double(2.0)*(x(i+1)-x(i-1)))
deriv1[0]=deriv1[1]
deriv1[n_elements(f)-1]=deriv1[n_elements(f)-2]
return,deriv1
end

function b0z,z
b0z=1000.0-800.0*z/1.16861e8  ;zmax
return,b0z
end

function fff,xi
fff=exp(-(xi^(double(2.0)))/(0.5e10^2.0))
if (abs(fff) lt 1e-5) then fff=0.0
return,fff
end

function bbz,x,z
bbz=fff((x-double(1.25e8))*b0z(z))*b0z(z)
return,bbz
end

nptx=100
nptz=100
dz=1168614.4 ;cm
dx=2500000.0 ;cm

bz=dblarr(nptx,nptz)
bx=dblarr(nptx,nptz)

b0zz=dblarr(nptz)
zax=dindgen(nptz)*dz
xax=dindgen(nptx)*dx





for i=0,nptz-1 do b0zz(i)=b0z(i*dx)

b0zd=deriv1(zax,b0zz)

for i=0,nptx-1 do begin
 for j=0,nptz-1 do begin
  x=i*dx
  z=j*dz


  bz(i,j)=bbz(x,z)
  bx(i,j)=-b0zd(j)*bbz(x,z)/b0z(z)*(x-double(1.25e8))

 endfor
endfor

close,10
openw,10,'bbb_gen.dat',/f77_unformatted
for i=0,nptx-1 do begin
 for j=0,nptz-1 do begin
   writeu,10,bz(i,j),bx(i,j)
 endfor
endfor
close,10




dbxdx=dblarr(nptx,nptz)
dbxdz=dblarr(nptx,nptz)
dbzdx=dblarr(nptx,nptz)
dbzdz=dblarr(nptx,nptz)

for i=0,nptx-1 do dbxdx(*,i)=deriv1(xax,reform(bx(*,i)))
for i=0,nptz-1 do dbzdz(i,*)=deriv1(zax,reform(bz(i,*)))

divb=dbxdx+dbzdz

bz=rotate(bz,2)
bx=rotate(bx,2)

loadct,4
tvframe,bz,/bar,xrange=[0,max(xax)],yrange=[0,max(zax)]

;stop

;bz=bz+1.0
;bx=bx+1.0

bcstep=200.0/(bz(*,0)+10.0)
bcstep1=bcstep
bcax=findgen(100)
for i=3,99 do bcstep1(i)=int_tabulated(bcax[0:i],bcstep[0:i])
bcstep1=100.0*bcstep1/max(bcstep1)

bcstep1=findgen(100)
tek_color

kkk=4

for i=0,nptx/kkk-1 do begin

ddz=dz
z=0.0
x=bcstep1(i*kkk)*dx

print,x
if (bz(x/dx,z/dz) ne 0) then begin
repeat begin
ddx=ddz*bx(x/dx,z/dz)/bz(x/dx,z/dz)

oplot,[x,x+ddx],[z,z+ddz];,psym=10

x=x+ddx
z=z+ddz




print,x/dx,z/dz

endrep until not((x ge 0) and (x lt dx*nptx) and (z ge 0) and (z lt dz*nptz))
endif

endfor



end
