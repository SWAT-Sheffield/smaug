pro rhoptoe,rho,pre,tem,ene,n1

close,1
openr,1,'./ttbl.dat'
readf,1,np,nr
pax=dblarr(np)
rax=dblarr(nr)
for i=0,np-1 do begin
 readf,1,paxv
 pax(i)=paxv
endfor
for i=0,nr-1 do begin
 readf,1,raxv
 rax(i)=raxv
endfor

ttbl=dblarr(np,nr)
etbl=dblarr(np,nr)


for i=0,np-1 do begin
 for j=0,nr-1 do begin
   readf,1,temv,enev
   ttbl(i,j)=temv
   etbl(i,j)=enev
 endfor
endfor


for i=0,n1-1 do begin


ip=interpol(indgen(np),alog10(pax),alog10(pre(i)))
ir=interpol(indgen(nr),alog10(rax),alog10(rho(i)))

tem(i)=bilinear(ttbl,ip,ir)
ene(i)=bilinear(etbl,ip,ir)

endfor

;stop


end


; z rho pre T pGas/pTot

npt=52

zzz=dblarr(npt)
rho=dblarr(npt)
pre=dblarr(npt)
tem=dblarr(npt)
pgt=dblarr(npt)

close,1
openr,1,'VALIIIC.dat'
for i=0,npt-1 do begin
 readf,1,zzzi,rhoi,prei,temi,pgti
 
 zzz(i)=zzzi
 rho(i)=rhoi
 pre(i)=prei
 tem(i)=temi
 pgt(i)=pgti
 
endfor

pgas=pre*pgt

babs=sqrt((pre-pgas)*8.d0*!pi)

nz=200

zax=dindgen(nz)*(max(zzz)-min(zzz))/(nz*1.d0)+min(zzz)

rho=interpol(rho,zzz,zax) ;/lsquadratic, /quadratic, /spline
pre=interpol(pre,zzz,zax)
tem=interpol(tem,zzz,zax)
pgt=interpol(pgt,zzz,zax)
pgas=interpol(pgas,zzz,zax)
babs=interpol(babs,zzz,zax)

zax=zax*1d5

dz=(zax(1)-zax(0))

ggg=27400.d0

p_1=dblarr(nz)


p_1(nz-1)=pgas(nz-1)

for i=nz-2,0,-1 do p_1(i)=p_1(i+1)+0.5d0*(rho(i+1)+rho(i))*ggg*(zax(i+1)-zax(i))


r_1=-deriv(zax,p_1)/ggg

r_2=smooth(r_1,10,/edge_truncate)

p_2=dblarr(nz)

p_2(nz-1)=pgas(nz-1)

for i=nz-2,0,-1 do p_2(i)=p_1(i+1)+0.5d0*(r_2(i)+r_2(i+1))*ggg*(zax(i+1)-zax(i))

r_2=-deriv(zax,p_2)/ggg



r_3=smooth(r_2,10,/edge_truncate)

p_3=dblarr(nz)

p_3(nz-1)=pgas(nz-1)

for i=nz-2,0,-1 do p_3(i)=p_3(i+1)+0.5d0*(r_3(i)+r_3(i+1))*ggg*(zax(i+1)-zax(i))

r_3=-deriv(zax,p_3)/ggg


t_3=dblarr(nz)
e_3=dblarr(nz)

t_4=dblarr(nz)
e_4=dblarr(nz)

rhoptoe,rho,pgas,t_3,e_3,nz

rhoptoe,r_3,p_3,t_4,e_4,nz

tek_color

plot,zax,t_4,/ylog
oplot,zax,t_3,color=2
oplot,zax,tem,color=3



end
