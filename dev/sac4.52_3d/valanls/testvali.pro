
function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
;for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2)
res(1)=res(2)
res(nel-1)=res(nel-3)
res(nel-2)=res(nel-3)
return,res
end

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

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

; z rho pre T pGas/pTot

restore, 'Temp_VALIIIC_old.sav'
restore, 'z_VALIIIC_old.sav'

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

zax=zax*1.0d5

dz=(zax(1)-zax(0))

ggg=27400.d0

p_1=dblarr(nz)


p_1(nz-1)=pgas(nz-1)

for i=nz-2,0,-1 do p_1(i)=p_1(i+1)+0.5d0*(rho(i+1)+rho(i))*ggg*(zax(i+1)-zax(i))

for i=0,30 do p_1 = SMOOTH( p_1, 2)

r_1=-deriv1(p_1,zax)/ggg


p_2=dblarr(nz)

p_2(nz-1)=pgas(nz-1)

for i=nz-2,0,-1 do p_2(i)=p_2(i+1)+0.5d0*(r_1(i)+r_1(i+1))*ggg*(zax(i+1)-zax(i))

for i=0,10 do p_2 = SMOOTH( p_2,2);,/edge_truncate)
r_2=-deriv1(p_2,zax)/ggg


r_3=r_2
p_3=p_2

;r_3=smooth(r_2,10,/edge_truncate)

;p_3=dblarr(nz)

;p_3(nz-1)=pgas(nz-1)

;for i=nz-2,0,-1 do p_3(i)=p_3(i+1)+0.5d0*(r_3(i)+r_3(i+1))*ggg*(zax(i+1)-zax(i))

;r_3=-deriv(zax,p_3)/ggg


t_3=dblarr(nz)
e_3=dblarr(nz)

t_4=dblarr(nz)
e_4=dblarr(nz)

rhoptoe,rho,pgas,t_3,e_3,nz

rhoptoe,r_3,p_3,t_4,e_4,nz
filename='testres/'
;filename='../../VAC_3D/vac4.52/mag_field/'
save, filename=filename+'z_rho_p_e_t_3.sav', zax,rho,pgas,e_3,t_3

tek_color

window,0

plot,zax,t_4,/ylog
oplot,zax,t_3,color=2
oplot,zax,tem,color=3
oplot,zz*100.d0,T(*,1),color=4

set_plot, 'ps' 

file_name='temp.eps'
device, filename=file_name, /ENCAPSULATED, /color

plot,zax,t_4,/ylog
oplot,zax,t_3,color=2
oplot,zax,tem,color=3

device, /close
set_plot, 'x'

;t_4, r_3, p_3




dz=(zax(1)-zax(0))



filemc=2


close,filemc
openw,filemc,'./atm4mc.dat'

dnz1=5
dnz2=0

theta=0.0
chi=0.0


printf,filemc,nz-dnz1-dnz2,' 1.'
    
for k=nz-1-dnz1,dnz2,-1 do begin
    
printf,filemc,format= $
'(F7.4,E14.5,F9.1,4E14.3,F10.2,2E14.3,2F10.3)',$
0.0,dz*(k+1),t_4(k),p_3(k),0.0,0.0,r_3(k),babs(k),2.0d5,45000.0,theta,chi 
    
endfor

close,filemc 

spawn,' ./modcon < ./modcon_5000.inp > /dev/null'

spawn,'./inv_fei' ; > /dev/null'; Fe 6301, 6302

filesp=3

wlax=dblarr(1)
sti=dblarr(1)
stv=dblarr(1)
stq=dblarr(1)
stu=dblarr(1)
ico=dblarr(1)
    
    
openr,filesp,'./sto4li.dat'
    
readf,filesp,nbl
for nbli=1,nbl do begin
readf,filesp,dum,wlct,nptwt,ict
ico=[ico,ict]
for nwl=1,nptwt do begin
readf,filesp,wlt,stit,stvt,stqt,stut       
  wlax=[wlax,wlt+wlct]
  sti=[sti,stit]
  stv=[stv,stvt]
  stq=[stq,stqt]
  stu=[stu,stut]            
endfor    
endfor

close,filesp
nel=n_elements(wlax)-1
wlax=wlax[1:nel]
sti=sti[1:nel]
stv=stv[1:nel]
stq=stq[1:nel]
stu=stu[1:nel]
ico=ico[1:n_elements(ico)-1]
    
tek_color    

 
window,1    
;openps,'linecomp.ps'    
    
plot,wlax,sti

restore,filename='obs6300.sav'

oplot,wl1,sp1-(max(sp1)-1.0), color=3

set_plot, 'ps' 

file_name='prof.eps'
device, filename=file_name, /ENCAPSULATED, /color

plot,wlax,sti
oplot,wl1,sp1-(max(sp1)-1.0), color=3

device, /close
set_plot, 'x'


end

