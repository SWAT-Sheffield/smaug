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

rho=interpol(rho,zzz,zax)
pre=interpol(pre,zzz,zax)
tem=interpol(tem,zzz,zax)
pgt=interpol(pgt,zzz,zax)
pgas=interpol(pgas,zzz,zax)
babs=interpol(babs,zzz,zax)

zax=zax*1d5

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
0.0,dz*(k+1),tem(k),pgas(k),0.0,0.0,rho(k),babs(k),2d5,45000.0,theta,chi 
    
endfor

close,filemc 

spawn,' ./modcon < ./modcon_5000.inp > /dev/null'

spawn,'./inv_fei' ; > /dev/null'

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
    
openps,'linecomp.ps'    
    
plot,wlax,sti

restore,filename='obs6300.sav'

oplot,wl1,sp1-(max(sp1)-1.0), color=3

wlax2bl=wlax[200:399]

restore,filename='~/Work/VIKTOR/RADI/spd_2d_test.sav'

oplot,wlax2bl,spd(1,*,0),color=5

restore,filename='~/Work/VIKTOR/RADI/spd_2d_test_mu155.sav'

;oplot,wlax,spd(1,*,0),color=6

restore,filename='~/Work/VIKTOR/RADI/spd_test_6302.sav'

oplot,wlax2bl,spd(1,*,0,0),color=6

closeps

end






