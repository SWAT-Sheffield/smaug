For use with vac4.52
Test 

setvac -d=23 -phi=0 -z=0 -g=260,260 -p=mhd -u=nul \
       -on=cd,tvdlf,tvd,poisson,resist -off=mc,fct,impl,ct,gencoord,rk,mpi


setvac -d=23 -phi=0 -z=0 -g=260,260 -p=hdadiab -u=example \
       -on=cd,poisson,rk -off=mc,fct,tvdlf,tvd,impl,ct,gencoord,resist,mpi
