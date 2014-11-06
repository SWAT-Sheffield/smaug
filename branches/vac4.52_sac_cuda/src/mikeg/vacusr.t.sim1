!##############################################################################
! module vacusr - sim1 ! setvac -d=22 -g=204,204 -p=hdadiab -u=sim1


INCLUDE:vacusr.gravity.t
!INCLUDE:vacnul.specialini.t
!INCLUDE:vacnul.specialbound.t
! INCLUDE:vacnul.specialsource.t
!INCLUDE:vacnul.specialio.t
INCLUDE:vacusr.viscosity.t
!INCLUDE:vacusr.diffusion.t
!=============================================================================
subroutine specialini(ix^L,w)

include 'vacdef.f'

integer:: ix^L
integer:: ix_1,ix_2,ix_3
double precision:: w(ixG^T,1:nw)

double precision:: rhoin,xcent^D,radius
double precision:: inirho,iniene
double precision:: onemor,inix,ddx
double precision:: p_1,p_2

integer:: iii_,iix_1,info,i,j
double precision:: pi,comi,eneu,sum,mode,bmax,l
character*79 atmfilename

double precision:: p1,p2,rho1,rho2,v1,v2,T1,T2, b1_1,b1_2,b2_1,b2_2
double precision:: vx,vz,bx,bz

!-----------------------------------------------------------------------------

p1=5.d0/12.d0/Pi
rho1=25.d0/36.d0/Pi

vx=1.d0
vz=1.d0

bx=1.d0/2.0/sqrt(Pi)
bz=1.d0/2.0/sqrt(Pi)

  w(ix^S,rho_)=0.d0
  w(ix^S,rhob_)=rho1
  
  w(ix^S,m2_)=-vx*rho1*sin(2.d0*Pi*x(ix^S,1)/x(ixmax1,10,1))
  
  w(ix^S,m1_)=vz*rho1*sin(2.d0*Pi*x(ix^S,2)/x(10,ixmax2,2))
  
  w(ix^S,b2_)=-bx*sin(2.d0*Pi*x(ix^S,1)/x(ixmax1,10,1))
  w(ix^S,b1_)=bz*sin(4.d0*Pi*x(ix^S,2)/x(10,ixmax2,2))
  
  w(ix^S,e_)=0.d0 
  w(ix^S,eb_)=p1/(eqpar(gamma_)-1.d0)


!  w(ix^S,e_)=half*(^C&w(ix^S,m^C_)**2.d0+)/(w(ix^S,rho_)+w(ix^S,rhob_))

  w(ix^S,eb_)=w(ix^S,eb_)+half*((^C&w(ix^S,b^C_)**2.d0+))+half*(^C&w(ix^S,m^C_)**2.d0+)/(w(ix^S,rho_)+w(ix^S,rhob_))


  w(ix^S,rho_)=w(ix^S,rhob_)
  w(ix^S,rhob_)=0.d0

  w(ix^S,e_)=w(ix^S,eb_)
  w(ix^S,eb_)=0.d0 

  w(ix^S,bg2_)=0.d0
  w(ix^S,bg1_)=0.d0
  
!  w(ix^S,bg2_)=w(ix^S,b2_)
!  w(ix^S,b2_)=0.d0           
!  w(ix^S,bg1_)=w(ix^S,b1_)
!  w(ix^S,b1_)=0.d0


return
end


!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)


include 'vacdef.f'

integer:: ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
double precision:: fdt,fdthalf2

double precision:: pre(ixG^T),tem(ixG^T),kapr(ixG^T),so(ixG^T),flux(ixG^T)
double precision:: tau(ixG^T),ine(ixG^T)

double precision:: preg(ixG^T),pret(ixG^T)

integer:: rix_1,i,j
double precision:: mol_0, rrr_

double precision:: fsokr,avgflux

integer:: iw,iiw,iix_1



!point rast source**********
double precision:: rad,sigma1,sigma2,Q0,qt0,xc1,xc2,zzs,qmin,qra
double precision:: rfc,tdep,sdep
double precision:: xs(100),zs(100),ts(100),qs(100)
double precision:: tlast,rn
integer:: ns
!logical:: filexist
integer:: singl
!***************

integer:: ix_1,ix_2,idim, ix^L

!*****************
double precision:: t01,t02,a1,a2,s1,s2,sf,xc1,xc2,rad,rfc,sdep,tdep,sigma2


!-----------------------------------------------------------------------------



eqpar(nu_)=1.d0
!eqpar(nu_)=0.d0



if(abs(eqpar(nu_))>smalldouble)&
   call addsource_visc(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
 
!call addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)


write(*,*) '***time=',qt



end



!=============================================================================
subroutine specialbound(qt,ix^L,iw,iB,w)
include 'vacdef.f'



integer:: iw^LIM,idim^LIM
double precision:: qt,w(ixG^T,1:nw)
integer:: ix,ix^D,ixe,ixf,ix^L,ixpair^L,idim,iw,iB
integer:: iwv,jdim

integer:: Ns,i,j
double precision:: ki

integer:: ix_1,ix_2

double precision:: tmpp1,tmpp2



return
end

!=============================================================================
subroutine getdt_special(w,ix^L)

! If the Coriolis force is made very strong it may require time step limiting,
! but this is not implemented here.

include 'vacdef.f'
double precision:: w(ixG^T,nw)
integer:: ix^L
!-----------------------------------------------------------------------------

!call getdt_diff(w,ix^L)


if(abs(eqpar(nu_))>smalldouble)&
   call getdt_visc(w,ix^L)

!call getdt_res(w,ix^L)

!call getdt_grav(w,ix^L)

return
end


subroutine specialeta(w,ix^L,idirmin)
 
include 'vacdef.f'
 
double precision:: w(ixG^T,nw)
integer:: ix^L,idirmin
!-----------------------------------------------------------------------------
 
stop 'specialeta is not defined'

end

