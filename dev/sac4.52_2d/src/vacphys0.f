!##############################################################################
! module vacphys0 - mhd

!=============================================================================
subroutine conserve(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Transform primitive variables into conservative ones

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------


! Calculate total energy from pressure, kinetic and magnetic energy

w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=w(ixmin1:ixmax1,ixmin2:ixmax2,p_)&
   /(eqpar(gamma_)-1)+half*((w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))*(w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)&
   **2+w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)**2)+((w(ixmin1:ixmax1,ixmin2:ixmax2,&
   b1_))**2+(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_))**2))+((w(ixmin1:ixmax1,&
   ixmin2:ixmax2,b1_)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)))


! Convert velocity to momentum
w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)=(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)
w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)=(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,v2_);


return
end

!=============================================================================
subroutine primitive(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Transform conservative variables into primitive ones

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------


! Calculate pressure

call getpthermal(w,ixmin1,ixmin2,ixmax1,ixmax2,tmp)

w(ixmin1:ixmax1,ixmin2:ixmax2,p_)=tmp(ixmin1:ixmax1,ixmin2:ixmax2)

! Convert momentum to velocity
w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)=w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)&
   /(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   rhob_))
w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)=w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)&
   /(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   rhob_));

return
end

!=============================================================================
subroutine getv(w,ixmin1,ixmin2,ixmax1,ixmax2,idim,v)

! Calculate v_idim=m_idim/rho within ix

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),v(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getv')>=1
if(oktest)write(*,*)'GetV w:',w(ixtest1,ixtest2,iwtest)

v(ixmin1:ixmax1,ixmin2:ixmax2)=w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
   +idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   rhob_))

if(oktest)write(*,*)'GetV v:',v(ixtest1,ixtest2)

return 
end


!=============================================================================
subroutine getcmax(new_cmax,w,ixmin1,ixmin2,ixmax1,ixmax2,idim,cmax)

! Calculate cmax_idim=cfast_i+abs(v_idim) within ix^L
! where cfast_i=sqrt(0.5*(cf**2+sqrt(cf**4-4*cs**2*b_i**2/rho)))
! and cf**2=b**2/rho+cs**2/rho is the square of the speed of the fast wave 
! perpendicular to the magnetic field, and cs is the sound speed.

include 'vacdef.f'

logical:: new_cmax
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),cmax(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision:: csound2(ixGlo1:ixGhi1,ixGlo2:ixGhi2),cfast2(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
save csound2,cfast2
!-----------------------------------------------------------------------------
oktest=index(teststr,'getcmax')>=1

!Direction independent part of getcmax:
if(new_cmax)then
   new_cmax=.false.
   call getcsound2(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)
   if(oktest)write(*,*)'csound2:',csound2(ixtest1,ixtest2)
   cfast2(ixmin1:ixmax1,ixmin2:ixmax2)=((w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)&
      +w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))**2+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
      b2_)+w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_))**2 )/(w(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))&
      +csound2(ixmin1:ixmax1,ixmin2:ixmax2)
end if
if(oktest)write(*,*)'cfast2:',cfast2(ixtest1,ixtest2)

cmax(ixmin1:ixmax1,ixmin2:ixmax2)=sqrt(half*(cfast2(ixmin1:ixmax1,&
   ixmin2:ixmax2)+ sqrt(cfast2(ixmin1:ixmax1,ixmin2:ixmax2)**2&
   -4*csound2(ixmin1:ixmax1,ixmin2:ixmax2)* ((w(ixmin1:ixmax1,ixmin2:ixmax2,&
   b0_+idim)+w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim))**2)/(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))))) &
   +abs(w(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)))

if(oktest) write(*,*)'cmax:',cmax(ixtest1,ixtest2)


return 
end

!=============================================================================
subroutine getcsound2prim(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L
! from the primitive variables in w.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),csound2(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct Getcsound2prim for NONIDEAL gas in vacphys.t.mhd')

csound2(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(gamma_)*(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,p_)+(eqpar(gamma_)-one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)&
   -half*( (w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))**2+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,bg2_))**2 )))/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))

return 
end

!=============================================================================
subroutine getcsound2(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),csound2(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct Getcsound2 for NONIDEAL gas in vacphys.t.mhd')

oktest=index(teststr,'getcsound2')>=1
if(oktest) write(*,*)'Getcsound2'

call getpthermal(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)
if(oktest) write(*,*)'p(ixtest)=',csound2(ixtest1,ixtest2)
csound2(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(gamma_)*(csound2(ixmin1:ixmax1,&
   ixmin2:ixmax2)+(eqpar(gamma_)-one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)&
   -half*( (w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))**2+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,bg2_))**2 )))/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))

return 
end

!=============================================================================
subroutine getpthermal(w,ixmin1,ixmin2,ixmax1,ixmax2,p)

!!! This subroutine should not use tmp,tmp2


include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),p(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------


p(ixmin1:ixmax1,ixmin2:ixmax2)=half*( w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)**2&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)**2 )/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))

p(ixmin1:ixmax1,ixmin2:ixmax2)=p(ixmin1:ixmax1,ixmin2:ixmax2)&
   + half*( (w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2)+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,b2_)**2) )+( (w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)&
   *w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)&
   *w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)) )

p(ixmin1:ixmax1,ixmin2:ixmax2)=(eqpar(gamma_)-one)*(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,e_)-p(ixmin1:ixmax1,ixmin2:ixmax2))


return 
end

!=============================================================================
subroutine getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,p)

include 'vacdef.f'

double precision::  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),p(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),gamma
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

gamma=eqpar(gamma_)

p(ixmin1:ixmax1,ixmin2:ixmax2)=(gamma-two)*(( (w(ixmin1:ixmax1,ixmin2:ixmax2,&
   b1_)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)) )+ half*( (w(ixmin1:ixmax1,&
   ixmin2:ixmax2,b1_))**2.d0+(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_))**2.d0 ))

p(ixmin1:ixmax1,ixmin2:ixmax2)=(gamma-one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,e_)&
   -half*( w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)**2.d0+w(ixmin1:ixmax1,&
   ixmin2:ixmax2,m2_)**2.d0 )/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)))-p(ixmin1:ixmax1,ixmin2:ixmax2)


return 
end

!=============================================================================
subroutine getptotal_bg(w,ixmin1,ixmin2,ixmax1,ixmax2,p)

include 'vacdef.f'

double precision::  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),p(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),gamma
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

gamma=eqpar(gamma_)

p(ixmin1:ixmax1,ixmin2:ixmax2)=(eqpar(gamma_)-one)*w(ixmin1:ixmax1,&
   ixmin2:ixmax2,eb_)-half*(eqpar(gamma_)-two)*( (w(ixmin1:ixmax1,&
   ixmin2:ixmax2,bg1_)**2.d0)+(w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)&
   **2.d0) )    

return 
end

