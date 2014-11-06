*##############################################################################
* module vacphys0 - mhd

*=============================================================================
      subroutine conserve(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Transform primitive variables into conservative ones

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'conserv').ge.1
      if(oktest)write(*,*)'Conserve w:',w(ixtest1,ixtest2,iwtest)

* Calculate total energy from pressure, kinetic and magnetic energy
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,e_)=w(ix_1,ix_2,p_)/(eqpar(gamma_)-1)+half*
     &      (w(ix_1,ix_2,rho_)*(w(ix_1,ix_2,v1_)**2+w(ix_1,ix_2,v2_)**
     &      2)+(w(ix_1,ix_2,b1_)**2+w(ix_1,ix_2,b2_)**2))
      enddo
      enddo

* Convert velocity to momentum
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,m1_)=w(ix_1,ix_2,rho_)*w(ix_1,ix_2,v1_)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,m2_)=w(ix_1,ix_2,rho_)*w(ix_1,ix_2,v2_)
         
      enddo
      enddo

      if(oktest)write(*,*)'New w:',w(ixtest1,ixtest2,iwtest)

      return
      end

*=============================================================================
      subroutine primitive(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Transform conservative variables into primitive ones

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'primitive').ge.1
      if(oktest)write(*,*)'Primitive w:',w(ixtest1,ixtest2,iwtest)

* Calculate pressure
      call getpthermal(.true.,w,ixmin1,ixmin2,ixmax1,ixmax2,tmp)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,p_)=tmp(ix_1,ix_2)
      enddo
      enddo

* Convert momentum to velocity
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,v1_)=w(ix_1,ix_2,m1_)/w(ix_1,ix_2,rho_)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,v2_)=w(ix_1,ix_2,m2_)/w(ix_1,ix_2,rho_)
         
      enddo
      enddo

      if(oktest)write(*,*)'New w:',w(ixtest1,ixtest2,iwtest)

      return
      end

*=============================================================================
      subroutine getv(w,ixmin1,ixmin2,ixmax1,ixmax2,idim,v)

* Calculate v_idim=m_idim/rho within ix

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),v(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getv').ge.1
      if(oktest)write(*,*)'GetV w:',w(ixtest1,ixtest2,iwtest)

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         v(ix_1,ix_2)=w(ix_1,ix_2,m0_+idim)/w(ix_1,ix_2,rho_)
      enddo
      enddo

      if(oktest)write(*,*)'GetV v:',v(ixtest1,ixtest2)

      return 
      end


*=============================================================================
      subroutine getcmax(new_cmax,w,ixmin1,ixmin2,ixmax1,ixmax2,idim,cmax)

* Calculate cmax_idim=cfast_i+abs(v_idim) within ix^L
* where cfast_i=sqrt(0.5*(cf**2+sqrt(cf**4-4*cs**2*b_i**2/rho)))
* and cf**2=b**2/rho+cs**2/rho is the square of the speed of the fast wave 
* perpendicular to the magnetic field, and cs is the sound speed.

      include 'vacdef.f'

      logical  new_cmax
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),cmax(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      double precision  csound2(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   cfast2(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      save csound2,cfast2
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getcmax').ge.1

*Direction independent part of getcmax:
      if(new_cmax)then
         new_cmax=.false.
         call getcsound2(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)
         if(oktest)write(*,*)'csound2:',csound2(ixtest1,ixtest2)
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            cfast2(ix_1,ix_2)=( w(ix_1,ix_2,b1_)**2+w(ix_1,ix_2,b2_)**
     &         2 )/w(ix_1,ix_2,rho_)+csound2(ix_1,ix_2)
         enddo
         enddo
      end if
      if(oktest)write(*,*)'cfast2:',cfast2(ixtest1,ixtest2)
*if(oktest.and.index(teststr,'discriminant')>=1)write(*,*)'Discriminant:',&
*   cfast2(ix^S)**2-4*csound2(ix^S)*(w(ix^S,b0_+idim)**2)/w(ix^S,rho_)

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cmax(ix_1,ix_2)=sqrt(half*(cfast2(ix_1,ix_2)+ sqrt(cfast2(ix_1,ix_2)
     &      **2-4*csound2(ix_1,ix_2)* (w(ix_1,ix_2,b0_+idim)**
     &      2)/w(ix_1,ix_2,rho_)))) +abs(w(ix_1,ix_2,m0_+idim)/
     &      w(ix_1,ix_2,rho_))
      enddo
      enddo

      if(oktest) write(*,*)'cmax:',cmax(ixtest1,ixtest2)

      return 
      end

*=============================================================================
      subroutine getcsound2prim(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)

* Calculate the square of the thermal sound speed csound2 within ix^L
* from the primitive variables in w.
* csound2=gamma*p/rho

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),csound2(
     &   ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ixmin1,ixmin2,ixmax1,ixmax2
*-----------------------------------------------------------------------------

      if(eqpar(gamma_).le.zero)call die(
     &   'Correct Getcsound2prim for NONIDEAL gas in vacphys.t.mhd')

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         csound2(ix_1,ix_2)=eqpar(gamma_)*w(ix_1,ix_2,p_)/w(ix_1,ix_2,rho_)
      enddo
      enddo

      return 
      end

*=============================================================================
      subroutine getcsound2(w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)

* Calculate the square of the thermal sound speed csound2 within ix^L.
* csound2=gamma*p/rho

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),csound2(
     &   ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ixmin1,ixmin2,ixmax1,ixmax2
*-----------------------------------------------------------------------------

      if(eqpar(gamma_).le.zero)call die(
     &   'Correct Getcsound2 for NONIDEAL gas in vacphys.t.mhd')

      oktest=index(teststr,'getcsound2').ge.1
      if(oktest) write(*,*)'Getcsound2'

      call getpthermal(.true.,w,ixmin1,ixmin2,ixmax1,ixmax2,csound2)
      if(oktest) write(*,*)'p(ixtest)=',csound2(ixtest1,ixtest2)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         csound2(ix_1,ix_2)=eqpar(gamma_)*csound2(ix_1,ix_2)/w(ix_1,ix_2,rho_)
      enddo
      enddo

      return 
      end

*=============================================================================
      subroutine getpthermal(clipping,w,ixmin1,ixmin2,ixmax1,ixmax2,p)

*!! This subroutine should not use tmp,tmp2

* Calculate thermal pressure within ix^L
* p=(g-1)*(e-0.5*(m**2/r+b**2))
* where g is the adiabatic index gamma, e is the total energy density,
* m is the momentum density and b is the magnetic field
* If clipping is .true. clip off negative pressures

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),p(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  ixmin1,ixmin2,ixmax1,ixmax2
      logical  clipping
*-----------------------------------------------------------------------------

      if(eqpar(gamma_).le.zero)call die(
     &   'Correct GetPthermal for NONIDEAL gas in vacphys.t.mhd')

      oktest=index(teststr,'getpthermal').ge.1

      if(oktest) write(*,*)'GetPthermal'

* First calculate kinetic energy*2=m**2/rho

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         p(ix_1,ix_2)=( w(ix_1,ix_2,m1_)**2+w(ix_1,ix_2,m2_)**
     &      2 )/w(ix_1,ix_2,rho_)
      enddo
      enddo
      if(oktest) write(*,*)'p(ixtest)=',p(ixtest1,ixtest2)

* Add magnetic energy*2=b**2
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         p(ix_1,ix_2)=p(ix_1,ix_2)+ w(ix_1,ix_2,b1_)**2+w(ix_1,ix_2,b2_)**2
      enddo
      enddo
      if(oktest) write(*,*)'p(ixtest)=',p(ixtest1,ixtest2)

* Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         p(ix_1,ix_2)=(eqpar(gamma_)-1)*(w(ix_1,ix_2,e_)-half*p(ix_1,ix_2))
      enddo
      enddo
      if(oktest) write(*,*)'p(ixtest)=',p(ixtest1,ixtest2)

      if(smallp.lt.zero)then
         minval_1=bigdouble
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            minval_1=min(minval_1,p(ix_1,ix_2))
         enddo
         enddo
         smallp=minval_1*smallpcoeff
         
         if(oktest) write(*,*)'smallp, smallpcoeff:',smallp,smallpcoeff
         minval_1=bigdouble
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            minval_1=min(minval_1,p(ix_1,ix_2))
         enddo
         enddo
         if(oktest)then
            write(*,*)'minval(p):',minval_1
         endif
         if(oktest) write(*,*)'ix-limits:',ixmin1,ixmin2,ixmax1,ixmax2
         if(smallp.lt.zero)then
            write(*,*)'Initial condition contains negative pressure!'
            maxval_1=-bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               maxval_1=max(maxval_1,p(ix_1,ix_2))
            enddo
            enddo
            smallp=maxval_1*smallpcoeff
            
            write(*,*)'Using smallp=',smallp
         endif
      endif

* Clip off negative pressure if clipping is set
      if(clipping)then
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            p(ix_1,ix_2)=max(p(ix_1,ix_2),smallp)
         enddo
         enddo
      endif

*!!if(clipping)p(ix^S)=abs(p(ix^S)) Hans-s solution for negative pressure

      return 
      end

*=============================================================================
      subroutine getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,p)

* Calculate total pressure within ix^L including magnetic pressure
* p=(g-1)*e-0.5*(g-1)*m**2/rho+(1-0.5*g)*b**2

*!! This subroutine should not use tmp

      include 'vacdef.f'

      double precision   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),p(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),gamma
      integer  ixmin1,ixmin2,ixmax1,ixmax2
*-----------------------------------------------------------------------------

      if(eqpar(gamma_).le.zero)call die(
     &   'Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

      oktest=index(teststr,'getptotal').ge.1
      if(oktest) write(*,*)'GetPtotal'

      gamma=eqpar(gamma_)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         p(ix_1,ix_2)=(one-half*gamma)*( w(ix_1,ix_2,b1_)**
     &      2+w(ix_1,ix_2,b2_)**2 )
      enddo
      enddo

* Add contribution of total energy=(g-1)*e
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         p(ix_1,ix_2)=p(ix_1,ix_2)+(gamma-one)*w(ix_1,ix_2,e_)
      enddo
      enddo
      if(oktest) write(*,*)'p(ixtest)=',p(ixtest1,ixtest2)

* Subtract contribution of kinetic energy=0.5*(g-1)*m**2/rho
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         p(ix_1,ix_2)=p(ix_1,ix_2)-half*(gamma-one)*(w(ix_1,ix_2,m1_)**
     &      2+w(ix_1,ix_2,m2_)**2)/w(ix_1,ix_2,rho_)
      enddo
      enddo

      if(oktest) write(*,*)'p(ixtest)=',p(ixtest1,ixtest2)

      return 
      end

