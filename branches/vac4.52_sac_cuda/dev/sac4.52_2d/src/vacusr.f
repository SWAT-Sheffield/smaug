!##############################################################################
! module vacusr - sim1 ! setvac -d=22 -g=204,204 -p=hdadiab -u=sim1


!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD GRAVITATIONAL SOURCE TERMS, SET GRAVITY
!
!------------------------------------------------------------------------------
!    See vacusr.t.gravity and vacusrpar.t.gravity for an example of usage
!
!    Gravitational force is added to the momentum equation:
!
!    d m_i/dt += rho*eqpar(grav0_+i)
!
!    Gravitational work is added to the energy equation (if present):
!
!    de/dt += Sum_i m_i*eqpar(grav0_+i)
!
!    The eqpar(grav1_),eqpar(grav2_),... coefficients are the components of 
!    the gravitational acceleration in each dimension. Set them to 0 for no
!    gravity in that direction. 
!    The !!! comments show how a grav array could be used for a spatially
!    (and maybe temporally) varying gravitational field.
!    The setgrav subroutine has to be completed then.
!
!============================================================================
subroutine addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

! Add gravity source calculated from w to wnew within ixO for all variables 
! in iws. w is at time qtC, wnew is advanced from qt to qt+qdt.

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: iiw,iw,idim
!!! ! For a spatially varying gravity define the common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common /gravity/ grav

!-----------------------------------------------------------------------------

!!! ! If grav needs to be calculated only once do it for the whole grid
!!! if(it==itmin)call setgrav(w,ixG^L,ixG^L,grav)
!!! ! Otherwise call setgrav in every time step
!!! call setgrav(w,ixI^L,ixO^L,grav)

! add sources from gravity
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(m1_,m2_)
      ! dm_i/dt= +rho*g_i
      idim=iw-m0_
      if(abs(eqpar(grav0_+idim))>smalldouble) wnew(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,m0_+idim)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
         +idim)+ qdt*eqpar(grav0_+idim)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_))

!          wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
!              qdt*eqpar(grav0_+idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

      !!! ! For a spatially varying gravity use instead of the above lines
      !!! wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
      !!!    qdt*grav(ixO^S,idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

   case(e_)
      ! de/dt= +g_i*m_i
      do idim=1,ndim
         if(abs(eqpar(grav0_+idim))>smalldouble) wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ee_)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ee_)&
            + qdt*eqpar(grav0_+idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_+idim)/(w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rhob_))

!            wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
!               qdt*eqpar(grav0_+idim)*w(ixO^S,m0_+idim)

         !!! ! For a spatially varying gravity use instead of the above lines
         !!! wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
         !!!    qdt*grav(ixO^S,idim)*w(ixO^S,m0_+idim)

      end do
   end select ! iw
end do        ! iiw

return
end
!=============================================================================
!!! subroutine setgrav(w,ixI^L,ixO^L,grav)

! Set the gravitational acceleration within ixO based on x(ixI,ndim) 
! and/or w(ixI,nw)

!!! include 'vacdef.f'

!!! double precision:: w(ixG^T,nw),grav(ixG^T,ndim)
!!! integer:: ixI^L,ixO^L
!----------------------------------------------------------------------------
!!! return
!!! end
!=============================================================================

subroutine getdt_grav(w,ixmin1,ixmin2,ixmax1,ixmax2)

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: dtgrav
save dtgrav

!!! ! For spatially varying gravity you need a common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common/gravity/grav

!----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

if(it==itmin)then
   ! If gravity is descibed by the equation parameters, use this:
   dtgrav=bigdouble
   do idim=1,ndim
      if(abs(eqpar(grav0_+idim))>zero)dtgrav=min(dtgrav,one&
         /sqrt(maxval(abs(eqpar(grav0_+idim))/dx(ixMmin1:ixMmax1,&
         ixMmin2:ixMmax2,1:ndim))))
   enddo
   !!! ! For spatially varying gravity use this instead of the lines above:
   !!! call setgrav(w,ixG^L,ixM^L,grav)
   !!! ! If gravity does not change with time, calculate dtgrav here:
   !!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))
endif

!!! ! If gravity changes with time, calculate dtgrav here:
!!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))

      call mpiallreduce(dtgrav,MPI_MIN)

! limit the time step
dt=min(dt,dtgrav)
if(oktest)write(*,*)'Gravity limit for dt:',dtgrav

return
end

!=============================================================================

!==============================================================================
subroutine addsource_visc(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

! Add viscosity source to wnew within ixO 

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer:: ix,ixmin1,ixmin2,ixmax1,ixmax2,idim,idir,jdir,iiw,iw
double precision:: tmp2(ixGlo1:ixGhi1,ixGlo2:ixGhi2),nushk(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ndim)


double precision:: tmprhoL(ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmprhoR&
   (ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmprhoC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision:: tmpVL(ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmpVR(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2), tmpVC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision:: tmpBL(ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmpBR(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2), tmpBC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision:: tmpL(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tmpR(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2), tmpC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision:: nuL(ixGlo1:ixGhi1,ixGlo2:ixGhi2),nuR(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

integer:: jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,hxmin2,hxmax1,hxmax2, hxOmin1,&
   hxOmin2,hxOmax1,hxOmax2

double precision:: c_ene,c_shk

integer:: i,j,k,l,m,ii0,ii1,t00

double precision:: sB

!-----------------------------------------------------------------------------

! Calculating viscosity sources 
! involves second derivatives, two extra layers
call ensurebound(2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,qtC,w)
ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

!sehr wichtig
call setnushk(w,ixmin1,ixmin2,ixmax1,ixmax2,nushk)

do idim=1,ndim
      tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_)
      call setnu(w,rho_,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)      
      call gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
      tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
         ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)         
      call gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
      tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
         ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
      wnew(ixImin1:ixImax1,ixImin2:ixImax2,rho_)=wnew(ixImin1:ixImax1,&
         ixImin2:ixImax2,rho_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
         -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
         ixImin2:ixImax2,idim)*qdt
enddo
   

do idim=1,ndim
      tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)-half*((w(ixImin1:ixImax1,ixImin2:ixImax2,b1_)**2&
         +w(ixImin1:ixImax1,ixImin2:ixImax2,b2_)**2)+(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,m1_)**2+w(ixImin1:ixImax1,ixImin2:ixImax2,m2_)**2)&
         /(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)+w(ixImin1:ixImax1,&
         ixImin2:ixImax2,rhob_)))
      call setnu(w,173,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)      
      call gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
      tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
         ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)      
      call gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
      tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
         ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
      wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
         ixImin2:ixImax2,e_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
         -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
         ixImin2:ixImax2,idim)*qdt
enddo




      tmprhoC(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,rho_)+w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_)



do k=1,ndim
        jxmin1=ixmin1+kr(k,1);jxmin2=ixmin2+kr(k,2);jxmax1=ixmax1+kr(k,1)
        jxmax2=ixmax2+kr(k,2); 
        hxmin1=ixmin1-kr(k,1);hxmin2=ixmin2-kr(k,2);hxmax1=ixmax1-kr(k,1)
        hxmax2=ixmax2-kr(k,2);
        tmprhoL(ixmin1:ixmax1,ixmin2:ixmax2)=((w(ixmin1:ixmax1,ixmin2:ixmax2,&
           rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))+(w(hxmin1:hxmax1,&
           hxmin2:hxmax2,rho_)+w(hxmin1:hxmax1,hxmin2:hxmax2,rhob_)))/two
        tmprhoR(ixmin1:ixmax1,ixmin2:ixmax2)=((w(jxmin1:jxmax1,jxmin2:jxmax2,&
           rho_)+w(jxmin1:jxmax1,jxmin2:jxmax2,rhob_))+(w(ixmin1:ixmax1,&
           ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)))/two

   do l=1,ndim
        call setnu(w,l+m0_,k,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)      
        tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,m0_+l)/(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)&
           +w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_))


      do ii1=0,1
        if (ii1 .eq. 0) then
                           i=k
                           ii0=l
                        else
                           i=l
                           ii0=k
        endif



        if (i .eq. k) then 
           tmpVL(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,ixmin2:ixmax2,&
              m0_+ii0)+w(hxmin1:hxmax1,hxmin2:hxmax2,m0_+ii0))/two
           tmpVR(ixmin1:ixmax1,ixmin2:ixmax2)=(w(jxmin1:jxmax1,jxmin2:jxmax2,&
              m0_+ii0)+w(ixmin1:ixmax1,ixmin2:ixmax2,m0_+ii0))/two

           call gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
           tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
              ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,k))&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
           call gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
           tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
              ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,k))&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2) 

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=(tmprhoR(ixImin1:ixImax1,&
              ixImin2:ixImax2)*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
              -tmprhoL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpL(ixImin1:ixImax1,&
              ixImin2:ixImax2))/dx(ixImin1:ixImax1,ixImin2:ixImax2,k)/two
 
           wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,m0_+ii0)+tmp2(ixImin1:ixImax1,ixImin2:ixImax2)&
              *qdt

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=(tmpVR(ixImin1:ixImax1,&
              ixImin2:ixImax2)*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
              -tmpVL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpL(ixImin1:ixImax1,&
              ixImin2:ixImax2))/dx(ixImin1:ixImax1,ixImin2:ixImax2,k)/two

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,e_)+tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*qdt
        endif




        if (i .ne. k) then
           call gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
           tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
              ixImin2:ixImax2)*(nuL(ixImin1:ixImax1,ixImin2:ixImax2)&
              +nuR(ixImin1:ixImax1,ixImin2:ixImax2)+two*nushk(ixImin1:ixImax1,&
              ixImin2:ixImax2,k))/two/two

           tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmprhoC(ixImin1:ixImax1,&
              ixImin2:ixImax2)*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
           call gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,i,tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,m0_+ii0)+tmpC(ixImin1:ixImax1,ixImin2:ixImax2)&
              *qdt

           tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
              ixImin2:ixImax2,m0_+ii0)*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
           call gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,i,tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,e_)+tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt
        endif

     enddo
    enddo
   enddo





do k=1,ndim
 do l=1,ndim

   if (k .ne. l) then

    call setnu(w,b0_+l,k,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)

    do ii1=0,1

      if (ii1 .eq. 0) then
              ii0=k
              m=l
              sB=-1.d0
              j=k
      endif

      if (ii1 .eq. 1) then 
              ii0=l    !ii0 is index B
              m=k      !first derivative
              sB=1.d0  !sign B
              j=l      !first B in energy
      endif



!print*,'k,l,m,j,ii0,ii1=',k,l,m,j,ii0,ii1



      if (m .eq. k) then

           jxmin1=ixmin1+kr(m,1);jxmin2=ixmin2+kr(m,2);jxmax1=ixmax1+kr(m,1)
           jxmax2=ixmax2+kr(m,2); 
           hxmin1=ixmin1-kr(m,1);hxmin2=ixmin2-kr(m,2);hxmax1=ixmax1-kr(m,1)
           hxmax2=ixmax2-kr(m,2);
           tmpBL(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,ixmin2:ixmax2,&
              b0_+j)+w(hxmin1:hxmax1,hxmin2:hxmax2,b0_+j))/two
           tmpBR(ixmin1:ixmax1,ixmin2:ixmax2)=(w(jxmin1:jxmax1,jxmin2:jxmax2,&
              b0_+j)+w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+j))/two

           tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
              ixImin2:ixImax2,b0_+l)

           call gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
           tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
              ixImin2:ixImax2))*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
           call gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
           tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
              ixImin2:ixImax2))*tmp2(ixImin1:ixImax1,ixImin2:ixImax2) 

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,b0_+ii0)+sB*(tmpR(ixImin1:ixImax1,&
              ixImin2:ixImax2)-tmpL(ixImin1:ixImax1,ixImin2:ixImax2))&
              /dx(ixImin1:ixImax1,ixImin2:ixImax2,k)*qdt

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,e_)+sB*(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
              *tmpBR(ixImin1:ixImax1,ixImin2:ixImax2)-tmpL(ixImin1:ixImax1,&
              ixImin2:ixImax2)*tmpBL(ixImin1:ixImax1,ixImin2:ixImax2))&
              /dx(ixImin1:ixImax1,ixImin2:ixImax2,k)*qdt


      endif



      if (m .ne. k) then

           tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
              ixImin2:ixImax2,b0_+l)

           call gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
              ixImin2:ixImax2)*(nuL(ixImin1:ixImax1,ixImin2:ixImax2)&
              +nuR(ixImin1:ixImax1,ixImin2:ixImax2))/two

           call gradient1(tmp2,ixmin1,ixmin2,ixmax1,ixmax2,m,tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,b0_+ii0)+sB*tmpC(ixImin1:ixImax1,&
              ixImin2:ixImax2)*qdt

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
              ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,b0_+j)

           call gradient1(tmp2,ixmin1,ixmin2,ixmax1,ixmax2,m,tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
              ixImin2:ixImax2,e_)+sB*tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt

      endif


      enddo
   endif
 enddo
enddo




return
end

!=============================================================================
subroutine setnu(w,iw,idim,ixmin1,ixmin2,ixmax1,ixmax2,nuR,nuL)

! Set the viscosity coefficient nu within ixO based on w(ixI). 

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision:: d1R(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1),d1L(ixGlo1:ixGhi1&
   +1,ixGlo2:ixGhi2+1)
double precision:: d3R(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1),d3L(ixGlo1:ixGhi1&
   +1,ixGlo2:ixGhi2+1)
double precision:: md3R(ixGlo1:ixGhi1,ixGlo2:ixGhi2),md3L(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision:: md1R(ixGlo1:ixGhi1,ixGlo2:ixGhi2),md1L(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision:: nuR(ixGlo1:ixGhi1,ixGlo2:ixGhi2),nuL(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

double precision:: c_tot, c_hyp,cmax(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
    tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim, iw
integer:: kxmin1,kxmin2,kxmax1,kxmax2,jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,&
   hxmin2,hxmax1,hxmax2,gxmin1,gxmin2,gxmax1,gxmax2,ixFFmin1,ixFFmin2,&
   ixFFmax1,ixFFmax2,jxFFmin1,jxFFmin2,jxFFmax1,jxFFmax2,hxFFmin1,hxFFmin2,&
   hxFFmax1,hxFFmax2
integer:: ix_1,ix_2,ix_3

integer:: ixFlo1,ixFlo2,ixFhi1,ixFhi2,ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixYlo1,&
   ixYlo2,ixYhi1,ixYhi2

logical:: new_cmax

double precision:: tmp_nuI(ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2)

integer:: k,iwc

integer:: ix,ixe

      

integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus


integer:: hpe,jpe

double precision:: tgtbufferR1(1,ixGlo2:ixGhi2+2),tgtbufferR2(ixGlo1:ixGhi1&
   +2,1)
double precision:: tgtbufferL1(1,ixGlo2:ixGhi2+2),tgtbufferL2(ixGlo1:ixGhi1&
   +2,1)
double precision:: srcbufferR1(1,ixGlo2:ixGhi2+2),srcbufferR2(ixGlo1:ixGhi1&
   +2,1)
double precision:: srcbufferL1(1,ixGlo2:ixGhi2+2),srcbufferL2(ixGlo1:ixGhi1&
   +2,1)

integer:: n



!----------------------------------------------------------------------------

new_cmax=.true.

call getcmax(new_cmax,w,ixmin1,ixmin2,ixmax1,ixmax2,idim,cmax)
c_tot=maxval(cmax(ixmin1:ixmax1,ixmin2:ixmax2))

      call mpiallreduce(c_tot,MPI_MAX)



c_hyp=0.4d0 ! 0.6

!if (iw.eq.b^D_|.or.) c_hyp=0.02d0

!if (iw .eq. rho_) c_hyp=0.045d0

!if (iw .eq. 173) c_hyp=0.02d0



if (iw.eq.b1_.or.iw.eq.b2_) c_hyp=0.02d0

if (iw .eq. rho_) c_hyp=0.02d0

if (iw .eq. 173) c_hyp=0.02d0

        
if (iw .ne. 173) then     
        tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw)
        if (iw.eq.m1_.or.iw.eq.m2_) tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)&
           =w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw)/(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
           rho_)+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,rhob_))
endif

if (iw .eq. 173) tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=w(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,e_)-half*((w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,b1_)**2&
   +w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,b2_)**2)+(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   m1_)**2+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,m2_)**2)/(w(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,rho_)+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,rhob_)))


ixYlo1=ixmin1-2;ixYlo2=ixmin2-2;ixYhi1=ixmax1+2;ixYhi2=ixmax2+2;

ixFlo1=ixYlo1+1;ixFlo2=ixYlo2+1;ixFhi1=ixYhi1+1;ixFhi2=ixYhi2+1;

tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2)=tmp_nu(ixYlo1:ixYhi1,ixYlo2:ixYhi2)


     

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

n = (ixFhi1-ixFlo1+1)*(ixFhi2-ixFlo2+1)   

select case(idim)
  case(1)

 n=n/(ixFhi1-ixFlo1+1)

  case(2)

 n=n/(ixFhi2-ixFlo2+1)


end select



select case(idim)
   case(1)

 if(npe1>1)then

   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL


!source
   srcbufferL1(1,ixFlo2:ixFhi2)=tmp_nuI(ixFlo1+4,ixFlo2:ixFhi2) !left, lower

   srcbufferR1(1,ixFlo2:ixFhi2)=tmp_nuI(ixFhi1-4,ixFlo2:ixFhi2) !right, upper

   call mpineighbors(1,hpe,jpe)


   if (mpiupperB(1)) nmpirequest=nmpirequest+1
   if (mpiupperB(1)) call MPI_IRECV(tgtbufferR1(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   if (mpilowerB(1)) nmpirequest=nmpirequest+1
   if (mpilowerB(1)) call MPI_IRECV(tgtbufferL1(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

   if (mpiupperB(1)) call MPI_RSEND(srcbufferR1(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

   if (mpilowerB(1)) call MPI_RSEND(srcbufferL1(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

!target
   tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2)=tgtbufferR1(1,ixFlo2:ixFhi2) !right, upper R

   tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2)=tgtbufferL1(1,ixFlo2:ixFhi2) !left, lower  L


  endif
    case(2)

 if(npe2>1)then

   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL


!source
   srcbufferL2(ixFlo1:ixFhi1,1)=tmp_nuI(ixFlo1:ixFhi1,ixFlo2+4) !left, lower

   srcbufferR2(ixFlo1:ixFhi1,1)=tmp_nuI(ixFlo1:ixFhi1,ixFhi2-4) !right, upper

   call mpineighbors(2,hpe,jpe)


   if (mpiupperB(2)) nmpirequest=nmpirequest+1
   if (mpiupperB(2)) call MPI_IRECV(tgtbufferR2(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   if (mpilowerB(2)) nmpirequest=nmpirequest+1
   if (mpilowerB(2)) call MPI_IRECV(tgtbufferL2(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

   if (mpiupperB(2)) call MPI_RSEND(srcbufferR2(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

   if (mpilowerB(2)) call MPI_RSEND(srcbufferL2(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

!target
   tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1)=tgtbufferR2(ixFlo1:ixFhi1,1) !right, upper R

   tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1)=tgtbufferL2(ixFlo1:ixFhi1,1) !left, lower  L


  endif

end select

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)



if (iw .eq. 173) then 
   iwc=e_ 
else 
   iwc=iw
endif

do k=0,1  !left-right bc

if (typeB(iwc,2*idim-1+k) .ne. 'mpi') then
      if (upperB(2*idim-1+k)) then

          select case(idim)
               case(1)
                     tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2)=tmp_nuI(ixFhi1&
                        -5,ixFlo2:ixFhi2)
                case(2)
                     tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1)=tmp_nuI(ixFlo1:ixFhi1,&
                        ixFhi2-5)
            
          end select

      else

          select case(idim)
               case(1)
                     tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2)=tmp_nuI(ixFlo1&
                        +5,ixFlo2:ixFhi2)
                case(2)
                     tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1)=tmp_nuI(ixFlo1:ixFhi1,&
                        ixFlo2+5)
            
          end select

      endif
endif

enddo 

        ixFmin1=ixFlo1+1;ixFmin2=ixFlo2+1;ixFmax1=ixFhi1-1;ixFmax2=ixFhi2-1; 

        kxmin1=ixFmin1+2*kr(idim,1);kxmin2=ixFmin2+2*kr(idim,2)
        kxmax1=ixFmax1+2*kr(idim,1);kxmax2=ixFmax2+2*kr(idim,2); !5:66
        jxmin1=ixFmin1+kr(idim,1);jxmin2=ixFmin2+kr(idim,2)
        jxmax1=ixFmax1+kr(idim,1);jxmax2=ixFmax2+kr(idim,2); !4:65
        hxmin1=ixFmin1-kr(idim,1);hxmin2=ixFmin2-kr(idim,2)
        hxmax1=ixFmax1-kr(idim,1);hxmax2=ixFmax2-kr(idim,2); !2:63
        gxmin1=ixFmin1-2*kr(idim,1);gxmin2=ixFmin2-2*kr(idim,2)
        gxmax1=ixFmax1-2*kr(idim,1);gxmax2=ixFmax2-2*kr(idim,2); !1:62

        ixFFmin1=ixFlo1;ixFFmin2=ixFlo2;ixFFmax1=ixFhi1;ixFFmax2=ixFhi2; !2:65
        jxFFmin1=ixFlo1+kr(idim,1);jxFFmin2=ixFlo2+kr(idim,2)
        jxFFmax1=ixFhi1+kr(idim,1);jxFFmax2=ixFhi2+kr(idim,2); !3:66
        hxFFmin1=ixFlo1-kr(idim,1);hxFFmin2=ixFlo2-kr(idim,2)
        hxFFmax1=ixFhi1-kr(idim,1);hxFFmax2=ixFhi2-kr(idim,2); !1:64
        
        d3R(ixFmin1:ixFmax1,ixFmin2:ixFmax2)=abs(3.d0*(tmp_nuI(jxmin1:jxmax1,&
           jxmin2:jxmax2)-tmp_nuI(ixFmin1:ixFmax1,ixFmin2:ixFmax2))&
           -(tmp_nuI(kxmin1:kxmax1,kxmin2:kxmax2)-tmp_nuI(hxmin1:hxmax1,&
           hxmin2:hxmax2))) !3:64
        d1R(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2)=abs(tmp_nuI&
           (jxFFmin1:jxFFmax1,jxFFmin2:jxFFmax2)-tmp_nuI(ixFFmin1:ixFFmax1,&
           ixFFmin2:ixFFmax2)) !2:65

        do ix_1=ixmin1,ixmax1
        do ix_2=ixmin2,ixmax2    !3:62  +1=4:63

          md3R(ix_1,ix_2)=maxval(d3R(ix_1+1-kr(idim,1):ix_1+1+kr(idim,1),ix_2&
             +1-kr(idim,2):ix_2+1+kr(idim,2)))
          md1R(ix_1,ix_2)=maxval(d1R(ix_1+1-2*kr(idim,1):ix_1+1&
             +2*kr(idim,1),ix_2+1-2*kr(idim,2):ix_2+1+2*kr(idim,2)))

        enddo
        enddo

        WHERE (md1R(ixmin1:ixmax1,ixmin2:ixmax2).gt.0.d0)
          nuR(ixmin1:ixmax1,ixmin2:ixmax2)=c_tot*c_hyp*md3R(ixmin1:ixmax1,&
             ixmin2:ixmax2)/md1R(ixmin1:ixmax1,ixmin2:ixmax2)&
             *dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)
        ELSEWHERE 
          nuR(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0
        END WHERE
        
        maxviscoef=max(maxval(nuR(ixmin1:ixmax1,ixmin2:ixmax2)), maxviscoef)


!************

        d3L(ixFmin1:ixFmax1,ixFmin2:ixFmax2)=abs(3.d0*(tmp_nuI&
           (ixFmin1:ixFmax1,ixFmin2:ixFmax2)-tmp_nuI(hxmin1:hxmax1,&
           hxmin2:hxmax2))-(tmp_nuI(jxmin1:jxmax1,jxmin2:jxmax2)&
           -tmp_nuI(gxmin1:gxmax1,gxmin2:gxmax2)))
        d1L(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2)=abs(tmp_nuI&
           (ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2)-tmp_nuI(hxFFmin1:hxFFmax1,&
           hxFFmin2:hxFFmax2))    

        do ix_1=ixmin1,ixmax1
        do ix_2=ixmin2,ixmax2

          md3L(ix_1,ix_2)=maxval(d3L(ix_1+1-kr(idim,1):ix_1+1+kr(idim,1),ix_2&
             +1-kr(idim,2):ix_2+1+kr(idim,2)))
          md1L(ix_1,ix_2)=maxval(d1L(ix_1+1-2*kr(idim,1):ix_1+1&
             +2*kr(idim,1),ix_2+1-2*kr(idim,2):ix_2+1+2*kr(idim,2)))

        enddo
        enddo

        WHERE (md1L(ixmin1:ixmax1,ixmin2:ixmax2).gt.0.d0)
          nuL(ixmin1:ixmax1,ixmin2:ixmax2)=c_tot*c_hyp*md3L(ixmin1:ixmax1,&
             ixmin2:ixmax2)/md1L(ixmin1:ixmax1,ixmin2:ixmax2)&
             *dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)
        ELSEWHERE 
          nuL(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0  
        END WHERE

        maxviscoef=max(maxval(nuL(ixmin1:ixmax1,ixmin2:ixmax2)), maxviscoef)

              call mpiallreduce(maxviscoef,MPI_MAX)

return
end


!=============================================================================
!=============================================================================
subroutine setnushk(w,ixmin1,ixmin2,ixmax1,ixmax2,nushk)

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),tmp2(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),nushk(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ndim)
double precision:: c_shk

double precision:: tmp3(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim, iw,i

integer:: ix_1,ix_2

do idim=1,ndim
nushk(ixmin1:ixmax1,ixmin2:ixmax2,idim)=0.d0
enddo


go to 100
c_shk=0.5d0

tmp3(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0

!**************************BEGIN shock viscosity*******************************
      do idim=1,ndim
         tmp(ixmin1:ixmax1,ixmin2:ixmax2)=w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
            +idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,&
            ixmin2:ixmax2,rhob_))
         call gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
         tmp3(ixmin1:ixmax1,ixmin2:ixmax2)=tmp3(ixmin1:ixmax1,ixmin2:ixmax2)&
            +tmp2(ixmin1:ixmax1,ixmin2:ixmax2)
       enddo
      do idim=1,ndim
        nushk(ixmin1:ixmax1,ixmin2:ixmax2,idim)=tmp3(ixmin1:ixmax1,&
           ixmin2:ixmax2)*(dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)**2.d0)*c_shk
        WHERE (tmp3(ixmin1:ixmax1,ixmin2:ixmax2) .ge. 0.d0)
!	  nushk(ix^S,idim)=0.d0
        END WHERE
        nushk(ixmin1:ixmax1,ixmin2:ixmax2,idim)=abs(nushk(ixmin1:ixmax1,&
           ixmin2:ixmax2,idim))
      enddo
!****************************END shock viscosity*******************************

100 continue


return
end



!=============================================================================
subroutine getdt_visc(w,ixmin1,ixmin2,ixmax1,ixmax2)

! Check diffusion time limit for dt < dtdiffpar * dx**2 / (nu/rho)

! Based on Hirsch volume 2, p.631, eq.23.2.17

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),dtdiff_visc
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim, ix_1,ix_2

integer:: aa

! For spatially varying nu you need a common nu array
 double precision::tmpdt(ixGlo1:ixGhi1,ixGlo2:ixGhi2), nuL(ixGlo1:ixGhi1,&
    ixGlo2:ixGhi2),nuR(ixGlo1:ixGhi1,ixGlo2:ixGhi2), nushk(ixGlo1:ixGhi1,&
    ixGlo2:ixGhi2,ndim)
 common/visc/nuL
 common/visc/nuR
!-----------------------------------------------------------------------------

 call setnushk(w,ixmin1,ixmin2,ixmax1,ixmax2,nushk)

 dtdiffpar=0.25d0

 do idim=1,ndim
   tmpdt(ixmin1:ixmax1,ixmin2:ixmax2)=(maxviscoef+nushk(ixmin1:ixmax1,&
      ixmin2:ixmax2,idim)) !/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))   ! 1/dt
   dtdiff_visc=dtdiffpar/maxval(tmpdt(ixmin1:ixmax1,ixmin2:ixmax2)&
      /(dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)**2))
         call mpiallreduce(dtdiff_visc,MPI_MIN)
   dt=min(dt,dtdiff_visc)
 end do
 
 maxviscoef=0.d0

return
end


!***** 2-point central finite difference gradient******

subroutine gradient1(q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)
include 'vacdef.f'
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
integer:: hxmin1,hxmin2,hxmax1,hxmax2,kxmin1,kxmin2,kxmax1,kxmax2
integer:: minx11,minx12,maxx11,maxx12,k
!-----------------------------------------------------------------------------

hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmax1=ixmax1-kr(idim,1)
hxmax2=ixmax2-kr(idim,2);
kxmin1=ixmin1+kr(idim,1);kxmin2=ixmin2+kr(idim,2);kxmax1=ixmax1+kr(idim,1)
kxmax2=ixmax2+kr(idim,2);
gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(q(kxmin1:kxmax1,kxmin2:kxmax2)&
   -q(hxmin1:hxmax1,hxmin2:hxmax2))/dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)/two

 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2)=0.d0
 gradq(maxx11,ixmin2:ixmax2)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2)=0.d0
 gradq(ixmin1:ixmax1,maxx12)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2)=0.d0
 gradq(minx11,ixmin2:ixmax2)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2)=0.d0
 gradq(ixmin1:ixmax1,minx12)=0.d0

 end select
 endif
 endif
 enddo


return
end

!=============================================================================


!*****left upwind forward 2-point non-central finite difference gradient******

subroutine gradient1L(q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)
include 'vacdef.f'
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
integer:: hxmin1,hxmin2,hxmax1,hxmax2
integer:: minx11,minx12,maxx11,maxx12,k
!-----------------------------------------------------------------------------

hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmax1=ixmax1-kr(idim,1)
hxmax2=ixmax2-kr(idim,2);
gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(q(ixmin1:ixmax1,ixmin2:ixmax2)&
   -q(hxmin1:hxmax1,hxmin2:hxmax2))/dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)

 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2)=0.d0
 gradq(maxx11,ixmin2:ixmax2)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2)=0.d0
 gradq(ixmin1:ixmax1,maxx12)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2)=0.d0
 gradq(minx11,ixmin2:ixmax2)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2)=0.d0
 gradq(ixmin1:ixmax1,minx12)=0.d0

 end select
 endif
 endif
 enddo


return
end

!=============================================================================

!*****right upwind forward 2-point non-central finite difference gradient*****

subroutine gradient1R(q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)
include 'vacdef.f'
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
integer:: hxmin1,hxmin2,hxmax1,hxmax2
integer:: minx11,minx12,maxx11,maxx12,k
!-----------------------------------------------------------------------------

hxmin1=ixmin1+kr(idim,1);hxmin2=ixmin2+kr(idim,2);hxmax1=ixmax1+kr(idim,1)
hxmax2=ixmax2+kr(idim,2);
gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(q(hxmin1:hxmax1,hxmin2:hxmax2)&
   -q(ixmin1:ixmax1,ixmin2:ixmax2))/dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)

 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2)=0.d0
 gradq(maxx11,ixmin2:ixmax2)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2)=0.d0
 gradq(ixmin1:ixmax1,maxx12)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2)=0.d0
 gradq(minx11,ixmin2:ixmax2)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2)=0.d0
 gradq(ixmin1:ixmax1,minx12)=0.d0

 end select
 endif
 endif
 enddo


return
end

!=============================================================================
subroutine specialini(ixmin1,ixmin2,ixmax1,ixmax2,w)

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)

double precision:: rhoin,xcent1,xcent2,radius
double precision:: inirho,iniene
double precision:: onemor,inix,ddx
double precision:: p_1,p_2

integer:: iii_,iix_1,info,i,j
double precision:: pi,comi,eneu,sum,mode,bmax,l
character*79 atmfilename

integer:: ix_1,ix_2

double precision:: ixc_1,ixc_2
double precision:: rfc,a1,a2


double precision:: bz,xc, bmin,bmax

!-----------------------------------------------------------------------------

      if (ipe .ne. 0) read(unitpar,'(a)') atmfilename

      if (ipe .eq. 0) then
write(*,*)'param load'
read(unitpar,'(a)') atmfilename
write(*,*) 'from file ',atmfilename
open(33,file=atmfilename,status='old')

!iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)

!for VALMc
!iniene=7.89e5*8.31e3*(9.0412855e-13)/0.6d0/(eqpar(gamma_)-1.0)
      
!for Phot
!iniene=7160.d0*8.31e3*(2.7479785e-10)/0.6d0/(eqpar(gamma_)-1.0)


!for Phot 727.2 tube case with 1050 Gauss
!iniene=5500.d0*8.31e3*(5.7413280e-07)/0.6d0/(eqpar(gamma_)-1.0)

!for Photosphere -- TR 2.6 Mm 248 
!iniene=770000.d0*8.31e3*(9.5690001e-12)/1.2d0/(eqpar(gamma_)-1.0)

!for Photosphere -- TR 2.6 Mm 505 
!iniene=770000.d0*8.31e3*(9.5815098e-12)/1.2d0/(eqpar(gamma_)-1.0)

!for Photosphere -- TR 2.6 Mm 252 
!iniene=770000.d0*8.31e3*(9.5815098e-12)/1.2d0/(eqpar(gamma_)-1.0)

!for Photosphere -- TR 2.6 Mm 1023
iniene=770000.d0*8.31e3*(9.5815101e-12)/1.2d0/(eqpar(gamma_)-1.0)

!for Photosphere -- TR 2.6 Mm 1976 
iniene=770000.d0*8.31e3*(7.1818159e-12)/1.2d0/(eqpar(gamma_)-1.0)

      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
      call MPI_BCAST(iniene,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)

do ix_1=ixGhi1,ixGlo1,-1
      if (ipe .eq. 0) read(33,*) inix,inirho

      if (ipe .eq. 0) print*,'inix,inirho=',inix,inirho


      call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
      call MPI_BCAST(inix,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)
      call MPI_BCAST(inirho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)

 do ix_2=ixGlo2,ixGhi2

   x(ix_1,ix_2,1)=inix !*1000.d0
   w(ix_1,ix_2,rho_)=inirho
   w(ix_1,ix_2,e_)=iniene
   w(ix_1,ix_2,m1_)=0.0
   w(ix_1,ix_2,m2_)=0.0

 enddo


enddo
      if (ipe .eq. 0) close(33)

      if (ipe .eq. 0) print*,'grav=',eqpar(grav0_),eqpar(grav1_),eqpar(grav2_)



print*, '###################',ixmin1,ixmin2,ixmax1,ixmax2

call primitive(ixmin1,ixmin2,ixmax1,ixmax2,w)

do ix_2=ixGlo2,ixGhi2
 do ix_1=ixGhi1-1,ixGlo1,-1 

comi=-abs(x(ix_1+1,ix_2,1)-x(ix_1,ix_2,1))

w(ix_1,ix_2,p_)=w(ix_1+1,ix_2,p_)+w(ix_1,ix_2,rho_)*comi*1.d0*eqpar(grav1_)

if (ix_2 .eq. ixGlo2) print*,'eee=',w(ix_1,ix_2,rho_),w(ix_1,ix_2,p_)

 enddo
enddo


!goto 200

do ix_2=ixGlo2,ixGhi2
 do ix_1=ixGlo1+2,ixGhi1-2
       
       w(ix_1,ix_2,rho_)=-(1.D0/eqpar(grav1_))*(1.D0/(12.D0*(x(ix_1&
          +1,ix_2,1)-x(ix_1,ix_2,1))))*(w(ix_1+2,ix_2,p_)-8.D0*w(ix_1&
          +1,ix_2,p_)+8.D0*w(ix_1-1,ix_2,p_)-w(ix_1-2,ix_2,p_))
               
if (ix_2 .eq. ixGlo2) print*,'eeee=',w(ix_1,ix_2,rho_),w(ix_1,ix_2,p_)

     enddo
   enddo



!lower boundary
do ix_1=ixmin1+4,ixmin1+2,-1
  do ix_2=ixmin2,ixmax2
!    do ix_3=ixmin3,ixmax3
!        p_1=w(ix_1+2,ix_2,ix_3,p_)-8.d0*w(ix_1+1,ix_2,ix_3,p_)+8.d0*w(ix_1-1,ix_2,ix_3,p_)
!        p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
!        w(ix_1-2,ix_2,ix_3,p_) = 12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1-1,ix_2,ix_3,1))*p_2+p_1
         p_1=w(ix_1+2,ix_2,p_)-8.d0*w(ix_1+1,ix_2,p_)+8.d0*w(ix_1-1,ix_2,p_)
         p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
         w(ix_1-2,ix_2,p_) = 12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2+p_1
!       enddo
    enddo
 enddo


!upper boundary
do ix_1=ixmax1-4,ixmax1-2
   do ix_2=ixmin2,ixmax2
!      do ix_3=ixmin3,ixmax3
         
!          p_1=w(ix_1-2,ix_2,ix_3,p_)-8.d0*w(ix_1-1,ix_2,ix_3,p_)+8.d0*w(ix_1+1,ix_2,ix_3,p_)
!          p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
!          w(ix_1+2,ix_2,ix_3,p_) = -12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1-1,ix_2,ix_3,1))*p_2+p_1
           p_1=w(ix_1-2,ix_2,p_)-8.d0*w(ix_1-1,ix_2,p_)+8.d0*w(ix_1+1,ix_2,p_)
           p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
           w(ix_1+2,ix_2,p_) = -12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2&
              +p_1
      enddo
   enddo
!enddo


200 continue



call conserve(ixmin1,ixmin2,ixmax1,ixmax2,w)




w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)=w(ixmin1:ixmax1,ixmin2:ixmax2,e_)
w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=0.d0

w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)=w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=0.d0


!w(ix^S,b1_)=0.04d0

!magnetic field
goto 50
! ********************** smooth  Bx magnetic field***************
xc=2800000.d0 !x(ixmax1,10,1)/2.0


  do ix_1=ixmin1,ixmax1
   do ix_2=ixmin2,ixmax2

      w(ix_1,ix_2,b2_)=((atan((x(ix_1,10,1)-xc)/x(ixmax1,10,1)*100.d0))&
         +Pi/2.d0)/Pi

   enddo
  enddo

bmax=maxval(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_))
bmin=minval(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_))

    w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)=0.04d0*(w(ixmin1:ixmax1,ixmin2:ixmax2,&
       b2_)-bmin)/bmax
    
   
! **************************************************************
50 continue
goto 101
! ********************** Bz(z) magnetic field***************

  do ix_1=ixmin1,ixmax1
   do ix_2=ixmin2,ixmax2

      w(ix_1,ix_2,b1_)=0.4d0*(1.d0-x(ix_1,10,1)/x(ixmax1,10,1))

   enddo
  enddo
    
! **************************************************************
101 continue
goto 100
! ********************** smooth  Bz magnetic field***************
xc=x(10,ixmax2,2)/2.0


  do ix_1=ixmin1,ixmax1
   do ix_2=ixmin2,ixmax2

      w(ix_1,ix_2,b1_)=((atan((x(10,ix_2,2)-xc)/x(10,ixmax2,2)*100.d0))&
         +Pi/2.d0)/Pi

   enddo
  enddo

bmax=maxval(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))
bmin=minval(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))

    w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)=0.04*(w(ixmin1:ixmax1,ixmin2:ixmax2,&
       b1_)-bmin)/bmax
    
! **************************************************************
100 continue
goto 600
! ********************** smooth  opoz Bz magnetic field***************
xc=x(10,ixmax2,2)/2.0


  do ix_1=ixmin1,ixmax1
   do ix_2=ixmin2,ixmax2

      w(ix_1,ix_2,b1_)=((atan((x(10,ix_2,2)-xc)/x(10,ixmax2,2)*100.d0))&
         +Pi/2.d0)/Pi

   enddo
  enddo

bmax=maxval(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))
bmin=minval(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))

    w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)=0.8d0*((w(ixmin1:ixmax1,ixmin2:ixmax2,&
       b1_)-bmin)/bmax-0.5d0)
    
! **************************************************************
600 continue
goto 700
! ********************** smooth  tube Bz magnetic field***************
xc=x(10,ixmax2,2)/2.0


  do ix_1=ixmin1,ixmax1
   do ix_2=ixmin2,ixmax2

      tmp(ix_1,ix_2)=((atan((x(10,ix_2,2)-xc)/x(10,ixmax2,2)*100.d0))&
         +Pi/2.d0)/Pi

   enddo
  enddo

bmax=maxval(tmp(ixmin1:ixmax1,ixmin2:ixmax2))
bmin=minval(tmp(ixmin1:ixmax1,ixmin2:ixmax2))

    tmp(ixmin1:ixmax1,ixmin2:ixmax2)=(tmp(ixmin1:ixmax1,ixmin2:ixmax2)&
       -bmin)/bmax
       
xc=3.0d0*x(10,ixmax2,2)/4.0


  do ix_1=ixmin1,ixmax1
   do ix_2=ixmin2,ixmax2

      w(ix_1,ix_2,b1_)=((atan((x(10,ix_2,2)-xc)/x(10,ixmax2,2)*100.d0))&
         +Pi/2.d0)/Pi

   enddo
  enddo

bmax=maxval(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))
bmin=minval(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))

    w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)=0.4*(tmp(ixmin1:ixmax1,ixmin2:ixmax2)&
       -(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)-bmin)/bmax)
    
    
    
! **************************************************************
700 continue


w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)=w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)&
   -(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2.d0+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   b2_)**2.d0)*(1.d0-eqpar(gamma_)/2.d0)/(eqpar(gamma_)-1.d0)

w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_)=w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)
w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)=w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)
w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)=0.d0
w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)=0.d0


300 continue

return
end


!=============================================================================
subroutine specialsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)


include 'vacdef.f'

integer:: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision:: fdt,fdthalf2

double precision:: pre(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tem(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),kapr(ixGlo1:ixGhi1,ixGlo2:ixGhi2),so(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),flux(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision:: tau(ixGlo1:ixGhi1,ixGlo2:ixGhi2),ine(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

double precision:: preg(ixGlo1:ixGhi1,ixGlo2:ixGhi2),pret(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

integer:: rix_1,i,j
double precision:: mol_0, rrr_

double precision:: fsokr,avgflux

integer:: iw,iiw,iix_1

integer:: ix_1,ix_2


double precision:: s_period,xc1,xc2,s_rad1,s_rad2,vvv,r1,r2
 

!*****************
double precision:: t01,t02,a1,a2,s1,s2,sf,rad,rfc,sdep,tdep,sigma2
double precision:: xc21,xc22,xc23,xc24,xc25, r21,r22,r23,r24,r25

!-----------------------------------------------------------------------------


eqpar(eta_)=0.d0

!call addsource_diff(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

eqpar(nu_)=1.0d0
!eqpar(nu_)=0.d0

call addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)



if(abs(eqpar(nu_))>smalldouble)call addsource_visc(qdt,ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)
 



s_period=30.d0

!s_period=24.d0

!xc1=0.1d6 !m
!xc2=10.0d6  !m
!s_rad1=0.02d6 !m
!s_rad2=1.0d6 !m

!****for Phot 727.2 tube case with 1050 Gauss
!xc1=0.2d5 !m
!xc2=5.0d5  !m
!s_rad1=0.02d5 !m
!s_rad2=0.1d6 !m
!****


!s temp bump
!xc1=0.2d6 !m
!xc2=5.0d5  !m
!s_rad1=0.02d5 !m
!s_rad2=0.1d6 !m


! single driver
!xc2=2.0d6  !m
!xc1=0.1d6 !m

! single driver shift
xc2=1.0d6  !m
xc1=0.1d6 !m

! single driver
s_rad1=0.1d5 !m
s_rad2=0.2d6 !m


! multiple driver
!s_rad1=0.1d5 !m
!s_rad2=0.1d6 !m


! muliple driver
!xc21=0.66d6  !m
!xc22=1.32d6  !m
!xc23=2.0d6  !m
!xc24=2.66d6  !m
!xc25=3.32d6  !m

!xc1=0.1d6 !m


;print*, '**Imin', ixImin1, ixImax1
;print*, '**Imax', ixImin2, ixImax2

;print*, '**Omin', ixOmin1, ixOmax1
;print*, '**Omax', ixOmin2, ixOmax2

 
do ix_1=ixOmin1,ixOmax1
 do ix_2=ixOmin2,ixOmax2

! single driver
     r1=(x(ix_1,ix_2,1)-xc1)**2.d0
     r2=(x(ix_1,ix_2,2)-xc2)**2.d0

! multiple driver
!     r1=(x(ix_1,ix_2,1)-xc1)**2.d0
!     r21=(x(ix_1,ix_2,2)-xc21)**2.d0
!     r22=(x(ix_1,ix_2,2)-xc22)**2.d0
!     r23=(x(ix_1,ix_2,2)-xc23)**2.d0          
!     r24=(x(ix_1,ix_2,2)-xc24)**2.d0
!     r25=(x(ix_1,ix_2,2)-xc25)**2.d0          
   
     tdep=sin(qt*2.d0*pi/s_period)


!     vvv=50.d0*tdep*exp(-r1/s_rad1**2.d0-r2/s_rad2**2.d0)     

!****for Phot 727.2 tube case with 1050 Gauss
!     w(ix_1,ix_2,e_)=w(ix_1,ix_2,e_)+(w(ix_1,ix_2,rho_)+w(ix_1,ix_2,rhob_))*(vvv**2.d0)*qdt/2.d0
!     w(ix_1,ix_2,m2_)=w(ix_1,ix_2,m2_)+(w(ix_1,ix_2,rho_)+w(ix_1,ix_2,rhob_))*vvv*qdt
!****

! single driver
     vvv=100.d0*tdep*exp(-r1/s_rad1**2.d0-r2/s_rad2**2.d0)

! multiple driver
!     vvv=exp(-r1/s_rad1**2.d0-r21/s_rad2**2.d0)
!     vvv=vvv+exp(-r1/s_rad1**2.d0-r22/s_rad2**2.d0)
!     vvv=vvv+exp(-r1/s_rad1**2.d0-r23/s_rad2**2.d0)
!     vvv=vvv+exp(-r1/s_rad1**2.d0-r24/s_rad2**2.d0)
!     vvv=vvv+exp(-r1/s_rad1**2.d0-r25/s_rad2**2.d0)
!     vvv=100.d0*tdep*vvv

    
     w(ix_1,ix_2,e_)=w(ix_1,ix_2,e_)+(w(ix_1,ix_2,rho_)+w(ix_1,ix_2,rhob_))&
        *(vvv**2.d0)*qdt/2.d0
     w(ix_1,ix_2,m2_)=w(ix_1,ix_2,m2_)+(w(ix_1,ix_2,rho_)+w(ix_1,ix_2,rhob_))&
        *vvv*qdt

    
 enddo
enddo



      if (ipe.eq.0) write(*,*) '***time=',qt


end



!=============================================================================
subroutine specialbound(qt,ixmin1,ixmin2,ixmax1,ixmax2,iw,iB,w)
include 'vacdef.f'

integer:: ix_1,ix_2

integer:: iwmin,iwmax,idimmin,idimmax
double precision:: qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
integer:: ix,ix1,ix2,ixe,ixf,ixmin1,ixmin2,ixmax1,ixmax2,ixpairmin1,&
   ixpairmin2,ixpairmax1,ixpairmax2,idim,iw,iB
integer:: iwv,jdim

integer:: Ns,i,j
double precision:: ki


double precision:: tmpp1,tmpp2

Ns=1

select case(iB)

case(1)

case(2)


case default
stop 'error iB!'
end select

return
end

!=============================================================================
subroutine getdt_special(w,ixmin1,ixmin2,ixmax1,ixmax2)

! If the Coriolis force is made very strong it may require time step limiting,
! but this is not implemented here.

include 'vacdef.f'
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

!call getdt_diff(w,ix^L)

if(abs(eqpar(nu_))>smalldouble)call getdt_visc(w,ixmin1,ixmin2,ixmax1,ixmax2)


call getdt_grav(w,ixmin1,ixmin2,ixmax1,ixmax2)

return
end


subroutine specialeta(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)
 
include 'vacdef.f'
 
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idirmin
!-----------------------------------------------------------------------------
 
stop 'specialeta is not defined'
end

