*##############################################################################
* module vacusr - gravity

*==============================================================================
*
*    THE FOLLOWING SUBROUTINES ADD GRAVITATIONAL SOURCE TERMS, SET GRAVITY
*
*------------------------------------------------------------------------------
*    See vacusr.t.gravity and vacusrpar.t.gravity for an example of usage
*
*    Gravitational force is added to the momentum equation:
*
*    d m_i/dt += rho*eqpar(grav0_+i)
*
*    Gravitational work is added to the energy equation (if present):
*
*    de/dt += Sum_i m_i*eqpar(grav0_+i)
*
*    The eqpar(grav1_),eqpar(grav2_),... coefficients are the components of 
*    the gravitational acceleration in each dimension. Set them to 0 for no
*    gravity in that direction. 
*    The !!! comments show how a grav array could be used for a spatially
*    (and maybe temporally) varying gravitational field.
*    The setgrav subroutine has to be completed then.
*
*============================================================================
      subroutine addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

* Add gravity source calculated from w to wnew within ixO for all variables 
* in iws. w is at time qtC, wnew is advanced from qt to qt+qdt.

      include 'vacdef.f'

      integer           ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws(niw_)
      double precision  qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  iiw,iw,idim
*!! ! For a spatially varying gravity define the common grav array
*!! double precision:: grav(ixG^T,ndim)
*!! common /gravity/ grav

*-----------------------------------------------------------------------------

*!! ! If grav needs to be calculated only once do it for the whole grid
*!! if(it==itmin)call setgrav(w,ixG^L,ixG^L,grav)
*!! ! Otherwise call setgrav in every time step
*!! call setgrav(w,ixI^L,ixO^L,grav)

* add sources from gravity
      do iiw=1,iws(niw_)
       iw=iws(iiw)
      if(iw.eq.m1_.or.iw.eq.m2_)then
*           dm_i/dt= +rho*g_i
            idim=iw-m0_
            if(abs(eqpar(grav0_+idim)).gt.smalldouble)then
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  wnew(ix_1,ix_2,m0_+idim)=wnew(ix_1,ix_2,m0_+idim)+ qdt*
     &               eqpar(grav0_+idim)*w(ix_1,ix_2,rho_)
               enddo
               enddo
            endif

*           !! ! For a spatially varying gravity use instead of the above lines
*           !! wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
*           !!    qdt*grav(ixO^S,idim)*w(ixO^S,rho_)

      else if(iw.eq.e_)then
*           de/dt= +g_i*m_i
            do idim=1,ndim
               if(abs(eqpar(grav0_+idim)).gt.smalldouble)then
                  do ix_2=ixOmin2,ixOmax2
                  do ix_1=ixOmin1,ixOmax1
                     wnew(ix_1,ix_2,ee_)=wnew(ix_1,ix_2,ee_)+ qdt*
     &                  eqpar(grav0_+idim)*w(ix_1,ix_2,m0_+idim)
                  enddo
                  enddo
               endif

*              !! ! For a spatially varying gravity use instead of the above lines
*              !! wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
*              !!    qdt*grav(ixO^S,idim)*w(ixO^S,m0_+idim)

            end do
*        iw
         end if
*     iiw
      end do

      return
      end
*=============================================================================
*!! subroutine setgrav(w,ixI^L,ixO^L,grav)

* Set the gravitational acceleration within ixO based on x(ixI,ndim) 
* and/or w(ixI,nw)

*!! include 'vacdef.f'

*!! double precision:: w(ixG^T,nw),grav(ixG^T,ndim)
*!! integer:: ixI^L,ixO^L
*----------------------------------------------------------------------------
*!! return
*!! end
*=============================================================================

      subroutine getdt_grav(w,ixmin1,ixmin2,ixmax1,ixmax2)

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim
      double precision  dtgrav
      save dtgrav

*!! ! For spatially varying gravity you need a common grav array
*!! double precision:: grav(ixG^T,ndim)
*!! common/gravity/grav

*----------------------------------------------------------------------------

      oktest=index(teststr,'getdt').ge.1

      if(it.eq.itmin)then
*        If gravity is descibed by the equation parameters, use this:
         dtgrav=bigdouble
         do idim=1,ndim
            maxval_1=-bigdouble
            do idim_3=1,ndim
            do ix_2=ixMmin2,ixMmax2
            do ix_1=ixMmin1,ixMmax1
               maxval_1=max(maxval_1,abs(eqpar(grav0_+idim))/
     &            dx(ix_1,ix_2,idim_3))
            enddo
            enddo
            enddo
            if(abs(eqpar(grav0_+idim)).gt.zero)then
               dtgrav=min(dtgrav,one/sqrt(maxval_1))
            endif
         enddo
*        !! ! For spatially varying gravity use this instead of the lines above:
*        !! call setgrav(w,ixG^L,ixM^L,grav)
*        !! ! If gravity does not change with time, calculate dtgrav here:
*        !! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))
      endif

*!! ! If gravity changes with time, calculate dtgrav here:
*!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))



* limit the time step
      dt=min(dt,dtgrav)
      if(oktest)write(*,*)'Gravity limit for dt:',dtgrav

      return
      end

*=============================================================================

*=============================================================================
      subroutine specialbound(qt,ixmin1,ixmin2,ixmax1,ixmax2,iw,iB,w)

* Calculates the boundary values in the iB-th boundary segment, user-defined

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,iw,iB
      double precision  qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      call die('Special boundary is not defined')
      end
*=============================================================================
*=============================================================================
      subroutine specialini(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Initialize w for VACINI, user-defined

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      call die('Special initial condition is not defined')
      end
*=============================================================================
*INCLUDE:vacnul.specialsource.t
*=============================================================================
      subroutine readfileini_special(w)

* Reads from unitini,filenameini in user-defined format.
* Check readfileini_asc and readfileini_bin in vacio.t on what should be done.

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      call die('Special readfileini is not defined')
      end
*=============================================================================
      subroutine savefileout_special(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)

* Save current results into filenameout in user-defined format.
* Check savefileout_asc and savefileout_bin in vacio.t on what should be done.

      include 'vacdef.f'

      integer  qunit,ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      call die('Special savefileout is not defined')
      end
*=============================================================================
      subroutine savefilelog_special(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)

* Save user-defined log data into filename(filelog_) in user-defined format.
* Check savefilelog_default on opening the file etc.

      include 'vacdef.f'

      integer  qunit,ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      call die('Special savefilelog is not defined')
      end
*=============================================================================

*=============================================================================
* There are 3 subroutines in this file:
*
* specialsource -- for sources other than resistivity
* getdt_special -- for time step conditions other than CFL or resistivity
* specialeta    -- for non-constant resistivity with eqpar(eta_)<zero
*
*=============================================================================
      subroutine specialsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

      include 'vacdef.f'

      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_)
      double precision  qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      call addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

      return
      end

*=============================================================================
      subroutine getdt_special(w,ixmin1,ixmin2,ixmax1,ixmax2)

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixmin1,ixmin2,ixmax1,ixmax2
*-----------------------------------------------------------------------------

      call getdt_grav(w,ixmin1,ixmin2,ixmax1,ixmax2)

      return
      end

*=============================================================================
      subroutine specialeta(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idirmin
*-----------------------------------------------------------------------------

      call die('vacusr.gravity: specialeta is not defined')
      end
*=============================================================================

* end module vacusr - gravity
*##############################################################################
