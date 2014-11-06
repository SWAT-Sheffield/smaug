*##############################################################################
* module vactvdlf 
* Subroutines for TVDLF with Hancock predictor. Also needed for TVD-MUSCL.
*=============================================================================
      subroutine hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,wnew)

* Predictor for TVDLF and TVD-MUSCL schemes due to Hancock 
* with flux limiting over conservative variables.

      include 'vacdef.f'

      double precision  qdt,qtC,qt
      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idimmin,idimmax
      double precision wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision fLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),fRC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),vLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),vRC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  ixmin1,ixmin2,ixmax1,ixmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,
     &   iiw,iw,idim
      logical  transport
*-----------------------------------------------------------------------------

      oktest=index(teststr,'hancock').ge.1
      if(oktest)write(*,*)'Hancock: wCT,wnew', wCT(ixtest1,ixtest2,iwtest),
     &   wnew(ixtest1,ixtest2,iwtest)

* Expand limits in each idim direction in which fluxes are added
      do idim= idimmin,idimmax
         ixmin1=ixOmin1-kr(idim,1)
         ixmin2=ixOmin2-kr(idim,2)
         ixmax1=ixOmax1+kr(idim,1)
         ixmax2=ixOmax2+kr(idim,2)
      enddo
      if(ixImin1.gt.ixmin1.or.ixImin2.gt.ixmin2.or.ixImax1.lt.ixmax1.or.
     &   ixImax2.lt.ixmax2) call die( 
     &   'Error in Hancock: Nonconforming input limits')

      do idim= idimmin,idimmax
*        Calculate w_j+g_j/2 and w_j-g_j/2
*        First copy all variables, then upwind wLC and wRC for the iws.
*        wLC is to the left of ixO, wRC is to the right of wCT.
*        SHIFT
         hxOmin1=ixOmin1-kr(idim,1)
         hxOmin2=ixOmin2-kr(idim,2)
         hxOmax1=ixOmax1-kr(idim,1)
         hxOmax2=ixOmax2-kr(idim,2)
*        SHIFT BEGIN
         do iw_3=1,nw
         do i_2=hxOmin2,hxOmax2
         do i_1=hxOmin1,hxOmax1
            wRC(i_1,i_2,iw_3)=wCT(i_1+(ixOmin1-hxOmin1),i_2+
     &         (ixOmin2-hxOmin2),iw_3)
         enddo
         enddo
         enddo
*        SHIFT END
         do iw_3=1,nw
         do ix_2=ixOmin2,ixOmax2
         do ix_1=ixOmin1,ixOmax1
            wLC(ix_1,ix_2,iw_3)=wCT(ix_1,ix_2,iw_3)
         enddo
         enddo
         enddo
         call upwindLR(ixOmin1,ixOmin2,ixOmax1,ixOmax2,hxOmin1,hxOmin2,
     &      hxOmax1,hxOmax2,iws,idim,wCT,wLC,wRC)

         

         if(oktest)write(*,*)'idim,wLC,wRC:', idim,wLC(ixtest1,ixtest2,
     &      iwtest),wRC(ixtest1,ixtest2,iwtest)

*        Calculate fLC=f(w_j+g_j/2) and fRC=f(w_j-g_j/2) 
*        Add fluxes (eq 4.38c with the choice W~=U, P=I, qdt is not halved here)

*        Calculate vLC and vRC velocities
         call getv(wRC,hxOmin1,hxOmin2,hxOmax1,hxOmax2,idim,vRC)
         call getv(wLC,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,vLC)
         if(oktest)write(*,*)'vLC,vRC',vLC(ixtest1,ixtest2),vRC(ixtest1,
     &      ixtest2)

*        Advect w(iws)
         do iiw=1,iws(niw_)
          iw=iws(iiw)
*           Calculate the fLC and fRC fluxes
            call getflux(wRC,hxOmin1,hxOmin2,hxOmax1,hxOmax2,iw,idim,fRC,
     &         transport)
            call getflux(wLC,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,idim,fLC,
     &         transport)
            if(transport)then
               do i_2=hxOmin2,hxOmax2
               do i_1=hxOmin1,hxOmax1
                  fRC(i_1,i_2)=fRC(i_1,i_2)+vRC(i_1,i_2)*wRC(i_1,i_2,iw)
               enddo
               enddo
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  fLC(ix_1,ix_2)=fLC(ix_1,ix_2)+vLC(ix_1,ix_2)*
     &               wLC(ix_1,ix_2,iw)
               enddo
               enddo
            endif
*           Advect w(iw)
            
            if(gencoord.and.vectoriw(iw).ge.0)then
            
                   call die('Error: Gencoord is off')
            else
                call addflux(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,idim,fLC,
     &             ixOmin1,ixOmin2,ixOmax1,ixOmax2,fRC,hxOmin1,hxOmin2,
     &             hxOmax1,hxOmax2,wnew)
            endif
         end do
         if(oktest)write(*,*)'wnew:',wnew(ixtest1,ixtest2,iwtest)
*     next idim
      end do

      if(typeaxial.ne.'slab'.and.idimmin.eq.1) call addgeometry(qdt,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,wCT,wnew)
      if(sourceunsplit)call addsource2(qdt*(idimmax-idimmin+one)/
     &   ndim,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,
     &   iws,qtC,wCT,qt,wnew)

      return
      end

*=============================================================================
      subroutine tvdmusclf(addfluxok,method,qdt,ixIImin1,ixIImin2,ixIImax1,
     &   ixIImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,
     &   qt,wnew)

* method=='tvdmu'  --> 2nd oreder (may be 3rd order in 1D) TVD-MUSCL scheme. 
* method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.).
* method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
* method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.
*
* addfluxok=.false. is used for dissipative filter
*
* Flux limiting is over primitive or conservative variables. 
* Based on section 4.4.2. 
*
* Equation 4.39b should read j+1/2 !
* Equation BARMIN(10) should read uL_i+1/2=u_n+1/2_i+...

      include 'vacdef.f'

      logical  addfluxok
      character*10   method
      double precision  qdt,qtC,qt
      integer  ixIImin1,ixIImin2,ixIImax1,ixIImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idimmin,idimmax,ix
      double precision wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision fLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),fRC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),vLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),vRC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),cmaxC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      double precision  courantmax
      integer  ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixImin1,ixImin2,ixImax1,
     &   ixImax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,jxCmin1,jxCmin2,jxCmax1,
     &   jxCmax2,iiw,iw,idim
      logical  transport,new_cmax
*-----------------------------------------------------------------------------

      oktest=index(teststr,'tvd').ge.1
      if(oktest)write(*,*)'TVDMUSCLF: wCT,wnew', wCT(ixtest1,ixtest2,iwtest),
     &   wnew(ixtest1,ixtest2,iwtest)

      if(idimmax.gt.idimmin .and. typelimited.eq.'original' .and.
     &    method.ne.'tvdlf1' .and. method.ne.'tvdmu1') call die( 
     &   'Error in TVDMUSCLF: Unsplit dim. and original is limited')

* The flux calculation contracts by one in the idim direction it is applied.
* The limiter contracts the same directions by one more, so expand ixO by 2.
      ixImin1=ixOmin1
      ixImin2=ixOmin2
      ixImax1=ixOmax1
      ixImax2=ixOmax2
      do idim= idimmin,idimmax
         ixImin1=ixImin1-2*kr(idim,1)
         ixImin2=ixImin2-2*kr(idim,2)
         ixImax1=ixImax1+2*kr(idim,1)
         ixImax2=ixImax2+2*kr(idim,2)
      enddo
      if(ixIImin1.gt.ixImin1.or.ixIImin2.gt.ixImin2.or.ixIImax1.lt.ixImax1.or.
     &   ixIImax2.lt.ixImax2)  call die(
     &   'Error in TVDMUSCLF: Nonconforming input limits')

      do idim= idimmin,idimmax
         hxOmin1=ixOmin1-kr(idim,1)
         hxOmin2=ixOmin2-kr(idim,2)
         hxOmax1=ixOmax1-kr(idim,1)
         hxOmax2=ixOmax2-kr(idim,2)
*        ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1
         ixCmax2=ixOmax2
          ixCmin1=hxOmin1
         ixCmin2=hxOmin2

*        Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 (eq.4.38a,b)
*        SHIFT
         jxCmin1=ixCmin1+kr(idim,1)
         jxCmin2=ixCmin2+kr(idim,2)
         jxCmax1=ixCmax1+kr(idim,1)
         jxCmax2=ixCmax2+kr(idim,2)
*        SHIFT BEGIN
         do iw_3=1,nw
         do ix_2=ixCmin2,ixCmax2
         do ix_1=ixCmin1,ixCmax1
            wRC(ix_1,ix_2,iw_3)=wCT(ix_1+(jxCmin1-ixCmin1),ix_2+
     &         (jxCmin2-ixCmin2),iw_3)
         enddo
         enddo
         enddo
*        SHIFT END
         do iw_3=1,nw
         do ix_2=ixCmin2,ixCmax2
         do ix_1=ixCmin1,ixCmax1
            wLC(ix_1,ix_2,iw_3)=wCT(ix_1,ix_2,iw_3)
         enddo
         enddo
         enddo
         if(method.ne.'tvdlf1'.and.method.ne.'tvdmu1')then
         if(typelimited.eq.'previous')then
               call upwindLR(ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCmin1,ixCmin2,
     &            ixCmax1,ixCmax2,iws,idim,wold,wLC,wRC)
         else if(typelimited.eq.'predictor')then
               call upwindLR(ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCmin1,ixCmin2,
     &            ixCmax1,ixCmax2,iws,idim,wCT,wLC,wRC)
         else if(typelimited.eq.'original')then
               call upwindLR(ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCmin1,ixCmin2,
     &            ixCmax1,ixCmax2,iws,idim,wnew,wLC,wRC)
         else
               call die('Error in TVDMUSCLF: No such base for limiter:'//
     &            typelimited)
            end if
         endif

         

         if(oktest)write(*,*)'idim,wLC,wRC:', idim,wLC(ixtest1,ixtest2,
     &      iwtest),wRC(ixtest1,ixtest2,iwtest)

         if(method.eq.'tvdlf'.or.method.eq.'tvdlf1')then
*           For the high order Lax-Friedrich TVDLF scheme the limiter is based on
*           the maximum eigenvalue, it is calculated in advance.
*           The spectral radius of the Jacobian BARMIN(eq.12)=cmaxC.
*           Bernardo&Lin's local Lax-Fri (J.Comp.Phys 1989, 84, 90), though
*           they propose max(cmaxL,cmaxR) instead of cmax((uR+uL)/2) that I use.
*           To save memory we use wLC temporarily to store the mean 0.5*(uR+uL)
            do iw_3=1,nw
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               wLC(ix_1,ix_2,iw_3)=half*(wLC(ix_1,ix_2,iw_3)+
     &            wRC(ix_1,ix_2,iw_3))
            enddo
            enddo
            enddo
            new_cmax=.true.
            call getcmax(new_cmax,wLC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,
     &         cmaxC)
            if(courantpar.gt.zero.and.istep.eq.nstep.and.implpar.lt.zero)then
*              Calculate dtcourant(idim) for next time step
               if(gencoord)then
                  maxval_1=-bigdouble
                  do ix_2=ixCmin2,ixCmax2
                  do ix_1=ixCmin1,ixCmax1
                     maxval_1=max(maxval_1,cmaxC(ix_1,ix_2)*
     &                  surfaceC(ix_1,ix_2,idim)*two/(dvolume(ix_1,ix_2)+
     &                  dvolume(ix_1+(jxCmin1-ixCmin1),ix_2+
     &                  (jxCmin2-ixCmin2))))
                  enddo
                  enddo
                  courantmax=maxval_1
               else
                  maxval_1=-bigdouble
                  do ix_2=ixCmin2,ixCmax2
                  do ix_1=ixCmin1,ixCmax1
                     maxval_1=max(maxval_1,cmaxC(ix_1,ix_2)*
     &                  two/(dx(ix_1,ix_2,idim)+dx(ix_1+(jxCmin1-
     &                  ixCmin1),ix_2+(jxCmin2-ixCmin2),idim)))
                  enddo
                  enddo
                  courantmax=maxval_1
               endif
               
               if(qdt*courantmax.gt.one)then
                  nerror(couranterr_)=nerror(couranterr_)+1
                  if(nerror(couranterr_).eq.1)write(*,*)
     &               'Courant condition error (code=',couranterr_,') at it=',
     &               it,' for direction ',idim,' CFL=',qdt*courantmax
               endif
               if(courantmax.gt.smalldouble)dtcourant(idim)=
     &            min(dtcourant(idim),courantpar/courantmax)
            endif
            if(oktest)write(*,*)',cmaxC,wLR:',cmaxC(ixtest1,ixtest2),
     &         wLC(ixtest1,ixtest2,iwtest)
*           We regain wLC for further use
            do iw_3=1,nw
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               wLC(ix_1,ix_2,iw_3)=2*wLC(ix_1,ix_2,iw_3)-wRC(ix_1,ix_2,iw_3)
            enddo
            enddo
            enddo
*        if TVDLF
         endif

         if(addfluxok)then
*           Calculate velocities for transport fluxes
            call getv(wLC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,vLC)
            call getv(wRC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,vRC)
         endif

         if(addfluxok.or.method.eq.'tvdlf'.or.method.eq.'tvdlf1')then
*           Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iws
            do iiw=1,iws(niw_)
             iw=iws(iiw)
               if(addfluxok)then
                  call getflux(wLC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idim,
     &               fLC,transport)
                  call getflux(wRC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idim,
     &               fRC,transport)
                  if(transport)then
                     do ix_2=ixCmin2,ixCmax2
                     do ix_1=ixCmin1,ixCmax1
                        fLC(ix_1,ix_2)=fLC(ix_1,ix_2)+vLC(ix_1,ix_2)*
     &                     wLC(ix_1,ix_2,iw)
                     enddo
                     enddo
                     do ix_2=ixCmin2,ixCmax2
                     do ix_1=ixCmin1,ixCmax1
                        fRC(ix_1,ix_2)=fRC(ix_1,ix_2)+vRC(ix_1,ix_2)*
     &                     wRC(ix_1,ix_2,iw)
                     enddo
                     enddo
                  endif
                  if(oktest.and.iw.eq.iwtest)write(*,*)'fLC,fRC:',
     &                fLC(ixtest1,ixtest2),fRC(ixtest1,ixtest2)
*                 To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
                  do ix_2=ixCmin2,ixCmax2
                  do ix_1=ixCmin1,ixCmax1
                     fLC(ix_1,ix_2)=half*(fLC(ix_1,ix_2)+fRC(ix_1,ix_2))
                  enddo
                  enddo
               endif

*              Add TVDLF dissipation to the flux
               if(method.eq.'tvdlf'.or.method.eq.'tvdlf1')then
*                 To save memory we use fRC to store -cmax*half*(w_R-w_L)
                  do ix_2=ixCmin2,ixCmax2
                  do ix_1=ixCmin1,ixCmax1
                     fRC(ix_1,ix_2)=-cmaxC(ix_1,ix_2)*half*
     &                  (wRC(ix_1,ix_2,iw)-wLC(ix_1,ix_2,iw))
                  enddo
                  enddo

*                 Reduce dissipation by acmcoef if set
                  if(acmcoef(iw).ge.zero)then
                     do ix_2=ixCmin2,ixCmax2
                     do ix_1=ixCmin1,ixCmax1
                        fRC(ix_1,ix_2)=acmcoef(iw)*fRC(ix_1,ix_2)
                     enddo
                     enddo
                  endif

*                 Reduce dissipation by ACM switch if exponent is set
                  if(acmexpo.gt.zero)then
*                     We use tmp for the jump.
                      if(acmnolim)then
*                        SHIFT BEGIN
                         do ix_2=ixCmin2,ixCmax2
                         do ix_1=ixCmin1,ixCmax1
                            tmp(ix_1,ix_2)=abs(wCT(ix_1+(jxCmin1-
     &                         ixCmin1),ix_2+(jxCmin2-ixCmin2),iw)-
     &                         wCT(ix_1,ix_2,iw))
                         enddo
                         enddo
*                        SHIFT END
                      else
                         do ix_2=ixCmin2,ixCmax2
                         do ix_1=ixCmin1,ixCmax1
                            tmp(ix_1,ix_2)=abs(wRC(ix_1,ix_2,iw)-
     &                         wLC(ix_1,ix_2,iw))
                         enddo
                         enddo
                      endif
*                     acmswitch will overwrite tmp but that is OK.
                      call acmswitch(tmp,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,
     &                   fRC)
                  endif

                  if(addfluxok)then
*                    fLC contains physical+dissipative fluxes
                     do ix_2=ixCmin2,ixCmax2
                     do ix_1=ixCmin1,ixCmax1
                        fLC(ix_1,ix_2)=fLC(ix_1,ix_2)+fRC(ix_1,ix_2)
                     enddo
                     enddo
                  else
*                    fLC contains dissipative flux only
                     do ix_2=ixCmin2,ixCmax2
                     do ix_1=ixCmin1,ixCmax1
                        fLC(ix_1,ix_2)=fRC(ix_1,ix_2)
                     enddo
                     enddo
                  endif
               endif

*              Add high order flux_idim to wnew=U_n+1 (eq.4.9)
               if(gencoord.and.vectoriw(iw).ge.0)then
                  
                         call die('Error: gencoord is off')
               else
                  call addflux(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,idim,
     &               fLC,ixOmin1,ixOmin2,ixOmax1,ixOmax2,fLC,hxOmin1,hxOmin2,
     &               hxOmax1,hxOmax2,wnew)
               endif

               if(oktest.and.iw.eq.iwtest)write(*,*)'fLR(ix),fLR(hx):',
     &             fLC(ixtest1,ixtest2),fLC(ixtest1-kr(idim,1),ixtest2-
     &            kr(idim,2))
*           Next iw
            end do

            if(oktest)write(*,*)'wnew:',wnew(ixtest1,ixtest2,iwtest)

         endif

*        For the MUSCL scheme apply the characteristic based limiter
         if(method.eq.'tvdmu'.or.method.eq.'tvdmu1')then
                  
            call tvdlimit2(method,qdt,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,
     &         ixOmin2,ixOmax1,ixOmax2,iws,idim,wLC,wRC,wnew)
            if(.false.)call die('TVD-MUSCL requires the TVD module to be on!')
         endif
*     Next idim
      enddo

      if(oktest)write(*,*)'wnew with fluxes:',wnew(ixtest1,ixtest2,iwtest)

* Add geometrical and physical sources if physical fluxes were added
      if(addfluxok)then
         if(typeaxial.ne.'slab'.and.idimmin.eq.1) call addgeometry(qdt,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,wCT,wnew)
         if(sourceunsplit)call addsource2(qdt*(idimmax-idimmin+one)/
     &      ndim,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &      ixOmax2,iws,qtC,wCT,qt,wnew)

         if(oktest)write(*,*)'wnew with source:',wnew(ixtest1,ixtest2,iwtest)
      endif

      return
      end

*=============================================================================
      subroutine upwindLR(ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,
     &   ixRmax1,ixRmax2,iws,idim,w,wLC,wRC)

* Determine the upwinded wLC(ixL) and wRC(ixR) from w. For musclomega=1 
* any typelimiter can be used, for musclomega>1 the 'muscl1','muscl2' limiters,
* minmod(x,omega*y) and minmod(omega*x,y) respectively, are used.

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wLC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision  ldw(ixGlo1:ixGhi1,ixGlo2:ixGhi2),dwC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,
     &   ixRmax2,jxRmin1,jxRmin2,jxRmax1,jxRmax2,ixCmin1,ixCmin2,ixCmax1,
     &   ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,iws(niw_),iiw,iw,idim
*-----------------------------------------------------------------------------

      oktest=index(teststr,'upwindlr').ge.1
      if(oktest)write(*,*)'UpwindLR omega:',musclomega

* Eqs.3.38a and 3.38b

      if(oktest.and.idim.eq.idimtest)write(*,*)'wLC,wRC:',wLC(ixtest1,ixtest2,
     &   iwtest),wRC(ixtest1,ixtest2,iwtest)
      if(oktest.and.idim.eq.idimtest)write(*,*)'wg,wh,wi,wj,wk:',w(ixtest1-2*
     &   kr(idim,1),ixtest2-2*kr(idim,2),iwtest),w(ixtest1-kr(idim,1),ixtest2-
     &   kr(idim,2),iwtest),w(ixtest1,ixtest2,iwtest),w(ixtest1+
     &   kr(idim,1),ixtest2+kr(idim,2),iwtest),w(ixtest1+2*
     &   kr(idim,1),ixtest2+2*kr(idim,2),iwtest)

* Transform w,wL,wR to primitive variables 
      if(useprimitive)then
         call primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,w)
         call primitive(ixLmin1,ixLmin2,ixLmax1,ixLmax2,wLC)
         call primitive(ixRmin1,ixRmin2,ixRmax1,ixRmax2,wRC)
      endif

*SHIFT
      jxRmin1=ixRmin1+kr(idim,1)
      jxRmin2=ixRmin2+kr(idim,2)
      jxRmax1=ixRmax1+kr(idim,1)
      jxRmax2=ixRmax2+kr(idim,2)
      ixCmax1=jxRmax1
      ixCmax2=jxRmax2
       ixCmin1=ixLmin1-kr(idim,1)
      ixCmin2=ixLmin2-kr(idim,2)
       
*SHIFT MORE
      jxCmin1=ixCmin1+kr(idim,1)
      jxCmin2=ixCmin2+kr(idim,2)
      jxCmax1=ixCmax1+kr(idim,1)
      jxCmax2=ixCmax2+kr(idim,2)
*SHIFT BEGIN
      do iiw=1,iws(niw_)
       iw=iws(iiw)
         do ix_2=ixCmin2,ixCmax2
         do ix_1=ixCmin1,ixCmax1
            dwC(ix_1,ix_2)=w(ix_1+(jxCmin1-ixCmin1),ix_2+(jxCmin2-
     &         ixCmin2),iw)-w(ix_1,ix_2,iw)
         enddo
         enddo
         if(musclomega.gt.one)then
            typelimiter(iw)='muscl1'
            call dwlimiter2(dwC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idim,ldw)
            do ix_2=ixRmin2,ixRmax2
            do ix_1=ixRmin1,ixRmax1
               wRC(ix_1,ix_2,iw)=wRC(ix_1,ix_2,iw)-muscleta1*
     &            ldw(ix_1+(jxRmin1-ixRmin1),ix_2+(jxRmin2-ixRmin2))
            enddo
            enddo
            do ix_2=ixLmin2,ixLmax2
            do ix_1=ixLmin1,ixLmax1
               wLC(ix_1,ix_2,iw)=wLC(ix_1,ix_2,iw)+muscleta2*ldw(ix_1,ix_2)
            enddo
            enddo
            typelimiter(iw)='muscl2'
            call dwlimiter2(dwC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idim,ldw)
            do ix_2=ixRmin2,ixRmax2
            do ix_1=ixRmin1,ixRmax1
               wRC(ix_1,ix_2,iw)=wRC(ix_1,ix_2,iw)-muscleta2*
     &            ldw(ix_1+(jxRmin1-ixRmin1),ix_2+(jxRmin2-ixRmin2))
            enddo
            enddo
            do ix_2=ixLmin2,ixLmax2
            do ix_1=ixLmin1,ixLmax1
               wLC(ix_1,ix_2,iw)=wLC(ix_1,ix_2,iw)+muscleta1*ldw(ix_1,ix_2)
            enddo
            enddo
         else
            call dwlimiter2(dwC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idim,ldw)
            do ix_2=ixRmin2,ixRmax2
            do ix_1=ixRmin1,ixRmax1
               wRC(ix_1,ix_2,iw)=wRC(ix_1,ix_2,iw)-half*ldw(ix_1+
     &            (jxRmin1-ixRmin1),ix_2+(jxRmin2-ixRmin2))
            enddo
            enddo
            do ix_2=ixLmin2,ixLmax2
            do ix_1=ixLmin1,ixLmax1
               wLC(ix_1,ix_2,iw)=wLC(ix_1,ix_2,iw)+half*ldw(ix_1,ix_2)
            enddo
            enddo
         endif
      end do
*SHIFT END

* Transform w,wL,wR back to conservative variables 
      if(useprimitive)then
         call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,w)
         call conserve(ixLmin1,ixLmin2,ixLmax1,ixLmax2,wLC)
         call conserve(ixRmin1,ixRmin2,ixRmax1,ixRmax2,wRC)
      endif

      return
      end

*=============================================================================
* end module vactvdlf
*##############################################################################
