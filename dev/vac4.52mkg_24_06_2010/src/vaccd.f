*#############################################################################
* module vaccd
* Centered difference scheme
*=============================================================================
      subroutine centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)

* Advance the iws flow variables from t to t+qdt within ixO^L by centered 
* differencing in space the dw/dt+dF_i(w)/dx_i=S type equation. 
* wCT contains the time centered variables at time qtC for flux and source.
* w is the old value at qt on input and the new value at qt+qdt on output.

      include 'vacdef.f'

      double precision  qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idimmin,idimmax
      logical   transport,tvdpreproc

      double precision  v(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ndim),f(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  iiw,iw,ixmin1,ixmin2,ixmax1,ixmax2,hxOmin1,hxOmin2,hxOmax1,
     &   hxOmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,
     &   jxCmax2,idim,idir

      logical  firstsweep,lastsweep
      common/first/firstsweep,lastsweep
*-----------------------------------------------------------------------------

      oktest=index(teststr,'centdiff').ge.1
      if(oktest)write(*,*)'CentDiff wCT,w:',wCT(ixtest1,ixtest2,iwtest),
     &   w(ixtest1,ixtest2,iwtest)

* An extra layer is needed in each direction for which fluxes are added.
      ixmin1=ixOmin1
      ixmin2=ixOmin2
      ixmax1=ixOmax1
      ixmax2=ixOmax2
      do idim= idimmin,idimmax
         ixmin1=ixmin1-kr(idim,1)
         ixmin2=ixmin2-kr(idim,2)
         ixmax1=ixmax1+kr(idim,1)
         ixmax2=ixmax2+kr(idim,2)
      enddo
      if(ixImin1.gt.ixmin1.or.ixImin2.gt.ixmin2.or.ixImax1.lt.ixmax1.or.
     &   ixImax2.lt.ixmax2) call die( 
     &   'Error in CentDiff: Non-conforming input limits')

      tvdpreproc=(typefull1.eq.'tvd').and.(typeadvance.eq.'onestep').and.
     &   (typeaxial.eq.'cylinder'.or.(typeaxial.eq.'slab'.and.sourceunsplit))

      if(index(teststr,'noryuprep').ge.1)tvdpreproc=.false.

      if(oktest.and.tvdpreproc)write(*,*)'TVD preproc is done:',firstsweep

      if(tvdpreproc.and.firstsweep)then
*        Add half of curveture and source terms following Ryu, ApJ 452,364, eq.(A9)
         if(typeaxial.eq.'cylinder')then
*           Apply curvature correction
*           wCT=w-(dt/2)*F_r(w)/r
*           In case of angmomfix
*           wCT=w-dt*F_r(w)/r
*           We use fC for velocity in the r_ direction
            call getv(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,r_,fC)
            do iiw=1,iws(niw_)
             iw=iws(iiw)
               call getflux(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,iw,r_,f,
     &            transport)
               if(transport)then
                  do ix_2=ixMmin2,ixMmax2
                  do ix_1=ixMmin1,ixMmax1
                     f(ix_1,ix_2)=f(ix_1,ix_2)+fC(ix_1,ix_2)*w(ix_1,ix_2,iw)
                  enddo
                  enddo
               endif
               if(angmomfix.and.iw.eq.mphi_)then
                  do ix_2=ixMmin2,ixMmax2
                  do ix_1=ixMmin1,ixMmax1
                     wCT(ix_1,ix_2,iw)=wCT(ix_1,ix_2,iw)-qdt*
     &                  f(ix_1,ix_2)/x(ix_1,ix_2,r_)
                  enddo
                  enddo
               else
                  do ix_2=ixMmin2,ixMmax2
                  do ix_1=ixMmin1,ixMmax1
                     wCT(ix_1,ix_2,iw)=wCT(ix_1,ix_2,iw)-(qdt/2)*
     &                  f(ix_1,ix_2)/x(ix_1,ix_2,r_)
                  enddo
                  enddo
               endif
            enddo
*           Add half of geometrical sources as well
            call addgeometry(qdt/2,ixMmin1,ixMmin2,ixMmax1,ixMmax2,iws,w,wCT)
         endif
*        Add half of unsplit sources
         if(sourceunsplit)call addsource2(qdt/2,ixGmin1,ixGmin2,ixGmax1,
     &      ixGmax2,ixMmin1,ixMmin2,ixMmax1,ixMmax2,iws,qt,w,qt,wCT)

         call getboundary(qt+qdt/2,1,nw,1,ndim,wCT)
*!! This is just a partial step, so processing is probably not useful here
*   if(nproc(2)/=0.and.implpar<=zero)then
*      if(nproc(2)<0.or.(it-itmin==((it-itmin)/nproc(2))*nproc(2)))then
*         call process(2,idim^LIM,wCT)
*      endif
*!! endif
      endif

* Add fluxes to w
      do idim= idimmin,idimmax
         ixmin1=ixOmin1-kr(idim,1)
         ixmin2=ixOmin2-kr(idim,2)
         ixmax1=ixOmax1+kr(idim,1)
         ixmax2=ixOmax2+kr(idim,2)
          ixCmin1=ixmin1
         ixCmin2=ixmin2
          ixCmax1=ixOmax1
         ixCmax2=ixOmax2
*        SHIFT
         jxCmin1=ixCmin1+kr(idim,1)
         jxCmin2=ixCmin2+kr(idim,2)
         jxCmax1=ixCmax1+kr(idim,1)
         jxCmax2=ixCmax2+kr(idim,2)
          
         hxOmin1=ixOmin1-kr(idim,1)
         hxOmin2=ixOmin2-kr(idim,2)
         hxOmax1=ixOmax1-kr(idim,1)
         hxOmax2=ixOmax2-kr(idim,2)

         if(gencoord)then
            do idir=1,ndim
               call getv(wCT,ixmin1,ixmin2,ixmax1,ixmax2,idir,f)
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  v(ix_1,ix_2,idir)=f(ix_1,ix_2)
               enddo
               enddo
            enddo
         else
            call getv(wCT,ixmin1,ixmin2,ixmax1,ixmax2,idim,f)
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               v(ix_1,ix_2,idim)=f(ix_1,ix_2)
            enddo
            enddo
         endif
         do iiw=1,iws(niw_)
          iw=iws(iiw)
            if(gencoord)then
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  fC(ix_1,ix_2)=zero
               enddo
               enddo
               do idir=1,ndim
                  call getflux(wCT,ixmin1,ixmin2,ixmax1,ixmax2,iw,idir,f,
     &               transport)
*                 Add transport flux
                  if(transport)then
                     do ix_2=ixmin2,ixmax2
                     do ix_1=ixmin1,ixmax1
                        f(ix_1,ix_2)=f(ix_1,ix_2)+v(ix_1,ix_2,idir)*
     &                     wCT(ix_1,ix_2,iw)
                     enddo
                     enddo
                  endif
*                 Center flux to interface 
*                 SHIFT BEGIN
                  do ix_2=ixCmin2,ixCmax2
                  do ix_1=ixCmin1,ixCmax1
                     tmp(ix_1,ix_2)=half*(f(ix_1,ix_2)+f(ix_1+
     &                  (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2)))
                  enddo
                  enddo
*                 SHIFT END

*                 Store flux for fluxCT or fluxCD schemes
*                 Check idir==idim since we do not project on the normal vector
                  if((typeconstrain.eq.'fluxCT'.or.typeconstrain.eq.
     &               'fluxCD').and.idim.eq.idir.and.iw.ne.b0_+
     &               idim .and.iw.gt.b0_.and.iw.le.b0_+ndim.and.
     &               istep.eq.nstep) call storeflux(qdt,tmp,ixCmin1,ixCmin2,
     &               ixCmax1,ixCmax2,idim,iw)

*                 Project to normal vector
                  do ix_2=ixCmin2,ixCmax2
                  do ix_1=ixCmin1,ixCmax1
                     fC(ix_1,ix_2)=fC(ix_1,ix_2)+normalC(ix_1,ix_2,idim,idir)*
     &                  tmp(ix_1,ix_2)
                  enddo
                  enddo
               enddo
            else
*              Get non-transported flux
               call getflux(wCT,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,f,
     &            transport)
               if(oktest.and.iw.eq.iwtest)write(*,*)'  fj,fi,fh:',f(ixtest1+
     &            kr(idim,1),ixtest2+kr(idim,2)),f(ixtest1,ixtest2),f(ixtest1-
     &            kr(idim,1),ixtest2-kr(idim,2))
*              Add transport flux
               if(transport)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     f(ix_1,ix_2)=f(ix_1,ix_2)+v(ix_1,ix_2,idim)*
     &                  wCT(ix_1,ix_2,iw)
                  enddo
                  enddo
               endif
               if(oktest.and.iw.eq.iwtest)write(*,*)'v+fj,fi,fh:',f(ixtest1+
     &            kr(idim,1),ixtest2+kr(idim,2)),f(ixtest1,ixtest2),f(ixtest1-
     &            kr(idim,1),ixtest2-kr(idim,2))
*              Center flux to interface
*              SHIFT BEGIN
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  fC(ix_1,ix_2)=half*(f(ix_1,ix_2)+f(ix_1+(jxCmin1-
     &               ixCmin1),ix_2+(jxCmin2-ixCmin2)))
               enddo
               enddo
*              SHIFT END
            endif

            call addflux(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,idim,fC,
     &         ixOmin1,ixOmin2,ixOmax1,ixOmax2,fC,hxOmin1,hxOmin2,hxOmax1,
     &         hxOmax2,w)

            if(oktest.and.iw.eq.iwtest)write(*,*)'fCi,fCh:',fC(ixtest1,
     &         ixtest2),fC(ixtest1-kr(idim,1),ixtest2-kr(idim,2))
            if(oktest.and.iw.eq.iwtest)write(*,*)'idim,wnew:',idim,w(ixtest1,
     &         ixtest2,iw)
*        next iw
         end do
*     next idim
      end do

      if(.not.tvdpreproc)then
         if(typeaxial.ne.'slab'.and.idimmin.eq.r_) call addgeometry(qdt,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,wCT,w)
         if(sourceunsplit) call addsource2(qdt*(idimmax-idimmin+one)/
     &      ndim, ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &      ixOmax2,iws,qtC,wCT,qt,w)
      endif

      return
      end

*=============================================================================
      subroutine centdiff4(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)

* Advance the iws flow variables from t to t+qdt within ixO^L by 
* fourth order centered  differencing in space the dw/dt+dF_i(w)/dx_i=S 
* type equation. 
* wCT contains the time centered variables at time qtC for flux and source.
* w is the old value at qt on input and the new value at qt+qdt on output.

      include 'vacdef.f'

      double precision  qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idimmin,idimmax
      logical   transport

      double precision  v(ixGlo1:ixGhi1,ixGlo2:ixGhi2),f(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  iiw,iw,ixmin1,ixmin2,ixmax1,ixmax2,idim,idir
*-----------------------------------------------------------------------------

      oktest=index(teststr,'centdiff').ge.1
      if(oktest)write(*,*)'CentDiff4 wCT,w:',wCT(ixtest1,ixtest2,iwtest),
     &   w(ixtest1,ixtest2,iwtest)

      if(gencoord)call die('CentDiff4 is not implemented for gen.coords yet!')

* Two extra layers are needed in each direction for which fluxes are added.
      ixmin1=ixOmin1
      ixmin2=ixOmin2
      ixmax1=ixOmax1
      ixmax2=ixOmax2
      do idim= idimmin,idimmax
         ixmin1=ixmin1-2*kr(idim,1)
         ixmin2=ixmin2-2*kr(idim,2)
         ixmax1=ixmax1+2*kr(idim,1)
         ixmax2=ixmax2+2*kr(idim,2)
      enddo
      if(ixImin1.gt.ixmin1.or.ixImin2.gt.ixmin2.or.ixImax1.lt.ixmax1.or.
     &   ixImax2.lt.ixmax2) call die( 
     &   'Error in CentDiff4: Non-conforming input limits')

* Add fluxes to w
      do idim= idimmin,idimmax
         ixmin1=ixOmin1-2*kr(idim,1)
         ixmin2=ixOmin2-2*kr(idim,2)
         ixmax1=ixOmax1+2*kr(idim,1)
         ixmax2=ixOmax2+2*kr(idim,2)

         call getv(wCT,ixmin1,ixmin2,ixmax1,ixmax2,idim,v)

         do iiw=1,iws(niw_)
          iw=iws(iiw)

*           Get non-transported flux
            call getflux(wCT,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,f,transport)

            if(oktest.and.iw.eq.iwtest)write(*,*)'  fj,fi,fh:',f(ixtest1+
     &         kr(idim,1),ixtest2+kr(idim,2)),f(ixtest1,ixtest2),f(ixtest1-
     &         kr(idim,1),ixtest2-kr(idim,2))

*           Add transport flux
            if(transport)then
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)=f(ix_1,ix_2)+v(ix_1,ix_2)*wCT(ix_1,ix_2,iw)
               enddo
               enddo
            endif
            if(oktest.and.iw.eq.iwtest)write(*,*)'v+fj,fi,fh:',f(ixtest1+
     &         kr(idim,1),ixtest2+kr(idim,2)),f(ixtest1,ixtest2),f(ixtest1-
     &         kr(idim,1),ixtest2-kr(idim,2))

*           Add divergence of flux
            call gradient4(.false.,f,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,tmp)
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               w(ix_1,ix_2,iw)=w(ix_1,ix_2,iw)-qdt*tmp(ix_1,ix_2)
            enddo
            enddo

            if(oktest.and.iw.eq.iwtest)write(*,*)'fk,fj,fh,fg:',f(ixtest1+2*
     &         kr(idim,1),ixtest2+2*kr(idim,2)),f(ixtest1+kr(idim,1),ixtest2+
     &         kr(idim,2)),f(ixtest1-kr(idim,1),ixtest2-kr(idim,2)),f(ixtest1-
     &         2*kr(idim,1),ixtest2-2*kr(idim,2))
            if(oktest.and.iw.eq.iwtest)write(*,*)'idim,wnew:',idim,w(ixtest1,
     &         ixtest2,iw)

*        next iw
         end do
*     next idim
      end do

* Add sources
      if(typeaxial.ne.'slab'.and.idimmin.eq.r_) call addgeometry(qdt,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,wCT,w)
      if(sourceunsplit) call addsource2(qdt*(idimmax-idimmin+one)/
     &   ndim, ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws,qtC,wCT,qt,w)

      return
      end

*=============================================================================
* end module vaccd
*#############################################################################
