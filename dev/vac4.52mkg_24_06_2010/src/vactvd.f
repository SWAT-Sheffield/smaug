*##############################################################################
* module vactvd
* Subroutines for TVD, TVD-MacCormack and TVD-MUSCL schemes
*=============================================================================
      subroutine tvdlimit(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,idimmin,idimmax,w,qt,wnew)

* Limit the iws flow variables in wnew. 
* Limiting is based on w or wnew depending of "typelimited".
* Create the left and right values by shifting then call tvdlimit2.
* "method" can be temporally second order 'tvd' or first order 'tvd1'
* For the 2nd order 'tvd' source terms are also added to second order accuracy

      include 'vacdef.f'

      character*10   method
      double precision  qdt,qt
      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idimmin,idimmax
      double precision w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*!! wL
      integer  ixICmin1,ixICmin2,ixICmax1,ixICmax2,jxICmin1,jxICmin2,jxICmax1,
     &   jxICmax2,iiw,iw,idim

      logical  firstsweep,lastsweep
      common/first/firstsweep,lastsweep
*-----------------------------------------------------------------------------

      oktest=index(teststr,'tvd').ge.1
      if(oktest)write(*,*)'TVDLimit w,wnew,method:',w(ixtest1,ixtest2,iwtest),
     &   wnew(ixtest1,ixtest2,iwtest),method

      if(typelimited.eq.'predictor')then
*        !!Side effect, w is overwritten, but after TVDlimit w is not needed anymore
         call getboundary(qt,1,nw,idimmin,idimmax,wnew)
         do iiw=1,iws(niw_)
          iw=iws(iiw)
            do ix_2=ixImin2,ixImax2
            do ix_1=ixImin1,ixImax1
               w(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)
            enddo
            enddo
         end do
      end if

*!!Side effect, w is overwritten
*!!if(typelimited=='previous')w(ixI^S,1:nw)=wold(ixI^S,1:nw)

      do idim= idimmin,idimmax
         ixICmax1=ixOmax1+kr(idim,1)
         ixICmax2=ixOmax2+kr(idim,2)
         ixICmin1=ixOmin1-2*kr(idim,1)
         ixICmin2=ixOmin2-2*kr(idim,2)

*        !! wL(ixIC^S,1:nw)=w(ixIC^S,1:nw)
*        SHIFT
         jxICmin1=ixICmin1+kr(idim,1)
         jxICmin2=ixICmin2+kr(idim,2)
         jxICmax1=ixICmax1+kr(idim,1)
         jxICmax2=ixICmax2+kr(idim,2)
*        SHIFT BEGIN
         do iw_3=1,nw
         do ix_2=ixICmin2,ixICmax2
         do ix_1=ixICmin1,ixICmax1
            wR(ix_1,ix_2,iw_3)=w(ix_1+(jxICmin1-ixICmin1),ix_2+
     &         (jxICmin2-ixICmin2),iw_3)
         enddo
         enddo
         enddo
*        SHIFT END
         
*           !!w-->wL
         call tvdlimit2(method,qdt,ixICmin1,ixICmin2,ixICmax1,ixICmax2,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,idim,w,wR,wnew)
      enddo

      if(method.eq.'tvd'.and.lastsweep.and.index(teststr,'noryuprep').lt.
     &   1.and.(typeaxial.eq.'cylinder'.or.(typeaxial.eq.'slab'.and.
     &   sourceunsplit)))then

         if(oktest)write(*,*)'Adding Ryu style (geometrical) sources'

*        Following Ryu et al (APJ 452,364) (geometrical) sources are added 
*        with second order accuracy:
*        whalf=(wold+whydro+S)/2     eq.A22 corrected by including S
*        wnew=whydro+S(whalf)        eq.A21

*        !!Side effect: We overwrite w, but that is OK, it is not needed anymore
         do iw_3=1,nw
         do ix_2=ixOmin2,ixOmax2
         do ix_1=ixOmin1,ixOmax1
            w(ix_1,ix_2,iw_3)=wnew(ix_1,ix_2,iw_3)
         enddo
         enddo
         enddo

         if(typeaxial.ne.'slab')call addgeometry(qdt,ixOmin1,ixOmin2,ixOmax1,
     &      ixOmax2,iws,wold,w)
         if(sourceunsplit)call addsource2(qdt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,t+qdt,wold,t+qdt,w)

         do iw_3=1,nw
         do ix_2=ixOmin2,ixOmax2
         do ix_1=ixOmin1,ixOmax1
            w(ix_1,ix_2,iw_3)=half*(wold(ix_1,ix_2,iw_3)+w(ix_1,ix_2,iw_3))
         enddo
         enddo
         enddo

         if(typeaxial.ne.'slab')call addgeometry(qdt,ixOmin1,ixOmin2,ixOmax1,
     &      ixOmax2,iws,w,wnew)
         if(sourceunsplit)call addsource2(qdt,ixOmin1,ixOmin2,ixOmax1,ixOmax2,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,t+qdt,w,t+qdt,wnew)
      endif

      return
      end

*=============================================================================
      subroutine tvdlimit2(method,qdt,ixICmin1,ixICmin2,ixICmax1,ixICmax2,
     &   ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,idim,wL,wR,wnew)

* Limit the iws flow variables in wnew according to typetvd. 
* wroeC is based on wL and wR.
* If method=='tvd' an extra adtdx**2*jumpC is added to phiC for 2nd order
* accuracy in time.

      include 'vacdef.f'

      character*10   method
      double precision  qdt
      integer  ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idim
      double precision wL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wR(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wroeC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision phiC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rphiC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),jumpC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   adtdxC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),smallaC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      double precision  courantmax,wtest(nw)
      integer  hxOmin1,hxOmin2,hxOmax1,hxOmax2,ixCmin1,ixCmin2,ixCmax1,
     &   ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,jxICmin1,jxICmin2,jxICmax1,
     &   jxICmax2,iiw,iw,il
*-----------------------------------------------------------------------------

      oktest=index(teststr,'tvdlimit').ge.1
      if(oktest)write(*,*)'TVDLimit2 wL,wR,wnew:',wL(ixtest1,ixtest2,iwtest),
     &   wR(ixtest1,ixtest2,iwtest),wnew(ixtest1,ixtest2,iwtest)
      if(typetvd.eq.'testing')then
         do iw_1=1,nw
            wtest(iw_1)=0.
         enddo
      endif

      hxOmin1=ixOmin1-kr(idim,1)
      hxOmin2=ixOmin2-kr(idim,2)
      hxOmax1=ixOmax1-kr(idim,1)
      hxOmax2=ixOmax2-kr(idim,2)
      ixCmax1=ixOmax1
      ixCmax2=ixOmax2
       ixCmin1=hxOmin1
      ixCmin2=hxOmin2
       

      jxCmin1=ixCmin1+kr(idim,1)
      jxCmin2=ixCmin2+kr(idim,2)
      jxCmax1=ixCmax1+kr(idim,1)
      jxCmax2=ixCmax2+kr(idim,2)
      jxICmin1=ixICmin1+kr(idim,1)
      jxICmin2=ixICmin2+kr(idim,2)
      jxICmax1=ixICmax1+kr(idim,1)
      jxICmax2=ixICmax2+kr(idim,2)

      call average(wL,wR,ixICmin1,ixICmin2,ixICmax1,ixICmax2,iws,idim,wroeC)

* A loop on characteristic variables to calculate the dissipative flux phiC.
      do il=1,nw
*        Calculate the jump in the il-th characteristic variable: L(wroe)*dw
         call geteigenjump(wL,wR,wroeC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,
     &      idim,smallaC,adtdxC,jumpC)

         if(oktest)write(*,*)'wL(ix,il),wR(ix,il):',wL(ixtest1,ixtest2,il),
     &      wR(ixtest1,ixtest2,il)
         if(oktest)write(*,*)'L*dwC,a:',jumpC(ixtest1,ixtest2),adtdxC(ixtest1,
     &      ixtest2)

*        Normalize the eigenvalue "a" (and its limit "smalla" if needed):
         if(gencoord)then
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               adtdxC(ix_1,ix_2)=adtdxC(ix_1,ix_2)*qdt*surfaceC(ix_1,ix_2,
     &            idim)*two/(dvolume(ix_1,ix_2)+dvolume(ix_1+
     &            (jxICmin1-ixICmin1),ix_2+(jxICmin2-ixICmin2)))
            enddo
            enddo
            if(typeentropy(il).eq.'harten' .or. typeentropy(il).eq.
     &         'powell')then
               do ix_2=ixICmin2,ixICmax2
               do ix_1=ixICmin1,ixICmax1
                  smallaC(ix_1,ix_2)=smallaC(ix_1,ix_2)*qdt*
     &               surfaceC(ix_1,ix_2,idim)*two/(dvolume(ix_1,ix_2)+
     &               dvolume(ix_1+(jxICmin1-ixICmin1),ix_2+
     &               (jxICmin2-ixICmin2)))
               enddo
               enddo
            endif
         else
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               adtdxC(ix_1,ix_2)=adtdxC(ix_1,ix_2)*qdt*two/
     &            (dx(ix_1,ix_2,idim)+dx(ix_1+(jxICmin1-ixICmin1),ix_2+
     &            (jxICmin2-ixICmin2),idim))
            enddo
            enddo
            if(typeentropy(il).eq.'harten' .or. typeentropy(il).eq.
     &         'powell')then
               do ix_2=ixICmin2,ixICmax2
               do ix_1=ixICmin1,ixICmax1
                  smallaC(ix_1,ix_2)=smallaC(ix_1,ix_2)*qdt*
     &               two/(dx(ix_1,ix_2,idim)+dx(ix_1+(jxICmin1-ixICmin1),ix_2+
     &               (jxICmin2-ixICmin2),idim))
               enddo
               enddo
            endif
         endif

         if(oktest)write(*,*)'qdt,adtdxC:',qdt,adtdxC(ixtest1,ixtest2)

*        For potentially extreme eigenvalues calculate dtcourant for next step
         if((il.eq.extremeLW_.or.il.eq.extremeRW_).and.courantpar.gt.zero.and.
     &      istep.eq.nstep.and.implpar.lt.zero)then
            maxval_1=-bigdouble
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               maxval_1=max(maxval_1,abs(adtdxC(ix_1,ix_2)))
            enddo
            enddo
            courantmax=maxval_1
            
            if(courantmax.gt.one) then
               nerror(couranterr_)=nerror(couranterr_)+1
               if(nerror(couranterr_).eq.1)write(*,*)
     &            'Courant condition error (code=',couranterr_,') at it=',it,
     &            ' for direction',idim,' CFL=',courantmax
            endif
            if(courantmax.gt.smalldouble)dtcourant(idim)=min(dtcourant(idim),
     &         qdt*courantpar/courantmax)
         endif

*        Calculate the flux limiter function phi
         call getphi(method,jumpC,adtdxC,smallaC,ixICmin1,ixICmin2,ixICmax1,
     &      ixICmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,il,idim,phiC)

         if(oktest)write(*,*)'il,phiC:',il,phiC(ixtest1,ixtest2)

*        Multiply phiC by 2>acmcoef>0 (kappa in eq.2.21 Yee, Sandham, Djomehri)
         if(acmcoef(il).ge.zero)then
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               phiC(ix_1,ix_2)=acmcoef(il)*phiC(ix_1,ix_2)
            enddo
            enddo
         endif

*        Multiply by a term of Artificial Compression Method (ACM)
*        (theta_j+1/2 in eq.2.21 Yee, Sandham, Djomehri) if axmexpo>0 is set

         if(acmexpo.gt.zero)call acmswitch(jumpC,ixCmin1,ixCmin2,ixCmax1,
     &      ixCmax2,idim,phiC)

         if(oktest)write(*,*)'phiCstar:',phiC(ixtest1,ixtest2)

*        Add R(iw,il)*phiC(il) to each variable iw in wnew
         do iiw=1,iws(niw_)
          iw=iws(iiw)
            call rtimes(phiC,wroeC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,il,idim,
     &         rphiC)

*           Based on Eqs 5.6a,b in Yee
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               rphiC(ix_1,ix_2)=rphiC(ix_1,ix_2)*quarter*(dvolume(ix_1+
     &            (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2))+
     &            dvolume(ix_1,ix_2))
            enddo
            enddo
            if(gencoord.and.vectoriw(iw).ge.0)then
*              The dt=-one tells addflux_rotate that this is an artificial flux
               
                      call die('Error: gencoord is off')
            else
*              The dt=-one tells addflux that this is an artificial flux
               call addflux(-one,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,idim,
     &            rphiC,ixOmin1,ixOmin2,ixOmax1,ixOmax2,rphiC,hxOmin1,hxOmin2,
     &            hxOmax1,hxOmax2,wnew)
            endif

            if(oktest)write(*,*)'rphiCi,rphiCh:',rphiC(ixtest1,ixtest2)/
     &         dvolume(ixtest1,ixtest2),rphiC(ixtest1-kr(idim,1),ixtest2-
     &         kr(idim,2))/dvolume(ixtest1,ixtest2)
            if(typetvd.eq.'testing')wtest(iw)=wtest(iw)+rphiC(ixtest1,ixtest2)
*        iw
         end do
         if(oktest)write(*,*)'wnew:',wnew(ixtest1,ixtest2,iwtest)
*     il
      end do

      if(typetvd.eq.'testing')then
         write(*,*)'wtest                   wR-wL'
         do iiw=1,iws(niw_)
          iw=iws(iiw)
            write(*,*)wtest(iw),wR(ixtest1,ixtest2,iw)-wL(ixtest1,ixtest2,iw)
         end do
      endif

      return
      end

*=============================================================================
      subroutine getphi(method,jumpC,adtdxC,smallaC,ixICmin1,ixICmin2,
     &   ixICmax1,ixICmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,il,idim,phiC)

* Calculate the dissipative flux from jumpC=L*dw and adtdx=eigenvalue*dt/dx.
* Add Lax-Wendroff type correction if method=='tvd'.
* Limit according to method and typetvd.

      include 'vacdef.f'

      character*10   method
      integer  ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixCmin1,ixCmin2,ixCmax1,
     &   ixCmax2
      double precision jumpC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),ljumpC(
     &   ixGlo1:ixGhi1,ixGlo2:ixGhi2),adtdxC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   smallaC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),phiC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  jxCmin1,jxCmin2,jxCmax1,jxCmax2,ixmin1,ixmin2,ixmax1,ixmax2,
     &   hxmin1,hxmin2,hxmax1,hxmax2,il,idim
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getphi').ge.1
      if(oktest)write(*,*)'Getphi jumpC,adtdxC:',jumpC(ixtest1,ixtest2),
     &   adtdxC(ixtest1,ixtest2)

      if(method.eq.'tvdmu'.or.method.eq.'tvdmu1')then
*        In the MUSCL scheme phi=|a|*jump, apply entropy fix to it
         if(typeentropy(il).eq.'nul'.or.typeentropy(il).eq.'ratio')then
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               phiC(ix_1,ix_2)=abs(adtdxC(ix_1,ix_2))*jumpC(ix_1,ix_2)
            enddo
            enddo
         else
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               if(abs(adtdxC(ix_1,ix_2)).ge.smallaC(ix_1,ix_2))then
                        phiC(ix_1,ix_2)=abs(adtdxC(ix_1,ix_2))*
     &                     jumpC(ix_1,ix_2)
                     else
                        phiC(ix_1,ix_2)=half*(smallaC(ix_1,ix_2)+
     &                     adtdxC(ix_1,ix_2)**2/smallaC(ix_1,ix_2))*
     &                     jumpC(ix_1,ix_2)
                     endif
            enddo
            enddo
         endif
*        That's all for the MUSCL scheme
         if(oktest)write(*,*)'GetPhi Final phiC:',phiC(ixtest1,ixtest2)
         return
      endif

      if(method.eq.'tvd')then
*        Entropy fix to |a|-a**2
         if(typeentropy(il).eq.'nul'.or.typeentropy(il).eq.'ratio')then
           do ix_2=ixICmin2,ixICmax2
           do ix_1=ixICmin1,ixICmax1
              phiC(ix_1,ix_2)=abs(adtdxC(ix_1,ix_2))-adtdxC(ix_1,ix_2)**2
           enddo
           enddo
         else
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               if(abs(adtdxC(ix_1,ix_2)).ge.smallaC(ix_1,ix_2))then
                        phiC(ix_1,ix_2)=abs(adtdxC(ix_1,ix_2))-
     &                     adtdxC(ix_1,ix_2)**2
                     else
                        phiC(ix_1,ix_2)=half*smallaC(ix_1,ix_2)+(half/
     &                     smallaC(ix_1,ix_2)-one)*adtdxC(ix_1,ix_2)**2
                     endif
            enddo
            enddo
         endif
         if(oktest)write(*,*)'abs(nu)-nu**2:',phiC(ixtest1,ixtest2)
      else
*        Entropy fix to |a|
         if(typeentropy(il).eq.'nul'.or.typeentropy(il).eq.'ratio')then
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               phiC(ix_1,ix_2)=abs(adtdxC(ix_1,ix_2))
            enddo
            enddo
         else
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               if(abs(adtdxC(ix_1,ix_2)).ge.smallaC(ix_1,ix_2))then
                        phiC(ix_1,ix_2)=abs(adtdxC(ix_1,ix_2))
                     else
                        phiC(ix_1,ix_2)=half*smallaC(ix_1,ix_2)+half/
     &                     smallaC(ix_1,ix_2)*adtdxC(ix_1,ix_2)**2
                     endif
            enddo
            enddo
         endif
         if(oktest)write(*,*)'abs(nu)-nu**2:',phiC(ixtest1,ixtest2)
      endif

*SHIFT
      jxCmin1=ixCmin1+kr(idim,1)
      jxCmin2=ixCmin2+kr(idim,2)
      jxCmax1=ixCmax1+kr(idim,1)
      jxCmax2=ixCmax2+kr(idim,2)

      hxmin1=ixICmin1
      hxmin2=ixICmin2
       hxmax1=ixICmax1-kr(idim,1)
      hxmax2=ixICmax2-kr(idim,2)
*SHIFT MORE
      ixmin1=hxmin1+kr(idim,1)
      ixmin2=hxmin2+kr(idim,2)
      ixmax1=hxmax1+kr(idim,1)
      ixmax2=hxmax2+kr(idim,2)

*SHIFT BEGIN
      if(typetvd.eq.'symmetric')then
*           eq.3.53 and eq.3.69d
            call dwlimiter3(jumpC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,
     &         ljumpC)
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               phiC(ix_1,ix_2)=phiC(ix_1,ix_2)*(jumpC(ix_1,ix_2)-
     &            ljumpC(ix_1,ix_2))
            enddo
            enddo
*           extra (a*lambda)**2*delta
            if(method.eq.'tvd')then
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=phiC(ix_1,ix_2)+adtdxC(ix_1,ix_2)**
     &               2*jumpC(ix_1,ix_2)
               enddo
               enddo
            endif
      else if(typetvd.eq.'roe')then
*           eq.3.52 with correction of sign of a
         if(typelimiter(il).eq.'roe'.or.typelimiter(il).eq.'superroe')then
               call dwlimiterroe(adtdxC,jumpC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,
     &            il,idim,ljumpC)
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=phiC(ix_1,ix_2)*(jumpC(ix_1,ix_2)-
     &               ljumpC(ix_1,ix_2))
               enddo
               enddo
         else
               call dwlimiter2(jumpC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,
     &            idim,ljumpC)
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  if(adtdxC(ix_1,ix_2).le.0)then
                              phiC(ix_1,ix_2)=phiC(ix_1,ix_2)*
     &                           (jumpC(ix_1,ix_2)-ljumpC(ix_1+
     &                           (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2)))
                           else
                              phiC(ix_1,ix_2)=phiC(ix_1,ix_2)*
     &                           (jumpC(ix_1,ix_2)-ljumpC(ix_1,ix_2))
                           endif
               enddo
               enddo
            end if
*           extra (a*lambda)**2*delta
            if(method.eq.'tvd')then
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=phiC(ix_1,ix_2)+adtdxC(ix_1,ix_2)**
     &               2*jumpC(ix_1,ix_2)
               enddo
               enddo
            endif
      else if(typetvd.eq.'sweby')then
*           Sweby eqs.4.11-4.15, but no 0.5 ?!
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               phiC(ix_1,ix_2)=phiC(ix_1,ix_2)*jumpC(ix_1,ix_2)
            enddo
            enddo
         if(typelimiter(il).eq.'roe'.or.typelimiter(il).eq.'superroe')then
               call dwlimiterroe(adtdxC,phiC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,
     &            il,idim,ljumpC)
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=phiC(ix_1,ix_2)-ljumpC(ix_1,ix_2)
               enddo
               enddo
         else
               call dwlimiter2(phiC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,
     &            idim,ljumpC)
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  if(adtdxC(ix_1,ix_2).le.0)then
                              phiC(ix_1,ix_2)=phiC(ix_1,ix_2)-ljumpC(ix_1+
     &                           (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2))
                           else
                              phiC(ix_1,ix_2)=phiC(ix_1,ix_2)-
     &                           ljumpC(ix_1,ix_2)
                           endif
               enddo
               enddo
            end if
*           extra (a*lambda)**2*delta
            if(method.eq.'tvd')then
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=phiC(ix_1,ix_2)+adtdxC(ix_1,ix_2)**
     &               2*jumpC(ix_1,ix_2)
               enddo
               enddo
            endif
      else if(typetvd.eq.'yee')then
*           eq.3.51 with correction
            call dwlimiter2(jumpC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,
     &         ljumpC)

            if(oktest)write(*,*)'Limiter jumpC,ljumpC(ix,jx):',jumpC(ixtest1,
     &         ixtest2),ljumpC(ixtest1,ixtest2),ljumpC(ixtest1+
     &         kr(idim,1),ixtest2+kr(idim,2))

*           Use phiC as 0.5*(|nu|-nu**2) eq.3.45e for tvd otherwise 0.5*|nu|
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               phiC(ix_1,ix_2)=half*phiC(ix_1,ix_2)
            enddo
            enddo
*           gamma*lambda eq.3.51d, use tmp to store agdtdxC
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               if(abs(jumpC(ix_1,ix_2)).gt.smalldouble)then
                        tmp(ix_1,ix_2)=adtdxC(ix_1,ix_2)+phiC(ix_1,ix_2)*
     &                     (ljumpC(ix_1+(jxCmin1-ixCmin1),ix_2+
     &                     (jxCmin2-ixCmin2))-ljumpC(ix_1,ix_2))/
     &                     jumpC(ix_1,ix_2)
                     else
                        tmp(ix_1,ix_2)=adtdxC(ix_1,ix_2)
                     endif
            enddo
            enddo

            if(oktest)write(*,*)'agdtdxC:',tmp(ixtest1,ixtest2)

*           eq.3.51a
            if(typeentropy(il).eq.'nul'.or.typeentropy(il).eq.'ratio')then
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=-phiC(ix_1,ix_2)*(ljumpC(ix_1+
     &               (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2))+
     &               ljumpC(ix_1,ix_2))+abs(tmp(ix_1,ix_2))*jumpC(ix_1,ix_2)
               enddo
               enddo
            else
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  if(abs(tmp(ix_1,ix_2)).ge.smallaC(ix_1,ix_2))then
                              phiC(ix_1,ix_2)=-phiC(ix_1,ix_2)*
     &                           (ljumpC(ix_1+(jxCmin1-ixCmin1),ix_2+
     &                           (jxCmin2-ixCmin2))+ljumpC(ix_1,ix_2))+
     &                           abs(tmp(ix_1,ix_2))*jumpC(ix_1,ix_2)
                           else
                              phiC(ix_1,ix_2)=-phiC(ix_1,ix_2)*
     &                           (ljumpC(ix_1+(jxCmin1-ixCmin1),ix_2+
     &                           (jxCmin2-ixCmin2))+ljumpC(ix_1,ix_2))+(half*
     &                           smallaC(ix_1,ix_2)+half/smallaC(ix_1,ix_2)*
     &                           tmp(ix_1,ix_2)**2)*jumpC(ix_1,ix_2)
                           endif
               enddo
               enddo
            endif
      else if(typetvd.eq.'harten')then
*           See Ryu, section 2.3
*           Use phiC as 0.5*(|nu|-nu**2)*jumpC eq.3.45b,e
            do ix_2=ixICmin2,ixICmax2
            do ix_1=ixICmin1,ixICmax1
               phiC(ix_1,ix_2)=half*phiC(ix_1,ix_2)*jumpC(ix_1,ix_2)
            enddo
            enddo
            call dwlimiter2(phiC,ixICmin1,ixICmin2,ixICmax1,ixICmax2,il,idim,
     &         ljumpC)

            if(artcomp(il))then
*              ARTIFICIAL COMPRESSION, ORIGINALLY IMPLEMENTED BY M. NAUTA
*              
*              Equations 5.8b, 5.7c, 5.8a, 5.7a of Harten change ljumpC,
*              or Ryu's  2.98, 2.97, 2.96:
*              
*              sigma=(1-a*|dt/dx|)/2
*              gbar_i=minmod(sigma*alfa_i-1/2,sigma*alfa_i+1/2)
*              theta_i=(|alfa_i+1/2 - alfa_i-1/2|)/(|alfa_i+1/2|+|alfa_i-1/2|)
*              ljumpC=ljumpC+theta*gbar
*              
*              To save memory tmp is used for sigma & theta, tmp2 is used for gbar

               do ix_2=ixICmin2,ixICmax2
               do ix_1=ixICmin1,ixICmax1
                  tmp(ix_1,ix_2)=half*(one-abs(adtdxC(ix_1,ix_2)))*
     &               jumpC(ix_1,ix_2)
               enddo
               enddo
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  tmp2(ix_1,ix_2)=sign(one,jumpC(ix_1,ix_2))*
     &               max(zero,min(sign(one,jumpC(ix_1,ix_2))*
     &               tmp(ix_1+(hxmin1-ixmin1),ix_2+(hxmin2-
     &               ixmin2)),abs(tmp(ix_1,ix_2))))
               enddo
               enddo
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  if(abs(jumpC(ix_1,ix_2))+abs(jumpC(ix_1+(hxmin1-
     &               ixmin1),ix_2+(hxmin2-ixmin2))).gt.smalldouble)then
                              tmp(ix_1,ix_2)=abs(jumpC(ix_1,ix_2)-jumpC(ix_1+
     &                           (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2)))/
     &                           (abs(jumpC(ix_1,ix_2))+abs(jumpC(ix_1+
     &                           (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2))))
                           else
                              tmp(ix_1,ix_2)=zero
                           endif
               enddo
               enddo
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  ljumpC(ix_1,ix_2)=ljumpC(ix_1,ix_2)+tmp(ix_1,ix_2)*
     &               tmp2(ix_1,ix_2)
               enddo
               enddo
            endif

*           gamma*lambda eq.3.45d, use tmp as agdtdxC
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               if(abs(jumpC(ix_1,ix_2)).gt.smalldouble)then
                        tmp(ix_1,ix_2)=adtdxC(ix_1,ix_2)+(ljumpC(ix_1+
     &                     (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2))-
     &                     ljumpC(ix_1,ix_2))/jumpC(ix_1,ix_2)
                     else
                        tmp(ix_1,ix_2)=adtdxC(ix_1,ix_2)
                     endif
            enddo
            enddo
            if(oktest)write(*,*)'agdtdxC',tmp(ixtest1,ixtest2)
*           eq.3.45a with correction
            if(typeentropy(il).eq.'nul'.or.typeentropy(il).eq.'ratio')then
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  phiC(ix_1,ix_2)=-ljumpC(ix_1+(jxCmin1-ixCmin1),ix_2+
     &               (jxCmin2-ixCmin2))-ljumpC(ix_1,ix_2)+jumpC(ix_1,ix_2)*
     &               abs(tmp(ix_1,ix_2))
               enddo
               enddo
            else
               do ix_2=ixCmin2,ixCmax2
               do ix_1=ixCmin1,ixCmax1
                  if(abs(tmp(ix_1,ix_2)).ge.smallaC(ix_1,ix_2))then
                              phiC(ix_1,ix_2)=-ljumpC(ix_1+
     &                           (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2))-
     &                           ljumpC(ix_1,ix_2)+jumpC(ix_1,ix_2)*
     &                           abs(tmp(ix_1,ix_2))
                           else
                              phiC(ix_1,ix_2)=-ljumpC(ix_1+
     &                           (jxCmin1-ixCmin1),ix_2+(jxCmin2-ixCmin2))-
     &                           ljumpC(ix_1,ix_2)+jumpC(ix_1,ix_2)*(half*
     &                           smallaC(ix_1,ix_2)+half/smallaC(ix_1,ix_2)*
     &                           tmp(ix_1,ix_2)**2)
                           endif
               enddo
               enddo
            endif
*           extra -(a*lambda)**2*delta
      else if(typetvd.eq.'testing')then
*           phiC:=jumpC, thus R*L*(w(jx)-w(ix))==w(jx)-w(ix) can be tested
*           Use tvdmc, otherwise the second order correction is applied
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               phiC(ix_1,ix_2)=jumpC(ix_1,ix_2)
            enddo
            enddo
      else
            call die('Error in TVDLimit: Unknown TVD type='//typetvd)
      end if
*SHIFT END

      if(oktest)write(*,*)'GetPhi Final phiC:',phiC(ixtest1,ixtest2)

      return
      end

*=============================================================================
      subroutine entropyfix(ixmin1,ixmin2,ixmax1,ixmax2,il,aL,aR,a,smalla)

* Apply entropyfix based on typeentropy(il),aL,aR, and a
* Calculate "smalla" (Harten,Powell) or modify "a" (ratio)
*!! tmp and tmp2 are not to be used in this subroutine

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,il
      double precision aL(ixGlo1:ixGhi1,ixGlo2:ixGhi2),aR(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),a(ixGlo1:ixGhi1,ixGlo2:ixGhi2),smalla(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

      if(typeentropy(il).eq.'harten')then
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            smalla(ix_1,ix_2)=max(zero,a(ix_1,ix_2)-aL(ix_1,ix_2),aR(ix_1,
     &         ix_2)-a(ix_1,ix_2))
         enddo
         enddo
      else if(typeentropy(il).eq.'powell')then
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            smalla(ix_1,ix_2)=max(zero,two*(aR(ix_1,ix_2)-aL(ix_1,ix_2)))
         enddo
         enddo
      else if(typeentropy(il).eq.'ratio')then
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            if(aL(ix_1,ix_2).lt.zero .and. aR(ix_1,ix_2).gt.zero)a(ix_1,ix_2)=
     &         a(ix_1,ix_2)-2*aR(ix_1,ix_2)*aL(ix_1,ix_2)/(aR(ix_1,ix_2)-
     &         aL(ix_1,ix_2))
         enddo
         enddo
      else if(typeentropy(il).eq.'yee')then
*        This has been done in geteigenjump already
      else if(typeentropy(il).eq.'nul')then
*        No entropyfix is applied
      else
         call die('No such type of entropy fix:'//typeentropy(il))
      end if

      return
      end

*=============================================================================
* end module vactvd
*##############################################################################
