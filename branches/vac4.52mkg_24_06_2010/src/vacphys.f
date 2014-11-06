*##############################################################################
* module vacphys - mhd

*##############################################################################
* module vacphys.mhd0 - common subroutines for mhd and mhdiso

*##############################################################################
* module vacphys.mhdroe - subroutines for Roe-type Riemann solver for MHD
*=============================================================================
      subroutine average(wL,wR,ixmin1,ixmin2,ixmax1,ixmax2,iws,idim,wroe)

* Eight-wave MHD Riemann solver. See Powell, Notes on the eigeinsystem, Gombosi
* Calculate the wroe average of primitive variables in wL and wR, assignment:
* rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
* Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
*
* wL,wR,wroe are all interface centered quantities

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim,iws(niw_),idir,jdir,iw
      double precision wL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wR(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision cfast(ixGlo1:ixGhi1,ixGlo2:ixGhi2),cslow(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),afast(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   aslow(ixGlo1:ixGhi1,ixGlo2:ixGhi2),csound2(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rhodv(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv

*-----------------------------------------------------------------------------

      if(ndir.eq.1)call die('MHD with d=11 is the same as HD')

      oktest=index(teststr,'average').ge.1
      if(oktest)write(*,*)'Average wL,wR:',wL(ixtest1,ixtest2,iwtest),
     &   wR(ixtest1,ixtest2,iwtest)

*Averaging primitive variables
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wroe(ix_1,ix_2,rho_)=half*(wL(ix_1,ix_2,rho_)+wR(ix_1,ix_2,rho_))
      enddo
      enddo

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wroe(ix_1,ix_2,v1_)=half*(wL(ix_1,ix_2,m1_)/wL(ix_1,ix_2,rho_)+
     &      wR(ix_1,ix_2,m1_)/wR(ix_1,ix_2,rho_))
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wroe(ix_1,ix_2,b1_)=half*(wL(ix_1,ix_2,b1_)+wR(ix_1,ix_2,b1_))
      enddo
      enddo


      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wroe(ix_1,ix_2,v2_)=half*(wL(ix_1,ix_2,m2_)/wL(ix_1,ix_2,rho_)+
     &      wR(ix_1,ix_2,m2_)/wR(ix_1,ix_2,rho_))
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wroe(ix_1,ix_2,b2_)=half*(wL(ix_1,ix_2,b2_)+wR(ix_1,ix_2,b2_))
      enddo
      enddo


* Use afast and aslow for pressures pL and pR
      call getpthermal(.true.,wL,ixmin1,ixmin2,ixmax1,ixmax2,afast)
      call getpthermal(.true.,wR,ixmin1,ixmin2,ixmax1,ixmax2,aslow)
      if(p_.gt.0)then
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wroe(ix_1,ix_2,pp_)=half*(afast(ix_1,ix_2)+aslow(ix_1,ix_2))
         enddo
         enddo
      endif

      if(oktest)write(*,*)'Calculate saved variables'

      if(useprimitive.or.p_.lt.1.or.eqpar(gamma_).le.zero)then
*        dp=pR-pL
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            dp(ix_1,ix_2)=aslow(ix_1,ix_2)-afast(ix_1,ix_2)
         enddo
         enddo
      else
*        CONSERVATIVE dp=(g-1)*(de-v*dm+0.5*v**2*drho-0.5*d(B**2))
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            dp(ix_1,ix_2)=(eqpar(gamma_)-1)*(wR(ix_1,ix_2,ee_)-
     &         wL(ix_1,ix_2,ee_)-(wroe(ix_1,ix_2,m1_)*(wR(ix_1,ix_2,m1_)-
     &         wL(ix_1,ix_2,m1_))+wroe(ix_1,ix_2,m2_)*(wR(ix_1,ix_2,m2_)-
     &         wL(ix_1,ix_2,m2_)))+half*(wroe(ix_1,ix_2,m1_)**
     &         2+wroe(ix_1,ix_2,m2_)**2)*(wR(ix_1,ix_2,rho_)-
     &         wL(ix_1,ix_2,rho_))-half*(wR(ix_1,ix_2,b1_)**
     &         2-wL(ix_1,ix_2,b1_)**2+wR(ix_1,ix_2,b2_)**2-wL(ix_1,ix_2,b2_)**
     &         2))
         enddo
         enddo
      endif

*CONSERVATIVE rho*dv_idim=dm_idim-v_idim*drho
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         rhodv(ix_1,ix_2)=wR(ix_1,ix_2,m0_+idim)-wL(ix_1,ix_2,m0_+
     &      idim)-wroe(ix_1,ix_2,m0_+idim)*(wR(ix_1,ix_2,rho_)-
     &      wL(ix_1,ix_2,rho_))
      enddo
      enddo

*Calculate csound2,cfast,cslow,alphafast and alphaslow

* get csound**2
      call getcsound2prim(wroe,ixmin1,ixmin2,ixmax1,ixmax2,csound2)

* aa=B**2/rho+a**2
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cfast(ix_1,ix_2)=( wroe(ix_1,ix_2,b1_)**2+wroe(ix_1,ix_2,b2_)**
     &      2 )/wroe(ix_1,ix_2,rho_)+csound2(ix_1,ix_2)
      enddo
      enddo

* cs**2=0.5*(aa+sqrt(aa**2-4*a**2*(b_i**2/rho)))
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cslow(ix_1,ix_2)=max(zero, half*(cfast(ix_1,ix_2)-
     &      sqrt(cfast(ix_1,ix_2)**2-4d0*csound2(ix_1,ix_2)*
     &      wroe(ix_1,ix_2,b0_+idim)**2/wroe(ix_1,ix_2,rho_))))
      enddo
      enddo

* cf**2=aa-cs**2
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cfast(ix_1,ix_2)=cfast(ix_1,ix_2)-cslow(ix_1,ix_2)
      enddo
      enddo

* alpha_f**2=(a**2-cs**2)/(cf**2-cs**2)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         afast(ix_1,ix_2)=(csound2(ix_1,ix_2)-cslow(ix_1,ix_2))/
     &      (cfast(ix_1,ix_2)-cslow(ix_1,ix_2))
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         afast(ix_1,ix_2)=min(one,max(afast(ix_1,ix_2),zero))
      enddo
      enddo

* alpha_s=sqrt(1-alpha_f**2)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         aslow(ix_1,ix_2)=sqrt(one-afast(ix_1,ix_2))
      enddo
      enddo

* alpha_f=sqrt(alpha_f**2)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         afast(ix_1,ix_2)=sqrt(afast(ix_1,ix_2))
      enddo
      enddo

* cf=sqrt(cf**2)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cfast(ix_1,ix_2)=sqrt(cfast(ix_1,ix_2))
      enddo
      enddo

* cs=sqrt(cs**2)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cslow(ix_1,ix_2)=sqrt(cslow(ix_1,ix_2))
      enddo
      enddo

      if(oktest)write(*,*)'Average:rho,csound2,dp,rhodv',wroe(ixtest1,ixtest2,
     &   rho_),csound2(ixtest1,ixtest2),dp(ixtest1,ixtest2),rhodv(ixtest1,
     &   ixtest2)
      if(oktest)write(*,*)'Average:cf,cs,af,as',cfast(ixtest1,ixtest2),
     &   cslow(ixtest1,ixtest2),afast(ixtest1,ixtest2),aslow(ixtest1,ixtest2)

*Replace the primitive variables with more useful quantities:
* rho -> sqrt(rho)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wroe(ix_1,ix_2,rho_)=sqrt(wroe(ix_1,ix_2,rho_))
      enddo
      enddo

* Avoid sgn(b_idim)==0
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         if(abs(wroe(ix_1,ix_2,b0_+idim)).lt.smalldouble)wroe(ix_1,ix_2,b0_+
     &      idim)=smalldouble
      enddo
      enddo
* B_idir,jdir -> beta_idir,jdir
      idir=idim+1-ndir*(idim/ndir)
      if(ndir.eq.2)then
          do ix_2=ixmin2,ixmax2
          do ix_1=ixmin1,ixmax1
             if(wroe(ix_1,ix_2,b0_+idir).ge.zero)then
                    wroe(ix_1,ix_2,b0_+idir)=one
                 else
                    wroe(ix_1,ix_2,b0_+idir)=-one
                 endif
          enddo
          enddo
      else
*         beta_j=B_j/sqrt(B_i**2+B_j**2); beta_i=B_i/sqrt(B_i**2+B_j**2)
          jdir=idir+1-ndir*(idir/ndir)
          do ix_2=ixmin2,ixmax2
          do ix_1=ixmin1,ixmax1
             tmp(ix_1,ix_2)=sqrt(wroe(ix_1,ix_2,b0_+idir)**
     &          2+wroe(ix_1,ix_2,b0_+jdir)**2)
          enddo
          enddo
          do ix_2=ixmin2,ixmax2
          do ix_1=ixmin1,ixmax1
             if(tmp(ix_1,ix_2).gt.smalldouble)then
                    wroe(ix_1,ix_2,b0_+idir)=wroe(ix_1,ix_2,b0_+idir)/
     &                 tmp(ix_1,ix_2)
                    wroe(ix_1,ix_2,b0_+jdir)=wroe(ix_1,ix_2,b0_+jdir)/
     &                 tmp(ix_1,ix_2)
                 else
                    wroe(ix_1,ix_2,b0_+idir)=sqrt(half)
                    wroe(ix_1,ix_2,b0_+jdir)=sqrt(half)
                 endif
          enddo
          enddo
      endif

      return
      end

*=============================================================================
      subroutine geteigenjump(wL,wR,wroe,ixmin1,ixmin2,ixmax1,ixmax2,il,idim,
     &   smalla,a,jump)

* Calculate the il-th characteristic speed and the jump in the il-th
* characteristic variable in the idim direction within ixL.
* The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
* jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
* variables. However part of the summation is done in advance and saved into
* bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
* used in the entropy fix.
*
* All the variables are centered on the cell interface, thus the 
* "*C" notation is omitted for sake of brevity.

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,il,idim,idir,jdir
      double precision wL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wR(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision smalla(ixGlo1:ixGhi1,ixGlo2:ixGhi2),a(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),jump(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      double precisioncfast(ixGlo1:ixGhi1,ixGlo2:ixGhi2),cslow(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),afast(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   aslow(ixGlo1:ixGhi1,ixGlo2:ixGhi2),csound2(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rhodv(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
      double precision bdv(ixGlo1:ixGhi1,ixGlo2:ixGhi2),bdb(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      save bdv,bdb
      double precision aL(ixGlo1:ixGhi1,ixGlo2:ixGhi2),aR(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),cs2L(ixGlo1:ixGhi1,ixGlo2:ixGhi2),cs2R(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),cs2ca2L(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   cs2ca2R(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      save cs2L,cs2R,cs2ca2L,cs2ca2R
*-----------------------------------------------------------------------------

      oktest=index(teststr,'geteigenjump').ge.1

      idir=idim+1-ndir*(idim/ndir)
      jdir=idir+1-ndir*(idir/ndir)

      if(il.eq.fastRW_)then
*        Fast and slow waves use bdv=sqrho**2*sign(bx)*(betay*dvy+betaz*dvz)
*        bdb=sqrho*a*          (betay*dBy+betaz*dBz)
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdv(ix_1,ix_2)=wroe(ix_1,ix_2,b0_+idir)* (wR(ix_1,ix_2,m0_+idir)/
     &         wR(ix_1,ix_2,rho_)-wL(ix_1,ix_2,m0_+idir)/wL(ix_1,ix_2,rho_))
         enddo
         enddo
         if(ndir.eq.3)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               bdv(ix_1,ix_2)=bdv(ix_1,ix_2)+wroe(ix_1,ix_2,b0_+jdir)*
     &             (wR(ix_1,ix_2,m0_+jdir)/wR(ix_1,ix_2,rho_)-
     &            wL(ix_1,ix_2,m0_+jdir)/wL(ix_1,ix_2,rho_))
            enddo
            enddo
         endif
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdv(ix_1,ix_2)=bdv(ix_1,ix_2)*sign(wroe(ix_1,ix_2,rho_)**
     &         2,wroe(ix_1,ix_2,b0_+idim))
         enddo
         enddo

         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdb(ix_1,ix_2)=wroe(ix_1,ix_2,b0_+idir)*(wR(ix_1,ix_2,b0_+
     &         idir)-wL(ix_1,ix_2,b0_+idir))
         enddo
         enddo
         if(ndir.eq.3)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               bdb(ix_1,ix_2)=bdb(ix_1,ix_2)+wroe(ix_1,ix_2,b0_+jdir)*
     &            (wR(ix_1,ix_2,b0_+jdir)-wL(ix_1,ix_2,b0_+jdir))
            enddo
            enddo
         endif
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdb(ix_1,ix_2)=bdb(ix_1,ix_2)*sqrt(csound2(ix_1,ix_2))*
     &         wroe(ix_1,ix_2,rho_)
         enddo
         enddo
         if(oktest)write(*,*)'rhobetadv,sqrhoabetadb:',bdv(ixtest1,ixtest2),
     &      bdb(ixtest1,ixtest2)
      endif

      if(il.eq.alfvRW_)then
*        Alfven waves use      bdv=0.5*sqrho**2*      (betaz*dvy-betay*dvz)
*        bdb=0.5*sqrho*sign(bx)*(betaz*dBy-betay*dBz)
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdv(ix_1,ix_2)=wroe(ix_1,ix_2,b0_+jdir)* (wR(ix_1,ix_2,m0_+idir)/
     &         wR(ix_1,ix_2,rho_)-wL(ix_1,ix_2,m0_+idir)/wL(ix_1,ix_2,rho_)) -
     &         wroe(ix_1,ix_2,b0_+idir)* (wR(ix_1,ix_2,m0_+jdir)/
     &         wR(ix_1,ix_2,rho_)-wL(ix_1,ix_2,m0_+jdir)/wL(ix_1,ix_2,rho_))
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdb(ix_1,ix_2)=wroe(ix_1,ix_2,b0_+jdir)*(wR(ix_1,ix_2,b0_+
     &         idir)-wL(ix_1,ix_2,b0_+idir)) -wroe(ix_1,ix_2,b0_+idir)*
     &         (wR(ix_1,ix_2,b0_+jdir)-wL(ix_1,ix_2,b0_+jdir))
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdv(ix_1,ix_2)=bdv(ix_1,ix_2)*half*wroe(ix_1,ix_2,rho_)**2
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            bdb(ix_1,ix_2)=bdb(ix_1,ix_2)*half*sign(wroe(ix_1,ix_2,rho_),
     &         wroe(ix_1,ix_2,b0_+idim))
         enddo
         enddo
         if(oktest)write(*,*)'rhobetaXdv/2,sqrhobetaXdb/2:',bdv(ixtest1,
     &      ixtest2),bdb(ixtest1,ixtest2)
      endif

      if(il.eq.fastRW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)+cfast(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=half/csound2(ix_1,ix_2)*(afast(ix_1,ix_2)*
     &            (+cfast(ix_1,ix_2)*rhodv(ix_1,ix_2)+dp(ix_1,ix_2))+
     &            aslow(ix_1,ix_2)*(-cslow(ix_1,ix_2)*bdv(ix_1,ix_2)+
     &            bdb(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.fastLW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)-cfast(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=half/csound2(ix_1,ix_2)*(afast(ix_1,ix_2)*
     &            (-cfast(ix_1,ix_2)*rhodv(ix_1,ix_2)+dp(ix_1,ix_2))+
     &            aslow(ix_1,ix_2)*(+cslow(ix_1,ix_2)*bdv(ix_1,ix_2)+
     &            bdb(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.slowRW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)+cslow(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=half/csound2(ix_1,ix_2)*(aslow(ix_1,ix_2)*
     &            (+cslow(ix_1,ix_2)*rhodv(ix_1,ix_2)+dp(ix_1,ix_2))+
     &            afast(ix_1,ix_2)*(+cfast(ix_1,ix_2)*bdv(ix_1,ix_2)-
     &            bdb(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.slowLW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)-cslow(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=half/csound2(ix_1,ix_2)*(aslow(ix_1,ix_2)*
     &            (-cslow(ix_1,ix_2)*rhodv(ix_1,ix_2)+dp(ix_1,ix_2))+
     &            afast(ix_1,ix_2)*(-cfast(ix_1,ix_2)*bdv(ix_1,ix_2)-
     &            bdb(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.entroW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=wR(ix_1,ix_2,rho_)-wL(ix_1,ix_2,rho_)-
     &            dp(ix_1,ix_2)/csound2(ix_1,ix_2)
            enddo
            enddo
      else if(il.eq.diverW_)then
            if(divbwave)then
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)
               enddo
               enddo
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  jump(ix_1,ix_2)=wR(ix_1,ix_2,b0_+idim)-wL(ix_1,ix_2,b0_+
     &               idim)
               enddo
               enddo
            else
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  a(ix_1,ix_2)=zero
               enddo
               enddo
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  jump(ix_1,ix_2)=zero
               enddo
               enddo
            endif
      else if(il.eq.alfvRW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)+abs(wroe(ix_1,ix_2,b0_+
     &            idim))/wroe(ix_1,ix_2,rho_)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=+bdv(ix_1,ix_2)-bdb(ix_1,ix_2)
            enddo
            enddo
      else if(il.eq.alfvLW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               a(ix_1,ix_2)=wroe(ix_1,ix_2,m0_+idim)-abs(wroe(ix_1,ix_2,b0_+
     &            idim))/wroe(ix_1,ix_2,rho_)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               jump(ix_1,ix_2)=-bdv(ix_1,ix_2)-bdb(ix_1,ix_2)
            enddo
            enddo
      end if

* Calculate "smalla" or modify "a" based on the "typeentropy" switch

      if(typeentropy(il).eq.'yee')then
*        Based on Yee JCP 68,151 eq 3.23
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            smalla(ix_1,ix_2)=entropycoef(il)
         enddo
         enddo
      else if(typeentropy(il).eq.'harten'.or.typeentropy(il).eq.'powell'.or.
     &   typeentropy(il).eq. 'ratio')then
*        Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
*        Initialize left and right eigenvalues by velocities
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            aL(ix_1,ix_2)= wL(ix_1,ix_2,m0_+idim)/wL(ix_1,ix_2,rho_)
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            aR(ix_1,ix_2)= wR(ix_1,ix_2,m0_+idim)/wR(ix_1,ix_2,rho_)
         enddo
         enddo
*        Calculate the final "aL" and "aR"
      if(il.eq.fastRW_)then
*           These quantities will be used for all the fast and slow waves
*           Calculate soundspeed**2 and cs**2+ca**2.
            call getcsound2(wL,ixmin1,ixmin2,ixmax1,ixmax2,cs2L)
            call getcsound2(wR,ixmin1,ixmin2,ixmax1,ixmax2,cs2R)
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               cs2ca2L(ix_1,ix_2)=cs2L(ix_1,ix_2)+(wL(ix_1,ix_2,b1_)**
     &            2+wL(ix_1,ix_2,b2_)**2)/wL(ix_1,ix_2,rho_)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               cs2ca2R(ix_1,ix_2)=cs2R(ix_1,ix_2)+(wR(ix_1,ix_2,b1_)**
     &            2+wR(ix_1,ix_2,b2_)**2)/wR(ix_1,ix_2,rho_)
            enddo
            enddo
*           Save the discriminants into cs2L and cs2R
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               cs2L(ix_1,ix_2)=sqrt(cs2ca2L(ix_1,ix_2)**2-4*cs2L(ix_1,ix_2)*
     &            wL(ix_1,ix_2,b0_+idim)**2/wL(ix_1,ix_2,rho_))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               cs2R(ix_1,ix_2)=sqrt(cs2ca2R(ix_1,ix_2)**2-4*cs2R(ix_1,ix_2)*
     &            wR(ix_1,ix_2,b0_+idim)**2/wR(ix_1,ix_2,rho_))
            enddo
            enddo

*           The left and right eigenvalues for the fast wave going to right
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aL(ix_1,ix_2)=aL(ix_1,ix_2) + sqrt(half*(cs2ca2L(ix_1,ix_2) +
     &             cs2L(ix_1,ix_2)))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aR(ix_1,ix_2)=aR(ix_1,ix_2) + sqrt(half*(cs2ca2R(ix_1,ix_2) +
     &             cs2R(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.fastLW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aL(ix_1,ix_2)=aL(ix_1,ix_2) - sqrt(half*(cs2ca2L(ix_1,ix_2) +
     &             cs2L(ix_1,ix_2)))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aR(ix_1,ix_2)=aR(ix_1,ix_2) - sqrt(half*(cs2ca2R(ix_1,ix_2) +
     &             cs2R(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.slowRW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aL(ix_1,ix_2)=aL(ix_1,ix_2) + sqrt(half*(cs2ca2L(ix_1,ix_2) -
     &             cs2L(ix_1,ix_2)))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aR(ix_1,ix_2)=aR(ix_1,ix_2) + sqrt(half*(cs2ca2R(ix_1,ix_2) -
     &             cs2R(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.slowLW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aL(ix_1,ix_2)=aL(ix_1,ix_2) - sqrt(half*(cs2ca2L(ix_1,ix_2) -
     &             cs2L(ix_1,ix_2)))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aR(ix_1,ix_2)=aR(ix_1,ix_2) - sqrt(half*(cs2ca2R(ix_1,ix_2) -
     &             cs2R(ix_1,ix_2)))
            enddo
            enddo
      else if(il.eq.entroW_.or.il.eq.diverW_)then
*           These propagate by the velocity
      else if(il.eq.alfvRW_)then
*           Store the Alfven speeds into cs2ca2L and cs2ca2R
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               cs2ca2L(ix_1,ix_2)=abs(wL(ix_1,ix_2,b0_+idim))/
     &            sqrt(wL(ix_1,ix_2,rho_))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               cs2ca2R(ix_1,ix_2)=abs(wR(ix_1,ix_2,b0_+idim))/
     &            sqrt(wR(ix_1,ix_2,rho_))
            enddo
            enddo

            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aL(ix_1,ix_2)=aL(ix_1,ix_2) + cs2ca2L(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aR(ix_1,ix_2)=aR(ix_1,ix_2) + cs2ca2R(ix_1,ix_2)
            enddo
            enddo
      else if(il.eq.alfvLW_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aL(ix_1,ix_2)=aL(ix_1,ix_2) - cs2ca2L(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               aR(ix_1,ix_2)=aR(ix_1,ix_2) - cs2ca2R(ix_1,ix_2)
            enddo
            enddo
         end if
      end if

      call entropyfix(ixmin1,ixmin2,ixmax1,ixmax2,il,aL,aR,a,smalla)

      return
      end

*=============================================================================
      subroutine rtimes(q,wroe,ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,rq)

* Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

      include 'vacdef.f'

      integer           ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,idir,jdir
      double precision  wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rq(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      double precision cfast(ixGlo1:ixGhi1,ixGlo2:ixGhi2),cslow(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),afast(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   aslow(ixGlo1:ixGhi1,ixGlo2:ixGhi2),csound2(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rhodv(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
      double precision bv(ixGlo1:ixGhi1,ixGlo2:ixGhi2),v2a2(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      save  bv,v2a2
*-----------------------------------------------------------------------------

      oktest=index(teststr,'rtimes').ge.1

      idir=idim+1-ndir*(idim/ndir)
      jdir=idir+1-ndir*(idir/ndir)

        if(iw.eq.rho_)then
         if(il.eq.fastRW_.or.il.eq.fastLW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)*afast(ix_1,ix_2)
              enddo
              enddo
         else if(il.eq.slowRW_.or.il.eq.slowLW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)*aslow(ix_1,ix_2)
              enddo
              enddo
         else if(il.eq.entroW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)
              enddo
              enddo
         else if(il.eq.diverW_.or.il.eq.alfvRW_.or.il.eq.alfvLW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=zero
              enddo
              enddo
          end if
        else if(iw.eq.e_)then
          if(il.eq.fastRW_)then
            if(eqpar(gamma_).gt.zero)then
*              !!IDEAL GAS
*              Store 0.5*v**2+(2-gamma)/(gamma-1)*a**2
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  v2a2(ix_1,ix_2)=half*(wroe(ix_1,ix_2,m1_)**
     &               2+wroe(ix_1,ix_2,m2_)**2)+ (2-eqpar(gamma_))/
     &               (eqpar(gamma_)-1)*csound2(ix_1,ix_2)
               enddo
               enddo
            else
               call die(
     &            'Correct rTimes for NONIDEAL gas in src/vacphys.mhdroe.t')
*              Express the partial derivative de/dp using wroe
*              v2a2(ix^S)=half*(^C&wroe(ix^S,m^C_)**2+)+ &
*              (??dedp(ix^S)??-1)*csound2(ix^S)
            endif
*           Store sgn(bx)*(betay*vy+betaz*vz) in bv
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               bv(ix_1,ix_2)=wroe(ix_1,ix_2,b0_+idir)*wroe(ix_1,ix_2,m0_+idir)
            enddo
            enddo
            if(ndir.eq.3)then
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  bv(ix_1,ix_2)=bv(ix_1,ix_2)+wroe(ix_1,ix_2,b0_+jdir)*
     &               wroe(ix_1,ix_2,m0_+jdir)
               enddo
               enddo
            endif
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               bv(ix_1,ix_2)=bv(ix_1,ix_2)*sign(one,wroe(ix_1,ix_2,b0_+idim))
            enddo
            enddo
            if(oktest)write(*,*)'v2/2+(2-g)/(g-1)a2,betav:',v2a2(ixtest1,
     &         ixtest2),bv(ixtest1,ixtest2)
          else if(il.eq.alfvRW_)then
*           Store betaz*vy-betay*vz in bv
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               bv(ix_1,ix_2)=(wroe(ix_1,ix_2,b0_+jdir)*wroe(ix_1,ix_2,m0_+
     &            idir)-wroe(ix_1,ix_2,b0_+idir)*wroe(ix_1,ix_2,m0_+jdir))
            enddo
            enddo
            if(oktest)write(*,*)'betaXv:',bv(ixtest1,ixtest2)
          endif
         if(il.eq.fastRW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)*(-aslow(ix_1,ix_2)*
     &              cslow(ix_1,ix_2)*bv(ix_1,ix_2)+afast(ix_1,ix_2)*
     &              (v2a2(ix_1,ix_2)+cfast(ix_1,ix_2)*(cfast(ix_1,ix_2)+
     &              wroe(ix_1,ix_2,m0_+idim))))
              enddo
              enddo
         else if(il.eq.fastLW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)*(+aslow(ix_1,ix_2)*
     &              cslow(ix_1,ix_2)*bv(ix_1,ix_2)+afast(ix_1,ix_2)*
     &              (v2a2(ix_1,ix_2)+cfast(ix_1,ix_2)*(cfast(ix_1,ix_2)-
     &              wroe(ix_1,ix_2,m0_+idim))))
              enddo
              enddo
         else if(il.eq.slowRW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)*(+afast(ix_1,ix_2)*
     &              cfast(ix_1,ix_2)*bv(ix_1,ix_2)+aslow(ix_1,ix_2)*
     &              (v2a2(ix_1,ix_2)+cslow(ix_1,ix_2)*(cslow(ix_1,ix_2)+
     &              wroe(ix_1,ix_2,m0_+idim))))
              enddo
              enddo
         else if(il.eq.slowLW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)*(-afast(ix_1,ix_2)*
     &              cfast(ix_1,ix_2)*bv(ix_1,ix_2)+aslow(ix_1,ix_2)*
     &              (v2a2(ix_1,ix_2)+cslow(ix_1,ix_2)*(cslow(ix_1,ix_2)-
     &              wroe(ix_1,ix_2,m0_+idim))))
              enddo
              enddo
         else if(il.eq.entroW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)= q(ix_1,ix_2)*half*(wroe(ix_1,ix_2,m1_)**
     &              2+wroe(ix_1,ix_2,m2_)**2)
              enddo
              enddo
         else if(il.eq.diverW_)then
              if(divbwave)then
                 do ix_2=ixmin2,ixmax2
                 do ix_1=ixmin1,ixmax1
                    rq(ix_1,ix_2)= q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+idim)
                 enddo
                 enddo
              else
                 do ix_2=ixmin2,ixmax2
                 do ix_1=ixmin1,ixmax1
                    rq(ix_1,ix_2)= zero
                 enddo
                 enddo
              endif
         else if(il.eq.alfvRW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=+q(ix_1,ix_2)*bv(ix_1,ix_2)
              enddo
              enddo
         else if(il.eq.alfvLW_)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=-q(ix_1,ix_2)*bv(ix_1,ix_2)
              enddo
              enddo
          end if
        else if(iw.eq.m1_.or.iw.eq.m2_)then
          if(iw.eq.m0_+idim)then
           if(il.eq.fastRW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*afast(ix_1,ix_2)*
     &                (wroe(ix_1,ix_2,iw)+cfast(ix_1,ix_2))
                enddo
                enddo
           else if(il.eq.fastLW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*afast(ix_1,ix_2)*
     &                (wroe(ix_1,ix_2,iw)-cfast(ix_1,ix_2))
                enddo
                enddo
           else if(il.eq.slowRW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*aslow(ix_1,ix_2)*
     &                (wroe(ix_1,ix_2,iw)+cslow(ix_1,ix_2))
                enddo
                enddo
           else if(il.eq.slowLW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*aslow(ix_1,ix_2)*
     &                (wroe(ix_1,ix_2,iw)-cslow(ix_1,ix_2))
                enddo
                enddo
           else if(il.eq.entroW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*wroe(ix_1,ix_2,iw)
                enddo
                enddo
           else if(il.eq.diverW_.or.il.eq.alfvLW_.or.il.eq.alfvRW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=zero
                enddo
                enddo
            end if
          else
           if(il.eq.fastRW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*(afast(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,iw)-aslow(ix_1,ix_2)*cslow(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,b0_-m0_+iw)*sign(one,wroe(ix_1,ix_2,b0_+
     &                idim)))
                enddo
                enddo
           else if(il.eq.fastLW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*(afast(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,iw)+aslow(ix_1,ix_2)*cslow(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,b0_-m0_+iw)*sign(one,wroe(ix_1,ix_2,b0_+
     &                idim)))
                enddo
                enddo
           else if(il.eq.slowRW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*(aslow(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,iw)+afast(ix_1,ix_2)*cfast(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,b0_-m0_+iw)*sign(one,wroe(ix_1,ix_2,b0_+
     &                idim)))
                enddo
                enddo
           else if(il.eq.slowLW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*(aslow(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,iw)-afast(ix_1,ix_2)*cfast(ix_1,ix_2)*
     &                wroe(ix_1,ix_2,b0_-m0_+iw)*sign(one,wroe(ix_1,ix_2,b0_+
     &                idim)))
                enddo
                enddo
           else if(il.eq.entroW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=q(ix_1,ix_2)*wroe(ix_1,ix_2,iw)
                enddo
                enddo
           else if(il.eq.diverW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=zero
                enddo
                enddo
           else if(il.eq.alfvRW_)then
                if(iw.eq.m0_+idir)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     rq(ix_1,ix_2)=+q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+jdir)
                  enddo
                  enddo
                else
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     rq(ix_1,ix_2)=-q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+idir)
                  enddo
                  enddo
                endif
           else if(il.eq.alfvLW_)then
                if(iw.eq.m0_+idir)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     rq(ix_1,ix_2)=-q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+jdir)
                  enddo
                  enddo
                else
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     rq(ix_1,ix_2)=+q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+idir)
                  enddo
                  enddo
                endif        
            end if
*         iw=m_idir,m_jdir
          end if
        else if(iw.eq.b1_.or.iw.eq.b2_)then
          if(iw.eq.b0_+idim)then
            if(il.eq.diverW_ .and. divbwave)then
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=q(ix_1,ix_2)
              enddo
              enddo
            else
              do ix_2=ixmin2,ixmax2
              do ix_1=ixmin1,ixmax1
                 rq(ix_1,ix_2)=zero
              enddo
              enddo
            endif
          else
           if(il.eq.fastRW_.or.il.eq.fastLW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=+q(ix_1,ix_2)*aslow(ix_1,ix_2)*
     &                sqrt(csound2(ix_1,ix_2))*wroe(ix_1,ix_2,iw)/
     &                wroe(ix_1,ix_2,rho_)
                enddo
                enddo
           else if(il.eq.slowRW_.or.il.eq.slowLW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=-q(ix_1,ix_2)*afast(ix_1,ix_2)*
     &                sqrt(csound2(ix_1,ix_2))*wroe(ix_1,ix_2,iw)/
     &                wroe(ix_1,ix_2,rho_)
                enddo
                enddo
           else if(il.eq.entroW_.or.il.eq.diverW_)then
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   rq(ix_1,ix_2)=zero
                enddo
                enddo
           else if(il.eq.alfvRW_.or.il.eq.alfvLW_)then
                if(iw.eq.b0_+idir)then
                   do ix_2=ixmin2,ixmax2
                   do ix_1=ixmin1,ixmax1
                      rq(ix_1,ix_2)=-q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+jdir)/
     &                   sign(wroe(ix_1,ix_2,rho_),wroe(ix_1,ix_2,b0_+idim))
                   enddo
                   enddo
                else
                   do ix_2=ixmin2,ixmax2
                   do ix_1=ixmin1,ixmax1
                      rq(ix_1,ix_2)=+q(ix_1,ix_2)*wroe(ix_1,ix_2,b0_+idir)/
     &                   sign(wroe(ix_1,ix_2,rho_),wroe(ix_1,ix_2,b0_+idim))
                   enddo
                   enddo
                end if
            end if
*         iw=b_idir,b_jdir
          end if
      end if

      return
      end
*=============================================================================
* end module vacphys.mhdroe
*##############################################################################
*##############################################################################
* module vacproc.projectb - Projection of B for mhd(iso) in 2 or 3D

*=============================================================================
      subroutine projectb(w)

* Project B according to B'=B-grad phi, where Laplace phi=div B within ixM

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      double precision  divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2),phi(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),divbmax,qdivbmax,qdivbmin
      integer  ixmin1,ixmin2,ixmax1,ixmax2,iB,idim,info,matvecmax
      character*3  typestop

      save divbmax,matvecmax,typestop
      data divbmax /0.D0/
*-----------------------------------------------------------------------------

      oktest=index(teststr,'projectb').ge.1
      if(oktest)write(*,*)'ProjectB w:',w(ixtest1,ixtest2,iwtest)

* Determine the div B error in ix=ixM+1
      call getdivb(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,divb)

* Calculate and save divbmax, the maximum div B allowed
      if(divbmax.lt.smalldouble)then
         maxval_1=-bigdouble
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            maxval_1=max(maxval_1,abs(divb(ix_1,ix_2)))
         enddo
         enddo
         divbmax=maxval_1
         
         if(verbose)write(*,*)'ProjectB it,divbmax:',it,divbmax
*        Default values if procpar was not set
         if(procpar(divbcoeff_).eq.-one)procpar(divbcoeff_)=one/10
         if(procpar(divbconst_).eq.-one)procpar(divbconst_)=zero

*        Determine relative or absolute limit for div B
         if(procpar(divbcoeff_).lt.zero.and.procpar(divbcoeff_).gt.-one)then
*           1 > -divbcoeff > 0 means relative stopping criterion
            divbmax=-procpar(divbcoeff_)
            typestop='rel'
            if(verbose)write(*,*)'Reduction factor for div B:',divbmax
         else
            divbmax=divbmax*max(zero,procpar(divbcoeff_))+max(zero,
     &         procpar(divbconst_))
            typestop='max'
            if(verbose)write(*,*)'Allowed maximum for div B:',divbmax
            if(divbmax.lt.smalldouble) call die(
     &         'Error in ProjectB: Too small value for divbmax')
         endif

*        -divbcoeff>1 or -divbconst>1 gives maximum number of iterations
         matvecmax=1000
         if(procpar(divbconst_).lt.-one)matvecmax=nint(-procpar(divbconst_))
         if(procpar(divbcoeff_).lt.-one)matvecmax=nint(-procpar(divbcoeff_))

         if(verbose)write(*,*)'Maximum number of matvecs:',matvecmax

      endif

      if(oktest)then
         call getdivb(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,divb)
         maxval_1=-bigdouble
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            maxval_1=max(maxval_1,divb(ix_1,ix_2))
         enddo
         enddo
         qdivbmax=maxval_1
         minval_1=bigdouble
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            minval_1=min(minval_1,divb(ix_1,ix_2))
         enddo
         enddo
         qdivbmin=minval_1
         if(verbose)write(*,*)'Max and min of divb:',qdivbmax,qdivbmin
*        if(oktest)write(*,*)'Indices for max:',maxloc(divb(ix^S))
      endif

* Determine boundary condition for the Poisson solver based on 
* the bpundary condition for the normal component of the magnetic field
      do iB=1,nB
      if(typeB(b0_+idimB(iB),iB).eq.'periodic')then
           typeBscalar(iB)='periodic'
      else if(typeB(b0_+idimB(iB),iB).eq.'symm'.or.typeB(b0_+
     &   idimB(iB),iB).eq.'symm0')then
           typeBscalar(iB)='asymm'
      else if(typeB(b0_+idimB(iB),iB).eq.'asymm')then
           typeBscalar(iB)='symm'
      else if(typeB(b0_+idimB(iB),iB).eq.'fixed'.or.typeB(b0_+
     &   idimB(iB),iB).eq.'fixed1')then
           typeBscalar(iB)='grad0'
         
      else
           typeBscalar(iB)='nul'
         end if
         if(it.eq.itmin.and.oktest)write(*,*)'iB,idim,typeB,typeBscalar:',iB,
     &      idimB(iB),typeB(b0_+idimB(iB),iB),typeBscalar(iB)
      enddo 


* Initial guess for phi is zero for the iterative solvers (nonzero='false')
      do ix_2=ixMmin2,ixMmax2
      do ix_1=ixMmin1,ixMmax1
         phi(ix_1,ix_2)=zero
      enddo
      enddo

* Solve the Poisson equation
              
      call poisson('project B ',divb,divbmax,typestop,matvecmax,info,.false.,
     &   phi)
* call die('Error: Poisson solver is OFF! setvac -on=poisson; make vac')

      if(oktest)write(*,*)'Poisson solver info:',info

* Do not do anything if the initial guess satisfied the stopping criterion
      if(info.eq.3)return

* Do not subtract grad(phi) from B if the iterations did not reduce the error
      if(info.lt.0)return

* Subtract tmp=grad(phi) from the first ndim components of the B field in ixM+1

* First get the ghost cell values for the solution
      call boundscalar(phi)

* For the full MHD equations correct total energy according to the value of
* |procpar(3)|=1, 2, or 3
      if(fourthorder)then
         ixmin1=ixMmin1
         ixmin2=ixMmin2
         ixmax1=ixMmax1
         ixmax2=ixMmax2
      else
         ixmin1=ixMmin1-1
         ixmin2=ixMmin2-1
         ixmax1=ixMmax1+1
         ixmax2=ixMmax2+1
      endif
      do idim=1,ndim
         if(fourthorder)then
            call gradient4(.true.,phi,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp)
         else
            call gradient(.true.,phi,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp)
         endif
         if(typephys.eq.'mhd')then
         if(nint(abs(procpar(divbbound_))).eq.2)then
*              Correct total energy so that thermal pressure is kept constant
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  w(ix_1,ix_2,ee_)=w(ix_1,ix_2,ee_)-w(ix_1,ix_2,b0_+idim)*
     &               tmp(ix_1,ix_2)+tmp(ix_1,ix_2)**2/2
               enddo
               enddo
         else if(nint(abs(procpar(divbbound_))).eq.3)then
*              Correct total energy so that total pressure is kept constant
               if(eqpar(gamma_).ne.two.and.eqpar(gamma_).ne.one)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     w(ix_1,ix_2,ee_)=w(ix_1,ix_2,ee_)+(eqpar(gamma_)-two)/
     &                  (eqpar(gamma_)-one)*(-w(ix_1,ix_2,b0_+idim)*
     &                  tmp(ix_1,ix_2)+tmp(ix_1,ix_2)**2/2)
                  enddo
                  enddo
               endif
            end if
         endif
*        B'=B-grad(Phi)
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            w(ix_1,ix_2,b0_+idim)=w(ix_1,ix_2,b0_+idim)-tmp(ix_1,ix_2)
         enddo
         enddo
      end do

* Recalculate boundaries for the first ndim components of B if 
* procpar(divbbound_) is positive
      if(procpar(divbbound_).gt.zero) call getboundary(t,b1_,b2_,1,ndim,w)

      if(oktest)then
         call getdivb(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,divb)
         maxval_1=-bigdouble
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            maxval_1=max(maxval_1,divb(ix_1,ix_2))
         enddo
         enddo
         qdivbmax=maxval_1
         minval_1=bigdouble
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            minval_1=min(minval_1,divb(ix_1,ix_2))
         enddo
         enddo
         qdivbmin=minval_1
         if(verbose) write(*,*)'New   extrema of divb:',qdivbmax,qdivbmin
      endif

      return
      end

*=============================================================================
* end module vacproc.projectb
*##############################################################################

*##############################################################################
* module vacphys.mhdres - subroutines for resistive mhd and mhdiso

*=============================================================================
      subroutine getdt_res(w,ixmin1,ixmin2,ixmax1,ixmax2)

* If resistivity is  not zero, check diffusion time limit for dt

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),dtdiff
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim,idirmin
      save dtdiff

      double precision  current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*
     &   ndir:3),eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradeta(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,ndim)
      common/resist/current,eta,gradeta
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getdt').ge.1
      if(oktest)write(*,*)'GetDt_Res'

      if(eqpar(eta_).eq.zero)return

      if(eqpar(eta_).gt.zero)then
         minval_1=bigdouble
         do idim_3=1,ndim
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            minval_1=min(minval_1,dx(ix_1,ix_2,idim_3))
         enddo
         enddo
         enddo
         dtdiff=dtdiffpar*minval_1**2/eqpar(eta_)
      else if(eqpar(eta_).lt.zero)then
         if(it.eq.itmin)then
            call getcurrent(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,idirmin)
            call specialeta(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2,idirmin)
         endif
         dtdiff=bigdouble
         do idim=1,ndim
            maxval_1=-bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               maxval_1=max(maxval_1,eta(ix_1,ix_2)/dx(ix_1,ix_2,idim)**2)
            enddo
            enddo
            dtdiff=min(dtdiff,dtdiffpar/(smalldouble+maxval_1))
         enddo
      endif


      dt=min(dt,dtdiff)

      if(oktest) write(*,*)'GetDt dtdiff:',dtdiff
      if(oktest) write(*,*)'GetDt dt    :',dt

      return
      end

*=============================================================================
      subroutine getcurrent(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)

* Calculate idirmin and the idirmin:3 components of the common current array

      include 'vacdef.f'

      integer  idirmin0
      PARAMETER(  idirmin0=7-2*ndir)
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer           ixmin1,ixmin2,ixmax1,ixmax2,idirmin

      integer  ixImin1,ixImin2,ixImax1,ixImax2,idir,jdir,kdir

* For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
      double precision  current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*
     &   ndir:3),eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradeta(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,ndim)
      common/resist/current,eta,gradeta
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getcurrent').ge.1
      if(oktest)write(*,*)'GetCurrent'

      ixImin1=ixmin1-1
      ixImin2=ixmin2-1
      ixImax1=ixmax1+1
      ixImax2=ixmax2+1

* Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
* Current can have components (idirmin0:3)
* Determine exact value of idirmin while doing the loop.

      idirmin=4
      do i_3=idirmin0,3
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         current(ix_1,ix_2,i_3)=zero
      enddo
      enddo
      enddo

      do idir=idirmin0,3
       do jdir=1,ndim
       do kdir=1,ndir
         if(lvc(idir,jdir,kdir).ne.0)then
            do ix_2=ixImin2,ixImax2
            do ix_1=ixImin1,ixImax1
               tmp(ix_1,ix_2)=w(ix_1,ix_2,b0_+kdir)
            enddo
            enddo
            call gradient(.true.,tmp,ixmin1,ixmin2,ixmax1,ixmax2,jdir,tmp2)
            if(lvc(idir,jdir,kdir).eq.1)then
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  current(ix_1,ix_2,idir)=current(ix_1,ix_2,idir)+
     &               tmp2(ix_1,ix_2)
               enddo
               enddo
            else
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  current(ix_1,ix_2,idir)=current(ix_1,ix_2,idir)-
     &               tmp2(ix_1,ix_2)
               enddo
               enddo
            endif
            if(idir.lt.idirmin)idirmin=idir
         endif
      enddo
       enddo
       enddo

      if(oktest)then
         write(*,*)'idirmin,J(idirmin:3):',idirmin,(current(ixtest1,ixtest2,
     &      i_1),i_1=idirmin,3)
      endif

      return
      end

*=============================================================================
      subroutine addsource_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,
     &   ixmin2,ixmax1,ixmax2,iws,qtC,w,qt,wnew)

* Add resistive source to wnew within ixL if possible, otherwise shrink ixL
* Uses 3 point stencil (1 neighbour) in each direction, non-conservative

      include 'vacdef.f'

      integer           ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,
     &   ixmax2,iws(niw_)
      double precision  qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      integer  ix,jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,hxmin2,hxmax1,hxmax2,
     &   idim,idir,jdir,kdir,idirmin,iiw,iw

* Resistivity "eta" may or may not vary in time and/or space
* For ndir=2 only 3rd component of J can exist, ndir=1 is not possible for MHD
      double precision  current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*
     &   ndir:3),eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradeta(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,ndim)
      common/resist/current,eta,gradeta
*-----------------------------------------------------------------------------

      oktest=index(teststr,'addsource_res').ge.1
      if(oktest)write(*,*)'AddSource_Res1'

* Compact resistive sources involve one extra layer only
      call ensurebound(1,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,
     &   ixmax2,qtC,w)

* Calculate current density within ixL and determine idirmin
      call getcurrent(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)

* Calculate and save eta for the first time
* for eqpar(eta_)<0 call specialeta, this will also set the common gradeta
      if(eqpar(eta_).gt.zero)then
         if(it.eq.itmin)then
            do ix_2=ixGmin2,ixGmax2
            do ix_1=ixGmin1,ixGmax1
               eta(ix_1,ix_2)=eqpar(eta_)
            enddo
            enddo
         endif
      else
         call specialeta(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)
      endif
      if(oktest)write(*,*)'eta    :',eta(ixtest1,ixtest2)
      if(oktest)then
         write(*,*)'gradeta:',(gradeta(ixtest1,ixtest2,idim_1),idim_1=1,ndim)
      endif

      do idir=1,ndir

*        Put B_idir into tmp2 and eta*Laplace B_idir into tmp
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=zero
         enddo
         enddo
         do ix_2=ixGmin2,ixGmax2
         do ix_1=ixGmin1,ixGmax1
            tmp2(ix_1,ix_2)=w(ix_1,ix_2,b0_+idir)
         enddo
         enddo
         if(gencoord)then 
*           Use contour integral of Grad(B_idir) along cell edges
*           !! Assumes that connected cell centers are orthogonal to interfaces
            do idim=1,ndim
*              SHIFT
               jxmin1=ixmin1+kr(idim,1)
               jxmin2=ixmin2+kr(idim,2)
               jxmax1=ixmax1+kr(idim,1)
               jxmax2=ixmax2+kr(idim,2)
*              SHIFT MORE
               hxmin1=ixmin1-kr(idim,1)
               hxmin2=ixmin2-kr(idim,2)
               hxmax1=ixmax1-kr(idim,1)
               hxmax2=ixmax2-kr(idim,2)
*              SHIFT BEGIN
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  tmp(ix_1,ix_2)=tmp(ix_1,ix_2) +surfaceC(ix_1,ix_2,idim)*
     &               (tmp2(ix_1+(jxmin1-ixmin1),ix_2+(jxmin2-ixmin2))-
     &               tmp2(ix_1,ix_2)) /sqrt((x(ix_1+(jxmin1-ixmin1),ix_2+
     &               (jxmin2-ixmin2),1)-x(ix_1,ix_2,1))**2+(x(ix_1+
     &               (jxmin1-ixmin1),ix_2+(jxmin2-ixmin2),2)-x(ix_1,ix_2,2))**
     &               2) +surfaceC(ix_1+(hxmin1-ixmin1),ix_2+
     &               (hxmin2-ixmin2),idim)*(tmp2(ix_1+(hxmin1-ixmin1),ix_2+
     &               (hxmin2-ixmin2))-tmp2(ix_1,ix_2)) /sqrt((x(ix_1+
     &               (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2),1)-x(ix_1,ix_2,1))**
     &               2+(x(ix_1+(hxmin1-ixmin1),ix_2+(hxmin2-ixmin2),2)-
     &               x(ix_1,ix_2,2))**2)
               enddo
               enddo
*              SHIFT END
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               tmp(ix_1,ix_2)=tmp(ix_1,ix_2)/dvolume(ix_1,ix_2)
            enddo
            enddo
         else
            do idim=1,ndim
               if(typeaxial.eq.'cylinder'.and.idim.eq.r_)then
*                 Calculate 1/r d(r d(B_idir)/dr)/dr
                  do ix=ixmin1,ixmax1
                     do ix_1=ixmin2,ixmax2
                        tmp(ix,ix_1)=tmp(ix,ix_1) +(areaC(ix)*
     &                     (tmp2(ix+1,ix_1)-tmp2(ix,ix_1))/
     &                     (x(ix+1,ix_1,r_)-x(ix,ix_1,r_)) +areaC(ix-1)*
     &                     (tmp2(ix-1,ix_1)-tmp2(ix,ix_1))/
     &                     (x(ix,ix_1,r_)-x(ix-1,ix_1,r_))) /x(ix,ix_1,r_)/
     &                     dx(ix,ix_1,r_)
                     enddo
                  enddo
               else
*                 SHIFT
                  jxmin1=ixmin1+kr(idim,1)
                  jxmin2=ixmin2+kr(idim,2)
                  jxmax1=ixmax1+kr(idim,1)
                  jxmax2=ixmax2+kr(idim,2)
                   
*                 SHIFT MORE
                  hxmin1=ixmin1-kr(idim,1)
                  hxmin2=ixmin2-kr(idim,2)
                  hxmax1=ixmax1-kr(idim,1)
                  hxmax2=ixmax2-kr(idim,2)
*                 SHIFT BEGIN
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     tmp(ix_1,ix_2)=tmp(ix_1,ix_2)+(tmp2(ix_1+
     &                  (jxmin1-ixmin1),ix_2+(jxmin2-ixmin2))-2*
     &                  tmp2(ix_1,ix_2)+tmp2(ix_1+(hxmin1-ixmin1),ix_2+
     &                  (hxmin2-ixmin2)))/dx(ix_1,ix_2,idim)**2 
                  enddo
                  enddo

                  
*                 SHIFT END
               endif
            enddo
         endif
*        Axial symmetry : Add -B_r/r**2 or -B_phi/r**2
         if(typeaxial.eq.'cylinder'.and.(idir.eq.r_.or.idir.eq.phi_))then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               tmp(ix_1,ix_2)=tmp(ix_1,ix_2)-tmp2(ix_1,ix_2)/x(ix_1,ix_2,r_)**
     &            2
            enddo
            enddo
         endif

*        Multiply by eta
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=tmp(ix_1,ix_2)*eta(ix_1,ix_2)
         enddo
         enddo

*        Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
         if(eqpar(eta_).lt.zero)then
            do jdir=1,ndim
             do kdir=idirmin,3
               if(lvc(idir,jdir,kdir).ne.0)then
                  if(lvc(idir,jdir,kdir).eq.1)then
                     do ix_2=ixmin2,ixmax2
                     do ix_1=ixmin1,ixmax1
                        tmp(ix_1,ix_2)=tmp(ix_1,ix_2)-gradeta(ix_1,ix_2,jdir)*
     &                     current(ix_1,ix_2,kdir)
                     enddo
                     enddo
                  else
                     do ix_2=ixmin2,ixmax2
                     do ix_1=ixmin1,ixmax1
                        tmp(ix_1,ix_2)=tmp(ix_1,ix_2)+gradeta(ix_1,ix_2,jdir)*
     &                     current(ix_1,ix_2,kdir)
                     enddo
                     enddo
                  endif
               endif
            enddo
             enddo
         endif

*        Add sources related to eta*laplB-grad(eta) x J to B and e
         do iiw=1,iws(niw_)
          iw=iws(iiw)
            if(iw.eq.b0_+idir)then
*              dB_idir/dt+=tmp
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+qdt*tmp(ix_1,ix_2)
               enddo
               enddo
            else if(iw.eq.e_)then
*              de/dt+=B.tmp
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+qdt*tmp(ix_1,ix_2)*
     &               w(ix_1,ix_2,b0_+idir)
               enddo
               enddo
            endif
*        iiw
         end do
*     idir
      enddo

* de/dt+=eta*J**2
      do iiw=1,iws(niw_)
       iw=iws(iiw)
         if(iw.eq.e_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               tmp(ix_1,ix_2)=zero
            enddo
            enddo
            do idir=idirmin,3
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  tmp(ix_1,ix_2)=tmp(ix_1,ix_2)+current(ix_1,ix_2,idir)**2
               enddo
               enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+qdt*eta(ix_1,ix_2)*
     &            tmp(ix_1,ix_2)
            enddo
            enddo
         endif
      enddo

      return
      end

*=============================================================================
      subroutine addsource_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

* Add resistive source to wnew within ixO if possible, otherwise shrink ixO
* Uses 5 point stencil (2 neighbours) in each direction, conservative

      include 'vacdef.f'

      integer           ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws(niw_)
      double precision  qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      integer  ixmin1,ixmin2,ixmax1,ixmax2,idir,jdir,kdir,idirmin,iiw,iw

* Resistivity may or may not vary in time and/or space
* For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
      double precision  current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*
     &   ndir:3),eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradeta(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,ndim)
      common/resist/current,eta,gradeta
*-----------------------------------------------------------------------------

      oktest=index(teststr,'addsource_res').ge.1
      if(oktest)write(*,*)'AddSource_Res2'

* Calculating resistive sources involves second derivatives, two extra layers
      if(oktest)write(*,*)'calling ensurebound in Addsource_Res2'
      call ensurebound(2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,qtC,w)
      if(oktest)write(*,*)'end calling ensurebound in Addsource_Res2'
      ixmin1=ixOmin1-1
      ixmin2=ixOmin2-1
      ixmax1=ixOmax1+1
      ixmax2=ixOmax2+1

* Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
* Determine exact value of idirmin while doing the loop.
      call getcurrent(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)

* Calculate and save eta for the first time
* for eqpar(eta_)<0 call specialeta
      if(eqpar(eta_).gt.zero)then
         if(it.eq.itmin)then
            do ix_2=ixGmin2,ixGmax2
            do ix_1=ixGmin1,ixGmax1
               eta(ix_1,ix_2)=eqpar(eta_)
            enddo
            enddo
         endif
      else
         call specialeta(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)
      endif
      if(oktest)write(*,*)'eta:',eta(ixtest1,ixtest2)

* Calculate sources from resistivity
      do iiw=1,iws(niw_)
       iw=iws(iiw)
      if(iw.eq.b1_.or.iw.eq.b2_)then
*           dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
            idir=iw-b0_
            do jdir=1,ndim
             do kdir=idirmin,3
               if(lvc(idir,jdir,kdir).ne.0)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     tmp(ix_1,ix_2)=current(ix_1,ix_2,kdir)*eta(ix_1,ix_2)*qdt
                  enddo
                  enddo
                  call gradient(.true.,tmp,ixOmin1,ixOmin2,ixOmax1,ixOmax2,
     &               jdir,tmp2)
                  if(lvc(idir,jdir,kdir).eq.1)then
                     do ix_2=ixOmin2,ixOmax2
                     do ix_1=ixOmin1,ixOmax1
                        wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-tmp2(ix_1,ix_2)
                     enddo
                     enddo
                  else
                     do ix_2=ixOmin2,ixOmax2
                     do ix_1=ixOmin1,ixOmax1
                        wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+tmp2(ix_1,ix_2)
                     enddo
                     enddo
                  endif
              endif
            enddo
             enddo
      else if(iw.eq.e_)then
*           de/dt+= div(B x Jeta), thus e-= dt*eps_ijk d_i B_j Jeta_k
            do idir=1,ndim
             do jdir=1,ndir
             do kdir=idirmin,3
               if(lvc(idir,jdir,kdir).ne.0)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     tmp(ix_1,ix_2)=w(ix_1,ix_2,b0_+jdir)*current(ix_1,ix_2,
     &                  kdir)*eta(ix_1,ix_2)*qdt
                  enddo
                  enddo
                  call gradient(.false.,tmp,ixOmin1,ixOmin2,ixOmax1,ixOmax2,
     &               idir,tmp2)
                  if(lvc(idir,jdir,kdir).eq.1)then
                     do ix_2=ixOmin2,ixOmax2
                     do ix_1=ixOmin1,ixOmax1
                        wnew(ix_1,ix_2,ee_)=wnew(ix_1,ix_2,ee_)+
     &                     tmp2(ix_1,ix_2)
                     enddo
                     enddo
                  else
                     do ix_2=ixOmin2,ixOmax2
                     do ix_1=ixOmin1,ixOmax1
                        wnew(ix_1,ix_2,ee_)=wnew(ix_1,ix_2,ee_)-
     &                     tmp2(ix_1,ix_2)
                     enddo
                     enddo
                  endif
              endif
            enddo
             enddo
             enddo
*        iw
         end if
*     iiw
      end do

      return
      end

*=============================================================================
* end module vacphys.mhdres
*##############################################################################

*=============================================================================
      subroutine physini

* Tell VAC which variables are vectors, set default entropy coefficients

      include 'vacdef.f'
      integer  il
*-----------------------------------------------------------------------------

      iw_vector(1)=m0_
       iw_vector(2)=b0_

* The values of the constants are taken from Ryu & Jones ApJ 442, 228
      do il=1,nw
      if(il.eq.fastRW_.or.il.eq.fastLW_.or.il.eq.slowRW_.or.il.eq.slowLW_)then
            entropycoef(il)=0.2
      else if(il.eq.alfvRW_.or.il.eq.alfvLW_)then
            entropycoef(il)=0.4
      else
            entropycoef(il)= -one
         end if
      end do

      return
      end

*=============================================================================
      subroutine process(count,idimmin,idimmax,w)

* Process w before it is advected in directions idim^LIM, or before save
* count=1 and 2 for first and second (half step) processing during advection
* count=ifile+2 for saving results into the file indexed by ifile

      include 'vacdef.f'

      integer  count,idimmin,idimmax
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      logical  oktime
      double precision  cputime,time1,timeproc
      data timeproc /0.D0/

* The processing should eliminate divergence of B.
*-----------------------------------------------------------------------------

      oktest=index(teststr,'process').ge.1
      oktime=index(teststr,'timeproc').ge.1

      if(oktest)write(*,*)'Process it,idim^LIM,count',it,idimmin,idimmax,count

      if(oktime)time1=cputime()


      if(count.eq.0)then
         if(divbconstrain)then
            call die('CT module is OFF: setvac -on=ct; make vac')
         endif
      else
*        Use the projection scheme 
                  call projectb(w)
         if(.false.)call die(
     &      'Poisson module is OFF: setvac -on=poisson;make vac')
      endif


      if(oktime)then
         time1=cputime()-time1
         timeproc=timeproc+time1
         write(*,*)'Time.Proc:',time1,timeproc
      endif

      return
      end

*=============================================================================
      subroutine getdt(w,ixmin1,ixmin2,ixmax1,ixmax2)

* If resistivity is  not zero, check diffusion time limit for dt

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixmin1,ixmin2,ixmax1,ixmax2
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getdt').ge.1
      if(oktest)write(*,*)'GetDt'

      if(eqpar(eta_).eq.zero)return

             call getdt_res(w,ixmin1,ixmin2,ixmax1,ixmax2)


      return
      end

*=============================================================================
      subroutine getdivb(w,ixOmin1,ixOmin2,ixOmax1,ixOmax2,divb)

* Calculate div B within ixO

      include 'vacdef.f'

      integer           ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixmin1,ixmin2,ixmax1,
     &   ixmax2,idim
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),divb(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getdivb').ge.1
      if(oktest)write(*,*)'getdivb ixO=',ixOmin1,ixOmin2,ixOmax1,ixOmax2

      if(fourthorder)then
         ixmin1=ixOmin1-2
         ixmin2=ixOmin2-2
         ixmax1=ixOmax1+2
         ixmax2=ixOmax2+2
      else
         ixmin1=ixOmin1-1
         ixmin2=ixOmin2-1
         ixmax1=ixOmax1+1
         ixmax2=ixOmax2+1
      endif
      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         divb(ix_1,ix_2)=zero
      enddo
      enddo
      do idim=1,ndim
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=w(ix_1,ix_2,b0_+idim)
         enddo
         enddo
         if(fourthorder)then
            call gradient4(.false.,tmp,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,
     &         tmp2)
         else
            call gradient(.false.,tmp,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,
     &         tmp2)
         endif
         do ix_2=ixOmin2,ixOmax2
         do ix_1=ixOmin1,ixOmax1
            divb(ix_1,ix_2)=divb(ix_1,ix_2)+tmp2(ix_1,ix_2)
         enddo
         enddo
      enddo

      if(oktest)then
         write(*,*)'divb:',divb(ixtest1,ixtest2)
*   write(*,*)'bx=',w(ixtest1-1:ixtest1+1,ixtest2,b1_)
*   write(*,*)'by=',w(ixtest1,ixtest2-1:ixtest2+1,b2_)
*   write(*,*)'x =',x(ixtest1-1:ixtest1+1,ixtest2,1)
*   write(*,*)'y =',x(ixtest1,ixtest2-1:ixtest2+1,2)
*   write(*,*)'dx=',dx(ixtest1,ixtest2,1)
*   write(*,*)'dy=',dx(ixtest1,ixtest2,2)
      endif

      return
      end

*=============================================================================
      subroutine getflux(w,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,f,transport)

* Calculate non-transport flux f_idim[iw] within ix^L.
* Set transport=.true. if a transport flux should be added

      include 'vacdef.f'

      integer           ixmin1,ixmin2,ixmax1,ixmax2,iw,idim
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),f(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      logical           transport
*-----------------------------------------------------------------------------

      oktest= index(teststr,'getflux').ge.1
      if(oktest.and.iw.eq.iwtest)write(*,*)'Getflux idim,w:',idim,w(ixtest1,
     &   ixtest2,iwtest)

      transport=.true.

*        f_i[rho]=v_i*rho
      if(iw.eq.rho_)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               f(ix_1,ix_2)=zero
            enddo
            enddo

*        f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
      else if(iw.eq.m1_)then
            if(idim.eq.1)then
               call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,f)
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)=f(ix_1,ix_2) -w(ix_1,ix_2,b1_)*
     &               w(ix_1,ix_2,b0_+idim)
               enddo
               enddo
            else
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)= -w(ix_1,ix_2,b1_)*w(ix_1,ix_2,b0_+idim)
               enddo
               enddo
            endif 
      else if(iw.eq.m2_)then
            if(idim.eq.2)then
               call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,f)
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)=f(ix_1,ix_2) -w(ix_1,ix_2,b2_)*
     &               w(ix_1,ix_2,b0_+idim)
               enddo
               enddo
            else
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)= -w(ix_1,ix_2,b2_)*w(ix_1,ix_2,b0_+idim)
               enddo
               enddo
            endif 

*        f_i[e]=v_i*e+(m_i*ptotal-b_i*(b_k*m_k))/rho
      else if(iw.eq.e_)then
            call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,f)
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               f(ix_1,ix_2)=(w(ix_1,ix_2,m0_+idim)*f(ix_1,ix_2)-
     &             w(ix_1,ix_2,b0_+idim)*( w(ix_1,ix_2,b1_)*
     &            w(ix_1,ix_2,m1_)+w(ix_1,ix_2,b2_)*w(ix_1,ix_2,m2_) ))/
     &            w(ix_1,ix_2,rho_)
            enddo
            enddo

*        f_i[b_k]=v_i*b_k-m_k/rho*b_i
      else if(iw.eq.b1_)then
            if(idim.eq.1) then
*              f_i[b_i] should be exactly 0, so we do not use the transport flux
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)=zero
               enddo
               enddo
               transport=.false.
            else
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)= -w(ix_1,ix_2,m1_)/w(ix_1,ix_2,rho_)*
     &               w(ix_1,ix_2,b0_+idim)
               enddo
               enddo
            endif  
      else if(iw.eq.b2_)then
            if(idim.eq.2) then
*              f_i[b_i] should be exactly 0, so we do not use the transport flux
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)=zero
               enddo
               enddo
               transport=.false.
            else
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  f(ix_1,ix_2)= -w(ix_1,ix_2,m2_)/w(ix_1,ix_2,rho_)*
     &               w(ix_1,ix_2,b0_+idim)
               enddo
               enddo
            endif  

      else
            call die('Error in getflux: unknown flow variable')
      end if

      if(oktest.and.iw.eq.iwtest)write(*,*)'transport,f:',transport,f(ixtest1,
     &   ixtest2)

      return
      end

*=============================================================================
      subroutine addsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

* Add sources from resistivity and Powell solver

      include 'vacdef.f'

      integer           ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws(niw_)
      double precision  qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'addsource').ge.1
      if(oktest)write(*,*)'Addsource, compactres,divbfix:',compactres,divbfix
      if(oktest)write(*,*)'Before adding source:',wnew(ixtest1,ixtest2,iwtest)

* Sources for resistivity in eqs. for e, B1, B2 and B3
      if(abs(eqpar(eta_)).gt.smalldouble)then
               
         if(compactres)then
            call addsource_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &         ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)
         else
            call addsource_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &         ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)
         endif
         if(oktest)write(*,*)'With resistive source:',wnew(ixtest1,ixtest2,
     &      iwtest)
        
         
      endif


* Sources related to div B in the Powell solver
       if(divbfix) call addsource_divb(qdt,ixImin1,ixImin2,ixImax1,ixImax2,
     &    ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

      if(oktest)write(*,*)'After adding source:',wnew(ixtest1,ixtest2,iwtest)

      return
      end

*=============================================================================
      subroutine addsource_divb(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

* Add Powell's divB related sources to wnew within ixO if possible, 
* otherwise shrink ixO

      include 'vacdef.f'

      integer           ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws(niw_),iiw,iw
      double precision  qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      double precision  divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

* Calculating div B involves first derivatives
      call ensurebound(1,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,qtC,w)

* We calculate now div B
      call getdivb(w,ixOmin1,ixOmin2,ixOmax1,ixOmax2,divb)
      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         divb(ix_1,ix_2)=qdt*divb(ix_1,ix_2)
      enddo
      enddo

      do iiw=1,iws(niw_)
       iw=iws(iiw)
         if(iw.eq.m1_)then
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-w(ix_1,ix_2,b1_)*
     &               divb(ix_1,ix_2)
               enddo
               enddo
         else if(iw.eq.m2_)then
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-w(ix_1,ix_2,b2_)*
     &               divb(ix_1,ix_2)
               enddo
               enddo
           
         else if(iw.eq.b1_)then
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-w(ix_1,ix_2,m1_)/
     &               w(ix_1,ix_2,rho_)*divb(ix_1,ix_2)
               enddo
               enddo
         else if(iw.eq.b2_)then
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-w(ix_1,ix_2,m2_)/
     &               w(ix_1,ix_2,rho_)*divb(ix_1,ix_2)
               enddo
               enddo
           
         else if(iw.eq.e_)then
               do ix_2=ixOmin2,ixOmax2
               do ix_1=ixOmin1,ixOmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-(w(ix_1,ix_2,m1_)*
     &               w(ix_1,ix_2,b1_)+w(ix_1,ix_2,m2_)*w(ix_1,ix_2,b2_))/
     &               w(ix_1,ix_2,rho_)*divb(ix_1,ix_2)
               enddo
               enddo
         end if
      end do

      return
      end

*=============================================================================
      subroutine addgeometry(qdt,ixmin1,ixmin2,ixmax1,ixmax2,iws,w,wnew)

* Add geometrical source terms to wnew

      include 'vacdef.f'

      integer           ixmin1,ixmin2,ixmax1,ixmax2,iws(niw_),ix,iiw,iw,idir
      double precision  qdt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'addgeometry').ge.1

      if(typeaxial.eq.'slab'.or.typeaxial.eq.'test')then
*            No source terms in slab symmetry
      else if(typeaxial.eq.'nozzle')then
             call die('MHD for nozzle geometry is not implemented')
      else if(typeaxial.eq.'sphere')then
            do iiw=1,iws(niw_)
             iw=iws(iiw)
*              s[mr]=ptotal*2/r+(-Bphi**2-Btheta**2+mphi**2/rho+mtheta**2/rho)/r
            if(iw.eq.mr_)then
                  call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,tmp)
*                 This discretization maintains an exact equilibrium for p=const.
*                 areaside=(areaCi-areaCh)/areadx
                  do ix= ixmin1,ixmax1
                     do ix_1=ixmin2,ixmax2
                        wnew(ix,ix_1,iw)=wnew(ix,ix_1,iw) +qdt*tmp(ix,ix_1)*
     &                     areaside(ix)
                     enddo
                  enddo
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     tmp(ix_1,ix_2)=zero
                  enddo
                  enddo
                  do idir=2,ndir
                     do ix_2=ixmin2,ixmax2
                     do ix_1=ixmin1,ixmax1
                        tmp(ix_1,ix_2)=tmp(ix_1,ix_2)-w(ix_1,ix_2,b0_+idir)**
     &                     2 +w(ix_1,ix_2,m0_+idir)**2/w(ix_1,ix_2,rho_)
                     enddo
                     enddo
                  end do
*              s[mphi]=(-mphi*mr/rho+Bphi*Br)/radius and phi-->theta
            else if(iw.eq.m2_)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     tmp(ix_1,ix_2)= -w(ix_1,ix_2,m2_)*w(ix_1,ix_2,mr_)/
     &                  w(ix_1,ix_2,rho_) +w(ix_1,ix_2,b2_)*w(ix_1,ix_2,br_) 
                  enddo
                  enddo
*              s[Bphi]=((Bphi*mr-Br*mphi)/rho)/radius and phi-->theta
            else if(iw.eq.b2_)then
                 do ix_2=ixmin2,ixmax2
                 do ix_1=ixmin1,ixmax1
                    tmp(ix_1,ix_2)=(w(ix_1,ix_2,b2_)*w(ix_1,ix_2,mr_)-
     &                 w(ix_1,ix_2,br_)*w(ix_1,ix_2,m2_)) /w(ix_1,ix_2,rho_) 
                 enddo
                 enddo
               end if
*              Divide by radius and add to wnew for variables other than rho and e
               if(iw.ne.rho_ .and. iw.ne.e_)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+qdt*
     &                  tmp(ix_1,ix_2)/x(ix_1,ix_2,r_)
                  enddo
                  enddo
               endif
               if(oktest.and.iw.eq.iwtest)write(*,*)'Geometrical source:',
     &            tmp(ixtest1,ixtest2)
            end do
      else if(typeaxial.eq.'cylinder')then
            do iiw=1,iws(niw_)
             iw=iws(iiw)
*              s[mr]=(ptotal-Bphi**2+mphi**2/rho)/radius
            if(iw.eq.mr_)then
                  call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,tmp)
                  if(.not.gencoord)then
*                    For nonuniform Cartesian grid this provides hydrostatic equil.
                     do ix= ixmin1,ixmax1
                        do ix_1=ixmin2,ixmax2
                           wnew(ix,ix_1,iw)=wnew(ix,ix_1,iw) +qdt*
     &                        tmp(ix,ix_1)*areaside(ix)
                        enddo
                     enddo
                     do ix_2=ixmin2,ixmax2
                     do ix_1=ixmin1,ixmax1
                        tmp(ix_1,ix_2)=zero
                     enddo
                     enddo
                  endif

               end if
*              Divide by radius and add to wnew
               if(iw.eq.mr_.or.(iw.eq.mphi_.and..not.angmomfix).or.
     &            iw.eq.bphi_)then
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+qdt*
     &                  tmp(ix_1,ix_2)/x(ix_1,ix_2,r_)
                  enddo
                  enddo
               endif
               if(oktest.and.iw.eq.iwtest)write(*,*)'Geometrical source:',
     &            tmp(ixtest1,ixtest2)
            end do
      end if

      return
      end

*=============================================================================
* end module vacphys.mhd0
*##############################################################################

*=============================================================================
*=============================================================================
      subroutine keeppositive_rho(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Keep density positive

      include 'vacdef.f'

      integer           ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      logical  toosmallr
*-----------------------------------------------------------------------------

* Overwrite smallrho if it is not yet defined
      if(smallrho.eq. -one)then
          minval_1=bigdouble
          do ix_2=ixmin2,ixmax2
          do ix_1=ixmin1,ixmax1
             minval_1=min(minval_1,w(ix_1,ix_2,rho_))
          enddo
          enddo
          smallrho=minval_1*smallrhocoeff
          
      endif

* Check for very small density
      any_1=.false.
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         any_1=any_1 .or. w(ix_1,ix_2,rho_).lt.max(zero,smallrho)
      enddo
      enddo
      toosmallr=any_1

      if(toosmallr)then
         nerror(toosmallr_)=nerror(toosmallr_)+1
         if(nerror(toosmallr_).eq.1)then
            write(*,'(a,i2,a,i7)')'Too small density (code=',toosmallr_,
     &         ') at it=',it
            minval_1=bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               minval_1=min(minval_1,w(ix_1,ix_2,rho_))
            enddo
            enddo
            write(*,*)'Value < smallrho: ',minval_1,smallrho

            
         endif
         if(smallrho.gt.zero)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               w(ix_1,ix_2,rho_)=max(w(ix_1,ix_2,rho_),smallrho)
            enddo
            enddo
         endif
      endif

      return
      end
*=============================================================================

      subroutine keeppositive(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Keep pressure and density positive (following D.Ryu)

      include 'vacdef.f'

      integer           ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      logical  toosmallp
*-----------------------------------------------------------------------------

      if(vacuumrho.lt.zero)then
*        Keep density positive
         call keeppositive_rho(ixmin1,ixmin2,ixmax1,ixmax2,w)
      else
*        Where rho is small use vacuum state: rho=vacuumrho, v=0, p=smallp, same B
*!!      ^C&w(ix^S,m^C_)=w(ix^S,m^C_)/w(ix^S,rho_)*vacuumrho;
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            if(w(ix_1,ix_2,rho_).lt.smallrho)then
                  w(ix_1,ix_2,m1_)=zero
                  w(ix_1,ix_2,m2_)=zero
            
                  w(ix_1,ix_2,rho_)=vacuumrho
                  w(ix_1,ix_2,e_)=smallp/(eqpar(gamma_)-one)+half*
     &               (w(ix_1,ix_2,b1_)**2+w(ix_1,ix_2,b2_)**2)
               endif
         enddo
         enddo
      endif

* Calculate pressure without clipping toosmall values (.false.)
      call getpthermal(.false.,w,ixmin1,ixmin2,ixmax1,ixmax2,tmp)

      any_1=.false.
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         any_1=any_1 .or. tmp(ix_1,ix_2).lt.max(zero,smallp)
      enddo
      enddo
      toosmallp=any_1

      if(toosmallp)then
         nerror(toosmallp_)=nerror(toosmallp_)+1
         if(nerror(toosmallp_).eq.1)then
            write(*,'(a,i2,a,i7)')'Too small pressure (code=',toosmallp_,
     &         ') at it=',it
            minval_1=bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               minval_1=min(minval_1,tmp(ix_1,ix_2))
            enddo
            enddo
            write(*,*)'Value < smallp: ',minval_1,smallp
       
            
         endif
         if(smallp.gt.zero)then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               w(ix_1,ix_2,e_)=max(tmp(ix_1,ix_2),smallp)/(eqpar(gamma_)-1)+
     &            half*((w(ix_1,ix_2,m1_)**2+w(ix_1,ix_2,m2_)**
     &            2)/w(ix_1,ix_2,rho_)+(w(ix_1,ix_2,b1_)**2+w(ix_1,ix_2,b2_)**
     &            2))
            enddo
            enddo
         endif
      endif

      return
      end

*=============================================================================
* end module vacphys - mhd
*##############################################################################
