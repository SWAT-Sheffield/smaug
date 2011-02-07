*#############################################################################
* module vacgrid
* Subroutines for boundaries, grid, divergence of flux, gradients
* Also the limiter functions for TVD, TVDLF, TVDMU schemes


*=============================================================================
      subroutine setnozzle(r,rC)

* Set global geometric factors for a nozzle:
*
*   area(ixG^T)     cross section at cell centers.
*   areaC(ixG^T)    cross section at cell boundaries.
*   areadx(ixG^T)   cell volume (integrated cross section).
*   areaside(ixG^T) =(areaC_i+1/2 - areaC_i-1/2)/dr
*
* The input arrays 'r' and 'rC' are the locations of cell centers and
* cell boundaries, respectively.

* This is Yee's nozzle with a cross section 1.398+0.347*tanh(0.8*x-4) 

      include 'vacdef.f'

      double precision  r(IXGLO1-1:IXGHI1+1),rC(IXGLO1-1:IXGHI1+1)
      double precision  qareaC(IXGLO1-1:IXGHI1+1)
*----------------------------------------------------------------------------

      oktest=index(teststr,'setnozzle').ge.1

* ADAPTOR does not know tanh function, so tanh(x)=1-2/(exp(2*x)+1)
      do ix_1=ixGmin1,ixGmax1
         area(ix_1) =1.398+0.347*(1-2/(exp(1.6*r(ix_1) -8)+1))
      enddo
      do ix_1=ixGmin1-1,ixGmax1
         qareaC(ix_1)=1.398+0.347*(1-2/(exp(1.6*rC(ix_1)-8)+1))
      enddo
      do ix_1=ixGmin1,ixGmax1
         areaC(ix_1)=qareaC(ix_1)
      enddo
* We need to integrate areaC for exact cell volume, e.g. with Mathematica
      do ix_1=ixGmin1,ixGmax1
         areadx(ix_1)=1.398*(rC(ix_1)-rC(ix_1+(ixGmin1-1-ixGmin1)))+0.347/0.8*
     &       log(cosh(0.8*rC(ix_1)-4.)/cosh(0.8*rC(ix_1+(ixGmin1-1-ixGmin1))-
     &      4.))
      enddo
* We use qareaC to calculate areaside
      do ix_1=ixGmin1,ixGmax1
         areaside(ix_1)=(qareaC(ix_1)-qareaC(ix_1+(ixGmin1-1-ixGmin1)))/
     &      areadx(ix_1)
      enddo

      if(oktest)write(*,*)'SetNozzle area,areaC,areadx,areaside:',
     &   area(ixtest1),areaC(ixtest1),areadx(ixtest1),areaside(ixtest1)

      return 
      end
*=============================================================================
*=============================================================================
      subroutine boundsetup

* Variables describing boundaries
*
*    iB=1..nB                                    - boundary ID
*    ixBmin(idim,iB):ixBmax(idim,iD)             - location of boundary
*    idimB(iB)                                   - direction orthogonal to the 
*                                                  boundary layer
*    upperB(iB)                                  - .true. if boundary is at a
*                                                  high index of the phys. grid
*    typeB(iw,iB)                                - boundary type string

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,iB,jB,iw,idim,idm,ixGmin(ndim),
     &   ixGmax(ndim),ixMmin(ndim),ixMmax(ndim)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'boundsetup').ge.1
      if(oktest)write(*,*)'BoundSetup'

* If ixBmax(1,1)==0, the user did not specify boundary regions. Setup default.
* A schematic picture for ndim=2, where the numbers are iB for each region, and
* *-s are the extra layers for qx. 
*
*                     ************
*                     *1444444424*    ixGmax2
*                     *4144444442*
*                     *11      22*    ixMmax2
*                     *11      22*
*                     *11      22*    ixMmin2
*                     *1333333323*
*                     *3133333332*    ixGmin2
*                     ************
*                      i i    i i
*                      x x    x x
*                      G M    M G
*                      m m    m m
*                      i i    a a
*                      n n    x x
*                      1 1    1 1

      ixGmin(1)=ixGmin1
      ixGmax(1)=ixGmax1
      ixGmin(2)=ixGmin2
      ixGmax(2)=ixGmax2
      ixMmin(1)=ixMmin1
      ixMmax(1)=ixMmax1
      ixMmin(2)=ixMmin2
      ixMmax(2)=ixMmax2

      if(ixBmax(1,1).eq.0)then
         nB=2*ndim
         do iB=1,nB
            idim=(iB+1)/2
            idimB(iB) = idim
            upperB(iB)= 2*idim.eq.iB
            do idm=1,ndim
               ixBmin(idm,iB)=ixGmin(idm)
               ixBmax(idm,iB)=ixGmax(idm)
            end do
            if(upperB(iB))then
               ixBmin(idim,iB)=ixMmax(idim)+1
               ixBmax(idim,iB)=ixGmax(idim)
            else
               ixBmin(idim,iB)=ixGmin(idim)
               ixBmax(idim,iB)=ixMmin(idim)-1
            end if
         end do
      else
*        Check if boundary regions are inside grid and in between mesh and grid
         do iB=1,nB
            ixmin1=ixBmin(1,iB)
            ixmin2=ixBmin(2,iB)
            ixmax1=ixBmax(1,iB)
            ixmax2=ixBmax(2,iB)
            
            if(ixmin1.lt.ixGmin1.or.ixmax1.gt.ixGmax1.or.ixmin2.lt.ixGmin2.or.
     &         ixmax2.gt.ixGmax2) then
               write(*,*)'Error for boundary region iB,ixL=',iB,ixmin1,ixmin2,
     &            ixmax1,ixmax2
               call die(
     &            'Error in BoundSetup: Boundary region is outside grid')
            endif
         if(idimB(iB).eq.1)then
               if(upperB(iB))then
                  if(ixmin1-1.ne.ixMmax1.or.ixmax1.ne.ixGmax1)write(*,*
     &               )'Warning in BoundSetup: Boundary does not fit, iB:',iB
               else
                  if(ixmax1+1.ne.ixMmin1.or.ixmin1.ne.ixGmin1)write(*,*
     &               )'Warning in BoundSetup: Boundary does not fit, iB:',iB
               endif 
         else if(idimB(iB).eq.2)then
               if(upperB(iB))then
                  if(ixmin2-1.ne.ixMmax2.or.ixmax2.ne.ixGmax2)write(*,*
     &               )'Warning in BoundSetup: Boundary does not fit, iB:',iB
               else
                  if(ixmax2+1.ne.ixMmin2.or.ixmin2.ne.ixGmin2)write(*,*
     &               )'Warning in BoundSetup: Boundary does not fit, iB:',iB
               endif 
            end if
         end do
      end if

* Identify the periodic pairs if they are not defined in boundlist
* Check type, direction, orientation, and size before matching pairs
      do iB=1,nB
         if(typeB(1,iB).eq.'periodic'.and.ipairB(iB).eq.0)then
            do jB=iB+1,nB
               if(typeB(1,jB).eq.'periodic'.and.ipairB(jB).eq.0.and.
     &            idimB(iB).eq.idimB(jB).and.(upperB(iB).neqv.upperB(jB)).and.
     &            ixBmax(1,iB)-ixBmin(1,iB).eq.ixBmax(1,jB)-
     &            ixBmin(1,jB).and.ixBmax(2,iB)-ixBmin(2,iB).eq.ixBmax(2,jB)-
     &            ixBmin(2,jB))then
                  ipairB(iB)=jB
                   ipairB(jB)=iB
               endif
            end do
            if(ipairB(iB).eq.0)call die(
     &         'Error in BoundSetup: No periodic pair')
         end if
      end do




* symm0 means zero orthogonal flux via the boundary 
      do iw=1,nw
         do iB=1,nB
           if(typeB(iw,iB).eq.'symm0') nofluxB(iw,idimB(iB))=.true.
         enddo
      enddo

      if(oktest)then
         do iB=1,nB
            
            write(unitterm,*)'iB,idimB,upperB,typeB:',iB,idimB(iB),upperB(iB),
     &         ' ',typeB(1,iB)
            write(unitterm,*)'ixBmin:',(ixBmin(idm,iB),idm=1,ndim)
            write(unitterm,*)'ixBmax:',(ixBmax(idm,iB),idm=1,ndim)
         end do
         do iB=1,nB
            if(typeB(1,iB).eq.'periodic')write(*,*)'Pairs:',iB,ipairB(iB)
         end do
      end if

      return
      end

*=============================================================================
      subroutine ensurebound(dix,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,qt,w)

* Check if there is enough information for calculating derivatives.
* The requirement is that ixI is wider than ixO by dix. 
* Adjust ixI and ixO. Call getboundary if needed.

      include 'vacdef.f'

      integer  dix,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2
      double precision  qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'ensurebound').gt.0
      if(oktest)write(*,*)'EnsureBound dix,ixI,ixO:',dix,',',ixImin1,ixImin2,
     &   ixImax1,ixImax2,',',ixOmin1,ixOmin2,ixOmax1,ixOmax2

* Check wether ixO+dix is within the grid
      if(ixGmin1.gt.ixOmin1-dix.or.ixGmin2.gt.ixOmin2-dix.or.
     &   ixGmax1.lt.ixOmax1+dix.or.ixGmax2.lt.ixOmax2+dix)then
         ixOmin1=ixMmin1
         ixOmin2=ixMmin2
         ixOmax1=ixMmax1
         ixOmax2=ixMmax2
      endif
* Check whether ixI is wider than ixO by at least dix otherwise getboundary
      if(ixImin1.gt.ixOmin1-dix.or.ixImin2.gt.ixOmin2-dix.or.
     &   ixImax1.lt.ixOmax1+dix.or.ixImax2.lt.ixOmax2+dix)then
         ixImin1=ixGmin1
         ixImin2=ixGmin2
         ixImax1=ixGmax1
         ixImax2=ixGmax2
         call getboundary(qt,1,nw,1,ndim,w)
         if(oktest)write(*,*)'calling getboundary'
      end if

      if(oktest)write(*,*)'Final       dix,ixI,ixO:',dix,',',ixImin1,ixImin2,
     &   ixImax1,ixImax2,',',ixOmin1,ixOmin2,ixOmax1,ixOmax2

      return
      end

*=============================================================================
      subroutine getboundary(qt,iwmin,iwmax,idimmin,idimmax,w)

      include 'vacdef.f'

      integer  iwmin,iwmax,idimmin,idimmax
      double precision  qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
      integer  ix,ix1,ix2,ixe,ixf,ixmin1,ixmin2,ixmax1,ixmax2,ixpairmin1,
     &   ixpairmin2,ixpairmax1,ixpairmax2,idim,iw,iB
      integer  iwv,jdim
      double precision  coeffnormal,coefftransv
      logical  initialized
      data initialized/.false./
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getboundary').ge.1
      if(oktest)write(*,*)'GetBoundary it,step:',it,step
      if(oktest)write(*,*)'GetBoundary wold:',w(ixtest1,ixtest2,iwtest)

      if(extraB)call specialbound(qt,ixMmin1,ixMmin2,ixMmax1,ixMmax2,0,0,w)
      if(smallfix)call keeppositive(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)



      iB=0
7701  continue
         iB=iB+1
         if(oktest)write(*,*)'iB  :',iB
         if(iB.gt.nB) goto 7702  
         idim=idimB(iB)
*        Only boundary segments in the required direction(s) are filled in
         if(idimmin.gt.idim.or.idimmax.lt.idim) goto 7701  

         ixmin1=ixBmin(1,iB)
         ixmin2=ixBmin(2,iB)
         ixmax1=ixBmax(1,iB)
         ixmax2=ixBmax(2,iB)

*        The possibly shifted coordinates parallel to the boundary layer 
*        are defined by the PAIR of the periodic/overlapping boundary.
*        Put the location of the source of information into ixpairL.
         if(ipairB(iB).gt.0)then
            ixpairmin1=ixBmin(1,ipairB(iB))
            ixpairmin2=ixBmin(2,ipairB(iB))
            ixpairmax1=ixBmax(1,ipairB(iB))
            ixpairmax2=ixBmax(2,ipairB(iB))
         if(idim.eq.1)then
               if(upperB(iB))then
                  ixpairmin1=ixpairmin1+dixBmin1
                  ixpairmax1=ixpairmax1+dixBmax1
               else
                  ixpairmin1=ixpairmin1-dixBmin1
                  ixpairmax1=ixpairmax1-dixBmax1
               endif
            
         else if(idim.eq.2)then
               if(upperB(iB))then
                  ixpairmin2=ixpairmin2+dixBmin2
                  ixpairmax2=ixpairmax2+dixBmax2
               else
                  ixpairmin2=ixpairmin2-dixBmin2
                  ixpairmax2=ixpairmax2-dixBmax2
               endif
            
            end if
         endif

         do iw= iwmin,iwmax
            if(oktest)write(*,*)'  iw:',iw
            if(oktest)write(*,*)'typeB(iw,iB):',typeB(iw,iB)
         if(typeB(iw,iB).eq.'cont'.or.typeB(iw,iB).eq.'fixed')then
*              For 'cont' copy w at the edge into the whole boundary region.
*              Fot 'fixed' copy w first, later use values stored in fixB.
*              For fullgridini=T store the w values read from the file in fixB.
            if(idim.eq.1)then
                  if(upperB(iB))then
                     ixe=ixmin1-1
                  else
                     ixe=ixmax1+1
                  endif
                  if(fixedB(iw,iB))then
*                    HPF$ INDEPENDENT
                     do ix1=ixmin1,ixmax1
               do ix2=ixmin2,ixmax2
                        w(ix1,ix2,iw)=fixB1(ix1-ixe,ix2,iw)
                     enddo
               enddo 
                  else if(typeB(iw,iB).eq.'cont' .or. .not.fullgridini) then
*                    HPF$ INDEPENDENT
                     do ix= ixmin1,ixmax1
                        do ix_1=ixmin2,ixmax2
                           w(ix,ix_1,iw)=w(ixe,ix_1,iw)
                        enddo
                     end do
                  end if 
            else if(idim.eq.2)then
                  if(upperB(iB))then
                     ixe=ixmin2-1
                  else
                     ixe=ixmax2+1
                  endif
                  if(fixedB(iw,iB))then
*                    HPF$ INDEPENDENT
                     do ix1=ixmin1,ixmax1
               do ix2=ixmin2,ixmax2
                        w(ix1,ix2,iw)=fixB2(ix1,ix2-ixe,iw)
                     enddo
               enddo 
                  else if(typeB(iw,iB).eq.'cont' .or. .not.fullgridini) then
*                    HPF$ INDEPENDENT
                     do ix= ixmin2,ixmax2
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix,iw)=w(ix_1,ixe,iw)
                        enddo
                     end do
                  end if 
               end if
         else if (typeB(iw,iB).eq.'cont1'.or.typeB(iw,iB).eq.'fixed1'.or.
     &      typeB(iw,iB).eq.'grad1')then
*              First order extrapolation from edge and inner edge to the boundary
*              'cont1'  extrapolates in every time step, can be numericly unstable.
*              'fixed1' extrapolates first, stores VALUES into fixB, then restores.
*              'grad1'  extrapolates first, stores DIFFERENCES into fixB, then 
*              adds the stored differences to the current edge value.
            if(idim.eq.1)then
                  if(upperB(iB))then
                     ixe=ixmin1-1
                      ixf=ixe-1
                  else
                     ixe=ixmax1+1
                      ixf=ixe+1
                  endif
                  if(fixedB(iw,iB))then
                     if(typeB(iw,iB).eq.'grad1')then
*                       HPF$ INDEPENDENT
                        do ix1=ixmin1,ixmax1
               do ix2=ixmin2,ixmax2
                           w(ix1,ix2,iw)=fixB1(ix1-ixe,ix2,iw)+w(ixe,ix2,iw)
                        enddo
               enddo
                     else
*                       HPF$ INDEPENDENT
                        do ix1=ixmin1,ixmax1
               do ix2=ixmin2,ixmax2
                           w(ix1,ix2,iw)=fixB1(ix1-ixe,ix2,iw)
                        enddo
               enddo
                     endif
*                 HPF_
                  else if(typeB(iw,iB).eq.'cont1'.or. .not.fullgridini)then
*                 HPF_ endif
*                 HPF_ if(.not.fixedB(iw,iB).and.&
*                 HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
*                    HPF$ INDEPENDENT
                     do ix= ixmin1,ixmax1
                        do ix_1=ixmin2,ixmax2
                           w(ix,ix_1,iw)=(abs(ix-ixe)+1)*w(ixe,ix_1,iw)-
     &                         abs(ix-ixe)   *w(ixf,ix_1,iw)  
                        enddo
                     end do
                  end if 
            else if(idim.eq.2)then
                  if(upperB(iB))then
                     ixe=ixmin2-1
                      ixf=ixe-1
                  else
                     ixe=ixmax2+1
                      ixf=ixe+1
                  endif
                  if(fixedB(iw,iB))then
                     if(typeB(iw,iB).eq.'grad1')then
*                       HPF$ INDEPENDENT
                        do ix1=ixmin1,ixmax1
               do ix2=ixmin2,ixmax2
                           w(ix1,ix2,iw)=fixB2(ix1,ix2-ixe,iw)+w(ix1,ixe,iw)
                        enddo
               enddo
                     else
*                       HPF$ INDEPENDENT
                        do ix1=ixmin1,ixmax1
               do ix2=ixmin2,ixmax2
                           w(ix1,ix2,iw)=fixB2(ix1,ix2-ixe,iw)
                        enddo
               enddo
                     endif
*                 HPF_
                  else if(typeB(iw,iB).eq.'cont1'.or. .not.fullgridini)then
*                 HPF_ endif
*                 HPF_ if(.not.fixedB(iw,iB).and.&
*                 HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
*                    HPF$ INDEPENDENT
                     do ix= ixmin2,ixmax2
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix,iw)=(abs(ix-ixe)+1)*w(ix_1,ixe,iw)-
     &                         abs(ix-ixe)   *w(ix_1,ixf,iw)  
                        enddo
                     end do
                  end if 
               end if
         else if(typeB(iw,iB).eq.'periodic')then
*              Update boundary by translation of w by width of mesh (and shift)
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  w(ix_1,ix_2,iw)=w(ix_1+(ixpairmin1-ixmin1),ix_2+
     &               (ixpairmin2-ixmin2),iw)
               enddo
               enddo
         else if(typeB(iw,iB).eq.'symm'.or.typeB(iw,iB).eq.'symm0'.or.
     &      typeB(iw,iB).eq.'asymm')then
*              Reflect w into the boundary region, multiply by -1 for "asymm"
*              In generalized coordinates take into account the other vector
*              components for vector variables. The symmetry of the transverse 
*              component is based on typeB for the jdim=idim+1 -th component.
*              ixe is used for the reflection, normal vectors are taken at ixf+1/2
               if(gencoord.and.vectoriw(iw).ge.0)then
*                 Determine direction of vector component, symmetry coefficients
*                 for normal and transverse vector components
                  iwv=vectoriw(iw)
                   jdim=idim+1-(idim/ndim)*ndim 
                  coeffnormal=1
                   if(typeB(iwv+idim,iB).eq.'asymm') coeffnormal=-1
                  coefftransv=1
                   if(typeB(iwv+jdim,iB).eq.'asymm') coefftransv=-1
               endif
            if(idim.eq.1)then
                  if(upperB(iB))then
                     ixe=2*ixmin1-1
                      ixf=ixmin1-1
                  else
                     ixe=2*ixmax1+1
                      ixf=ixmax1
                  endif
                  if(gencoord.and.vectoriw(iw).ge.0)then
*                    HPF$ INDEPENDENT
                     do ix= ixmin1,ixmax1
                        do ix_1=ixmin2,ixmax2
                           w(ix,ix_1,iw)=zero
                        enddo
                        do jdim=1,ndim
                           do ix_1=ixmin2,ixmax2
                              w(ix,ix_1,iw)=w(ix,ix_1,iw)+normalC(ixf,ix_1,
     &                           idim,jdim)*w(ixe-ix,ix_1,iwv+jdim)
                           enddo
                        end do
                        do ix_1=ixmin2,ixmax2
                           w(ix,ix_1,iw)=w(ix,ix_1,iw)*normalC(ixf,ix_1,idim,
     &                        iw-iwv)*(coeffnormal-coefftransv)+
     &                        w(ixe-ix,ix_1,iw)*coefftransv
                        enddo
                     end do
                  else
*                    HPF$ INDEPENDENT
                     do ix= ixmin1,ixmax1
                        do ix_1=ixmin2,ixmax2
                           w(ix,ix_1,iw)=w(ixe-ix,ix_1,iw)
                        enddo
                     end do
                     if(typeB(iw,iB).eq.'asymm')then
                        do ix_2=ixmin2,ixmax2
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix_2,iw)=-w(ix_1,ix_2,iw)
                        enddo
                        enddo
                     endif
                  endif  
               
            else if(idim.eq.2)then
                  if(upperB(iB))then
                     ixe=2*ixmin2-1
                      ixf=ixmin2-1
                  else
                     ixe=2*ixmax2+1
                      ixf=ixmax2
                  endif
                  if(gencoord.and.vectoriw(iw).ge.0)then
*                    HPF$ INDEPENDENT
                     do ix= ixmin2,ixmax2
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix,iw)=zero
                        enddo
                        do jdim=1,ndim
                           do ix_1=ixmin1,ixmax1
                              w(ix_1,ix,iw)=w(ix_1,ix,iw)+normalC(ix_1,ixf,
     &                           idim,jdim)*w(ix_1,ixe-ix,iwv+jdim)
                           enddo
                        end do
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix,iw)=w(ix_1,ix,iw)*normalC(ix_1,ixf,idim,
     &                        iw-iwv)*(coeffnormal-coefftransv)+
     &                        w(ix_1,ixe-ix,iw)*coefftransv
                        enddo
                     end do
                  else
*                    HPF$ INDEPENDENT
                     do ix= ixmin2,ixmax2
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix,iw)=w(ix_1,ixe-ix,iw)
                        enddo
                     end do
                     if(typeB(iw,iB).eq.'asymm')then
                        do ix_2=ixmin2,ixmax2
                        do ix_1=ixmin1,ixmax1
                           w(ix_1,ix_2,iw)=-w(ix_1,ix_2,iw)
                        enddo
                        enddo
                     endif
                  endif  
               
               end if
         else if(typeB(iw,iB).eq.'special')then
*                 Skip special now, we do it after normal boundary type variables
*                 HPF_ if(.false.)write(*,*)'Avoiding xlhpf90 compiler bug'
            
         else
               write(uniterr,*)'Error in GetBoundary: boundary type(', iw,iB,
     &            ')=',typeB(iw,iB),' is not implemented!'
               call die(' ')
*           typeB(iw,iB)
            end if
*        next iw
         end do
*        Do special boundaries
         do iw= iwmin,iwmax
            if(oktest)write(*,*)'special, iw:',iw
            if(typeB(iw,iB).eq.'special')call specialbound(qt,ixmin1,ixmin2,
     &         ixmax1,ixmax2,iw,iB,w)
*        next iw
         end do
*     next iB
      goto 7701  
7702  continue

* Fixed boundaries (fixed,fixed1) or fixed gradients (grad1) are stored into
* fixB after all boundaries have been updated.
* This needs to be done in the very first time step only.
      if(.not.initialized)then
         initialized=.true.
         do iB= 1,nB
            ixmin1=ixBmin(1,iB)
            ixmin2=ixBmin(2,iB)
            ixmax1=ixBmax(1,iB)
            ixmax2=ixBmax(2,iB)
            do iw= iwmin,iwmax
               if( (typeB(iw,iB).eq.'fixed'.or.typeB(iw,iB).eq.'fixed1'.or.
     &            typeB(iw,iB).eq.'grad1') .and. .not.fixedB(iw,iB))then
                  fixedB(iw,iB)=.true.
                  if(idimB(iB).eq.1)then
                        if(upperB(iB))then
                           ixe=ixmin1-1
                        else
                           ixe=ixmax1+1
                        endif
                        if(typeB(iw,iB).eq.'grad1')then
*                          HPF$ INDEPENDENT
                           do ix1= ixmin1,ixmax1
                     do ix2= ixmin2,ixmax2
                            fixB1(ix1-ixe,ix2,iw)=w(ix1,ix2,iw)-w(ixe,ix2,iw)
                           enddo
                     enddo
                        else
*                          HPF$ INDEPENDENT
                           do ix1= ixmin1,ixmax1
                     do ix2= ixmin2,ixmax2
                            fixB1(ix1-ixe,ix2,iw)=w(ix1,ix2,iw)
                           enddo
                     enddo
                        endif
                     
                  else if(idimB(iB).eq.2)then
                        if(upperB(iB))then
                           ixe=ixmin2-1
                        else
                           ixe=ixmax2+1
                        endif
                        if(typeB(iw,iB).eq.'grad1')then
*                          HPF$ INDEPENDENT
                           do ix1= ixmin1,ixmax1
                     do ix2= ixmin2,ixmax2
                            fixB2(ix1,ix2-ixe,iw)=w(ix1,ix2,iw)-w(ix1,ixe,iw)
                           enddo
                     enddo
                        else
*                          HPF$ INDEPENDENT
                           do ix1= ixmin1,ixmax1
                     do ix2= ixmin2,ixmax2
                            fixB2(ix1,ix2-ixe,iw)=w(ix1,ix2,iw)
                           enddo
                     enddo
                        endif
                     
                  end if
               end if 
*           iw
            end do
*        iB
         end do
      end if

      if(oktest)write(*,*)'GetBoundary wnew:',w(ixtest1,ixtest2,iwtest)

      return 
      end

*=============================================================================
      subroutine setnoflux(iw,idim,ixmin1,ixmin2,ixmax1,ixmax2,fRC,ixRmin1,
     &   ixRmin2,ixRmax1,ixRmax2,fLC,ixLmin1,ixLmin2,ixLmax1,ixLmax2)

* Set flux in direction idim to zero for variable iw if typeB is 'symm0'
* in a boundary region

      include 'vacdef.f'

      integer  iw,idim,ixmin1,ixmin2,ixmax1,ixmax2,ixLmin1,ixLmin2,ixLmax1,
     &   ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2
      double precision  fRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2), fLC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  iB,ixe,ixBmin1,ixBmin2,ixBmax1,ixBmax2

*-----------------------------------------------------------------------------

      oktest=index(teststr,'setnoflux').ge.1

      do iB=1,nB
         if(typeB(iw,iB).eq.'symm0'.and.idimB(iB).eq.idim)then
            ixBmin1=ixBmin(1,iB)
            ixBmin2=ixBmin(2,iB)
            ixBmax1=ixBmax(1,iB)
            ixBmax2=ixBmax(2,iB)
*           Calculate edge index and set the appropriate flux to zero
         if(idim.eq.1)then
            if(upperB(iB))then
               ixe=ixBmin1-1+ixRmin1-ixmin1
               if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:', it,idim,iw,
     &            ixe,fRC(ixe,ixtest2)
               do ix_1=ixBmin2,ixBmax2
                  fRC(ixe,ix_1)=zero
               enddo
            else
               ixe=ixBmax1+1+ixLmin1-ixmin1
               if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:',it,idim,iw,ixe,
     &            fLC(ixe,ixtest2)
               do ix_1=ixBmin2,ixBmax2
                  fLC(ixe,ix_1)=zero
               enddo
            endif
            
         else if(idim.eq.2)then
            if(upperB(iB))then
               ixe=ixBmin2-1+ixRmin2-ixmin2
               if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:', it,idim,iw,
     &            ixe,fRC(ixtest1,ixe)
               do ix_1=ixBmin1,ixBmax1
                  fRC(ix_1,ixe)=zero
               enddo
            else
               ixe=ixBmax2+1+ixLmin2-ixmin2
               if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:',it,idim,iw,ixe,
     &            fLC(ixtest1,ixe)
               do ix_1=ixBmin1,ixBmax1
                  fLC(ix_1,ixe)=zero
               enddo
            endif
            
            end if
         endif
      enddo

      return
      end

*=============================================================================
      subroutine gridsetup1

* Cartesian or polar grid. Determine x at the boundaries.
* Determine often needed combinations of x, such as dx or dvolume.
* Determine variables for axial symmetry
*
* ixe          - edge coordinate of the grid touching the boundary region
* ixf          - coordinate inside of ixe
* qx           - x with an extended index range for calculation of dx

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,hxmin1,hxmin2,hxmax1,hxmax2,jxmin1,
     &   jxmin2,jxmax1,jxmax2
      integer  ix,ixe,ixf,idim,jdim
      double precision  qx(IXGlo1-1:IXGhi1+1,IXGlo2-1:IXGhi2+1,ndim)
      double precision  r(IXGLO1-1:IXGHI1+1),rC(IXGLO1-1:IXGHI1+1)

*-----------------------------------------------------------------------------

      oktest=index(teststr,'gridsetup').ge.1
      if(oktest)write(*,*)'GridSetup1'

* Calculate qx in the boundary layers by linearly extrapolating x from
* the touching edge of the grid (ixe) and the inner edge (ixf).


      do idim_3=1,ndim
      do ix_2=ixGlo2-1,ixGhi2+1
      do ix_1=ixGlo1-1,ixGhi1+1
         qx(ix_1,ix_2,idim_3)=zero
      enddo
      enddo
      enddo
      do idim_3=1,ndim
      do ix_2=ixGmin2,ixGmax2
      do ix_1=ixGmin1,ixGmax1
         qx(ix_1,ix_2,idim_3) = x(ix_1,ix_2,idim_3)
      enddo
      enddo
      enddo
      do idim=1,ndim
         ixmin1=ixGmin1-1
         ixmin2=ixGmin2-1
         ixmax1=ixGmax1+1
         ixmax2=ixGmax2+1
         if(idim.eq.1)then
*              First do the upper layers
               ixmax1=ixGmax1+1
                ixmin1=ixMmax1+1
               if(fullgridini ) ixmin1=ixGmax1+1
               ixe=ixmin1-1
                ixf=ixe-1
*              !!forall replaced by do for sake of sequentializing (Adaptor)
      	 do jdim=1,ndim
                  do ix= ixmin1,ixmax1
                     do ix_1=ixmin2,ixmax2
                        qx(ix,ix_1,jdim)=(1+abs(ixe-ix))*qx(ixe,ix_1,jdim)-
     &                      abs(ixe-ix) *qx(ixf,ix_1,jdim)
                     enddo
                  end do
               end do
*              Next do the lower layers
               ixmin1=ixGmin1-1
                ixmax1=ixMmin1-1
               if(fullgridini ) ixmax1=ixGmin1-1
               ixe=ixmax1+1
                ixf=ixe+1
      	 do jdim=1,ndim
                  do ix= ixmin1,ixmax1
                      do ix_1=ixmin2,ixmax2
                         qx(ix,ix_1,jdim)=(1+abs(ixe-ix))*qx(ixe,ix_1,jdim)-
     &                       abs(ixe-ix) *qx(ixf,ix_1,jdim)
                      enddo
                  end do
               end do 
         else if(idim.eq.2)then
*              First do the upper layers
               ixmax2=ixGmax2+1
                ixmin2=ixMmax2+1
               if(fullgridini ) ixmin2=ixGmax2+1
               ixe=ixmin2-1
                ixf=ixe-1
*              !!forall replaced by do for sake of sequentializing (Adaptor)
      	 do jdim=1,ndim
                  do ix= ixmin2,ixmax2
                     do ix_1=ixmin1,ixmax1
                        qx(ix_1,ix,jdim)=(1+abs(ixe-ix))*qx(ix_1,ixe,jdim)-
     &                      abs(ixe-ix) *qx(ix_1,ixf,jdim)
                     enddo
                  end do
               end do
*              Next do the lower layers
               ixmin2=ixGmin2-1
                ixmax2=ixMmin2-1
               if(fullgridini ) ixmax2=ixGmin2-1
               ixe=ixmax2+1
                ixf=ixe+1
      	 do jdim=1,ndim
                  do ix= ixmin2,ixmax2
                      do ix_1=ixmin1,ixmax1
                         qx(ix_1,ix,jdim)=(1+abs(ixe-ix))*qx(ix_1,ixe,jdim)-
     &                       abs(ixe-ix) *qx(ix_1,ixf,jdim)
                      enddo
                  end do
               end do 
         end if
      enddo

      do idim_3=1,ndim
      do ix_2=ixGmin2,ixGmax2
      do ix_1=ixGmin1,ixGmax1
         x(ix_1,ix_2,idim_3)=qx(ix_1,ix_2,idim_3)
      enddo
      enddo
      enddo

* Loop with ^D instead of idim to avoid an SP2 xlphf90 compiler error
* if qx is distributed. But it should not be
*{
*^D&jx^L=ixG^L+kr(^D,^DD); hx^L=ixG^L-kr(^D,^DD);
*   dx(ixG^S,^D)=half*(qx(jx^S,^D)-qx(hx^S,^D))
*   if(oktest)write(*,*)'dx,qxj,qxh:',dx(ixtest^DD,^D),&
*      qx(ixtest^DD+kr(^D,^DD),^D),qx(ixtest^DD-kr(^D,^DD),^D)
*\}

      do idim=1,ndim
         jxmin1=ixGmin1+kr(idim,1)
         jxmin2=ixGmin2+kr(idim,2)
         jxmax1=ixGmax1+kr(idim,1)
         jxmax2=ixGmax2+kr(idim,2)
         hxmin1=ixGmin1-kr(idim,1)
         hxmin2=ixGmin2-kr(idim,2)
         hxmax1=ixGmax1-kr(idim,1)
         hxmax2=ixGmax2-kr(idim,2)
         do ix_2=ixGmin2,ixGmax2
         do ix_1=ixGmin1,ixGmax1
            dx(ix_1,ix_2,idim)=half*(qx(ix_1+(jxmin1-ixGmin1),ix_2+
     &         (jxmin2-ixGmin2),idim)-qx(ix_1+(hxmin1-ixGmin1),ix_2+
     &         (hxmin2-ixGmin2),idim))
         enddo
         enddo
      end do

* Calculate geometrical factors for axial symmetry based on Boris FCT.
* Fluxes are multiplied by areaC. The cell volume is proportional to areadx.
* Gradient type source terms are multiplied by areaside = darea/areadx.

      if(oktest)write(*,*)'Start calculating geometrical factors'

      if(typeaxial.ne.'slab')then
*        r is the radial coordinate defined for ixGmin1..ixGmax1+1
         ixmin1=ixGmin1-1
         ixmax1=ixGmax1+1
         do ix_1=ixmin1,ixmax1
            r(ix_1)=qx(ix_1,ixGmin2,1)
         enddo
*        rC is the centered radial coordinate defined for ixGmin1-1/2..ixGmax1+1/2
         ixmin1=ixGmin1-1
         ixmax1=ixGmax1
         do ix_1=ixmin1,ixmax1
            rC(ix_1)=half*(qx(ix_1,ixGmin2,1)+qx(ix_1+(ixmin1+
     &         1-ixmin1),ixGmin2,1))
         enddo
      if(typeaxial.eq.'test')then
*           This is slab symmetry really, but uses the geometrical factors
            do ix_1=ixGmin1,ixGmax1
               area(ix_1)=one
            enddo
            do ix_1=ixGmin1,ixGmax1
               areaC(ix_1)=one
            enddo
            do ix_1=ixGmin1,ixGmax1
               areadx(ix_1)=rC(ix_1)-rC(ix_1+(ixGmin1-1-ixGmin1))
            enddo
            do ix_1=ixGmin1,ixGmax1
               areaside(ix_1)=zero
            enddo
      else if(typeaxial.eq.'cylinder')then
            do ix_1=ixGmin1,ixGmax1
               area(ix_1)=r(ix_1)
            enddo
            do ix_1=ixGmin1,ixGmax1
               areaC(ix_1)=rC(ix_1)
            enddo
            do ix_1=ixGmin1,ixGmax1
               areadx(ix_1)=(rC(ix_1)**2-rC(ix_1+(ixGmin1-1-ixGmin1))**2)/2
            enddo
            do ix_1=ixGmin1,ixGmax1
               areaside(ix_1)=(rC(ix_1)-rC(ix_1+(ixGmin1-1-ixGmin1)))/
     &            areadx(ix_1)
            enddo
      else if(typeaxial.eq.'sphere')then
            do ix_1=ixGmin1,ixGmax1
               area(ix_1)=r(ix_1)**2
            enddo
            do ix_1=ixGmin1,ixGmax1
               areaC(ix_1)=rC(ix_1)**2
            enddo
            do ix_1=ixGmin1,ixGmax1
               areadx(ix_1)=(rC(ix_1)**3-rC(ix_1+(ixGmin1-1-ixGmin1))**3)/3
            enddo
            do ix_1=ixGmin1,ixGmax1
               areaside(ix_1)=(rC(ix_1)**2-rC(ix_1+(ixGmin1-1-ixGmin1))**
     &            2)/areadx(ix_1)
            enddo
      else if(typeaxial.eq.'nozzle')then
            call setnozzle(r,rC)
      else
            call die('Error in GridSetup1: Unknown type of axial symmetry:'//
     &         typeaxial)
         end if
      endif

      if(oktest)write(*,*)'Start calculating cell volumes'

* Calculate volume of cells and total volume of mesh
      if(typeaxial.eq.'slab')then
         do ix_2=ixGmin2,ixGmax2
         do ix_1=ixGmin1,ixGmax1
            dvolume(ix_1,ix_2)= dx(ix_1,ix_2,1)*dx(ix_1,ix_2,2)
         enddo
         enddo
      else
         do ix= ixGmin1,ixGmax1
            do ix_1=ixGmin2,ixGmax2
               dvolume(ix,ix_1)= areadx(ix)*dx(ix,ix_1,2)
            enddo
         enddo
      endif
      sum_1=0.
      do ix_2=ixMmin2,ixMmax2
      do ix_1=ixMmin1,ixMmax1
         sum_1=sum_1+dvolume(ix_1,ix_2)
      enddo
      enddo
      volume=sum_1


* For polar grid dx_phi=r*dphi. 
      if(polargrid)then
         do ix_2=ixGmin2,ixGmax2
         do ix_1=ixGmin1,ixGmax1
            dx(ix_1,ix_2,pphi_)=x(ix_1,ix_2,r_)*dx(ix_1,ix_2,pphi_)
         enddo
         enddo
      endif

      if(oktest)write(*,*)'Finish GridSetup1'

      return 
      end

*=============================================================================
      subroutine addflux(qdt,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,fRC,ixRmin1,
     &   ixRmin2,ixRmax1,ixRmax2,fLC,ixLmin1,ixLmin2,ixLmax1,ixLmax2,wnew)

* If qdt>0 physical flux, otherwise dissipative flux
*
* For qdt>0:
*
* In generalized coordinates 
*    wnew=wnew-qdt*(fRC*surfaceRC_idim-fLC*surfaceLC_idim)/dvolume
* In Cartesian coordinates 
*    wnew=wnew-qdt*(fRC(ixR)-fLC(ixL))/dx_idim
* with possible geometrical terms in axial symmetry for idim==r_==1:
*    wnew=wnew-qdt*(areaRC(r)*fRC-areaLC(r)*fLC)/areadx(r)
* or to conserve angular momentum for the mphi_ variable
*    wnew=wnew-qdt*(areaRC(r)**2*fRC-areaLC(r)**2*fLC)/(area(r)*areadx(r))
* i.e. mphi*r is integrated instead of mphi. The geometrical source terms 
* are not needed in this case (for HD(ADIAB) and MHD(ISO) at least).
*
* For qdt<0:
*
* Both in Cartesian and generalized coordinates
*    wnew=wnew+(fRC(ixR)-fLC(ixL))/dvolume
* But for angular momentum conservation
*    wnew=wnew+(fRC(ixR)*(x(ixR,r_)+x(jxR,r_))-fLC(ixL)*(x(ixL,r_)+x(jxL,r_)))
*              *half/x(ix,r_)/dvolume

      include 'vacdef.f'

      double precision  qdt,fRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   fLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ixmin1,ixmin2,ixmax1,ixmax2,hxmin1,hxmin2,hxmax1,hxmax2,jxmin1,
     &   jxmin2,jxmax1,jxmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixLmin1,ixLmin2,
     &   ixLmax1,ixLmax2,iw,idim

      integer  ix
*-----------------------------------------------------------------------------

      oktest= index(teststr,'addflux').ge.1 .and. iw.eq.iwtest
      if(oktest)write(*,*)'AddFlux wold:',wnew(ixtest1,ixtest2,iwtest)
      if(oktest)write(*,*)'idim,ixR,ixL:',idim,ixRmin1,ixRmin2,ixRmax1,
     &   ixRmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2

*!! The SHIFT relations for ixL and ixR are not valid for the 
*!! (TVD) MacCormack scheme

* Set flux to zero if required by boundary condition
* In generalized coordinates addflux with a vector variable can only be 
* called from vaccd where the flux is zero anyway given the correct symmetries
      if(nofluxB(iw,idim).and. .not.(gencoord.and.vectoriw(iw).ge.
     &   0)) call setnoflux(iw,idim,ixmin1,ixmin2,ixmax1,ixmax2,fRC,ixRmin1,
     &   ixRmin2,ixRmax1,ixRmax2,fLC,ixLmin1,ixLmin2,ixLmax1,ixLmax2)

      if(qdt.ge.zero)then
*        Physical flux
         if(gencoord)then
*           SHIFT
*           S ixRmin1=ixmin1;
*           SHIFT MORE
*           S ixLmin1=ixmin1-1
*           SHIFT MORE
            hxmin1=ixmin1-kr(idim,1)
            hxmin2=ixmin2-kr(idim,2)
            hxmax1=ixmax1-kr(idim,1)
            hxmax2=ixmax2-kr(idim,2)
            if(iw.ne.mphi_.or..not.angmomfix)then
*              Add conservative fluxes in generalized coordinates:
*              d(surface*F)/dvolume
*              SHIFT BEGIN
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-qdt/
     &               dvolume(ix_1,ix_2) *(surfaceC(ix_1,ix_2,idim)*
     &               fRC(ix_1+(ixRmin1-ixmin1),ix_2+(ixRmin2-ixmin2))-
     &               surfaceC(ix_1+(hxmin1-ixmin1),ix_2+(hxmin2-ixmin2),idim)*
     &               fLC(ix_1+(ixLmin1-ixmin1),ix_2+(ixLmin2-ixmin2)))
               enddo
               enddo
*              SHIFT END
            else
*              Conserve angular momentum in generalized coordinates
*              (1/r)*d(r*surface*F)/dvolume
*              SHIFT MORE 
               jxmin1=ixmin1+kr(idim,1)
               jxmin2=ixmin2+kr(idim,2)
               jxmax1=ixmax1+kr(idim,1)
               jxmax2=ixmax2+kr(idim,2)
*              SHIFT BEGIN
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-half*
     &               qdt/dvolume(ix_1,ix_2)/x(ix_1,ix_2,r_) *
     &               (surfaceC(ix_1,ix_2,idim)*(x(ix_1+(jxmin1-ixmin1),ix_2+
     &               (jxmin2-ixmin2),r_)+x(ix_1,ix_2,r_))*fRC(ix_1+
     &               (ixRmin1-ixmin1),ix_2+(ixRmin2-ixmin2))-surfaceC(ix_1+
     &               (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2),idim)*
     &               (x(ix_1+(hxmin1-ixmin1),ix_2+(hxmin2-ixmin2),r_)+
     &               x(ix_1,ix_2,r_))*fLC(ix_1+(ixLmin1-ixmin1),ix_2+
     &               (ixLmin2-ixmin2)))
               enddo
               enddo
*              SHIFT END
            endif
         else if(typeaxial.eq.'slab'.or.idim.gt.1) then
*           Add conservative flux in Cartesian coordinates or Z direction:
*           SHIFT
*           S ixRmin1=ixmin1;
*           SHIFT MORE
*           S ixLmin1=ixmin1-1
*           SHIFT BEGIN
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)-qdt/dx(ix_1,ix_2,idim) *
     &            (fRC(ix_1+(ixRmin1-ixmin1),ix_2+(ixRmin2-ixmin2))-fLC(ix_1+
     &            (ixLmin1-ixmin1),ix_2+(ixLmin2-ixmin2)))
            enddo
            enddo
*           SHIFT END
         else if(iw.ne.mphi_.or..not.angmomfix) then
*           Add radial conservative flux in Cartesian coordinates with axial symm.:
*           (1/r)*d(r*F)/dr
*           SHIFT BEGIN
            do ix= ixmin1,ixmax1
               do ix_1=ixmin2,ixmax2
                  wnew(ix,ix_1,iw)=wnew(ix,ix_1,iw)-qdt/areadx(ix) *
     &               (areaC(ix)  *fRC(ix+ixRmin1-ixmin1,ix_1+
     &               (ixRmin2-ixmin2))-areaC(ix-1)*fLC(ix+ixLmin1-ixmin1,ix_1+
     &               (ixLmin2-ixmin2)))
               enddo
            enddo
*           SHIFT END
         else
*           Conserve angular momentum in Cartesian/polar coordinates:
*           (1/r**2)*d(r**2*F)/dr
            if(it.eq.itmin)write(*,*)'Ang.mom.conservation for variable',iw
*           SHIFT BEGIN
            do ix= ixmin1,ixmax1
               do ix_1=ixmin2,ixmax2
                  wnew(ix,ix_1,iw)=wnew(ix,ix_1,iw)-qdt/area(ix)/areadx(ix) *
     &               (areaC(ix)**2*  fRC(ix+ixRmin1-ixmin1,ix_1+
     &               (ixRmin2-ixmin2))-areaC(ix-1)**2*fLC(ix+
     &               ixLmin1-ixmin1,ix_1+(ixLmin2-ixmin2)))
               enddo
            enddo
*           SHIFT END
         endif

*        Store flux for fluxCT or fluxCD schemes
*        In generalized coordinates the fluxes of subroutine addflux are projected 
*        to the normal vector so they are not good estimates of the electric field 
*        and therefore the flux is collected elsewhere in the vaccd module
         if(.not.gencoord.and.index(typeconstrain,'flux').eq.1.and.
     &      iw.gt.b0_.and.iw.le.b0_+ndim.and.iw.ne.b0_+idim.and.
     &      istep.eq.nstep) call storeflux(qdt,fRC,ixLmin1,ixLmin2,ixRmax1,
     &      ixRmax2,idim,iw)
      else
*        Dissipative numerical flux
*        SHIFT
*        S ixRmin1=ixmin1;
*        SHIFT MORE
*        S ixLmin1=ixmin1-1
         if(iw.ne.mphi_.or.idim.gt.1.or..not.angmomfix) then
*           Normal case (no angular momentum fix)
*           SHIFT BEGIN
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+(fRC(ix_1+
     &            (ixRmin1-ixmin1),ix_2+(ixRmin2-ixmin2))-fLC(ix_1+
     &            (ixLmin1-ixmin1),ix_2+(ixLmin2-ixmin2)))/dvolume(ix_1,ix_2)
            enddo
            enddo
*           SHIFT END
            if((index(typeconstrain,'flux').gt.0.or.typeconstrain.eq.
     &         'ryuCT').and.iw.gt.b0_.and.iw.le.b0_+ndim.and.iw.ne.b0_+
     &         idim.and.istep.eq.nstep)call storeflux(-one,fRC,ixLmin1,
     &         ixLmin2,ixRmax1,ixRmax2,idim,iw)
         else
*           Ang.mom. conservation by changing the difference formula
*           (1/r)*d(r*R.phi)/dvolume
*           SHIFT MORE
            jxmin1=ixmin1+kr(1,r_)
            jxmin2=ixmin2+kr(2,r_)
            jxmax1=ixmax1+kr(1,r_)
            jxmax2=ixmax2+kr(2,r_)
*           SHIFT MORE
            hxmin1=ixmin1-kr(idim,1)
            hxmin2=ixmin2-kr(idim,2)
            hxmax1=ixmax1-kr(idim,1)
            hxmax2=ixmax2-kr(idim,2)
*           SHIFT BEGIN
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               wnew(ix_1,ix_2,iw)=wnew(ix_1,ix_2,iw)+(fRC(ix_1+
     &            (ixRmin1-ixmin1),ix_2+(ixRmin2-ixmin2))*(x(ix_1,ix_2,r_)+
     &            x(ix_1+(jxmin1-ixmin1),ix_2+(jxmin2-ixmin2),r_))-fLC(ix_1+
     &            (ixLmin1-ixmin1),ix_2+(ixLmin2-ixmin2))*(x(ix_1,ix_2,r_)+
     &            x(ix_1+(hxmin1-ixmin1),ix_2+(hxmin2-ixmin2),r_)))*
     &            half/x(ix_1,ix_2,r_)/dvolume(ix_1,ix_2)
            enddo
            enddo
*           SHIFT END
         endif
      endif

      if(oktest)write(*,*)'qdt,fRC,fLC:',qdt,fRC(ixRmin1-ixmin1+
     &   ixtest1,ixRmin2-ixmin2+ixtest2),fLC(ixLmin1-ixmin1+
     &   ixtest1,ixLmin2-ixmin2+ixtest2)
      if(oktest)write(*,*)'AddFlux wnew:',wnew(ixtest1,ixtest2,iwtest)

      return
      end

*=============================================================================
      subroutine storeflux(qdt,fC,ixmin1,ixmin2,ixmax1,ixmax2,idim,iw)

*!! This subroutine cannot use tmp or tmp2 !!!

* Store the fluxes added to the B field. 
* qdt<0 indicates a numerical flux.
* This is used by the flux-CT, transport-flux-CT, and flux-CD schemes

      include 'vacdef.f'

      double precision  qdt, fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim,iw

*-----------------------------------------------------------------------------

      oktest= index(teststr,'storeflux').ge.1
      if(oktest)write(*,*)'StoreFlux idim,iw,fC,qdt:',idim,iw,fC(ixtest1,
     &   ixtest2),qdt

      call die('Error: CT/CD module is off. Use setvac -on=ct and recompile!')

      return
      end

*=============================================================================
      subroutine gradient(realgrad,q,ixmin1,ixmin2,ixmax1,ixmax2,idir,gradq)

*!! This subroutine should not use tmp or tmp2

* Calculate gradient of q within ixL in Cartesian direction idir
* If realgrad is .true., add corrections from axial symmetry, otherwise 
* a term in divergence is calculated

      include 'vacdef.f'

      logical  realgrad
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idir
      double precision  q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)

      double precision  qC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ix,jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,hxmin2,hxmax1,hxmax2,
     &   ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,idim
*-----------------------------------------------------------------------------

      oktest= index(teststr,'gradient').ge.1
      if(oktest)write(*,*)'Gradient q(hx,ix,jx):',q(ixtest1-
     &   kr(idir,1),ixtest2-kr(idir,2)),q(ixtest1,ixtest2),q(ixtest1+
     &   kr(idir,1),ixtest2+kr(idir,2))

      if(gencoord)then
*        Integrate surface averaged q*normal_idim over all the surfaces
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            gradq(ix_1,ix_2)=zero
         enddo
         enddo
         do idim=1,ndim
*           SHIFT
            hxmin1=ixmin1-kr(idim,1)
            hxmin2=ixmin2-kr(idim,2)
            hxmax1=ixmax1-kr(idim,1)
            hxmax2=ixmax2-kr(idim,2)
            ixCmin1=hxmin1
            ixCmin2=hxmin2
            ixCmax1=ixmax1
            ixCmax2=ixmax2
*           SHIFT MORE
            jxCmin1=ixCmin1+kr(idim,1)
            jxCmin2=ixCmin2+kr(idim,2)
            jxCmax1=ixCmax1+kr(idim,1)
            jxCmax2=ixCmax2+kr(idim,2)
*           SHIFT BEGIN
            do ix_2=ixCmin2,ixCmax2
            do ix_1=ixCmin1,ixCmax1
               qC(ix_1,ix_2)=surfaceC(ix_1,ix_2,idim)*normalC(ix_1,ix_2,idim,
     &            idir)*half*(q(ix_1,ix_2)+q(ix_1+(jxCmin1-ixCmin1),ix_2+
     &            (jxCmin2-ixCmin2)))
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               gradq(ix_1,ix_2)=gradq(ix_1,ix_2)+(qC(ix_1,ix_2)-qC(ix_1+
     &            (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2)))/dvolume(ix_1,ix_2)
            enddo
            enddo
*           SHIFT END
         enddo
*        Contribution to the radial gradient from the ignored direction
         if(realgrad .and. idir.eq.r_ .and. typeaxial.eq.'cylinder')then
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               gradq(ix_1,ix_2)=gradq(ix_1,ix_2)-q(ix_1,ix_2)/x(ix_1,ix_2,r_)
            enddo
            enddo
         endif

      elseif(typeaxial.eq.'slab'.or.idir.gt.1)then
*        SHIFT
         jxmin1=ixmin1+kr(idir,1)
         jxmin2=ixmin2+kr(idir,2)
         jxmax1=ixmax1+kr(idir,1)
         jxmax2=ixmax2+kr(idir,2)
*        SHIFT MORE
         hxmin1=ixmin1-kr(idir,1)
         hxmin2=ixmin2-kr(idir,2)
         hxmax1=ixmax1-kr(idir,1)
         hxmax2=ixmax2-kr(idir,2)
*        SHIFT BEGIN
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            gradq(ix_1,ix_2)=half*(q(ix_1+(jxmin1-ixmin1),ix_2+
     &         (jxmin2-ixmin2))-q(ix_1+(hxmin1-ixmin1),ix_2+(hxmin2-ixmin2)))/
     &         dx(ix_1,ix_2,idir)
         enddo
         enddo
*        SHIFT END
      else
         do ix= ixmin1-1,ixmax1
            do ix_1=ixmin2,ixmax2
               qC(ix,ix_1)=areaC(ix)*half*(q(ix+1,ix_1)+q(ix,ix_1))
            enddo
         enddo
         do ix= ixmin1,ixmax1
            do ix_1=ixmin2,ixmax2
               gradq(ix,ix_1)=(qC(ix,ix_1)-qC(ix-1,ix_1))/areadx(ix)
            enddo
         enddo
         if(realgrad)then
            do ix= ixmin1,ixmax1
               do ix_1=ixmin2,ixmax2
                  gradq(ix,ix_1)=gradq(ix,ix_1)-q(ix,ix_1)*areaside(ix)
               enddo
            enddo
         endif
      endif
      if(oktest)write(*,*)'gradq:',gradq(ixtest1,ixtest2)

      return
      end

*=============================================================================
      subroutine gradient4(realgrad,q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)

*!! This subroutine should not use tmp or tmp2

* Calculate 4th order gradient of q within ixL in Cartesian direction idim

*!! We assume uniform Cartesian grid in slab symmetry for now

*!! If realgrad is .true., add corrections from axial symmetry, otherwise 
*!! a term in divergence is calculated

      include 'vacdef.f'

      logical  realgrad
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim
      double precision  q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)

      integer  kxmin1,kxmin2,kxmax1,kxmax2,jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,
     &   hxmin2,hxmax1,hxmax2,gxmin1,gxmin2,gxmax1,gxmax2
*-----------------------------------------------------------------------------

      if(gencoord)call die('Error: gradient4 does not work for gen.coords!')
      if(typeaxial.ne.'slab'.and..not.realgrad)call die(
     &   'Error: gradient4 does not work for axial symmetry yet!')

      oktest= index(teststr,'gradient').ge.1
      if(oktest)write(*,*)'Gradient4 q(kx,jx,hx,gx):',q(ixtest1+2*
     &   kr(idim,1),ixtest2+2*kr(idim,2)),q(ixtest1+kr(idim,1),ixtest2+
     &   kr(idim,2)),q(ixtest1-kr(idim,1),ixtest2-kr(idim,2)),q(ixtest1-2*
     &   kr(idim,1),ixtest2-2*kr(idim,2))

*SHIFT
      kxmin1=ixmin1+2*kr(idim,1)
      kxmin2=ixmin2+2*kr(idim,2)
      kxmax1=ixmax1+2*kr(idim,1)
      kxmax2=ixmax2+2*kr(idim,2)
       
*SHIFT MORE
      jxmin1=ixmin1+kr(idim,1)
      jxmin2=ixmin2+kr(idim,2)
      jxmax1=ixmax1+kr(idim,1)
      jxmax2=ixmax2+kr(idim,2)
       
*SHIFT MORE
      hxmin1=ixmin1-kr(idim,1)
      hxmin2=ixmin2-kr(idim,2)
      hxmax1=ixmax1-kr(idim,1)
      hxmax2=ixmax2-kr(idim,2)
*SHIFT MORE
      gxmin1=ixmin1-2*kr(idim,1)
      gxmin2=ixmin2-2*kr(idim,2)
      gxmax1=ixmax1-2*kr(idim,1)
      gxmax2=ixmax2-2*kr(idim,2)

*SHIFT BEGIN
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         gradq(ix_1,ix_2)=-(q(ix_1+(kxmin1-ixmin1),ix_2+(kxmin2-ixmin2))-8*
     &      (q(ix_1+(jxmin1-ixmin1),ix_2+(jxmin2-ixmin2))-q(ix_1+
     &      (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2)))-q(ix_1+(gxmin1-
     &      ixmin1),ix_2+(gxmin2-ixmin2)))/dx(ix_1,ix_2,idim)/12
      enddo
      enddo
*SHIFT END

      if(oktest)write(*,*)'gradq:',gradq(ixtest1,ixtest2)

      return
      end

*=============================================================================
      subroutine laplace4(q,ixmin1,ixmin2,ixmax1,ixmax2,laplaceq)

*!! This subroutine should not use tmp or tmp2

* Calculate 4th order laplace of q within ixL in Cartesian direction idim

*!! We assume uniform Cartesian grid in slab symmetry for now

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),laplaceq(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)

      integer  idim,kxmin1,kxmin2,kxmax1,kxmax2,jxmin1,jxmin2,jxmax1,jxmax2,
     &   hxmin1,hxmin2,hxmax1,hxmax2,gxmin1,gxmin2,gxmax1,gxmax2
*-----------------------------------------------------------------------------

      if(gencoord)call die('Error: laplace4 does not work for gen.coords!')
      if(typeaxial.ne.'slab')call die(
     &   'Error: laplace4 does not work for axial symmetry yet!')

      oktest= index(teststr,'lpalace').ge.1

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         laplaceq(ix_1,ix_2)=zero
      enddo
      enddo

      do idim=1,ndim
*        SHIFT
         kxmin1=ixmin1+2*kr(idim,1)
         kxmin2=ixmin2+2*kr(idim,2)
         kxmax1=ixmax1+2*kr(idim,1)
         kxmax2=ixmax2+2*kr(idim,2)
          
*        SHIFT MORE
         jxmin1=ixmin1+kr(idim,1)
         jxmin2=ixmin2+kr(idim,2)
         jxmax1=ixmax1+kr(idim,1)
         jxmax2=ixmax2+kr(idim,2)
          
*        SHIFT MORE
         hxmin1=ixmin1-kr(idim,1)
         hxmin2=ixmin2-kr(idim,2)
         hxmax1=ixmax1-kr(idim,1)
         hxmax2=ixmax2-kr(idim,2)
*        SHIFT MORE
         gxmin1=ixmin1-2*kr(idim,1)
         gxmin2=ixmin2-2*kr(idim,2)
         gxmax1=ixmax1-2*kr(idim,1)
         gxmax2=ixmax2-2*kr(idim,2)

*        SHIFT BEGIN
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            laplaceq(ix_1,ix_2)=laplaceq(ix_1,ix_2)+(q(ix_1+
     &         (kxmin1-ixmin1),ix_2+(kxmin2-ixmin2))+q(ix_1+
     &         (gxmin1-ixmin1),ix_2+(gxmin2-ixmin2))+30*q(ix_1,ix_2)-16*
     &         (q(ix_1+(jxmin1-ixmin1),ix_2+(jxmin2-ixmin2))+q(ix_1+
     &         (hxmin1-ixmin1),ix_2+(hxmin2-ixmin2))))/dx(ix_1,ix_2,idim)**
     &         2/12
         enddo
         enddo
*        SHIFT END

         if(oktest)write(*,*)'idim,q(kx,jx,ix,hx,gx):',idim,q(ixtest1+2*
     &      kr(idim,1),ixtest2+2*kr(idim,2)),q(ixtest1+kr(idim,1),ixtest2+
     &      kr(idim,2)),q(ixtest1,ixtest2),q(ixtest1-kr(idim,1),ixtest2-
     &      kr(idim,2)),q(ixtest1-2*kr(idim,1),ixtest2-2*kr(idim,2))
      enddo

      if(oktest)write(*,*)'laplaceq:',laplaceq(ixtest1,ixtest2)

      return
      end

*=============================================================================
      subroutine dwlimiter2(dwC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idim,ldw)

* Limit the centered dwC differences within ixC for iw in direction idim.
* The limiter is chosen according to typelimiter.

      include 'vacdef.f'

      double precision dwC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),ldw(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,iw,idim
      double precision  qsmall, qsmall2
      PARAMETER(  qsmall=1.D-12, qsmall2=2.D-12)
*-----------------------------------------------------------------------------

* Contract indices in idim for output.
      ixOmin1=ixCmin1+kr(idim,1)
      ixOmin2=ixCmin2+kr(idim,2)
       ixOmax1=ixCmax1
      ixOmax2=ixCmax2
*SHIFT
      hxOmin1=ixOmin1-kr(idim,1)
      hxOmin2=ixOmin2-kr(idim,2)
      hxOmax1=ixOmax1-kr(idim,1)
      hxOmax2=ixOmax2-kr(idim,2)

* Store the sign of dwC in tmp
      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         tmp(ix_1,ix_2)=sign(one,dwC(ix_1,ix_2))
      enddo
      enddo

*SHIFT BEGIN
      if(typelimiter(iw).eq.'muscl1')then
*           Minmod limiter eq.3.38c
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldw(ix_1,ix_2)=tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &            ix_2)),tmp(ix_1,ix_2)*musclomega*dwC(ix_1+
     &            (hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))))
            enddo
            enddo
      else if(typelimiter(iw).eq.'muscl2')then
*           Minmod limiter eq.3.38d
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldw(ix_1,ix_2)=tmp(ix_1,ix_2)* max(zero,min(musclomega*
     &            abs(dwC(ix_1,ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+
     &            (hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))))
            enddo
            enddo
      else if(typelimiter(iw).eq.'minmod')then
*           Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldw(ix_1,ix_2)=tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &            ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &            (hxOmin2-ixOmin2))))
            enddo
            enddo
      else if(typelimiter(iw).eq.'woodward')then
*           Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldw(ix_1,ix_2)=2*tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &            ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &            (hxOmin2-ixOmin2)),tmp(ix_1,ix_2)*0.25*(dwC(ix_1+
     &            (hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))+dwC(ix_1,ix_2))))
            enddo
            enddo
      else if(typelimiter(iw).eq.'superbee')then
*           Roe's superbee limiter (eq.3.51i)
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldw(ix_1,ix_2)=tmp(ix_1,ix_2)* max(zero,min(2*
     &            abs(dwC(ix_1,ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+
     &            (hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))),min(abs(dwC(ix_1,
     &            ix_2)),2*tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &            (hxOmin2-ixOmin2))))
            enddo
            enddo
      else if(typelimiter(iw).eq.'vanleer')then
*          van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
           do ix_2=ixOmin2,ixOmax2
           do ix_1=ixOmin1,ixOmax1
              ldw(ix_1,ix_2)=2*max(dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &           (hxOmin2-ixOmin2))*dwC(ix_1,ix_2),0d0)/(dwC(ix_1,ix_2)+
     &           dwC(ix_1+(hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))+qsmall)
           enddo
           enddo
      else if(typelimiter(iw).eq.'albada')then
*          Albada limiter (eq.3.51g) with delta2=1D.-12
           do ix_2=ixOmin2,ixOmax2
           do ix_1=ixOmin1,ixOmax1
              ldw(ix_1,ix_2)=(dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &           (hxOmin2-ixOmin2))*(dwC(ix_1,ix_2)**2+qsmall)+dwC(ix_1,ix_2)*
     &           (dwC(ix_1+(hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))**
     &           2+qsmall))/(dwC(ix_1,ix_2)**2+dwC(ix_1+(hxOmin1-
     &           ixOmin1),ix_2+(hxOmin2-ixOmin2))**2+qsmall2)
           enddo
           enddo
      else
            call die('Error in dwLimiter: No such TVD limiter:'//
     &         typelimiter(iw))
      end if
*SHIFT END
                              
      return
      end

*=============================================================================
      subroutine dwlimiter3(dwC,ixImin1,ixImin2,ixImax1,ixImax2,iw,idim,ldwC)

* Limit dwC differences within ixI for iw in direction idim.
* The limiter is chosen according to typelimiter.
* These limiters are used for Yee's symmetric TVD scheme.

      include 'vacdef.f'

      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,jxOmin1,jxOmin2,jxOmax1,jxOmax2,hxOmin1,hxOmin2,hxOmax1,
     &   hxOmax2,iw,idim
      double precision dwC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),ldwC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

* Contract indices in idim for output.
      ixOmin1=ixImin1+kr(idim,1)
      ixOmin2=ixImin2+kr(idim,2)
      ixOmax1=ixImax1-kr(idim,1)
      ixOmax2=ixImax2-kr(idim,2)
*SHIFT
      hxOmin1=ixOmin1-kr(idim,1)
      hxOmin2=ixOmin2-kr(idim,2)
      hxOmax1=ixOmax1-kr(idim,1)
      hxOmax2=ixOmax2-kr(idim,2)
*SHIFT MORE
      jxOmin1=ixOmin1+kr(idim,1)
      jxOmin2=ixOmin2+kr(idim,2)
      jxOmax1=ixOmax1+kr(idim,1)
      jxOmax2=ixOmax2+kr(idim,2)

* Store the sign of dwC in tmp
      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         tmp(ix_1,ix_2)=sign(one,dwC(ix_1,ix_2))
      enddo
      enddo

*SHIFT BEGIN
      if(typelimiter(iw).eq.'minmod')then
*           Minmod limiters according to eq(3.53b)
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldwC(ix_1,ix_2)=tmp(ix_1,ix_2)* (max(zero,min(abs(dwC(ix_1,
     &            ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &            (hxOmin2-ixOmin2))))+max(zero,min(abs(dwC(ix_1,ix_2)),
     &            tmp(ix_1,ix_2)*dwC(ix_1+(jxOmin1-ixOmin1),ix_2+
     &            (jxOmin2-ixOmin2)))))-dwC(ix_1,ix_2)
            enddo
            enddo
      else if(typelimiter(iw).eq.'woodward')then
*           Minmod limiter with 3 arguments according to eq(3.53c)
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldwC(ix_1,ix_2)=tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &            ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &            (hxOmin2-ixOmin2)),tmp(ix_1,ix_2)*dwC(ix_1+
     &            (jxOmin1-ixOmin1),ix_2+(jxOmin2-ixOmin2))))
            enddo
            enddo
      else if(typelimiter(iw).eq.'superbee')then
*           Minmod limiter with 4 arguments according to eq(3.53d)
            do ix_2=ixOmin2,ixOmax2
            do ix_1=ixOmin1,ixOmax1
               ldwC(ix_1,ix_2)=2*tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &            ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &            (hxOmin2-ixOmin2)),tmp(ix_1,ix_2)*dwC(ix_1+
     &            (jxOmin1-ixOmin1),ix_2+(jxOmin2-ixOmin2)),tmp(ix_1,ix_2)*
     &            0.25*(dwC(ix_1+(hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))+
     &            dwC(ix_1+(jxOmin1-ixOmin1),ix_2+(jxOmin2-ixOmin2)))))
            enddo
            enddo
      else
            call die('Error in dwCLimiter: No such TVD limiter:'//
     &         typelimiter(iw))
      end if
*SHIFT END
                  
      return
      end

*=============================================================================
      subroutine dwlimiterroe(adtdxC,dwC,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iw,
     &   idim,ldw)

* Limit the centered dwC differences within ixC for iw in direction idim.
* The limiter is chosen according to typelimiter from
* the "well-behaved" limiters of Arora and Roe JCP 132, 3

      include 'vacdef.f'

      double precision dwC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),adtdxC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),ldw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ixOmin1,ixOmin2,ixOmax1,ixOmax2,hxOmin1,hxOmin2,hxOmax1,
     &   hxOmax2,jxOmin1,jxOmin2,jxOmax1,jxOmax2,iw,idim
*-----------------------------------------------------------------------------

      oktest=index(teststr,'limiterroe').ge.1
      if(oktest)write(*,*)'dwLimiterRoe'

*SHIFT
      hxOmin1=ixOmin1-kr(idim,1)
      hxOmin2=ixOmin2-kr(idim,2)
      hxOmax1=ixOmax1-kr(idim,1)
      hxOmax2=ixOmax2-kr(idim,2)
*SHIFT MORE
      jxOmin1=ixOmin1+kr(idim,1)
      jxOmin2=ixOmin2+kr(idim,2)
      jxOmax1=ixOmax1+kr(idim,1)
      jxOmax2=ixOmax2+kr(idim,2)

* Store the sign of dwC in tmp
      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         tmp(ix_1,ix_2)=sign(one,dwC(ix_1,ix_2))
      enddo
      enddo

*SHIFT BEGIN
      if(typelimiter(iw).eq.'superroe')then
         do ix_2=ixOmin2,ixOmax2
         do ix_1=ixOmin1,ixOmax1
            if(adtdxC(ix_1,ix_2).gt.zero)then
                  ldw(ix_1,ix_2)=2*tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &               ix_2))/(1-adtdxC(ix_1,ix_2)),tmp(ix_1,ix_2)*
     &               dwC(ix_1+(hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))/
     &               (smalldouble+adtdxC(ix_1,ix_2)),tmp(ix_1,ix_2)*
     &               ((1+adtdxC(ix_1,ix_2))*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &               (hxOmin2-ixOmin2))+(2-adtdxC(ix_1,ix_2))*
     &               dwC(ix_1,ix_2))/6))
               else
                  ldw(ix_1,ix_2)=2*tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &               ix_2))/(1+adtdxC(ix_1,ix_2)),tmp(ix_1,ix_2)*
     &               dwC(ix_1+(jxOmin1-ixOmin1),ix_2+(jxOmin2-ixOmin2))/
     &               (smalldouble-adtdxC(ix_1,ix_2)),tmp(ix_1,ix_2)*
     &               ((1-adtdxC(ix_1,ix_2))*dwC(ix_1+(jxOmin1-ixOmin1),ix_2+
     &               (jxOmin2-ixOmin2))+(2+adtdxC(ix_1,ix_2))*
     &               dwC(ix_1,ix_2))/6))
               endif
         enddo
         enddo
      else if(typelimiter(iw).eq.'roe')then
*        Arora and Roe JCP 132, 3 eq. 2.0,2.1 with s1=PHImax=2
         do ix_2=ixOmin2,ixOmax2
         do ix_1=ixOmin1,ixOmax1
            if(adtdxC(ix_1,ix_2).gt.zero)then
                  ldw(ix_1,ix_2)=2*tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &               ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(hxOmin1-ixOmin1),ix_2+
     &               (hxOmin2-ixOmin2)),tmp(ix_1,ix_2)*((1+adtdxC(ix_1,ix_2))*
     &               dwC(ix_1+(hxOmin1-ixOmin1),ix_2+(hxOmin2-ixOmin2))+
     &               (2-adtdxC(ix_1,ix_2))*dwC(ix_1,ix_2))/6))
               else
                   ldw(ix_1,ix_2)=2*tmp(ix_1,ix_2)* max(zero,min(abs(dwC(ix_1,
     &                ix_2)),tmp(ix_1,ix_2)*dwC(ix_1+(jxOmin1-ixOmin1),ix_2+
     &                (jxOmin2-ixOmin2)),tmp(ix_1,ix_2)*((1-
     &                adtdxC(ix_1,ix_2))*dwC(ix_1+(jxOmin1-ixOmin1),ix_2+
     &                (jxOmin2-ixOmin2))+(2+adtdxC(ix_1,ix_2))*
     &                dwC(ix_1,ix_2))/6))
               endif
         enddo
         enddo
      end if
*SHIFT END

      return
      end

*=============================================================================
      subroutine acmswitch(jumpC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,phiC)

* Reduce phiC in ixC by the ACM switch which is calculated from jumpC at ixIC
* The ACM switch parameters are the exponent acmexpo and the width acmwidth

      include 'vacdef.f'

      double precision jumpC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),phiC(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      integer  ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim

      integer  jxCmin1,jxCmin2,jxCmax1,jxCmax2,hxCmin1,hxCmin2,hxCmax1,
     &   hxCmax2,ixmin1,ixmin2,ixmax1,ixmax2,hxmin1,hxmin2,hxmax1,hxmax2,
     &   iwidth
*-----------------------------------------------------------------------------

      oktest=index(teststr,'acmswitch').ge.1
      if(oktest)write(*,*)'ACMswitch phiC, ixC:',phiC(ixtest1,ixtest2),
     &   ixCmin1,ixCmin2,ixCmax1,ixCmax2

*SHIFT
      jxCmin1=ixCmin1+kr(idim,1)
      jxCmin2=ixCmin2+kr(idim,2)
      jxCmax1=ixCmax1+kr(idim,1)
      jxCmax2=ixCmax2+kr(idim,2)
*SHIFT MORE
      hxCmin1=ixCmin1-kr(idim,1)
      hxCmin2=ixCmin2-kr(idim,2)
      hxCmax1=ixCmax1-kr(idim,1)
      hxCmax2=ixCmax2-kr(idim,2)

      ixmin1=ixCmin1
      ixmin2=ixCmin2
       ixmax1=jxCmax1
      ixmax2=jxCmax2
*SHIFT MORE
      hxmin1=ixmin1-kr(idim,1)
      hxmin2=ixmin2-kr(idim,2)
      hxmax1=ixmax1-kr(idim,1)
      hxmax2=ixmax2-kr(idim,2)

      do ix_2=ixCmin2,ixCmax2
      do ix_1=ixCmin1,ixCmax1
         tmp(ix_1,ix_2)=abs(jumpC(ix_1,ix_2))
      enddo
      enddo

* Add two bondary values
*!! For periodic boundary, we should copy from other side !!!
      if(idim.eq.1)then
          do ix_1=ixCmin2,ixCmax2
             tmp(ixCmin1-1,ix_1)=tmp(ixCmin1,ix_1)
          enddo
          do ix_1=ixCmin2,ixCmax2
             tmp(ixCmax1+1,ix_1)=tmp(ixCmax1,ix_1) 
          enddo
      else if(idim.eq.2)then
          do ix_1=ixCmin1,ixCmax1
             tmp(ix_1,ixCmin2-1)=tmp(ix_1,ixCmin2)
          enddo
          do ix_1=ixCmin1,ixCmax1
             tmp(ix_1,ixCmax2+1)=tmp(ix_1,ixCmax2) 
          enddo
      end if

*SHIFT BEGIN
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         tmp2(ix_1,ix_2)=abs((tmp(ix_1,ix_2)-tmp(ix_1+(hxmin1-ixmin1),ix_2+
     &      (hxmin2-ixmin2)))/(tmp(ix_1,ix_2)+tmp(ix_1+(hxmin1-ixmin1),ix_2+
     &      (hxmin2-ixmin2))+smalldouble))**acmexpo
      enddo
      enddo

      do ix_2=ixCmin2,ixCmax2
      do ix_1=ixCmin1,ixCmax1
         tmp(ix_1,ix_2)=max(tmp2(ix_1,ix_2),tmp2(ix_1+(jxCmin1-ixCmin1),ix_2+
     &      (jxCmin2-ixCmin2)))
      enddo
      enddo
*SHIFT END

*Extend stencil of theta_j+1/2 by repeatedly taking maximum of neighbors
*This part is executed for acmwidth>1 only
      do iwidth=1,acmwidth-1
         do ix_2=ixCmin2,ixCmax2
         do ix_1=ixCmin1,ixCmax1
            tmp2(ix_1,ix_2)=tmp(ix_1,ix_2)
         enddo
         enddo
*        Fill in extra rows of tmp2 to avoid problems at boundaries.
*        !! For periodic boundary, we should copy from other side !!!
      if(idim.eq.1)then
            do ix_1=ixCmin2,ixCmax2
               tmp2(ixCmin1-1,ix_1)=tmp(ixCmin1,ix_1)
            enddo
            do ix_1=ixCmin2,ixCmax2
               tmp2(ixCmax1+1,ix_1)=tmp(ixCmax1,ix_1) 
            enddo
      else if(idim.eq.2)then
            do ix_1=ixCmin1,ixCmax1
               tmp2(ix_1,ixCmin2-1)=tmp(ix_1,ixCmin2)
            enddo
            do ix_1=ixCmin1,ixCmax1
               tmp2(ix_1,ixCmax2+1)=tmp(ix_1,ixCmax2) 
            enddo
         end if
*        Take maximum of neighbors
*        SHIFT BEGIN
         do ix_2=ixCmin2,ixCmax2
         do ix_1=ixCmin1,ixCmax1
            tmp(ix_1,ix_2)=max(tmp2(ix_1+(hxCmin1-ixCmin1),ix_2+
     &         (hxCmin2-ixCmin2)),tmp2(ix_1+(jxCmin1-ixCmin1),ix_2+
     &         (jxCmin2-ixCmin2)))
         enddo
         enddo
*        SHIFT END
      enddo

*Multiply by theta_j+1/2 coefficient
      do ix_2=ixCmin2,ixCmax2
      do ix_1=ixCmin1,ixCmax1
         phiC(ix_1,ix_2)=phiC(ix_1,ix_2)*tmp(ix_1,ix_2)
      enddo
      enddo

      if(oktest)write(*,*)'ACM phiC, thetai, thetaj:',phiC(ixtest1,ixtest2),
     &   tmp2(ixtest1,ixtest2),tmp2(ixtest1+kr(idim,1),ixtest2+kr(idim,2))

      return
      end

*=============================================================================
* end module vacgrid
*##############################################################################



