*##############################################################################
* module vacpoisson 
* Poisson solver and iterative schemes 
* for a scalar quantity defined on the mesh
*=============================================================================
      subroutine poisson(purpose,rhs,tolerance,typestop,matvecmax,info,
     &   nonzero,phi)

* Solve grad div phi = rhs for the purpose defined in that string
* The right hand side of the Poisson equation is in rhs
* The initial guess for the solution is given by phi and nonzero
*    nonzero=.true. : the initial guess is non-zero and it is in phi
*    nonzero=.false.: the initial guess is zero and phi=0
* The required accuracy of the solution is given by tolerance and typestop:
*    typestop='max': resid=max(abs(Laplace(phi)-rhs)) < tolerance
*    typestop='abs': resid=sum((Laplace(phi)-rhs)**2) < tolerance
*    typestop='rel': resid=sum((Laplace(phi)-rhs)**2) < tolerance*resid_init
* At most matvecmax matrix vector multiplications can be performed
* The success of the iteration is given by returned value of info:
*     abs(info)=  0 - solution found satisfying given tolerance.
*                 1 - iteration aborted due to division by very small value.
*                 2 - no convergence within maximum number of iterations.
*                 3 - initial guess satisfies the stopping criterion.
*    sign(info)=  + - residual decreased
*                 - - residual did not reduce

      include 'vacdef.f'

      double precision  rhs(ixGlo1:ixGhi1,ixGlo2:ixGhi2),phi(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),tolerance
      integer  matvecmax,info
      character*10   purpose
      character*3  typestop
      logical  nonzero

      double precision  resid
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim,iter
      logical  okprint

      external matvec_poisson
*-----------------------------------------------------------------------------

      oktest=index(teststr,'poisson').ge.1
      okprint=index(teststr,'printiter').ge.1
      if(oktest)write(*,*)'Poisson solver called for ',purpose

      if(typepoisson.eq.'default')then
*        Select method for iterative solver depending on the grid type
         if(typeaxial.eq.'slab'.and. .not.gencoord)then
            typepoisson='cg'
         else
            typepoisson='bicgstab'
         endif
      endif
      if(it.eq.itmin.and.verbose) write(*,*) 'Using ',typepoisson,
     &   ' for solving the Poisson equation for ',purpose

* Initialize parameters and phi for the iterative solvers
      iter=matvecmax
      resid=tolerance

* Solve the Poisson equation
      if(typepoisson.eq.'cg')then
         call cgscalar(okprint,rhs,ixMmin1,ixMmin2,ixMmax1,ixMmax2,nonzero,
     &      phi,matvec_poisson, iter,resid,typestop,info)
      else if(typepoisson.eq.'bicgstab')then
         call bicgstabscalar(okprint,rhs,ixMmin1,ixMmin2,ixMmax1,ixMmax2,
     &      nonzero,phi,matvec_poisson, iter,resid,typestop,info)
      else if(typepoisson.eq.'minres')then
         call minresscalar(okprint,rhs,ixMmin1,ixMmin2,ixMmax1,ixMmax2,
     &      nonzero,phi,matvec_poisson, iter,resid,typestop,info)
      else
         call die('Error in Poisson: Unknown type of iterative method:'//
     &      typepoisson)
      end if

      if(oktest)write(*,*)'Poisson info,nmatvec,resid',info,iter,resid

      if(info.ne.0.and.info.ne.3)then
         nerror(poissonerr_)=nerror(poissonerr_)+1
         if(nerror(poissonerr_).eq.1)then
            write(*,*)'No convergence for Poisson eq. for ',purpose
            write(*,"(a,i2,a,i5,a,i2,a,i5,a,g10.3)")'Error code=',poissonerr_,
     &         ' it=',it,' info=',info,' iter=',iter,' resid=',resid
         if(abs(info).eq.1)then
               write(*,*)'Breakdown due to division by a small value'
         else if(abs(info).eq.2)then
               write(*,*)'No convergence within maximum number of iterations'
            end if
            if(info.gt.0)write(*,*)'The residual decreased'
            if(info.lt.0)write(*,*)'The residual did not decrease'
         end if
      end if

      return
      end

*=============================================================================
      subroutine matvec_poisson(qx,qy)

* Calculate qy=laplace(qx) for the Poisson solvers

      include 'vacdef.f'

      double precision  qx(ixGlo1:ixGhi1,ixGlo2:ixGhi2),qy(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)

      double precision  qdx(ixGlo1:ixGhi1,ixGlo2:ixGhi2),qddx(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)

      integer  ix,ix1,ix2,ixmin1,ixmin2,ixmax1,ixmax2,idim
*-----------------------------------------------------------------------------

      oktest=index(teststr,'matvec').ge.1
      if(oktest)write(*,*)'Matvec_Poisson'

      if(oktest)write(*,*)'before bound qx:',qx(ixtest1,ixtest2)

      call boundscalar(qx)

      if(oktest)write(*,*)'after  bound qx:',qx(ixtest1,ixtest2)

      if(index(teststr,'uselap4').gt.0)then
         call laplace4(qx,ixMmin1,ixMmin2,ixMmax1,ixMmax2,qy)
         return
      endif

      if(fourthorder)then
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            qy(ix_1,ix_2)=zero
         enddo
         enddo
         do idim=1,ndim
            call gradient4(.true. ,qx ,ixMmin1,ixMmin2,ixMmax1,ixMmax2 ,idim,
     &         qdx )
            call boundgradient(idim,qdx)
            call gradient4(.false.,qdx,ixMmin1,ixMmin2,ixMmax1,ixMmax2,idim,
     &         qddx)
            do ix_2=ixMmin2,ixMmax2
            do ix_1=ixMmin1,ixMmax1
               qy(ix_1,ix_2)=qy(ix_1,ix_2)+qddx(ix_1,ix_2)
            enddo
            enddo
         enddo
         return
      endif

* Calculate y=Laplace(phi)
      do ix_2=ixMmin2,ixMmax2
      do ix_1=ixMmin1,ixMmax1
         qy(ix_1,ix_2)=zero
      enddo
      enddo

      ixmin1=ixMmin1-1
      ixmin2=ixMmin2-1
      ixmax1=ixMmax1+1
      ixmax2=ixMmax2+1
      do idim=1,ndim
*        qdx=d qx/dx_idim (gradient)
         call gradient(.true. ,qx ,ixmin1,ixmin2,ixmax1,ixmax2 ,idim,qdx )
*        qddx=d qdx/dx_idim (divergence of gradient)
         call gradient(.false.,qdx,ixMmin1,ixMmin2,ixMmax1,ixMmax2,idim,qddx)
*        qy=Laplace(qx)=Sum_idim(qddx)
         do ix_2=ixMmin2,ixMmax2
         do ix_1=ixMmin1,ixMmax1
            qy(ix_1,ix_2)=qy(ix_1,ix_2)+qddx(ix_1,ix_2)
         enddo
         enddo

         if(oktest)write(*,*)'idim,qdx,qddx,qy:',idim,qdx(ixtest1,ixtest2),
     &      qddx(ixtest1,ixtest2),qy(ixtest1,ixtest2)
      enddo

      return
      end

*=============================================================================
      subroutine boundscalar(phi)

* Calculate boundary for phi based on typeBscalar

      include 'vacdef.f'

      double precision  phi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ix,ixe,ixf,ixmin1,ixmin2,ixmax1,ixmax2,ixpairmin1,ixpairmin2,
     &   ixpairmax1,ixpairmax2,idim,iB
*-----------------------------------------------------------------------------

      oktest=index(teststr,'boundscalar').ge.1
      if(oktest)write(*,*)'BoundScalar phi:',phi(ixtest1,ixtest2)



      do iB=1,nB
         idim=idimB(iB)
         ixmin1=ixBmin(1,iB)
         ixmin2=ixBmin(2,iB)
         ixmax1=ixBmax(1,iB)
         ixmax2=ixBmax(2,iB)

         if(oktest)write(*,*)'iB,idim,ixL,typeB',iB,idim,ixmin1,ixmin2,ixmax1,
     &      ixmax2,typeBscalar(iB)

      if(typeBscalar(iB).eq.'periodic')then
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
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               phi(ix_1,ix_2)=phi(ix_1+(ixpairmin1-ixmin1),ix_2+
     &            (ixpairmin2-ixmin2))
            enddo
            enddo
      else if(typeBscalar(iB).eq.'cont')then
*           ghost cells = edge
         if(idim.eq.1)then
               if(upperB(iB))then
                 ixe=ixmin1-1
               else
                 ixe=ixmax1+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  do ix_1=ixmin2,ixmax2
                     phi(ix,ix_1)=phi(ixe,ix_1)
                  enddo
               end do 
            
         else if(idim.eq.2)then
               if(upperB(iB))then
                 ixe=ixmin2-1
               else
                 ixe=ixmax2+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ix)=phi(ix_1,ixe)
                  enddo
               end do 
            
            end if
      else if(typeBscalar(iB).eq.'cont1')then
*           ghost cells are extrapolated from edge
         if(idim.eq.1)then
               if(upperB(iB))then
                 ixe=ixmin1-1
                  ixf=ixe-1
               else
                 ixe=ixmax1+1
                  ixf=ixe+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  do ix_1=ixmin2,ixmax2
                     phi(ix,ix_1)=(abs(ix-ixe)+1)*phi(ixe,ix_1)-abs(ix-
     &                  ixe)   *phi(ixf,ix_1)
                  enddo
               end do 
            
         else if(idim.eq.2)then
               if(upperB(iB))then
                 ixe=ixmin2-1
                  ixf=ixe-1
               else
                 ixe=ixmax2+1
                  ixf=ixe+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ix)=(abs(ix-ixe)+1)*phi(ix_1,ixe)-abs(ix-
     &                  ixe)   *phi(ix_1,ixf)
                  enddo
               end do 
            
            end if
      else if(typeBscalar(iB).eq.'symm')then
         if(idim.eq.1)then
               if(upperB(iB))then
                  ixe=2*ixmin1-1
               else
                  ixe=2*ixmax1+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  do ix_1=ixmin2,ixmax2
                     phi(ix,ix_1)=+phi(ixe-ix,ix_1)
                  enddo
               end do
            
         else if(idim.eq.2)then
               if(upperB(iB))then
                  ixe=2*ixmin2-1
               else
                  ixe=2*ixmax2+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ix)=+phi(ix_1,ixe-ix)
                  enddo
               end do
            
            end if
      else if(typeBscalar(iB).eq.'asymm')then
         if(idim.eq.1)then
               if(upperB(iB))then
                  ixe=2*ixmin1-1
               else
                  ixe=2*ixmax1+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  do ix_1=ixmin2,ixmax2
                     phi(ix,ix_1)=-phi(ixe-ix,ix_1)
                  enddo
               end do
            
         else if(idim.eq.2)then
               if(upperB(iB))then
                  ixe=2*ixmin2-1
               else
                  ixe=2*ixmax2+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ix)=-phi(ix_1,ixe-ix)
                  enddo
               end do
            
            end if
      else if(typeBscalar(iB).eq.'grad0')then
*           phi should have 0 gradient in all directions at the first ghost cells,
*           therefore phi is 0 at 1st, and copy of mesh edge at 2nd ghost cells
*           In generalized coordinates and/or in axial symmetry, multiply by 
*           the ratio of surfaces on the two sides of the 1st ghost cell.
         if(idim.eq.1)then
               if(upperB(iB))then
                  do ix_1=ixmin2,ixmax2
                     phi(ixmin1,ix_1)=zero
                  enddo
                  if(gencoord)then
                     do ix_1=ixmin2,ixmax2
                        phi(ixmin1+1,ix_1)=phi(ixmin1-1,ix_1)*
     &                     surfaceC(ixmin1-1,ix_1,idim)/surfaceC(ixmin1,ix_1,
     &                     idim)
                     enddo
                  else if(typeaxial.ne.'slab'.and.idim.eq.r_)then
                     do ix_1=ixmin2,ixmax2
                        phi(ixmin1+1,ix_1)=phi(ixmin1-1,ix_1)*
     &                     areaC(ixmin1-1)/areaC(ixmin1)
                     enddo
                  else
                     do ix_1=ixmin2,ixmax2
                        phi(ixmin1+1,ix_1)=phi(ixmin1-1,ix_1)
                     enddo
                  endif
               else
                  do ix_1=ixmin2,ixmax2
                     phi(ixmax1,ix_1)=zero
                  enddo
                  if(gencoord)then
                     do ix_1=ixmin2,ixmax2
                        phi(ixmax1-1,ix_1)=phi(ixmax1+1,ix_1)*
     &                     surfaceC(ixmax1,ix_1,idim)/surfaceC(ixmax1-
     &                     1,ix_1,idim)
                     enddo
                  else if(typeaxial.ne.'slab'.and.idim.eq.1)then
                     do ix_1=ixmin2,ixmax2
                        phi(ixmax1-1,ix_1)=phi(ixmax1+1,ix_1)*
     &                     areaC(ixmax1)/areaC(ixmax1-1)
                     enddo
                  else
                     do ix_1=ixmin2,ixmax2
                        phi(ixmax1-1,ix_1)=phi(ixmax1+1,ix_1)
                     enddo
                  endif
               endif
            
         else if(idim.eq.2)then
               if(upperB(iB))then
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ixmin2)=zero
                  enddo
                  if(gencoord)then
                     do ix_1=ixmin1,ixmax1
                        phi(ix_1,ixmin2+1)=phi(ix_1,ixmin2-1)*
     &                     surfaceC(ix_1,ixmin2-1,idim)/surfaceC(ix_1,ixmin2,
     &                     idim)
                     enddo
                  else if(typeaxial.ne.'slab'.and.idim.eq.r_)then
                     do ix_1=ixmin1,ixmax1
                        phi(ix_1,ixmin2+1)=phi(ix_1,ixmin2-1)*
     &                     areaC(ixmin1-1)/areaC(ixmin1)
                     enddo
                  else
                     do ix_1=ixmin1,ixmax1
                        phi(ix_1,ixmin2+1)=phi(ix_1,ixmin2-1)
                     enddo
                  endif
               else
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ixmax2)=zero
                  enddo
                  if(gencoord)then
                     do ix_1=ixmin1,ixmax1
                        phi(ix_1,ixmax2-1)=phi(ix_1,ixmax2+1)*
     &                     surfaceC(ix_1,ixmax2,idim)/surfaceC(ix_1,ixmax2-
     &                     1,idim)
                     enddo
                  else if(typeaxial.ne.'slab'.and.idim.eq.1)then
                     do ix_1=ixmin1,ixmax1
                        phi(ix_1,ixmax2-1)=phi(ix_1,ixmax2+1)*
     &                     areaC(ixmax1)/areaC(ixmax1-1)
                     enddo
                  else
                     do ix_1=ixmin1,ixmax1
                        phi(ix_1,ixmax2-1)=phi(ix_1,ixmax2+1)
                     enddo
                  endif
               endif
            
            end if
      else if(typeBscalar(iB).eq.'nul')then
             do ix_2=ixmin2,ixmax2
             do ix_1=ixmin1,ixmax1
                phi(ix_1,ix_2)=zero
             enddo
             enddo
         
      else
             write(*,*)'Error in BoundScalar, unknown boundary type:',
     &          typeBscalar(iB),' iB=',iB
             call die('Correct parameter file')
         end if
      end do

      if(oktest)write(*,*)'final phi:',phi(ixtest1,ixtest2)

      return
      end

*=============================================================================
      subroutine boundgradient(idim,phi)

* Calculate boundary for gradient of phi based on typeBscalar 
* and the direction idir in which the gradient was taken
* This is only needed for fourth order scheme. Note that the symm and antisymm
* conditions are reversed for the gradient of phi. 

      include 'vacdef.f'

      double precision  phi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      integer  ix,ixe,ixmin1,ixmin2,ixmax1,ixmax2,ixpairmin1,ixpairmin2,
     &   ixpairmax1,ixpairmax2,idim,iB
*-----------------------------------------------------------------------------

      oktest=index(teststr,'boundgrad').ge.1
      if(oktest)write(*,*)'BoundGradient grad phi:',phi(ixtest1,ixtest2)



      do iB=1,nB
         if(idim.eq.idimB(iB))then
            ixmin1=ixBmin(1,iB)
            ixmin2=ixBmin(2,iB)
            ixmax1=ixBmax(1,iB)
            ixmax2=ixBmax(2,iB)

            if(oktest)write(*,*)'iB,idim,ixL,typeB',iB,idim,ixmin1,ixmin2,
     &         ixmax1,ixmax2,typeBscalar(iB)

         if(typeBscalar(iB).eq.'periodic')then
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
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  phi(ix_1,ix_2)=phi(ix_1+(ixpairmin1-ixmin1),ix_2+
     &               (ixpairmin2-ixmin2))
               enddo
               enddo
         else if(typeBscalar(iB).eq.'symm')then
            if(idim.eq.1)then
               if(upperB(iB))then
                  ixe=2*ixmin1-1
               else
                  ixe=2*ixmax1+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  do ix_1=ixmin2,ixmax2
                     phi(ix,ix_1)=-phi(ixe-ix,ix_1)
                  enddo
               end do
               
            else if(idim.eq.2)then
               if(upperB(iB))then
                  ixe=2*ixmin2-1
               else
                  ixe=2*ixmax2+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ix)=-phi(ix_1,ixe-ix)
                  enddo
               end do
               
               end if
         else if(typeBscalar(iB).eq.'asymm')then
            if(idim.eq.1)then
               if(upperB(iB))then
                  ixe=2*ixmin1-1
               else
                  ixe=2*ixmax1+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  do ix_1=ixmin2,ixmax2
                     phi(ix,ix_1)=+phi(ixe-ix,ix_1)
                  enddo
               end do
               
            else if(idim.eq.2)then
               if(upperB(iB))then
                  ixe=2*ixmin2-1
               else
                  ixe=2*ixmax2+1
               endif
*              HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     phi(ix_1,ix)=+phi(ix_1,ixe-ix)
                  enddo
               end do
               
               end if
         else if(typeBscalar(iB).eq.'cont'.or.typeBscalar(iB).eq.'nul')then
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  phi(ix_1,ix_2)=zero
               enddo
               enddo
            
         else
               write(*,*)'Error in BoundGradient, unknown boundary type:',
     &            typeBscalar(iB),' iB=',iB
               call die('Correct parameter file')
            end if
         endif
      end do

      if(oktest)write(*,*)'final grad phi:',phi(ixtest1,ixtest2)

      return
      end

*=============================================================================
      subroutine cgscalar(okprint,rhs,ixmin1,ixmin2,ixmax1,ixmax2,nonzero,qx,
     &   matvec,iter,tol,typestop,info)

* The CG-algorithm is implemented as shown on page 12 of the thesis
* "Preconditioning for sparse matrices with applications."
* Auke van der Ploeg, University of Groningen, 1994.
* Rewritten to F90 by G. Toth based on the F77 subroutine in src/conjgrad.f 

* This subroutine determines the solution of A.QX=RHS, where
* the matrix-vector multiplication with A is performed by 
* the subroutine 'matvec'.

*!! If the matrix is not symmetric positive definite, CG is likely to fail.

*     Description of arguments:

*     okprint: (input) (boolean) 
*        Determines whether of not output is printed.
*     matvec: external routine for matrix-vector multiplication.
*     rhs: (input/output)
*        on input:  right-hand side vector.
*        on output: residual vector.
*     ixL:
*        Region of unknowns within the computational grid ixG
*     nonzero: (input) (boolean)
*        Tells CG if initial guess in qx is zero or not. 
*        If nonzero is .FALSE., one MATVEC call is saved.
*     qx: (input/output)
*        on input:  initial guess for the solution vector.
*        on output: solution vector.
*     matvec: (subroutine)
*        performes the action of the A matrix
*     iter: (input/output)
*       on input:  maximum number of iterations to be performed.
*       on output: actual  number of iterations done.
*     tol: (input/output) 
*       on input:  required (relative) 2-norm or maximum norm of residual
*       on output: achieved (relative) 2-norm or maximum norm of residual
*     typestop (input) (character*3)
*       Determine stopping criterion (||.|| denotes the 2-norm):
*       typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
*       typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
*       typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol
*     info (output)
*         Gives reason for returning:
*     abs(info)=  0 - solution found satisfying given tolerance.
*                 1 - iteration aborted due to division by very small value.
*                 2 - no convergence within maximum number of iterations.
*                 3 - initial guess satisfies the stopping criterion.
*    sign(info)=  + - residual decreased
*                 - - residual did not reduce

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,iter,info
      double precision  qx(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rhs(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),tol
      logical  okprint,nonzero
      character*3  typestop

      integer  i,itr,matv
      double precision  rho,rhonew,res,res0,bet,alf,assumedzero

      external matvec
*----------------------------------------------------------------------------
      if(okprint)write(*,*)'CGscalar tol,mxmv,ixL:',tol,iter,ixmin1,ixmin2,
     &   ixmax1,ixmax2

      assumedzero=1.D-16
       itr=0
       matv=0

      if (typestop.ne.'rel'.and.typestop.ne.'abs'.and.typestop.ne.'max') then
         write(*,*)'Error in CG:'
         call die('Parameter typestop='//typestop//
     &      ' should be one of rel/abs/max')
      end if

      if(okprint) write(*,*)'n gives the number of CG-iterations.'

* Calculate the initial residual R:=RHS-A*X and its 2-norm.

      if (nonzero) then
         call matvec(qx,tmp)
         matv = matv + 1
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            rhs(ix_1,ix_2)=rhs(ix_1,ix_2)-tmp(ix_1,ix_2)
         enddo
         enddo
      endif

* rho=||rhs||
      sum_1=0.
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         sum_1=sum_1+rhs(ix_1,ix_2)**2
      enddo
      enddo
      rho=sum_1

      rho=sqrt(rho)

      res0=rho
      if (typestop.eq.'max')then
         maxval_1=-bigdouble
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            maxval_1=max(maxval_1,abs(rhs(ix_1,ix_2)))
         enddo
         enddo
         res0=maxval_1
         
      endif
      res=res0

      assumedzero = assumedzero*res0

      if (okprint) then
         IF (typestop.eq.'max') THEN
            write(*,*)'n:',itr,' Maximum norm initial residual:',res0
         ELSE
            write(*,*)'n:',itr,' 2-norm intial residual:',res0
         END IF
      end if

      if (res0.lt.smalldouble.or.(typestop.ne.'rel'.and.res0.le.tol)) then
         info = 3
      else
*        Initialize rho and tmp=Z
         rho=rho*rho
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=rhs(ix_1,ix_2)
         enddo
         enddo

*        Do iteration
7701     continue
*           AZ=A.Z
            call matvec(tmp,tmp2) 
            matv=matv+1

*           alf=A.AZ
            sum_1=0.
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               sum_1=sum_1+tmp(ix_1,ix_2)*tmp2(ix_1,ix_2)
            enddo
            enddo
            alf=sum_1
            
            if(okprint)write(*,*)'alf=',alf

            if (abs(alf).le.assumedzero**2) then
               info = 1
               goto 7702  
            end if
            alf=rho/alf
            if(okprint)write(*,*)'alf=',alf

            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               qx(ix_1,ix_2) =qx(ix_1,ix_2)  + alf*tmp(ix_1,ix_2)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               rhs(ix_1,ix_2)=rhs(ix_1,ix_2) - alf*tmp2(ix_1,ix_2)
            enddo
            enddo

*           rhonew=||rhs||
            sum_1=0.
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               sum_1=sum_1+rhs(ix_1,ix_2)**2
            enddo
            enddo
            rhonew=sum_1
            
            rhonew=sqrt(rhonew)
            if(okprint)write(*,*)'rhonew=',rhonew

         if(typestop.eq.'max')then
               maxval_1=-bigdouble
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  maxval_1=max(maxval_1,abs(rhs(ix_1,ix_2)))
               enddo
               enddo
               res=maxval_1
               
         else if(typestop.eq.'rel')then
               res=rhonew/res0
         else if(typestop.eq.'abs')then
               res=rhonew
            end if
            rhonew=rhonew*rhonew
            itr=itr+1
            IF (okprint) write(*,*)'n:',itr,' ',typestop,
     &         '. norm of residual:',res

            if (res.le.tol) then
               info = 0
               goto 7702  
            end if
            if (itr.ge.iter) then
               info = 2
               goto 7702  
            end if
            if (rho.le.assumedzero**2) then
               info = 1
               goto 7702  
            end if 

            bet=rhonew/rho
            if(okprint)write(*,*)'bet=',bet

            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               tmp(ix_1,ix_2)=rhs(ix_1,ix_2)+bet*tmp(ix_1,ix_2)
            enddo
            enddo

            IF (okprint) write(*,*)'alf,bet,rho,rhonew:',alf,bet,rho,rhonew
            rho=rhonew
         goto 7701  
7702     continue
      endif

* return number of iterations and achieved residual
      iter=itr
      tol =res
      if((typestop.eq.'rel'.and.res.gt.one).or.(typestop.ne.'rel'.and.
     &   res.gt.res0))info=-info

* report results
      if(okprint)then
         write(*,*)'Total Number of CG-iterations:',itr
         write(*,*)'Number of matrix-vector mult.:',matv
      if(abs(info).eq.0)then
            write(*,*)'Successful iteration, norm of res. is:',tol
      else if(abs(info).eq.1)then
            write(*,*)'Iteration aborted due to division by a'
            write(*,*)'very small value.' 
      else if(abs(info).eq.2)then
            write(*,*)'Stopping crit. not fulfilled within given'
            write(*,*)'maximum number of iterations.'
      else if(abs(info).eq.3)then
            write(*,*)'Initial guess for the solution satisfies'
            write(*,*)'given stopping criterion.' 
      else
            write(*,*)'Impossible value for info:',info
         end if
         if(info.lt.0)write(*,*)'The residual did not reduce'
      endif

      return
* -------------------------- end of cgscalar -----------------------------
      end

*=============================================================================
      subroutine bicgstabscalar(okprint,rhs,ixmin1,ixmin2,ixmax1,ixmax2,
     &   nonzero,qx,matvec,mxmv,tol,typestop,info)

* Simple BiCGstab(\ell=1) iterative method
* Modified by G.Toth from the \ell<=2 version written
* by M.A.Botchev, Jan.'98
*
* This is the "vanilla" version of BiCGstab(\ell) as described
* in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements 
* to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
* 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
*    properties of BiCGstab methods in finite precision arithmetic",
*    Numerical Algorithms, 10, 1995, pp.203-223
* 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
*    hybrid BiCG methods", Computing, 56, 1996, 141-163
*
* {{ This code is based on:
* subroutine bistbl v1.0 1995
*
* Copyright (c) 1995 by D.R. Fokkema.
* Permission to copy all or part of this work is granted,
* provided that the copies are not made or distributed
* for resale, and that the copyright notice and this
* notice are retained.  }}
*
* ix^L    == (input) INTEGER sizes of the system to solve 
* qx      == (input/output) DOUBLE PRECISION array dimension n
*            initial guess on input, solution on output
* rhs     == (input) DOUBLE PRECISION array dimension n
*            right-hand side (rhs) vector
*            (output) changed to initial residual=rhs-A.qx
* matvec  == (input) EXTERNAL name of matrix vector subroutine
*            to deliver y:=A*x by CALL matvec(n,x,y)
* nonzero == (input) LOGICAL tells
*            BiCGstab(\ell) if initial guess in x is zero or not. 
*            If nonzero is .FALSE., one MATVEC call is saved.
* tol     == (input/output) DOUBLE PRECISION tolerance for all possible
*            stopping criterions (see the 'typestop' parameter)
*            On output, if info=0 or 1, tol is actually achieved
*            residual reduction or whatever (see the 'typestop' parameter)
* typestop== (input) CHARACTER*3 stopping criterion (||.|| denotes 
*            the 2-norm):
*            typestop='rel' -- relative stopping crit.: ||res|| <= tol*||res0||
*            typestop='abs' -- absolute stopping crit.: ||res|| <= tol
*            typestop='max' -- maximum  stopping crit.: max(abs(res)) <= tol
* NOTE(for typestop='rel' and 'abs'): To save computational work, the value of 
*            residual norm used to check the convergence inside the main 
*            iterative loop is computed from 
*            projections, i.e. it can be smaller than the true residual norm
*            (it may happen when e.g. the 'matrix-free' approach is used).
*            Thus, it is possible that the true residual does NOT satisfy
*            the stopping criterion ('rel' or 'abs').
*            The true residual norm (or residual reduction) is reported on 
*            output in parameter TOL -- this can be changed to save 1 MATVEC
*            (see comments at the end of the subroutine)
* mxmv   ==  (input/output) INTEGER.  On input: maximum number of matrix 
*            vector multiplications allowed to be done.  On output: 
*            actual number of matrix vector multiplications done
* info   ==  (output) INTEGER. 
*     abs(info)=  0 - solution found satisfying given tolerance.
*                 1 - iteration aborted due to division by very small value.
*                 2 - no convergence within maximum number of iterations.
*                 3 - initial guess satisfies the stopping criterion.
*    sign(info)=  + - residual decreased
*                 - - residual did not reduce
*
* Memory requirement:
* nwrk (see vacdef.t) >= 2*lmax+5
* wrk    ==  (workspace) DOUBLE PRECISION array dimension (ixG^T,nwrk)
*            

* ----------------------------------------------------------
      include 'vacdef.f'
*
*     .. Work Aliases ..
*
      integer   lmax,l,r,u,xp
      PARAMETER(   lmax=1,l=1,r=1,u=r+(l+1),xp=u+(l+1))
      integer   z,zz,y0,yl,y
      PARAMETER(   z=1,zz=z+(l+1),y0=zz+(l+1),yl=y0+1,y=yl+1)
*
*     .. Parameters ..
* 
      double precision  qx(ixGlo1:ixGhi1,ixGlo2:ixGhi2),rhs(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),tol
      integer   ixmin1,ixmin2,ixmax1,ixmax2, mxmv, info
      logical   okprint,nonzero
      character*3   typestop
*
*     .. Matrix ..
*
      external matvec
*
*     .. Local ..
*
      double precision   rwork(lmax+1,2*lmax+5)

      logical GoOn, rcmp, xpdt
      integer i, j, k, nmv
      double precision alpha, beta, omega, rho0, rho1, sigma
      double precision varrho, hatgamma
      double precision assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
      double precision mxnrmx, mxnrmr, kappa0, kappal
      double precision qtmp

*----------------------------------------------------------------------------
      if(okprint)write(*,*)'BiCGSTABscalar tol,mxmv:',tol,mxmv

      if(nwrk.lt.2*lmax+3) then
        write(*,*) 'Error in BiCGSTABscalar: in vacdef.t adjust nwrk>=',2*
     &     lmax+3
        call die('Edit src/vacdef.t or the parameter file')
      endif

      info = 0

      if (tol.le.zero) call die('Error in BiCGSTABscalar: tolerance < 0')
      if (mxmv.le.1)   call die('Error in BiCGSTABscalar: maxmatvec < 1')

*
*     --- Initialize first residual
*
      assumedzero = 1.d-16
      if (nonzero) then
         call matvecxV ( qx, wrk,r, matvec )
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,r) = rhs(ix_1,ix_2) - wrk(ix_1,ix_2,r)
         enddo
         enddo
         nmv = 1
      else
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,r) = rhs(ix_1,ix_2)
         enddo
         enddo
         nmv = 0
      endif
*
*     --- Initialize iteration loop
*

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         rhs(ix_1,ix_2) = wrk(ix_1,ix_2,r)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wrk(ix_1,ix_2,xp) = qx(ix_1,ix_2)
      enddo
      enddo

      if (nonzero) then
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            qx(ix_1,ix_2) = zero
         enddo
         enddo
      endif

      sum_1=0.
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         sum_1=sum_1+ wrk(ix_1,ix_2,r)**2 
      enddo
      enddo
      rnrm0 = sum_1

      rnrm0 = sqrt(rnrm0)
      rnrm = rnrm0
      if(okprint) print *,'initial rnrm:',rnrm

      mxnrmx = rnrm0
      mxnrmr = rnrm0
      rcmp = .false.
      xpdt = .false.

      alpha = zero
      omega = one
      sigma = one
      rho0 =  one
*
*     --- Iterate
*
      if(typestop.eq.'rel')then
         GoOn = rnrm.gt.tol*rnrm0 .and. nmv.lt.mxmv
         assumedzero = assumedzero*rnrm0
         rnrmMax = 0
         rnrmMax0 = 0
      else if(typestop.eq.'abs')then
         GoOn = rnrm.gt.tol       .and. nmv.lt.mxmv
         assumedzero = assumedzero*rnrm0
         rnrmMax = 0
         rnrmMax0 = 0
      else if(typestop.eq.'max')then
         maxval_1=-bigdouble
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            maxval_1=max(maxval_1, abs( wrk(ix_1,ix_2,r) ) )
         enddo
         enddo
         rnrmMax0 = maxval_1
         
         rnrmMax  = rnrmMax0
         if(okprint) print *,'initial rnrmMax:',rnrmMax
         GoOn = rnrmMax.gt.tol    .and. nmv.lt.mxmv
         assumedzero = assumedzero*rnrmMax
      else
         call die('Error in BiCGSTABScalar: unknown typestop value:'//
     &      typestop)
      end if

      if (.not.GoOn) then
         if(okprint) print *,'BiCGSTABScalar: nothing to do. info = ',info
         mxmv = nmv
         info = 3
         return
      endif

      do while (GoOn)
*
*     =====================
*     --- The BiCG part ---
*     =====================
*
         rho0 = -omega*rho0
         do k=1,l
            sum_1=0.
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               sum_1=sum_1+ rhs(ix_1,ix_2)*wrk(ix_1,ix_2,r+k-1) 
            enddo
            enddo
            rho1 = sum_1
            

            if (abs(rho0).lt.assumedzero**2) then
               info = 1
               return
            endif
            beta = alpha*(rho1/rho0)
            rho0 = rho1
            do j=0,k-1
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wrk(ix_1,ix_2,u+j) = wrk(ix_1,ix_2,r+j) - beta*
     &               wrk(ix_1,ix_2,u+j)
               enddo
               enddo
            enddo
            
            call matvecVV (wrk,u+k-1,u+k, matvec)
            nmv = nmv+1
            
            sum_1=0.
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               sum_1=sum_1+ rhs(ix_1,ix_2)*wrk(ix_1,ix_2,u+k) 
            enddo
            enddo
            sigma = sum_1
            

            if (abs(sigma).lt.assumedzero**2) then
               info = 1
               return
            endif

            alpha = rho1/sigma
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               qx(ix_1,ix_2) = alpha*wrk(ix_1,ix_2,u) + qx(ix_1,ix_2)
            enddo
            enddo
            
            do j=0,k-1
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wrk(ix_1,ix_2,r+j) = -alpha*wrk(ix_1,ix_2,u+j+1)    +
     &                wrk(ix_1,ix_2,r+j)
               enddo
               enddo
            enddo

            call matvecVV(wrk,r+k-1,r+k, matvec)
            nmv = nmv+1
            
            sum_1=0.
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               sum_1=sum_1+ wrk(ix_1,ix_2,r)**2 
            enddo
            enddo
            rnrm = sum_1
            
            rnrm = sqrt( rnrm )
            
            mxnrmx = max (mxnrmx, rnrm)
            mxnrmr = max (mxnrmr, rnrm)
         enddo

*        
*        ==================================
*        --- The convex polynomial part ---
*        ================================== 
*        
*        --- Z = R'R
*        
         do i=1,l+1 
            do j=i-1,l
               sum_1=0.
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  sum_1=sum_1+ wrk(ix_1,ix_2,r+j)*wrk(ix_1,ix_2,r+i-1) 
               enddo
               enddo
               qtmp = sum_1
               
               rwork(j+1,z+i-1) = qtmp
               rwork(z+i-1,j+1) = rwork(j+1,z+i-1) 
            enddo
         enddo

*        
*        --- tilde r0 and tilde rl (small vectors)
*        
         do i_2=zz,zz+l
         do i_1=1,l+1
            rwork(i_1,i_2)   = rwork(i_1,i_2+(z-zz)) 
         enddo
         enddo
         rwork(1,y0) = -one
         rwork(2,y0) = rwork(2,z)
         rwork(2,y0) = rwork(2,y0) / rwork(2,zz+1)
         rwork(l+1,y0) = zero

         rwork(1,yl) = zero
         rwork(2,yl) = rwork(2,z+l)
         rwork(2,yl) = rwork(2,yl) / rwork(2,zz+1)
         rwork(l+1,yl) = -one
*        
*        --- Convex combination
*        
         do i_1=1,l+1
            rwork(i_1,y) = zero
         enddo
         do j=1,l+1
            do i_1=1,l+1
               rwork(i_1,y) = rwork(i_1,y) + rwork(j,yl)*rwork(i_1,z+j-1)
            enddo
         enddo
         sum_1=0.
         do i_1=1,l+1
            sum_1=sum_1+ rwork(i_1,yl)*rwork(i_1,y) 
         enddo
         kappal = sqrt( sum_1 )
         do i_1=1,l+1
            rwork(i_1,y) = zero
         enddo
         do j=1,l+1
            do i_1=1,l+1
               rwork(i_1,y) = rwork(i_1,y) + rwork(j,y0)*rwork(i_1,z+j-1)
            enddo
         enddo
         sum_1=0.
         do i_1=1,l+1
            sum_1=sum_1+ rwork(i_1,y0)*rwork(i_1,y)  
         enddo
         kappa0 = sqrt( sum_1 )

         sum_1=0.
         do i_1=1,l+1
            sum_1=sum_1+ rwork(i_1,yl)*rwork(i_1,y) 
         enddo
         varrho = sum_1  
         varrho = varrho / (kappa0*kappal)
         
         hatgamma = sign(1d0,varrho)*max(abs(varrho),7d-1) * (kappa0/kappal)

         do i_1=1,l+1
            rwork(i_1,y0) = -hatgamma*rwork(i_1,yl) + rwork(i_1,y0)
         enddo

*        
*        --- Update
*        
         omega = rwork(l+1,y0)

         do j=1,l
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               wrk(ix_1,ix_2,u) = wrk(ix_1,ix_2,u) - rwork(1+j,y0)*
     &            wrk(ix_1,ix_2,u+j)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               qx(ix_1,ix_2)    = qx(ix_1,ix_2)    + rwork(1+j,y0)*
     &            wrk(ix_1,ix_2,r+j-1)
            enddo
            enddo
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               wrk(ix_1,ix_2,r) = wrk(ix_1,ix_2,r) - rwork(1+j,y0)*
     &            wrk(ix_1,ix_2,r+j)
            enddo
            enddo
         enddo

         do i_1=1,l+1
            rwork(i_1,y) = zero
         enddo
         do j=1,l+1
            do i_1=1,l+1
               rwork(i_1,y) = rwork(i_1,y) + rwork(j,y0)*rwork(i_1,z+j-1)
            enddo
         enddo

         sum_1=0.
         do i_1=1,l+1
            sum_1=sum_1+ rwork(i_1,y0)*rwork(i_1,y) 
         enddo
         rnrm = sqrt( sum_1 )

      if(typestop.eq.'rel')then
            GoOn = rnrm.gt.tol*rnrm0 .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/
     &         rnrm0
      else if(typestop.eq.'abs')then
            GoOn = rnrm.gt.tol       .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' matvecs, ||rn|| =',rnrm
      else if(typestop.eq.'max')then
            maxval_1=-bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               maxval_1=max(maxval_1, abs( wrk(ix_1,ix_2,r) ) )
            enddo
            enddo
            rnrmMax = maxval_1
            
            GoOn = rnrmMax.gt.tol    .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' matvecs, max(rn) =',rnrmMax
         end if

      enddo
*
*     =========================
*     --- End of iterations ---
*     =========================
*

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         qx(ix_1,ix_2) = wrk(ix_1,ix_2,xp) + qx(ix_1,ix_2)
      enddo
      enddo

* --------------------- One matvec can be saved by commenting out this:
*
*     --- Check stopping criterion
*
*!call matvecxV (qx, wrk,r, matvec)
*!nmv = nmv + 1
*!wrk(ix^S,r) = rhs(ix^S) - wrk(ix^S,r)   
*!rnrm = sqrt( sum( wrk(ix^S,r)**2 ) )
* --------------------- One matvec can be saved by commenting out this

      if(typestop.eq.'rel')then
         if (rnrm.gt.tol*rnrm0) info = 2
         tol = rnrm/rnrm0
      else if(typestop.eq.'abs')then
         if (rnrm.gt.tol) info = 2
         tol = rnrm
      else if(typestop.eq.'max')then
         if (rnrmMax.gt.tol) info = 2
         tol = rnrmMax
      end if

      if((typestop.ne.'max'.and.rnrm.gt.rnrm0).or.(typestop.eq.'max'.and.
     &   rnrmMax.gt.rnrmMax0)) info=-info

      mxmv = nmv

      return
      end  
*=======================================================================
      subroutine minresscalar(okprint,rhs,ixmin1,ixmin2,ixmax1,ixmax2,nonzero,
     &   qx,matvec,mxmv,tol,  typestop,info)
* By M.A.Botchev, Apr.'98 
* Matlab implementation of MINRES by G.Sleijpen was used as a base 
*
* parameters are standard for VAC
* info = 2 -- no convergence achieved and residual decreased
* info =-2 -- no convergence achieved and residual did not decrease
* info = 3 -- nothing to do

      include 'vacdef.f'

      character*3  typestop
      integer  ixmin1,ixmin2,ixmax1,ixmax2,info,mxmv
      double precision  rhs(ixGlo1:ixGhi1,ixGlo2:ixGhi2),qx(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),tol
      logical  okprint,nonzero
      external matvec
*-----------------------------------------------
      integer  v,vold,vtld,qw,wtld,wtld2,r
      PARAMETER(  v=1,vold=2,vtld=3,qw=4,wtld=5,wtld2=6,r=7)
      integer  nmv
      double precision  beta,btld,alpha,atld,el0,el1,el2,ci,es, rnrm,rnrm0,
     &   rnrmMax,rnrmMax0
      logical  GoOn
*-----------------------------------------------------------------
      if(okprint)write(*,*)'MINRESscalar tol,mxmv,ixL:',tol,mxmv,ixmin1,
     &   ixmin2,ixmax1,ixmax2

      if (typestop.ne.'rel'.and.typestop.ne.'abs'.and.typestop.ne.'max') then
         write(*,*)'Error in MINRES:'
         write(*,*)'Parameter typestop=',typestop,
     &      ' should be one of rel/abs/max.'
         call die('Error in vacpoisson.t')
      end if

      if(nwrk.lt.7) write(*,*) 
     &   'Error in MINRESscalar: in vacdef.t adjust nwrk>=7'

      info=0

* initial residual:
      if (nonzero) then
         call matvecxV ( qx, wrk,r, matvec )
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,r) = rhs(ix_1,ix_2) - wrk(ix_1,ix_2,r)
         enddo
         enddo
         nmv = 1
      else
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,r) = rhs(ix_1,ix_2)
         enddo
         enddo
         nmv = 0
      endif
      sum_1=0.
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         sum_1=sum_1+ wrk(ix_1,ix_2,r)**2 
      enddo
      enddo
      rnrm0 = sum_1

      rnrm0 = sqrt(rnrm0)
      rnrm = rnrm0
      rnrmMax = zero
      rnrmMax0= zero

      if(typestop.eq.'rel')then
         GoOn = rnrm.gt.tol*rnrm0 .and. nmv.lt.mxmv
      else if(typestop.eq.'abs')then
         GoOn = rnrm.gt.tol       .and. nmv.lt.mxmv
      else if(typestop.eq.'max')then
         maxval_1=-bigdouble
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            maxval_1=max(maxval_1, abs( wrk(ix_1,ix_2,r) ) )
         enddo
         enddo
         rnrmMax = maxval_1
         
         rnrmMax0=rnrmMax
         GoOn = rnrmMax.gt.tol    .and. nmv.lt.mxmv
      end if

      if (.not.GoOn) then
         info = 3
         if(okprint) print *,'VACMINRES90: nothing to do'
         return
      endif

      beta = zero
      btld = zero
      ci = - one
      es = zero
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wrk(ix_1,ix_2,v)     = wrk(ix_1,ix_2,r)/rnrm0
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wrk(ix_1,ix_2,wtld2) = wrk(ix_1,ix_2,v)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wrk(ix_1,ix_2,vold)  = zero
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wrk(ix_1,ix_2,qw)     = zero
      enddo
      enddo
           
* main loop:
      do while (GoOn)

         call matvecVV (wrk,v,vtld, matvec)
         nmv = nmv + 1
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,vtld) = wrk(ix_1,ix_2,vtld) - beta*
     &         wrk(ix_1,ix_2,vold)
         enddo
         enddo

         sum_1=0.
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            sum_1=sum_1+ wrk(ix_1,ix_2,v)*wrk(ix_1,ix_2,vtld) 
         enddo
         enddo
         alpha = sum_1
         
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,vtld) = wrk(ix_1,ix_2,vtld) - alpha*wrk(ix_1,ix_2,v)
         enddo
         enddo

         sum_1=0.
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            sum_1=sum_1+ wrk(ix_1,ix_2,vtld)**2 
         enddo
         enddo
         beta = sum_1
         
         beta = sqrt(beta)
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,vold) = wrk(ix_1,ix_2,v)
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,v) = wrk(ix_1,ix_2,vtld)/beta
         enddo
         enddo

         el1 = es*alpha - ci*btld
         el2 = es*beta
         atld = -es*btld-ci*alpha
         btld = ci*beta
         el0 = sqrt(atld**2 + beta**2)
         ci = atld/el0
         es = beta/el0

         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,wtld)  = wrk(ix_1,ix_2,wtld2) - el1*
     &         wrk(ix_1,ix_2,qw)
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,wtld2) = wrk(ix_1,ix_2,v)     - el2*
     &         wrk(ix_1,ix_2,qw)
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            wrk(ix_1,ix_2,qw) = wrk(ix_1,ix_2,wtld)/el0
         enddo
         enddo

         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            qx(ix_1,ix_2) = qx(ix_1,ix_2) + (rnrm*ci)*wrk(ix_1,ix_2,qw)
         enddo
         enddo

         rnrm = rnrm*es

      if(typestop.eq.'rel')then
            GoOn = rnrm.gt.tol*rnrm0 .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/
     &         rnrm0
      else if(typestop.eq.'abs')then
            GoOn = rnrm.gt.tol       .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' matvecs, ||rn|| =',rnrm
      else if(typestop.eq.'max')then
            if(mod(nmv,3).eq.0) then
               call matvecxV ( qx, wrk,r, matvec )
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  wrk(ix_1,ix_2,r) = rhs(ix_1,ix_2) - wrk(ix_1,ix_2,r)
               enddo
               enddo
               maxval_1=-bigdouble
               do ix_2=ixmin2,ixmax2
               do ix_1=ixmin1,ixmax1
                  maxval_1=max(maxval_1, abs( wrk(ix_1,ix_2,r) ) )
               enddo
               enddo
               rnrmMax = maxval_1
               
               GoOn = rnrmMax.gt.tol    .and. nmv.lt.mxmv
               if(okprint) print *, nmv,' matvecs, max(rn) =',rnrmMax
            endif
         end if

      enddo

      if(typestop.eq.'rel')then
         if (rnrm.gt.tol*rnrm0) info = 2
         tol = rnrm/rnrm0
      else if(typestop.eq.'abs')then
         if (rnrm.gt.tol) info = 2
         tol = rnrm
      else if(typestop.eq.'max')then
         if (rnrmMax.gt.tol) info = 2
         tol = rnrmMax
      end if

      if((typestop.ne.'max'.and.rnrm.gt.rnrm0).or.(typestop.eq.'max'.and.
     &   rnrmMax.gt.rnrmMax0)) info=-info

      mxmv = nmv

      return
      end  
*=============================================================================
      subroutine matvecxV (q, workVv,k, matvec)
* Performs v_k := A . q for BICGSTAB
      include 'vacdef.f'
      integer   lmax
      PARAMETER(   lmax=1)
      integer  k
      double precision   q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),workVv(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,2*lmax+5)
      external matvec
*----------------------------------------------------------------------
      call matvec(q,tmp)
      do ix_2=ixMmin2,ixMmax2
      do ix_1=ixMmin1,ixMmax1
         workVv(ix_1,ix_2,k)=tmp(ix_1,ix_2)
      enddo
      enddo
      return
      end

*=======================================================================
      subroutine matvecVV (workVv,k1,k2, matvec)
* Performs v_k2 := A . v_k1
      include 'vacdef.f'
      integer   lmax
      PARAMETER(   lmax=1)
      integer  k1,k2
      double precision   workVv(ixGlo1:ixGhi1,ixGlo2:ixGhi2,2*lmax+5)
      external matvec
*-----------------------------------------------------------------------
      do ix_2=ixMmin2,ixMmax2
      do ix_1=ixMmin1,ixMmax1
         tmp2(ix_1,ix_2)=workVv(ix_1,ix_2,k1)
      enddo
      enddo
      call matvecxV(tmp2, workVv,k2, matvec)
      return
      end
*=======================================================================
* end module vacpoisson
*##############################################################################
