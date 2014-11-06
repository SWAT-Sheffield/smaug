*##############################################################################
* module vac

*=============================================================================
      program vac

* Versatile Advection Code, (c) Gabor Toth. Started on Nov 8, 1994. 

*     declare common variables and parameters
      include 'vacdef.f'

      integer  ifile,ierrcode,iw
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnrm2,dtold,time0,
     &   time1

* functions
      logical  timetofinish,timetosave
      double precision  cputime
*-----------------------------------------------------------------------------


      verbose=.true. 
      if(verbose)then
         write(*,'(a)')'VAC 4.52 configured to'
         write(*,'(a)')'  -d=22 -phi=0 -z=0 -g=260,260 -p=mhd -u=gravity'
         write(*,'(a)')'  -on=cd,tvdlf,tvd,poisson,resist'
         write(*,'(a)')'  -off=mc,fct,impl,ct,gencoord,rk,mpi'
         
      endif

        time0=cputime()

*     Initialize physics dependent variables
      call physini
*     Read filenames and parameters for advection
      call readparameters(w)
*                                   Read initial data, set ixM,ixG,gencoord

             
      if(gencoord)then
         write(*,*) 'Error: input file contains general grid'
         write(*,*) 'Recompile vac after setvac -on=gencoord is set.'
      endif

*     Initialize boundary data
      call boundsetup
* Initialize grid geometry
      if(gencoord)then
         
                call die('Error: gencoord module is off')
      else
         call gridsetup1
      endif

*     Initialize it, t, headers etc.
      call startup

      if(verbose)write(*,'(a,f10.3,a)')'Start Advance  ',cputime()-
     &   time0,' sec'

      call getboundary(t,1,nw,1,ndim,w)
7701  continue
         do ifile=1,nfile
            if(timetosave(ifile)) call savefile(ifile,w)
         end do

*        Determine time step
         if(dtpar.gt.zero)then
            dt=dtpar
         else
            if(courantpar.gt.zero)call getdt_courant(w,ixMmin1,ixMmin2,
     &         ixMmax1,ixMmax2)
            call getdt(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2)
            call getdt_special(w,ixMmin1,ixMmin2,ixMmax1,ixMmax2)

            if(dtcantgrow.and.it.gt.itmin)dt=min(dt,dtold)
            dtold=dt
         endif
         if(dtmrpc.gt.zero)dt=min(dt,dtmrpc)

         if (timetofinish(time0)) goto 7702  

*        For slowsteps == 1, use dtpar in the first time step ONLY
         if(slowsteps.eq.1.and.it.eq.itmin)dtpar=-one

*        For slowsteps > 1, reduce dt for the first few steps 
         if(slowsteps.gt.it-itmin+1) dt=dt*(one-(one-(it-itmin+one)/
     &      slowsteps)**2)

         if(tmaxexact)dt=min(dt,tmax-t+smalldouble)

*        Store w into wold for residual calculations and 
*        for TVD limiting based on the previous time step.
         do iw_3=1,nw
         do ix_2=ixGmin2,ixGmax2
         do ix_1=ixGmin1,ixGmax1
            wold(ix_1,ix_2,iw_3)=w(ix_1,ix_2,iw_3)
         enddo
         enddo
         enddo

*        Advance w (except variables with typefull='nul') by dt in the full grid
         call advance(iw_full,w)

         if(residmin.gt.zero .or. residmax.lt.bigdouble)then
*           calculate true residual ||w_n+1-w_n|| for steady-state calculations
            residual=zero
            do iw=1,nw
               sum_1=0.
               do ix_2=ixMmin2,ixMmax2
               do ix_1=ixMmin1,ixMmax1
                  sum_1=sum_1+w(ix_1,ix_2,iw)**2
               enddo
               enddo
               wnrm2=sum_1
               
               if(wnrm2.lt.smalldouble)wnrm2=one
               sum_1=0.
               do ix_2=ixMmin2,ixMmax2
               do ix_1=ixMmin1,ixMmax1
                  sum_1=sum_1+(w(ix_1,ix_2,iw)-wold(ix_1,ix_2,iw))**2
               enddo
               enddo
               residual = residual + sum_1/wnrm2
            enddo
            
            residual=sqrt(residual/nw)
         endif  

         it=it+1
         t=t+dt

      goto 7701  
7702  continue

      time1=cputime()-time0

      do ifile=1,nfile
         if(itsavelast(ifile).lt.it)call savefile(ifile,w)
         close(unitini+ifile)
      enddo

      if(verbose)write(*,'(a,f10.3,a)')'Finish Advance ',time1,' sec'



      if(dt.lt.dtmin)write(unitterm,*)'Warning: dt<dtmin !!!'
      if(time1.gt.cputimemax)write(unitterm,*)
     &   'Warning: cputimemax exceeded !!!'
      do ierrcode=1,nerrcode
         if(nerror(ierrcode).gt.0)then
            write(*,"(a,i2,a,i5,a)")'Error (code=',ierrcode,') occurred ',
     &         nerror(ierrcode),' times !!!'
         if(ierrcode.eq.toosmallp_)then
               write(*,"(a)")'Error description: Pressure below psmall'
         else if(ierrcode.eq.toosmallr_)then
               write(*,"(a)")'Error description: Density below rhosmall'
         else if(ierrcode.eq.couranterr_)then
               write(*,"(a)")'Error description: Courant number above 1'
         else if(ierrcode.eq.poissonerr_)then
               write(*,"(a)")'Error description: Poisson solver failed'
            end if
         endif
      end do


      if(verbose)then
         if(implpar.gt.zero)then
            write(*,*)'Number of explicit evaluations:',nexpl
            write(*,*)'Number of Newton iterations   :',nnewton
            write(*,*)'Number of linear iterations   :',niter
            write(*,*)'Number of MatVecs             :',nmatvec
         endif

         if(residmin.gt.zero .or. residmax.lt.bigdouble)then
            write(*,*)'Number of time steps          :',it-itmin
            write(*,*)'Residual and residmin         :',residual,residmin
         endif

         write(*,'(a,f10.3,a)')'Finished VAC   ',cputime()-time0,' sec'
      endif



      end

*=============================================================================
      subroutine startup

      include 'vacdef.f'
      integer  ifile,iw,ivector,idim,qnvector
*-----------------------------------------------------------------------------

* Initialize dtmrpc which will be calculated by MRPC
      dtmrpc=-one

* Initialize dtcourant, which will be calculated by TVD, TVD-MUSCL or TVDLF
      do idim=1,ndim
         dtcourant(idim)=bigdouble
      enddo

* If dtpar is set, and not only for the first time step (slowsteps/=1)
* then set courantpar<0, so that dtcourant is not calculated at all
      if(dtpar.gt.zero.and.slowsteps.ne.1)courantpar= -one

      itmin=it
      do ifile=1,nfile
         tsavelast(ifile)=t
         itsavelast(ifile)=it
      end do
      isaveout=0

* Initialize vectoriw based on iw_vector. vectoriw=-1 for scalar variables,
      do iw=1,nw
         vectoriw(iw)=-1
      end do
* It points to the 0-th component (iw_vector=m0_,b0_,...) for vector variables.
* Only the first ndim components of the vector variables are rotated in
* generalized coordinates. 
* qnvector is only used to avoid compiler warning when nvector=0

      qnvector=nvector
      do ivector=1,qnvector
         do idim=1,ndim
            vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
         end do
      end do



* Initial value for residual and counters
      residual=bigdouble
      nexpl=0
      nnewton=0
      nmatvec=0
      niter=0



      return
      end

*=============================================================================
      subroutine advance(iws,w)

* w(iws,t) -> w(iws,t+qdt) based on typedimsplit and typesourcesplit
*
* Add split sources and fluxes with unsplit sources

      include 'vacdef.f'

      integer  iws(niw_)
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw), w1(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

* Add split sources berforehand if this is required
      if(sourcesplit)then
         do iw_3=1,nw
         do ix_2=ixGmin2,ixGmax2
         do ix_1=ixGmin1,ixGmax1
            w1(ix_1,ix_2,iw_3)=w(ix_1,ix_2,iw_3)
         enddo
         enddo
         enddo
         if(typesourcesplit.eq.'sf')then
               call addsource2(dt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixMmin1,
     &            ixMmin2,ixMmax1,ixMmax2,iws,t,w1,t,w)
         else if(typesourcesplit.eq.'sfs')then
               call addsource2(dt/2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixMmin1,
     &            ixMmin2,ixMmax1,ixMmax2,iws,t,w1,t,w)
         else if(typesourcesplit.eq.'ssf')then
               call addsource2(dt/2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixGmin1,
     &            ixGmin2,ixGmax1,ixGmax2,iws,t,w,t,w1)
               call addsource2(dt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixMmin1,
     &            ixMmin2,ixMmax1,ixMmax2,iws,t,w1,t,w)
         else if(typesourcesplit.eq.'ssfss')then
               call addsource2(dt/4,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixGmin1,
     &            ixGmin2,ixGmax1,ixGmax2,iws,t,w,t,w1)
               call addsource2(dt/2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixMmin1,
     &            ixMmin2,ixMmax1,ixMmax2,iws,t,w1,t,w)
         else
               call die('Error: Unknown typesourcesplit='//typesourcesplit)
         end if
         call getboundary(t,1,nw,1,ndim,w)
      endif

* Add fluxes and unsplit sources explicitly or implicitly
      if(typeimpl1.eq.'nul')then
         call advance_expl(typefull1,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iws,w1,w)
      else
          call die('IMPL module is switched off')
      endif

* Add split sources afterwards if this is required
      if(sourcesplit)then
      if(typesourcesplit.eq.'sfs')then
            do iw_3=1,nw
            do ix_2=ixGmin2,ixGmax2
            do ix_1=ixGmin1,ixGmax1
               w1(ix_1,ix_2,iw_3)=w(ix_1,ix_2,iw_3)
            enddo
            enddo
            enddo
            call addsource2(dt/2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixMmin1,
     &         ixMmin2,ixMmax1,ixMmax2,iws,t+dt,w1,t+dt,w)
            call getboundary(t+dt,1,nw,1,ndim,w)
      else if(typesourcesplit.eq.'ssfss')then
            do iw_3=1,nw
            do ix_2=ixGmin2,ixGmax2
            do ix_1=ixGmin1,ixGmax1
               w1(ix_1,ix_2,iw_3)=w(ix_1,ix_2,iw_3)
            enddo
            enddo
            enddo
            call addsource2(dt/4,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixGmin1,
     &         ixGmin2,ixGmax1,ixGmax2,iws,t+dt,w ,t+dt,w1)
            call addsource2(dt/2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixMmin1,
     &         ixMmin2,ixMmax1,ixMmax2,iws,t+dt,w1,t+dt, w)
            call getboundary(t+dt,1,nw,1,ndim,w)
         end if
      endif

      return
      end

*=============================================================================
      subroutine advance_expl(method,ixmin1,ixmin2,ixmax1,ixmax2,iws,w1,w)

* w(t) -> w(t+qdt) within ix^L based on typedimsplit, typesourcesplit, nproc
*
* Add fluxes and unsplit sources, possibly with dimensional splitting
* Boundaries should be kept updated by addsource2 and advect
*
* w1 can be ised freely.

      include 'vacdef.f'

      character*10   method
      integer  ixmin1,ixmin2,ixmax1,ixmax2,iws(niw_)
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),w1(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw)

      logical  firstsweep,lastsweep
      common/first/firstsweep,lastsweep
*-----------------------------------------------------------------------------

      oktest=index(teststr,'advance').ge.1
      if(oktest)write(*,*)'Advance method,it,w:',method,' ',it,w(ixtest1,
     &   ixtest2,iwtest)

      if(ixmin1.ne.ixGmin1.or.ixmin2.ne.ixGmin2.or.ixmax1.ne.ixGmax1.or.
     &   ixmax2.ne.ixGmax2)call die(
     &   'Error in Advance: No subgrids implemented yet...')

      nexpl=nexpl+1
      firstsweep=.true.
      if(dimsplit)then
         if((it/2)*2.eq.it .or. typedimsplit.eq.'xy')then
*           If typedimsplit='xy', always do the sweeps in order of increasing idim,
*           otherwise for even parity of "it" only, and reverse order for odd. 
            do idimsplit=1,ndim
               lastsweep= idimsplit.eq.ndim
               call advect(method,ixmin1,ixmin2,ixmax1,ixmax2,iws,idimsplit,
     &            idimsplit,w1,w)
            enddo
         else
*           If the parity of "it" is odd and typedimsplit=xyyx, do sweeps backwards
            do idimsplit=ndim,1,-1
               lastsweep= idimsplit.eq.1
               call advect(method,ixmin1,ixmin2,ixmax1,ixmax2,iws,idimsplit,
     &            idimsplit,w1,w)
            enddo
         endif
      else
*        Add fluxes from all directions at once
         lastsweep= .true.
         call advect(method,ixmin1,ixmin2,ixmax1,ixmax2,iws,1,ndim,w1,w)
      endif

      if(typefilter1.ne.'nul')then
*        We use updated w for the filter fluxes
         do iw_3=1,nw
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            w1(ix_1,ix_2,iw_3)=w(ix_1,ix_2,iw_3)
         enddo
         enddo
         enddo

*        Filter according to typefilter1
               
      if(typefilter1.eq.'tvd1')then
*           Call tvdlimit with tvd1 (2nd order Lax-Wendroff terms off)
            call tvdlimit(typefilter1,dt,ixmin1,ixmin2,ixmax1,ixmax2,ixmin1+
     &         2,ixmin2+2,ixmax1-2,ixmax2-2,iw_filter,1,ndim,w1,t+dt,w)
        
                 
      else if(typefilter1.eq.'tvdlf'.or.typefilter1.eq.'tvdmu'.or.
     &   typefilter1.eq.'tvdlf1'.or.typefilter1.eq.'tvdmu1')then
*           Call tvdmusclf with filter method and physical fluxes off
            call tvdmusclf(.false.,typefilter1,dt,ixmin1,ixmin2,ixmax1,ixmax2,
     &         ixmin1+2,ixmin2+2,ixmax1-2,ixmax2-2,iw_filter,1,ndim,t+dt,w1,t+
     &         dt,w)
        
      else
            call die('Error in Advance: typefilter='//typefilter1//
     &         ' is unknown or module is switched off!')
         end if
         call getboundary(t+dt,1,nw,1,ndim,w)
      endif

      call process(0,1,ndim,w)

      if(oktest)write(*,*)'Advance new w:',w(ixtest1,ixtest2,iwtest)

      return
      end

*=============================================================================
      subroutine advect(method,ixmin1,ixmin2,ixmax1,ixmax2,iws,idimmin,
     &   idimmax,w1,w) 
* Process w if nproc/=0:   		call process
* Add fluxes and unsplit sources in 
* directions idim=idimmin..idimmax:	call advect1
*
* Depending on typeadvance and implpar call advect1 several times

      include 'vacdef.f'

      character*10   method
      integer  ixmin1,ixmin2,ixmax1,ixmax2,iws(niw_),idimmin,idimmax
      double precision  w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),w(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw)

* For most Runge-Kutta type schemes one more full array is needed
* For classical RK4 another array is needed


*!!MEMORY Needed for typeadvance='adams2' only


      logical  firstsweep,lastsweep
      common/first/firstsweep,lastsweep
*-----------------------------------------------------------------------------

      oktest=index(teststr,'advect').ge.1
      if(oktest)write(*,*)'Advect method w:',method,' ',w(ixtest1,ixtest2,
     &   iwtest)

* For negative "nproc(1)" call process, if positive check whether this is the
* first sweep and if "it-itmin" is an integer multiple of "nproc(1)" 
* (the frequency of processing before the whole timestep)
* Processing is done in advance_impl for implicit methods
      if(nproc(1).ne.0.and.implpar.le.zero)then
         if(nproc(1).lt.0.or.(firstsweep.and.it-itmin.eq.((it-itmin)/
     &      nproc(1))*nproc(1)))call process(1,idimmin,idimmax,w)
      end if

* Typically use "method" and at least one extra variable w1
      do iw_3=1,nw
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w1(ix_1,ix_2,iw_3)=w(ix_1,ix_2,iw_3)
      enddo
      enddo
      enddo

      istep=0
      if(typeadvance.eq.'onestep')then
         call advect1(method,dt,ixmin1,ixmin2,ixmax1,ixmax2,iws,idimmin,
     &      idimmax,t,w1,t,w)

      else if(typeadvance.eq.'twostep')then
*        do predictor step with typepred method to calculate w1 from w, then
*        full step with method. Fluxes and unsplit sources are taken at w1.
         call advect1(typepred1,dt/2,ixmin1,ixmin2,ixmax1,ixmax2,iws,idimmin,
     &      idimmax,t     ,w,t,w1)
         call advect1(method   ,dt  ,ixmin1,ixmin2,ixmax1,ixmax2,iws,idimmin,
     &      idimmax,t+dt/2,w1,t,w)




      else
         write(*,*)'typeadvance=',typeadvance
         write(*,*)
     &      'Error in Advect: Unknown time integration method or RK is off'
         call die('Correct typeadvance or: cd src; setvac -on=rk; make vac')
      end if

      if(oktest)write(*,*)'Advect final w:',w(ixtest1,ixtest2,iwtest)

      firstsweep=.false.

      return
      end

*=============================================================================
      subroutine advect1(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,iws,
     &   idimmin,idimmax,qtC,wCT,qt,w)

* Process if not first advection and nproc<0 is set
* Advect w to w+qdt*dF(wCT)_idim/dx_idim+qdt*((idimmax-idimmin+1)/ndim)*S(wCT)
* getboundaries

      include 'vacdef.f'

      character*10   method
      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,iws(niw_),idimmin,idimmax,idim
      double precision  qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      logical  firstsweep,lastsweep
      common/first/firstsweep,lastsweep
*-----------------------------------------------------------------------------

      istep=istep+1

      if(index(teststr,'saveadvect1').ge.1) call savefile(fileout_,wCT)

      oktest=index(teststr,'advect1').ge.1
      if(oktest)write(*,*)'Advect1 istep,wCT,w:',istep,wCT(ixtest1,ixtest2,
     &   iwtest),w(ixtest1,ixtest2,iwtest)

* In the first step wCT=w thus wCT is already processed if there is processing.
* Otherwise for negative "nproc(2)" call process, if positive check whether 
* this is the first sweep and if "it-itmin" is an integer multiple of 
* "nproc(2)" (the frequency of processing before intermediate steps)
* No processing here for implicit methods
      if(istep.gt.1.and.nproc(2).ne.0.and.implpar.le.zero)then
         if(nproc(2).lt.0.or.(firstsweep.and.it-itmin.eq.((it-itmin)/
     &      nproc(2))*nproc(2)))call process(2,idimmin,idimmax,w)
      end if

* Shrink ixO^L in all directions by 2
      ixOmin1=ixImin1+2
      ixOmin2=ixImin2+2
      ixOmax1=ixImax1-2
      ixOmax2=ixImax2-2

              
      if(method.eq.'cd')then
         call centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &      ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)
      else if(method.eq.'cd4')then
         call centdiff4(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &      ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)



              
      else if(method.eq.'hancock')then
         call hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &      ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)
      else if(method.eq.'tvdlf'.or.method.eq.'tvdlf1'.or.method.eq.'tvdmu'.or.
     &   method.eq.'tvdmu1')then
         call tvdmusclf(.true.,method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)

            
                    
      else if(method.eq.'tvd'.or.method.eq.'tvd1')then
         call centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &      ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)
         call tvdlimit(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &      ixOmin2,ixOmax1,ixOmax2,iws,idimmin,idimmax,wCT,qt+qdt,w)

      else if(method.eq.'source')then
         if(sourceunsplit)call addsource2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,
     &      ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)
      else if(method.eq.'nul')then
*        There is nothing to do
*        HPF_ if(.false.)write(*,*)'This avoids an xlhpf compiler bug'
      else
         write(*,*)'Error in Advect1:',method,' is unknown or switched off!'
         call die('Error in Advect1:'//method//' is unknown or switched off!')
      end if

      call getboundary(qt+qdt,1,nw,1,ndim,w)

      if(oktest)write(*,*)'Advect1 final w:',w(ixtest1,ixtest2,iwtest)

      return
      end

*=============================================================================
      subroutine addsource2(qdt,ixIImin1,ixIImin2,ixIImax1,ixIImax2,ixOOmin1,
     &   ixOOmin2,ixOOmax1,ixOOmax2,iws,qtC,wCT,qt,w)

* Add source within ixOO for iws: w=w+qdt*S[wCT]

      include 'vacdef.f'

      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,ixIImin1,ixIImin2,ixIImax1,ixIImax2,ixOOmin1,ixOOmin2,
     &   ixOOmax1,ixOOmax2,iws(niw_)
      double precision  qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),
     &   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      oktest=index(teststr,'addsource').ge.1
      if(oktest)write(*,*)'Add Source qdt,wCT,w:',qdt,wCT(ixtest1,ixtest2,
     &   iwtest),w(ixtest1,ixtest2,iwtest)

* AddSource and SpecialSource may shrink ixO or expand ixI for derivatives 
      ixImin1=ixIImin1
      ixImin2=ixIImin2
      ixImax1=ixIImax1
      ixImax2=ixIImax2
      ixOmin1=ixOOmin1
      ixOmin2=ixOOmin2
      ixOmax1=ixOOmax1
      ixOmax2=ixOOmax2

      call specialsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

      call     addsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

      if(oktest)write(*,*)'wnew:',w(ixtest1,ixtest2,iwtest)

* If AddSource/SpecialSource shrunk ixO, getboundary is needed.
      if(ixOmin1.gt.ixOOmin1.or.ixOmin2.gt.ixOOmin2.or.ixOmax1.lt.ixOOmax1.or.
     &   ixOmax2.lt.ixOOmax2)then
         call getboundary(qt+qdt,1,nw,1,ndim,w)
         if(oktest)write(*,*)'wnew after getboundary:',w(ixtest1,ixtest2,
     &      iwtest)
      end if

      return
      end

*=============================================================================
      logical function timetofinish(time0)

* Finish when it or t reached its maximum expected value, or dt is too small,
* or residual is small enough. Other conditions may be included.

      include 'vacdef.f'

      double precision  time0, cputime
      logical  okfinish
*-----------------------------------------------------------------------------

      okfinish = it.ge.itmax .or. t.ge.tmax .or. dt.lt.dtmin .or.
     &    (it.gt.itmin.and.(residual.lt.residmin .or. residual.gt.residmax))

      if(cputimemax .lt. bigdouble .and. .not.okfinish) okfinish=
     &    cputimemax .le. cputime()-time0

      timetofinish=okfinish

      return
      end

*=============================================================================
      logical function timetosave(ifile)

* Save times are defined by either tsave(isavet(ifile),ifile) or 
* itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
* Other conditions may be included.

      include 'vacdef.f'

      integer  ifile
      logical  oksave
*-----------------------------------------------------------------------------

      oksave=.false.
      if(t.ge.tsave(isavet(ifile),ifile))then
         oksave=.true.
         isavet(ifile)=isavet(ifile)+1
      end if
      if(it.eq.itsave(isaveit(ifile),ifile))then
         oksave=.true.
         isaveit(ifile)=isaveit(ifile)+1
      end if
      if(it.eq.itsavelast(ifile)+ditsave(ifile)) oksave=.true.
      if(t .ge.tsavelast(ifile) +dtsave(ifile) ) oksave=.true.
      if(oksave)then
         tsavelast(ifile) =t
         itsavelast(ifile)=it
      end if
      timetosave=oksave

      return
      end

*=============================================================================
      subroutine getdt_courant(w,ixmin1,ixmin2,ixmax1,ixmax2)

* Ensure that the courant conditions is met
* Calculate the time for the  maximum propagation speed cmax_i to cross dx_i
* in each i directions then take minimum for all grid points in the mesh and 
* for all i directions, finally multiply by courantpar.
*
* If TVD or TDLF provides dtcourant(idim) we take the minimum of those.
* In case of generalized coordinates dtcourant(idim) is correct due to the
* rotations while the value calculated here does not use a rotation.

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),cmax(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),courantmax,dtold
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim
      logical  new_cmax
*-----------------------------------------------------------------------------

      oktest=index(teststr,'getdt').ge.1

      if(oktest) write(*,*)'getdt_courant'

      dtold=dt
      dt=bigdouble
      courantmax=zero
      new_cmax=.true.
      do idim=1,ndim
         if(dtcourant(idim).lt.bigdouble)then
*           If dtcourant(idim) is calculated, use it
*!!      if(it==itmin+1)write(*,*)'second order correction in dt_courant!!!'
*!!      dt=min(dt,dtcourant(idim),dtcourant(idim)**2/dtold,1.1*dtold)
            dt=min(dt,dtcourant(idim))
            if(oktest) write(*,*)'idim,dtcourant(idim)',idim,dtcourant(idim)
         else
*           dx>0, but cmax>=0 may actually be 0, thus we calculate 
*           max(cmax/dx) rather than min(dx/cmax).

            call getcmax(new_cmax,w,ixmin1,ixmin2,ixmax1,ixmax2,idim,cmax)
            maxval_1=-bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               maxval_1=max(maxval_1,cmax(ix_1,ix_2)/dx(ix_1,ix_2,idim))
            enddo
            enddo
            courantmax=max(courantmax,maxval_1)

            if(gencoord.and.it.eq.itmin+1.and.verbose)write(*,*
     &         )'Warning in GetDtCourant: for gencoord approx. only',
     &         ', better use TVD-type methods'
            if(oktest) write(*,*)'idim,cmax:',idim,cmax(ixtest1,ixtest2)
            maxval_1=-bigdouble
            do ix_2=ixmin2,ixmax2
            do ix_1=ixmin1,ixmax1
               maxval_1=max(maxval_1,cmax(ix_1,ix_2)/dx(ix_1,ix_2,idim))
            enddo
            enddo
            if(oktest)then
               write(*,*)'max(c/dx)',maxval_1
            endif
         endif
      end do

      if(index(teststr,'dtdecline').lt.1)then
         do idim=1,ndim
            dtcourant(idim)=bigdouble
         enddo
      endif
      if(courantmax.gt.smalldouble) dt=min(dt,courantpar/courantmax)

      if(oktest) write(*,*)'GetDtCourant dt=',dt

      return 
      end

*=============================================================================
      double precision function cputime()

* Return cputime in seconds as a double precision number.
* For g77 compiler replace F77_ with F77_ everywhere in this function
* so that f90tof77 does not touch the system_clock function.


      real  etime,total,tarray(2)
      external etime
*HPF_ real:: timef
*-----------------------------------------------------------------------------

*     No timing
      cputime=-1.D0



      total = etime(tarray)
      cputime=tarray(1)
*HPF_ cputime=timef()/1.0D3

*cputime=second()   ! Cray CF77 or F90 (total CPU time for more CPU-s)
*cputime=secondr()  ! Cray CF77 or F90 (elapsed time for more CPU-s)

      return
      end
*=============================================================================
* end module vac
*##############################################################################
