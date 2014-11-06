*##############################################################################
* module vacio

*=============================================================================
      subroutine readparameters(w)

* This subroutine sets or reads all default parameters from par/DEFAULT,
* then it reads the par/PROBLEM parameter file through standard input,
* and the initial data from data/PROBLEM.ini as soon as the filename is read 
* from the parameter file.

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      character*10   typepred(nw),typefull(nw),typeimpl(nw),typefilter(nw)
      double precision  muscleta
      integer  i,j,k,iw,idim,iB,ifile,isave
      logical  implmrpc

* The use of NAMELIST is not recommended by Fortran 90. It could be replaced
* by some string manipulations, but that is difficult in Fortran 77/Adaptor.
* I use NAMELIST, since it is simple and convenient.

      namelist /testlist/   teststr,ixtest1,ixtest2,ixtest3,iwtest,idimtest
      namelist /filelist/   filenameini,filename,varnames, typefileini,
     &   typefileout,typefilelog,snapshotini,snapshotout,fullgridini,
     &   fullgridout,dixBmin1,dixBmin2,dixBmax1,dixBmax2
      namelist /savelist/   tsave,itsave,dtsave,ditsave
      namelist /stoplist/   itmax,tmax,tmaxexact,dtmin,residmin,residmax,t,it,
     &   cputimemax
      namelist /methodlist/ wnames,fileheadout,eqpar,typeadvance,typefull,
     &   typepred,typeimpl,typefilter,typelimiter,typeentropy,entropycoef,
     &   artcomp,typelimited,typefct,typetvd,typeaxial,typepoisson,
     &   typeconstrain,constraincoef,useprimitive,muscleta,musclomega,
     &   acmwidth,acmnolim,acmcoef,acmexpo,fourthorder,implmrpc,dimsplit,
     &   typedimsplit,sourcesplit,typesourcesplit,sourceunsplit,divbfix,
     &   divbwave,divbconstrain,angmomfix,compactres,smallfix,smallp,
     &   smallpcoeff,smallrho,smallrhocoeff,vacuumrho,nproc,procpar
      namelist /boundlist/  nB,typeB,ixBmin,ixBmax,idimB,upperB,extraB,
     &   typeBscalar,ipairB
      namelist /paramlist/  courantpar,dtpar,dtdiffpar,dtcantgrow,slowsteps,
     &   implmrpcpar,implpar,impldiffpar,implerror,implrelax,impldwlimit,
     &   implrestart,implrestart2,impliter,impliternr,typeimplinit,
     &   typeimpliter,typeimplmat,implnewton,implconserv,implcentered,
     &   implnewmat,implpred,impl3level,impljacfast,implsource
*-----------------------------------------------------------------------------

*!! Set new scalars for sake of unaltered par/DEFAULT
      constraincoef=one
      cputimemax=bigdouble
*!!

* Set default values for arrays (except the ones read from the inifile)

      do ifile=1,nfile
         do isave=1,nsavehi
*           t  of saves into the output files 
            tsave(isave,ifile)=bigdouble
*           it of saves into the output files 
            itsave(isave,ifile)=biginteger
         end do
*        time between saves
         dtsave(ifile)=bigdouble
*        timesteps between saves
         ditsave(ifile)=biginteger
*        index for saves by t
         isavet(ifile)=1
*        index for saves by it
         isaveit(ifile)=1
      end do

      do iw=1,nw
*        Predictor scheme (will be adjusted later)
         typepred(iw)='default'
*        Full step scheme
         typefull(iw)='tvdlf'
*        Implicit step scheme
         typeimpl(iw)='nul'
*        Filter scheme
         typefilter(iw)='nul'
*        Limiter type for flow variables/characteristics
         typelimiter(iw)='minmod'
*        Entropy fix type
         typeentropy(iw)='nul'
*        No artificial compression for Harten type TVD
         artcomp(iw)=.false.
      end do

      do iw=1,nw
*        Coefficients (0,1) for the dissipative fluxes
         acmcoef(iw)=-one
*     negative value means no coefficient is used
      enddo

*     Elements define processing for the fullstep,
      do i=1,nfile+2
*        halfstep and the nfile output files. If the value
         nproc(i)=0
*     is 0, no processing. For nproc(1) and nproc(2)
      end do
*                               the value defines the proc. frequency. Negative
*                               value results in a call at every sweep. Positive
*                               value N results in a call at every N-th step before
*                               the first sweep. For nproc(ifile+2) the nonzero 
*                               values cause processing for that file.
      do i=1,nprocpar
*        Parameters for processing
         procpar(i)=-one
      end do
*     Default small pressure, redefined in getpthermal
      smallp=-one
*     Default small density, redefined in keeppositive
      smallrho=-one
*     Density representing vacuum
      vacuumrho=-one

*     If nB is not specified by the user, gridsetup 
      nB=2*ndim
*     will create 2*ndim boundary regions by default.
      do iB=1,nhiB
         do iw=1,nw             
*           Default boundary type
            typeB(iw,iB) ='cont'
*           Fixed boundaries are not extrapolated into yet
            fixedB(iw,iB)=.false.
         end do
*        periodic pair is unknown, but can be set or guessed
         ipairB(iB)=0
      end do
      do iw=1,nw
         do idim=1,ndim
*           No zero flux condition for variables
            nofluxB(iw,idim)=.false.
         enddo
      end do
*     Default width of boundary regions
      dixBmin1=2
      dixBmin2=2
      dixBmax1=2
      dixBmax2=2
*     An impossible value if user specifies boundaries.
      ixBmax(1,1)=0

* Read scalar parameters from the par/DEFAULT file

      unitpar=unitini-1
      open(unitpar,file='par/DEFAULT',status='old')

      read(unitpar,testlist)
      read(unitpar,filelist)
      read(unitpar,savelist)
      read(unitpar,stoplist)
      read(unitpar,methodlist)
      read(unitpar,boundlist)
      read(unitpar,paramlist)

      close(unitpar)

* end defaults

* Initialize Kronecker delta, and Levi-Civita tensor
      do i=1,3
         do j=1,3
            if(i.eq.j)then
               kr(i,j)=1
            else
               kr(i,j)=0
            endif
            do k=1,3
               if(i.eq.j.or.j.eq.k.or.k.eq.i)then
                  lvc(i,j,k)=0
               else if(i+1.eq.j.or.i-2.eq.j)then
                  lvc(i,j,k)=1
               else
                  lvc(i,j,k)=-1
               endif
            enddo
         enddo
      enddo

* Initialize error conunters and equation parameters
      do i=1,nerrcode
         nerror(i)=0
      end do
      do i=1,neqpar+nspecialpar
         eqpar(i)=zero
      end do

* read from STDIN
      unitpar=unitstdin


* Start reading parameters from standard input, i.e. "< par/PROBLEM"

      read(unitpar,testlist)

      oktest=index(teststr,'readparameters').ge.1
      if(oktest) write(unitterm,*)'ReadParameters'
      if(oktest) write(unitterm,testlist)

      varnames='default'
      read(unitpar,filelist)


      if(oktest) then 
         
         write(unitterm,*)filenameini
         do ifile=1,nfile
            write(unitterm,*)filename(ifile)
         enddo
         write(unitterm,*)'Type of ini/out and log files:',typefileini,
     &      typefileout,typefilelog
         if(varnames.ne.'default')write(unitterm,*)'Varnames:',varnames
         if(snapshotini.gt.0)write(unitterm,*)'Snapshotini:',snapshotini
         if(snapshotout.gt.0)write(unitterm,*)'Snapshotout:',snapshotout
         write(unitterm,*)'Fullgridini,out:',fullgridini,fullgridout
      endif

      call readfileini(w)

* Default for output header line
      fileheadout=fileheadini



      read(unitpar,savelist)
      do ifile=1,nfile
         if(dtsave(ifile).lt.bigdouble/2.and.oktest) write(unitterm,
     &      '(" DTSAVE  for file",i2," =",g10.5)') ifile,dtsave(ifile)
         if(ditsave(ifile).lt.biginteger.and.oktest) write(unitterm,
     &      '(" DITSAVE for file",i2," =",i10)') ifile,ditsave(ifile)
         if(tsave(1,ifile).eq.bigdouble.and.itsave(1,ifile).eq.biginteger.and.
     &       dtsave(ifile).eq.bigdouble.and.ditsave(ifile).eq.biginteger.and.
     &      verbose) write(uniterr,*)'Warning in ReadParameters: ',
     &      'No save condition for file ',ifile
      enddo

      read(unitpar,stoplist)
      if(oktest)then
         if(tmax.lt.bigdouble)         write(unitterm,*) 'TMAX= ',tmax
         if(tmaxexact.and.oktest)   write(unitterm,*) 'TMAXEXACT=',tmaxexact
         if(itmax.lt.biginteger)       write(unitterm,*) 'ITMAX=',itmax
         if(dtmin.gt.smalldouble)      write(unitterm,*) 'DTMIN=',dtmin
         if(residmin.gt.smalldouble)   write(unitterm,*) 'RESIDMIN=',residmin
         if(residmax.lt.bigdouble)     write(unitterm,*) 'RESIDMAX=',residmax
         if(cputimemax.lt.bigdouble) write(unitterm,*) 'CPUTIMEMAX=',
     &      cputimemax
      endif

      if(tmax.eq.bigdouble.and.itmax.eq.biginteger) write(uniterr,*
     &   ) 'Warning in ReadParameters: Neither tmax nor itmax are given!'

      read(unitpar,methodlist)


      typefull1='nul'
      typepred1='nul'
      typeimpl1='nul'
      typefilter1='nul'
      iw_full(niw_)=0
      iw_impl(niw_)=0
      iw_semi(niw_)=0
      iw_filter(niw_)=0
      do iw=1,nw
          if(typefull(iw).ne.'nul')then
              typefull1=typefull(iw)
              iw_full(niw_)=iw_full(niw_)+1
              iw_full(iw_full(niw_))=iw
          end if
          if(typepred(iw).ne.'nul')then
              typepred1=typepred(iw)
          endif
          if(typeimpl(iw).ne.'nul')then
              typeimpl1=typeimpl(iw)
              iw_impl(niw_)=iw_impl(niw_)+1
              iw_impl(iw_impl(niw_))=iw
          end if
          if(typefull(iw).ne.'nul'.and.typeimpl(iw).eq.'nul')then
              iw_semi(niw_)=iw_semi(niw_)+1
              iw_semi(iw_semi(niw_))=iw
          end if
          if(typefilter(iw).ne.'nul')then
              typefilter1=typefilter(iw)
              iw_filter(niw_)=iw_filter(niw_)+1
              iw_filter(iw_filter(niw_))=iw
          end if
      end do

      if(typefull1.eq.'source'.and.dimsplit.and.sourceunsplit) write(*,*
     &   )'Warning: dimensional splitting for source terms!'

      if(typepred1.eq.'default')then
              if(typefull1.eq.'fct')then
                     typepred1='fct'
              else if(typefull1.eq.'tvdlf'.or.typefull1.eq.'tvdmu')then
                     typepred1='hancock'
              else if(typefull1.eq.'source')then
                     typepred1='source'
              else if(typefull1.eq.'tvdlf1'.or.typefull1.eq.'tvdmu1'.or.
     &           typefull1.eq.'tvd'.or.typefull1.eq.'tvd1'.or.
     &           typefull1.eq.'tvdmc'.or.typefull1.eq.'cd'.or.
     &           typefull1.eq.'cd4'.or.typefull1.eq.'mc'.or.
     &           typefull1.eq.'nul')then
                     typepred1='nul'
              else
                     call die('No default predictor for full step='//
     &                  typefull1)
           end if
      endif

      muscleta1=(one-muscleta)/4
      muscleta2=(one+muscleta)/4

      polargrid=.false.
      if(typeaxial.eq.'cylinder')then
         if(phi_.gt.ndir)call die(
     &      'Error: typeaxial=cylinder with phi>ndir is impossible!')
         if(phi_.le.ndim.and.phi_.ge.2)then
            polargrid=.true.
            write(*,*) 'Using polar coordinates...'
            if(ndir.eq.3.and.(z_.lt.2.or.z_.gt.ndir.or.z_.eq.
     &         phi_))call die( 
     &'z direction is not set correctly! Use setup -z=? and make vac')
            if(gencoord)write(*,*)
     &         'generalized coordinates in the r,z plane...'
         endif
      endif

      if(angmomfix)then
         write(*,*)
     &      'Angular momentum conservation is on (angmomfix=T) for iw=',mphi_
         if(typeaxial.ne.'cylinder') call die(
     &      'angmomfix works in cylindrical symmetry only!'//
     &      ' Set typeaxial in par file')
         if(phi_.lt.2.or.phi_.gt.ndir) call die(
     &      'phi direction does not satisfy 1<phi<=ndir! Use setvac -phi!')
         if(gencoord)then
            write(*,*)'angmomfix=T in generalized coordinates...'
            if(polargrid.and.phi_.ne.ndir) call die(
     &         'phi=ndir is required for angmomfix in gen. coordinates!')
         endif
      endif

* Artificial compression requires Harten type TVD
      do iw=1,nw
         if(artcomp(iw))typetvd='harten'
      enddo

* Harmonize the parameters for dimensional splitting and source splitting
      if(typedimsplit   .eq.'default'.and.     dimsplit)   typedimsplit='xyyx'
      if(typedimsplit   .eq.'default'.and..not.dimsplit)   typedimsplit=
     &   'unsplit'
      if(typesourcesplit.eq.'default'.and.     sourcesplit)typesourcesplit=
     &   'sfs'
      if(typesourcesplit.eq.'default'.and..not.sourcesplit)typesourcesplit=
     &   'unsplit'
      dimsplit   = typedimsplit   .ne.'unsplit'
      sourcesplit= typesourcesplit.ne.'unsplit'
      if(sourcesplit)sourceunsplit=.false.


      if(.not.divbconstrain.or.b0_.eq.0)typeconstrain='nul'
      if(typeconstrain.eq.'nul')divbconstrain=.false.

      if(typeconstrain.ne.'nul')then
         divbfix=.false.
*        ???divbwave=.false.???
         do i_1=1,nfile+2
            nproc(i_1)=0
         enddo
         if(typeaxial.ne.'slab')write(*,*)'constrainb=T keeps div B = 0 ',
     &      'only approximately in cylindrical symm.'
      endif

      if(divbfix.and.(typephys.eq.'mhd'.or.typephys.eq.'mhdiso').and.
     &   (.not.sourcesplit).and.(.not.sourceunsplit))call die(
     &   'divbfix=T requires unsplitsource=T or splitsource=T !')

      if(typefilter1.eq.'tvdmu'.or.typefilter1.eq.'tvdlf')typelimited=
     &   'predictor'
      if((.not.dimsplit.or.typeimpl1.eq.'tvdlf1'.or.typeimpl1.eq.
     &   'tvdmu1').and.typelimited.eq.'original')typelimited='previous'


      read(unitpar,boundlist)
      if(nB.gt.nhiB) call die(
     &   'Error in ReadParameters: Too many boundary regions')
      if(ixBmax(1,1).eq.0.and.oktest)then
         do iB=1,nB
            write(unitterm,'(" TYPEB(",i2,")       = ",100a10)') iB,(typeB(iw,
     &         iB),iw=1,nw)
         end do
         write(unitterm,*)' EXTRAB = ',extraB
      end if
      if(ixBmax(1,1).ne.0.and.oktest) write(unitterm,boundlist)
      if(dixBmin1.gt.dixBlo.or.dixBmin2.gt.dixBlo.or.dixBmax1.gt.dixBhi.or.
     &   dixBmax2.gt.dixBhi) call die( 
     &   'Error in ReadParameters: adjust dixBlo,dixBhi')

* Initialize implicit parameters as needed
      implsource=sourceunsplit
      if(typeimpl1.eq.'nul')then
*        Explicit time integration
         implpar=-one
         if(implmrpc)write(*,*)'Warning: implmrpc=T but typeimpl is not set!'
         typeimplinit = 'unused'
         typeimplmat  = 'unused'
      else if(implmrpc)then
*        MR-PC time integration
         implpar=one
         typeimpliter='vac_mrpc'
         typeimplmat ='free'
         impliter=       1
         implrestart=    5
         implerror=      1.D-3
         if(residmin.gt.0)then
*           MR-PC for steady state
            impl3level=     .false.
            typeimplinit=   'nul'
         else
*           MR-PC for time accurate
            impl3level=     .true.
            typeimplinit=   'explicit2'
         endif
      else if(typeimpl1.eq.'source')then
*        Semi-implicit time integration for sources

         if(.not.sourceunsplit)call die('Implicit sources must be unsplit!')

         if(residmin.gt.zero)then
            write(*,*)'Warning: Semi-implicit sources for steady state!'
            implerror=residmin
         else
            implerror=1.D-5
         endif

*        Explicit fluxes, implicit sources --> Trapezoidal scheme
         impl3level=.false.
         implpar=half
         typeimplinit='explicit'

*        Sources only --> Matrix free approach
         typeimplmat='free'

      else
*        Implicit time integration
         if(residmin.gt.zero)then
*           Steady state calculation --> Backward Euler scheme
            impl3level=.false.
            implpar=one
            typeimplinit='nul'
            implerror=residmin
            if(iw_semi(niw_).gt.0)write(*,*)
     &         'Warning: some variables are explicit for steady state!'
         else
*           Implicit fluxes (and sources) --> BDF2 scheme
            impl3level=.true.
*           First step is Backward Euler, but it could be trapezoidal !!!
            implpar=one
            typeimplinit='nul'
            implerror=1.D-3
         endif
         
         
             typeimpliter='vac_bicg'
             typeimplmat='prec'
        
         if((typeimpl1.eq.'tvdlf1'.or.typeimpl1.eq.'cd').and.
     &      .not.sourceunsplit)then
             impljacfast=.true.
         endif
      endif

      read(unitpar,paramlist)


      if(typeimpl1.ne.'nul')then
*        Check and/or set parameters for implicit calculations

         if(implpar.le.zero)call die(
     &      'Error: implpar<=0 for implicit calculation!')
         if(implpar.gt.one)call die(
     &      'Error: implpar>1 for implicit calculation!')

         if(implmrpc.and.typeimpliter.ne.'vac_mrpc')call die(
     &      'Error: implmrpc=T requires typeimpliter=vac_mrpc')

         if(implpred)then
            implconserv=.false.
            typeadvance='onestep'
            implpar=one
         endif

         if(impl3level)then
            implpred=.false.
            implconserv=.false.
         endif

         if(typeimpliter.eq.'tridiag'.and.ndim.eq.1)typeimplmat='with'
         if(typeimplmat.eq.'prec'.and.ndim.eq.1)typeimpliter='tridiag'

*        steady state warnings
         if(residmin.gt.zero)then
            if(implpar.ne.one)write(*,*)'Warning: implpar<1 for steady state!'
            if(impl3level)write(*,*)'Warning: 3 level for steady state!'
         endif

         if(impljacfast.and.typeimpl1.ne.'tvdlf1'.and.typeimpl1.ne.
     &      'cd')call die(
     &      'Error: impljacfast=T works with typeimpl=tvdlf1 or cd!')

         if(typeimpliter.eq.'vac_mrpc')then
            if(implnewton)call die(
     &         'Error: MRPC with Newton-Raphson iteration!')
            if(residmin.gt.zero.and.typeimplinit.ne.'nul')write(*,*
     &         )'Warning: MRPC for steady state',
     &         ' should not use typeimplinit=',typeimplinit
         endif

      endif

* If all predictor steps are 'nul' then 'twostep' method reduces to 'onestep'
      if(typeadvance.eq.'twostep'.and.typepred1.eq.'nul')typeadvance='onestep'

      if(typefull1.eq.'tvd'.and.typeadvance.ne.'onestep')call 
     &   die('tvd method should only be used as a onestep method')

      if(typefull1.eq.'tvd' .and. .not.dimsplit)write(*,*)
     &   'Warning: One step TVD without dimensional splitting !!!'

      if(typeadvance.eq.'onestep'.or.typeadvance.eq.'adams2')then
            nstep=1
      else if(typeadvance.eq.'twostep')then
            nstep=2
      else if(typeadvance.eq.'threestep')then
            nstep=3
      else if(typeadvance.eq.'fourstep'.or.typeadvance.eq.'sterck'.or.
     &   typeadvance.eq.'jameson')then
            nstep=4
      else
            call die('Unknown typeadvance='//typeadvance)
      end if

      if(oktest) write(unitterm,methodlist)

      if(oktest) write(unitterm,paramlist)

      call setheaderstrings

      return 
      end

*=============================================================================
      subroutine readfileini(w)

* Reads from unitini named filenameini in typefilini format.
*
* The file may contain more than one snapshots, in which case the last set is 
* read. The compatibility of initial data with internal parameters is checked.

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      logical  fileexist
      character*91  fhead
*-----------------------------------------------------------------------------

      if(typefileini.eq.'auto')then
         inquire(FILE=filenameini,EXIST=fileexist)
         if(.not.fileexist) call die(
     &      'Stop: file does not exist, filenameini='//filenameini)
         open(unitini,FILE=filenameini,STATUS='old')
         read(unitini,'(a91)')fhead
         close(unitini)

               if(ichar(fhead(1:1)).ne.0.and.ichar(fhead(2:2)).ne.0.and.
     &             ichar(fhead(3:3)).ne.0.and.ichar(fhead(4:4)).ne.0)then
            typefileini='ascii'
               else if(ichar(fhead(89:89)).eq.0.and.ichar(fhead(90:90)).eq.
     &            0.and. (ichar(fhead(88:88)).eq.24.or.ichar(fhead(91:91)).eq.
     &            24))then
            typefileini='binary'
         else
            typefileini='special'
         endif
         if(verbose)then
            write(*,*)'Auto typefileini=',typefileini
                     if(typefileini.eq.'special'.and. ichar(fhead(89:89)).eq.
     &                  0.and.ichar(fhead(90:90)).eq.0.and.
     &                   (ichar(fhead(88:88)).eq.20.or.ichar(fhead(91:91)).eq.
     &                  20)) write(*,*)'Looks like a real*4 file'
         endif
      endif
      if(typefileout.eq.'auto')then
         typefileout=typefileini
         if(verbose)write(*,*)'Auto typefileout=',typefileout
      endif

      if(typefileini.eq.'ascii')then
            call readfileini_asc(w)
      else if(typefileini.eq.'binary')then
            call readfileini_bin(w)
      else if(typefileini.eq.'special')then
            call readfileini_special(w)
      else
            call die('Error in VAC: Unknown typefileini='//typefileini)
      end if

      return
      end

*=============================================================================
      subroutine readfileini_asc(w)

* Reads from unitini, filenameini in ASCII format.
*
* The file may contain more than one snapshots, in which case the last set is 
* read. The compatibility of the initial data with internal parameters checked.
*
* Variables in the order they are read from the file:
*
*   fileheadini - a header identifying the input file
*   it,t        - the initial timestep and time
*   ndimini     - dimensionality of grid,   Test: ==ndim
*   neqparini   - number of eq. parameters, Test: <=neqpar+nspecialpar
*   nwini       - number of flow variables, Test: ==nw
*   nx          - the grid dimensions,      Test: <=ixGhi-dixBmax-ixMmin+1
*   eqpar       - equation parameters from filenameini
*   varnamesini - names of the coordinates, variables, equation parameters
*                 eg. 'x y rho mx my e bx by  gamma eta'
*   x           - the (initial) coordinates
*   w           - the initial flow variables

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      logical  fileexist
*     0 if not EOF, -1 if EOF, >0 if error
      integer  ios
*     values describing input data
      integer  ndimini,neqparini,neqparin,nwini,nwin
      integer  ixmin1,ixmin2,ixmax1,ixmax2,ix1,ix2,idim,iw,ieqpar,snapshot
      double precision  eqparextra,wextra
      character*79   varnamesini
*-----------------------------------------------------------------------------

      oktest=index(teststr,'readfileini').ge.1

      if(oktest) write(unitterm,*)'ReadFileIni'

      inquire(file=filenameini,exist=fileexist)
      if(.not.fileexist) call die('Stop: file does not exist, filenameini='//
     &   filenameini)
      open(unitini,file=filenameini,status='old')

      snapshot=0
7701  continue
          read(unitini,'(a)',iostat=ios,end=100)fileheadini
*             Cycle until the last recorded state
              if(ios.lt.0)goto 7702  
              if(oktest)then
                 write(unitterm,*)'fileheadini=',fileheadini(1:30)
              endif
          read(unitini,*,iostat=ios)it,t,ndimini,neqparini,nwini
              if(oktest) write(unitterm, 
     &           "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")
     &           it,t,ndimini,neqparini,nwini
              gencoord= ndimini.lt.0
              call checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)
          read(unitini,*,iostat=ios)nx
              if(oktest) write(unitterm,"('nx =',3i4)")nx
              call setixGixMix(ixmin1,ixmin2,ixmax1,ixmax2)
          read(unitini,*,iostat=ios)(eqpar(ieqpar),ieqpar=1,neqparin),
     &       (eqparextra,ieqpar=neqparin+1,neqparini)
              if(oktest) write(unitterm,*)eqpar
          read(unitini,'(a)',iostat=ios)varnamesini
          if(varnames.eq.'default')varnames=varnamesini
              if(oktest) write(unitterm,*)varnames

          do ix2=ixmin2,ixmax2
          do ix1=ixmin1,ixmax1
          read(unitini,*,iostat=ios)(x(ix1,ix2,idim),idim=1,ndim),(w(ix1,ix2,
     &       iw),iw=1,nwin),(wextra,iw=nwin+1,nwini)
          enddo
          enddo
          if(ios.ne.0)then
              write(uniterr,*)'Stop: iostat=',ios
              call die('Error in reading file')
          end if
          snapshot=snapshot+1
          if(snapshot.eq.snapshotini)goto 7702  
      goto 7701  
7702  continue

100       continue

      close(unitini)

      if(oktest) write(*,*)'x,w:',x(ixtest1,ixtest2,idimtest),w(ixtest1,
     &   ixtest2,iwtest)
      if(oktest)then
         write(*,*)'x,w:',x(ixtest1,ixtest2,idimtest),(w(ixtest1,ixtest2,
     &      iw_1),iw_1=1,nw)
      endif

      return
      end

*=============================================================================
      subroutine readfileini_bin(w)

* Reads from unitini,filenameini in binary format.
*
* The file may contain more than one snapshots, in which case the last set is 
* read unless snapshotini is set. 
* The compatibility of initial data with internal parameters is checked.

      include 'vacdef.f'

      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

      logical  fileexist
*     0 if not EOF, -1 if EOF, >0 if error
      integer  ios
*     values describing input data
      integer  ndimini,neqparini,neqparin,nwini,nwin
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim,iw,ieqpar,snapshot
      double precision  eqparextra
      character*79   varnamesini
*-----------------------------------------------------------------------------

      oktest=index(teststr,'readfileini').ge.1

      if(oktest) write(unitterm,*)'ReadFileIni'

      inquire(file=filenameini,exist=fileexist)
      if(.not.fileexist) call die('Stop: file does not exist, filenameini='//
     &   filenameini)
      open(unitini,file=filenameini,status='old',form='unformatted')

      snapshot=0
7703  continue
          read(unitini,iostat=ios,end=100)fileheadini
*             Cycle until the last recorded state
              if(ios.lt.0)goto 7704  
              if(oktest)then
                 write(unitterm,*)'fileheadini=',fileheadini(1:30)
              endif
          read(unitini,iostat=ios)it,t,ndimini,neqparini,nwini
              if(oktest) write(unitterm, 
     &           "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")
     &           it,t,ndimini,neqparini,nwini
              gencoord= ndimini.lt.0
              call checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)
          read(unitini,iostat=ios)nx
              if(oktest) write(unitterm,"('nx =',3i4)")nx
              call setixGixMix(ixmin1,ixmin2,ixmax1,ixmax2)
          read(unitini,iostat=ios)(eqpar(ieqpar),ieqpar=1,neqparin),
     &       (eqparextra,ieqpar=neqparin+1,neqparini)
              if(oktest) write(unitterm,*)eqpar
          read(unitini,iostat=ios)varnamesini
          if(varnames.eq.'default')varnames=varnamesini
              if(oktest) write(unitterm,*)varnames

          read(unitini,iostat=ios)(((x(ix_1,ix_2,idim),ix_1=
     &       ixmin1,ixmax1),ix_2=ixmin2,ixmax2),idim=1,ndim)
*         To conform savefileout_bin we use loop for iw
          do iw=1,nwin
             read(unitini,iostat=ios)((w(ix_1,ix_2,iw),ix_1=
     &          ixmin1,ixmax1),ix_2=ixmin2,ixmax2)
          end do
          if(ios.ne.0)then
              write(uniterr,*)'Error in ReadFileIni: iostat=',ios
              call die('Error in reading file')
          end if
          snapshot=snapshot+1
          if(snapshot.eq.snapshotini)goto 7704  
      goto 7703  
7704  continue

100       continue

      close(unitini)

      if(oktest) write(*,*)'x,w:',x(ixtest1,ixtest2,idimtest),w(ixtest1,
     &   ixtest2,iwtest)
      if(oktest)then
         write(*,*)'x,w:',x(ixtest1,ixtest2,idimtest),(w(ixtest1,ixtest2,
     &      iw_1),iw_1=1,nw)
      endif

      return
      end

*=============================================================================
      subroutine checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)

      include 'vacdef.f'
      integer  ndimini,neqparini,nwini,neqparin,nwin
*-----------------------------------------------------------------------------

      if(ndim.ne.abs(ndimini))then
         write(*,*)'Error in ReadFileini: ndimini=',ndimini
         call die('Incompatible dimensionalities')
      endif

      if(neqpar+nspecialpar.ne.neqparini)write(*,"(a,i3,a,i3)")
     &   'Warning in ReadFileini: number of eq.params=',neqpar,
     &   ' /= neqparini=',neqparini

      if(nw.ne.nwini)write(*,"(a,i3,a,i3)")
     &   'Warning in ReadFileini: number of variables nw=',nw,' /= nwini=',
     &   nwini

      if((neqpar+nspecialpar.ne.neqparini.or.nw.ne.nwini).and.
     &   varnames.eq.'default')call die(
     &   'Define varnames (in &filelist for VAC, in 3rd line for VACINI)!')

* The number of equation parameters and variables to be read
      neqparin=min(neqparini,neqpar+nspecialpar)
      nwin=min(nwini,nw)

      return
      end


*=============================================================================
      subroutine setixGixMix(ixmin1,ixmin2,ixmax1,ixmax2)

      include 'vacdef.f'
      integer  ixmin1,ixmin2,ixmax1,ixmax2
*-----------------------------------------------------------------------------

      ixGmin1=ixGlo1
      ixGmin2=ixGlo2
      ixMmin1=ixGmin1+dixBmin1
      ixMmin2=ixGmin2+dixBmin2

* Shave off ghost cells from nx
      if(fullgridini)then
         
                nx(1)=nx(1)-dixBmin1
          nx(1)=nx(1)-dixBmax1
         
         
                nx(2)=nx(2)-dixBmin2
          nx(2)=nx(2)-dixBmax2
         
      endif

* Calculate mesh and grid sizes
      ixMmax1=ixMmin1+nx(1)-1
      ixMmax2=ixMmin2+nx(2)-1
      ixGmax1=ixMmax1+dixBmax1
      ixGmax2=ixMmax2+dixBmax2

* Set the index range for this grid
      if(fullgridini)then
         ixmin1=ixGmin1
         ixmin2=ixGmin2
         ixmax1=ixGmax1
         ixmax2=ixGmax2
         
      else
         ixmin1=ixMmin1
         ixmin2=ixMmin2
         ixmax1=ixMmax1
         ixmax2=ixMmax2
      endif

      if(ixGmax1.gt.ixGhi1.or.ixGmax2.gt.ixGhi2)then
         write(uniterr,*)'Stop: nxhi=',ixGhi1-dixBmax1-ixMmin1+
     &      1,ixGhi2-dixBmax2-ixMmin2+1
         call die('Error in SetixGixMix')
      end if

      nx1=nx(1)
      nx2=nx(2)



      return
      end

*=============================================================================
      subroutine setheaderstrings

* Check and/or put physics and equation parameter names into file header

      include 'vacdef.f'

      integer  i
      character*10   physics
*-----------------------------------------------------------------------------

* Check physics or add _typephysNDIMNDIR

      write(physics,'(a,i1,i1)')typephys,2,2

      i=index(fileheadout,'_')
      if(i.ge.1)then
               if(physics.ne.fileheadout(i+1:i+10))then
            write(*,*)'This code is configured to ',physics
                     call die('Error:  physics in file is '//
     &                  fileheadout(i+1:i+10))
         endif
      else
         i=79-len(typephys)-3
7705     continue
                     if(fileheadout(i:i).ne.' ' .or. i.eq.1)goto 7706  
            i=i-1
         goto 7705  
7706     continue
               fileheadout=fileheadout(1:i)//'_'//physics
         write(*,*)'Warning: incomplete input headline.',
     &      ' Added to output headline _',physics
      endif

* Check for equation parameter names in varnames, add them if missing

      if(varnames.ne.'default' .and. index(varnames,eqparname).le.0)then
         i=79-len(eqparname)-3
7707     continue
                     if(varnames(i:i).ne.' ' .or. i.eq.1)goto 7708  
            i=i-1
         goto 7707  
7708     continue
               varnames=varnames(1:i)//'   '//eqparname
      endif

* Check for special equation parameter names in varnames, add them if missing

      if(varnames.ne.'default' .and. index(varnames,specialparname).le.0)then
         i=79-len(specialparname)-3
7709     continue
                     if(varnames(i:i).ne.' ' .or. i.eq.1)goto 7710  
            i=i-1
         goto 7709  
7710     continue
               varnames=varnames(1:i)//'   '//specialparname
      endif

      return
      end

*=============================================================================
      subroutine savefile(ifile,w)

      include 'vacdef.f'
      integer  ifile,ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      character*10  itstring
*-----------------------------------------------------------------------------

      if(nproc(ifile+2).gt.0) call process(ifile+2,1,ndim,w)

* In most cases the mesh should be saved
      ixmin1=ixMmin1
      ixmin2=ixMmin2
      ixmax1=ixMmax1
      ixmax2=ixMmax2

      if(ifile.eq.fileout_)then
*           Produce the output file name
            filenameout=filename(fileout_)
            if(snapshotout.gt.0.and.isaveout.gt.0)then
               if(isaveout.eq.snapshotout*(isaveout/snapshotout))then
                  close(unitini+ifile)
                  write(itstring,'(i10)')isaveout+1
                                 filenameout=filenameout(1:index(filenameout,
     &                              ' ')-1)//'_'// itstring(10-
     &                              int(alog10(isaveout+1.5)):10) 
               endif
            endif
            isaveout=isaveout+1
            if(fullgridout)then
               ixmin1=ixGmin1
               ixmin2=ixGmin2
               ixmax1=ixGmax1
               ixmax2=ixGmax2
               
            end if
            if(typefileout.eq.'ascii')then
                  call savefileout_asc(unitini+ifile,w,ixmin1,ixmin2,ixmax1,
     &               ixmax2)
            else if(typefileout.eq.'binary')then
                  call savefileout_bin(unitini+ifile,w,ixmin1,ixmin2,ixmax1,
     &               ixmax2)
            else if(typefileout.eq.'special')then
                  call savefileout_special(unitini+ifile,w,ixmin1,ixmin2,
     &               ixmax1,ixmax2)
            else
                  call die('Error in SaveFile: Unknown typefileout:'//
     &               typefileout)
            end if
      else if(ifile.eq.filelog_)then
            if(typefilelog.eq.'default')then
                  call savefilelog_default(unitini+ifile,w,ixmin1,ixmin2,
     &               ixmax1,ixmax2)
            else if(typefilelog.eq.'special')then
                  call savefilelog_special(unitini+ifile,w,ixmin1,ixmin2,
     &               ixmax1,ixmax2)
            else
                  call die('Error in SaveFile: Unknown typefilelog:'//
     &               typefilelog)
            end if
      else
            write(*,*) 'No save method is defined for ifile=',ifile
            call die(' ')
      end if

      return 
      end

*=============================================================================
      subroutine savefileout_asc(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)

* This version saves into filename(fileout_) ASCII data at every save time
* in full accordance with the ReadFileini subroutine, except that the first
* line is fileheadout and not fileheadini.

      include 'vacdef.f'

      integer  qunit,ixmin1,ixmin2,ixmax1,ixmax2,ix1,ix2,iw,idim,ndimout
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),qw(nw)
      logical  fileopen
*-----------------------------------------------------------------------------

      inquire(qunit,opened=fileopen)

      if(.not.fileopen)open(qunit,file=filenameout,status='unknown')

      if(gencoord)then
         ndimout= -ndim
      else
         ndimout= ndim
      endif

      write(qunit,"(a)")fileheadout
      write(qunit,"(i7,1pe13.5,3i3)")it,t,ndimout,neqpar+nspecialpar,nw
      write(qunit,"(3i4)") ixmax1-ixmin1+1,ixmax2-ixmin2+1
      write(qunit,"(100(1pe13.5))")eqpar
      write(qunit,"(a)")varnames
      do ix2= ixmin2,ixmax2 
      do ix1= ixmin1,ixmax1 
*        Values with magnitude less than smalldouble are written as 0d0
         do iw_1=1,nw
            if(abs(w(ix1,ix2,iw_1)).gt.5.0d-16)then
                  qw(iw_1)=w(ix1,ix2,iw_1)
               else
                  qw(iw_1)=0d0
               endif
         enddo
         write(qunit,"(100(1pe18.10))")(x(ix1,ix2,idim_1),idim_1=
     &      1,ndim),(qw(iw_1),iw_1=1,nw)
      enddo
      enddo

      call flushunit(qunit)

      return 
      end

*=============================================================================
      subroutine savefileout_bin(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)

* This version saves into filename(fileout_) binary data at every save time
* in full accordance with the ReadFileini subroutine, except that the first
* line is fileheadout and not fileheadini.

      include 'vacdef.f'

      integer  qunit,ixmin1,ixmin2,ixmax1,ixmax2,idim,iw,ndimout
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      logical  fileopen
*-----------------------------------------------------------------------------

      inquire(qunit,opened=fileopen)
      if(.not.fileopen)open(qunit,file=filenameout,status='unknown',form=
     &   'unformatted')

      if(gencoord)then
         ndimout= -ndim
      else
         ndimout= ndim
      endif

      write(qunit)fileheadout
      write(qunit)it,t,ndimout,neqpar+nspecialpar,nw
      write(qunit) ixmax1-ixmin1+1,ixmax2-ixmin2+1
      write(qunit)eqpar
      write(qunit)varnames
      write(qunit)(((x(ix_1,ix_2,idim),ix_1=ixmin1,ixmax1),ix_2=
     &   ixmin2,ixmax2),idim=1,ndim)

* write(qunit)w(ix^S,1:nw) produces segmentation fault on Alpha, thus loop

      do iw=1,nw
         write(qunit)((w(ix_1,ix_2,iw),ix_1=ixmin1,ixmax1),ix_2=ixmin2,ixmax2)
      end do

      call flushunit(qunit)

      return 
      end

*=============================================================================
      subroutine savefilelog_default(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)

* This version saves into filename(filelog_) the following formatted data:
*
*   fileheadout
*   STRING_DESCRIBING_COLUMNS
*   it t dt wmean(1) wmean(2) ... wmean(nw) residual
*   it t dt wmean(1) wmean(2) ... wmean(nw) residual
*   it t dt wmean(1) wmean(2) ... wmean(nw) residual
*   etc.
*
* at every save time. wmean is the volume averaged w, residual is saved 
* if residmin>0 is set in the parfile.

      include 'vacdef.f'

      integer  qunit,ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  iw
      logical  fileopen
      double precision  wmean(nw)
*-----------------------------------------------------------------------------



      inquire(qunit,opened=fileopen)
      if(.not.fileopen)then
         open(qunit,file=filename(filelog_),status='unknown')
         write(qunit,'(a)')fileheadout
         if(residmin.gt.zero.or.residmax.lt.bigdouble)then
            write(qunit,'(a15,a55,a9)')'it   t   dt    ',wnames,' residual'
         else
            write(qunit,'(a15,a64)')   'it   t   dt    ',wnames
         endif
      endif



      do iw=1,nw 
         sum_1=0.
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            sum_1=sum_1+dvolume(ix_1,ix_2)*w(ix_1,ix_2,iw)
         enddo
         enddo
         wmean(iw)=sum_1/volume
         
      end do



      if(residmin.gt.zero.or.residmax.lt.bigdouble)then
         write(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean,residual
      else
         write(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean
      endif
      call flushunit(qunit)



      return 
      end

*=============================================================================
      subroutine die(message)

      character(*)   message
*-----------------------------------------------------------------------------

      write(*,*)message

      stop
      end
*=============================================================================
      subroutine flushunit(qunit)

*use F90_UNIX_IO,only : flush ! F90=f95 (NAG)
      implicit none

      integer qunit, ierror

*call flush(qunit,ierror) ! F90=f95 (NAG)
*     OS=Linux, SunOS, UNICOS, T3E, Fujitsu
      call flush(qunit)
*call flush_(qunit)  ! OS=AIX, F90=xlf

* no flush on Linux IA64 with Intel compiler

      return
      end
*=============================================================================
* end module vacio
*##############################################################################



