*##############################################################################
* module vacini

*=============================================================================
      subroutine process(count,idimmin,idimmax,w)

* Process w before it is advected in directions idim^LIM, or before save
* count=1 and 2 for first and second (half step) processing during advection
* count=ifile+2 for saving results into the file indexed by ifile

      include 'vacdef.f'

      integer  count,idimmin,idimmax
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

*-----------------------------------------------------------------------------

      return
      end
*=============================================================================
*=============================================================================
      program vacini

      include 'vacdef.f'

      integer  iw,ieqpar,ixmin1,ixmin2,ixmax1,ixmax2
      character*20   typeini
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),wpar(nw)
      logical  lastiw
*-----------------------------------------------------------------------------


      verbose=.true. 
      if(verbose)then
         write(*,'(a)')'VACINI 4.52 configured to'
         write(*,'(a)')'  -d=22 -phi=0 -z=0 -g=260,260 -p=mhd -u=gravity'
         write(*,'(a)')'  -on=cd,tvdlf,tvd,poisson,resist'
         write(*,'(a)')'  -off=mc,fct,impl,ct,gencoord,rk,mpi'
         
      endif

* Some default values
      t=zero
       it=0
       
      typefileout='ascii'
       typefileini='auto'
      snapshotini=0
      fullgridini=.false.
       fullgridout=.false.
      gencoord=   .false.
* There are no ghost cells in VACINI except when "readmesh" is used.
      dixBmin1=0
      dixBmin2=0
      dixBmax1=0
      dixBmax2=0
       ixMmin1=ixGlo1
      ixMmin2=ixGlo2
       
* Test cell
      ixtest1=ixMmin1
      ixtest2=ixMmin2
* Read parameters from STDIN
      unitpar=unitstdin



      if(verbose)write(*,*)'Filename for new initial file:'
      read(unitpar,'(a)')filename(fileout_)

      if(verbose)write(*,*)'Fileheader:'
      read(unitpar,'(a)')fileheadout
      if(verbose)write(*,*)'Variable names, e.g. "x y rho m1 m2":'
      read(unitpar,'(a)')varnames

      call setheaderstrings

      if(verbose)then
         write(*,*)'Select action by typing one of the following words: '
         write(*,*)'   test'
         write(*,*)'   typefileini,snapshotini,read,readmesh,readnext,',
     &      'typefileout,write,save'
         write(*,*)'   domain,grid,sheargrid,shiftgrid,polargrid,',
     &      'roundgrid,rotategrid'
         write(*,*)'   transpose,regrid,stretchgrid,stretchgridcent'
         write(*,*)'   polarvar,spherevar,cartesian,rotatevar,',
     &      'setvar,perturbvar,conserve,primitive,multiply,divide'
         write(*,*)'   uniform,shocktube,wave,wave1,special'
         write(*,*)'   eqpar,gencoord'
      endif

7701  continue
         if(verbose)write(*,*)'Action:'
         read(unitpar,'(a)')typeini
         if(verbose)write(*,*)'> ',typeini
         if(typeini.eq.'verbose')then
                if(verbose)write(*,*)'Verbose:'
                read(unitpar,*)verbose
                
         else if(typeini.eq.'test')then
                if(verbose)write(*,*)'Teststring:'
                read(unitpar,'(a)')teststr
                if(verbose)write(*,*)'ixtest, idimtest and iwtest:'
                read(unitpar,*) ixtest1,ixtest2,idimtest,iwtest
                
         else if(typeini.eq.'typefileini')then
                if(verbose)write(*,*)
     &             'Type of old initial file: ascii/binary/special'
                read(unitpar,'(a)')typefileini
         else if(typeini.eq.'snapshotini')then
                if(verbose)write(*,*)'Number of snapshot to be read:'
                read(unitpar,*)snapshotini
         else if(typeini.eq.'read'.or.typeini.eq.'readmesh')then
                if(verbose)write(*,*)'Filename for old initial file:'
                read(unitpar,'(a)')filenameini
                
                if(typeini.eq.'readmesh')then
                    if(verbose)write(*,*)
     &                 'Specify boundary width for old initial file:'
                    read(unitpar,*) dixBmin1,dixBmin2,dixBmax1,dixBmax2
                    fullgridini=.true.
                endif
                call readfileini(w)
         else if(typeini.eq.'readnext')then
                snapshotini=snapshotini+1
                call readfileini(w)
         else if(typeini.eq.'typefileout')then
                if(verbose)write(*,*)
     &             'Type of new initial file: ascii/binary/special'
                read(unitpar,'(a)')typefileout
         else if(typeini.eq.'write')then
               call savefile(fileout_,w)
         else if(typeini.eq.'save')then
               call savefile(fileout_,w)
               close(unitini+fileout_)
               goto 7702  
         else if(typeini.eq.'grid')then
               call ini_grid(.true.,ixMmin1,ixMmin2,ixMmax1,ixMmax2)
         else if(typeini.eq.'domain')then
               call ini_grid(.false.,ixMmin1,ixMmin2,ixMmax1,ixMmax2)
         else if(typeini.eq.'shiftgrid')then
               
               call shiftgrid(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'sheargrid')then
               
               call sheargrid(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'polargrid')then
                call makepolargrid(ixMmin1,ixMmin2,ixMmax1,ixMmax2)
               gencoord=.true.
               if(.false.)call die('Polar grid is meaningless in 1D')
         else if(typeini.eq.'spheregrid')then
               call die('Spherical grid is meaningful in 3D only')
         else if(typeini.eq.'roundgrid')then
                call roundgrid(ixMmin1,ixMmin2,ixMmax1,ixMmax2)
               gencoord=.true. 
               if(.false.)call die('Round grid is meaningless in 1D')
         else if(typeini.eq.'rotategrid')then
               call rotatevar(ixMmin1,ixMmin2,ixMmax1,ixMmax2,1,ndim,x)
               gencoord=.true.
         else if(typeini.eq.'transpose')then
               
                call transposexy(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
               if(.false.)call die('Transpose is implemented in 2D only.')
         else if(typeini.eq.'regrid')then
               
               call regrid(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'stretchgrid')then
               
               call stretchgrid(.true.,ixMmin1,ixMmin2,ixMmax1,ixMmax2)
         else if(typeini.eq.'stretchgridcent')then
               
               call stretchgrid(.false.,ixMmin1,ixMmin2,ixMmax1,ixMmax2)
         else if(typeini.eq.'polarvar')then
                call polarvar(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
               if(.false.)call die('Polar variables are meaningless in 1D')
         else if(typeini.eq.'spherevar')then
               call die('Spherical variables are meaningful '// 'in 3D only')
         else if(typeini.eq.'cartesian')then
                call cartesian(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
               if(.false.) call die('Polar variables are meaningless in 1D')
         else if(typeini.eq.'rotatevar')then
               call rotatevar(ixMmin1,ixMmin2,ixMmax1,ixMmax2,1,nw,w)
         else if(typeini.eq.'setvar')then
               if(verbose)write(*,*)'Give ix limits:'
               read(unitpar,*) ixmin1,ixmin2,ixmax1,ixmax2
               
7703           continue
                  if(verbose)write(*,*) 
     &               'Variable index, variable value, lastiw (T/F)?'
                  read(unitpar,*)iw,wpar(iw),lastiw
                  do ix_2=ixmin2,ixmax2
                  do ix_1=ixmin1,ixmax1
                     w(ix_1,ix_2,iw)=wpar(iw)
                  enddo
                  enddo
                  if(lastiw)goto 7704  
               goto 7703  
7704           continue
         else if(typeini.eq.'perturbvar')then
               call perturbvar(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'multiply')then
               if(verbose)write(*,*)
     &            'Give multiplying factors for each variable:'
               read(unitpar,*) (wpar(iw_1),iw_1=1,nw)
               do iw=1,nw
                  do ix_2=ixMmin2,ixMmax2
                  do ix_1=ixMmin1,ixMmax1
                     w(ix_1,ix_2,iw)=w(ix_1,ix_2,iw)*wpar(iw)
                  enddo
                  enddo
               end do
         else if(typeini.eq.'divide')then
               if(verbose)write(*,*)'Give dividing factors for each variable:'
               read(unitpar,*) (wpar(iw_1),iw_1=1,nw)
               do iw=1,nw
                  do ix_2=ixMmin2,ixMmax2
                  do ix_1=ixMmin1,ixMmax1
                     w(ix_1,ix_2,iw)=w(ix_1,ix_2,iw)/wpar(iw)
                  enddo
                  enddo
               end do
         else if(typeini.eq.'conserv'.or.typeini.eq.'conserve')then
               call conserve(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'primitive')then
               call primitive(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'uniform')then
               if(verbose)write(*,*)'Give values for each variable:'
               read(unitpar,*) (wpar(iw_1),iw_1=1,nw)
               do iw=1,nw
                  do ix_2=ixMmin2,ixMmax2
                  do ix_1=ixMmin1,ixMmax1
                     w(ix_1,ix_2,iw)=wpar(iw)
                  enddo
                  enddo
               end do
         else if(typeini.eq.'shocktube')then
               call ini_shocktube(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'wave')then
               call ini_wave(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'wave1')then
               call wave1(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'special')then
               call specialini(ixMmin1,ixMmin2,ixMmax1,ixMmax2,w)
         else if(typeini.eq.'eqpar')then
               if(verbose)write(*,*)'Equation params:',neqpar+nspecialpar
               read(unitpar,*)(eqpar(ieqpar),ieqpar=1,neqpar+nspecialpar)
         else if(typeini.eq.'gencoord')then
               if(verbose)write(*,*)'Generalized coordinates (T/F):'
               read(unitpar,*)gencoord
         else
               call die('Error in VACIni: no such action')
         end if
      goto 7701  
7702  continue



      end

*=============================================================================
      subroutine ini_grid(coord,ixmin1,ixmin2,ixmax1,ixmax2)

* Setup a uniform grid. When coord is .true., the user provides the coordinates
* for the centers of the grid, otherwise the boundaries of the computational
* domaines, thus the centers start at xmin+dx/2, and end at xmax-dx/2

      include 'vacdef.f'

      logical  coord
      integer  ixmin1,ixmin2,ixmax1,ixmax2,ix1,ix2,idim
      double precision  dx1,dx2,xmax(ndim),xmin(ndim)
*-----------------------------------------------------------------------------

      if(verbose)write(*,'(a,3i6)')'Size of mesh. Max: ',ixGhi1,ixGhi2
      read(unitpar,*) nx1,nx2
      if(verbose)write(*,'(a,3i6)')'Size of mesh: ',nx1,nx2
      if(coord)then
         if(verbose)write(*,*)'Coordinates of cell centers at the edges'
      else
         if(verbose)write(*,*)'Boundaries of the computational domain'
      endif
      if(verbose)write(*,*)'xmin coordinates:'
      read(unitpar,*)(xmin(idim),idim=1,ndim)
      if(verbose)write(*,*)'xmax coordinates:'
      read(unitpar,*)(xmax(idim),idim=1,ndim)

* Calculate cell sizes and modify coordinates for 'domain' action
      if(coord)then
         dx1=(xmax(1)-xmin(1))/(nx1-1)
         dx2=(xmax(2)-xmin(2))/(nx2-1)
      else
         dx1=(xmax(1)-xmin(1))/nx1
         dx2=(xmax(2)-xmin(2))/nx2
         
         xmax(1)=xmax(1)-dx1/2
         xmin(1)=xmin(1)+dx1/2
         
         xmax(2)=xmax(2)-dx2/2
         xmin(2)=xmin(2)+dx2/2
        
      endif



      ixmax1=ixmin1+nx1-1
      ixmax2=ixmin2+nx2-1
      if(ixmax1.gt.ixGhi1.or.ixmax2.gt.ixGhi2) call die(
     &   'Error in IniGrid: Too big grid')



      do ix1=ixmin1,ixmax1
      do ix2=ixmin2,ixmax2
         x(ix1,ix2,1)= ((ix1-ixmin1)*xmax(1)+ (ixmax1-ix1)*
     &      xmin(1)) /(ixmax1-ixmin1) 
      enddo
      enddo
      do ix1=ixmin1,ixmax1
      do ix2=ixmin2,ixmax2
         x(ix1,ix2,2)= ((ix2-ixmin2)*xmax(2)+ (ixmax2-ix2)*
     &      xmin(2)) /(ixmax2-ixmin2) 
      enddo
      enddo

      return
      end

*=============================================================================

      subroutine makepolargrid(ixmin1,ixmin2,ixmax1,ixmax2)

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  pi2
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)
     &   'First coordinate is interpreted as radius, second as angle/2pi.'

      pi2=8*atan(one)

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         tmp(ix_1,ix_2)=x(ix_1,ix_2,1)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         x(ix_1,ix_2,1)=tmp(ix_1,ix_2)*cos(x(ix_1,ix_2,2)*pi2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         x(ix_1,ix_2,2)=tmp(ix_1,ix_2)*sin(x(ix_1,ix_2,2)*pi2)
      enddo
      enddo

      return
      end

*=============================================================================

*=============================================================================
       
      subroutine roundgrid(ixmin1,ixmin2,ixmax1,ixmax2)

*      Calculate the shrink factor to shrink a rectangle to an ellipse. For a sguare
*   1               in directions parallel to x and y
*   1-r+r*sqrt(0.5) in diagonal directions, where r is the radial distance
*                      normalized to 1. 

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  dist1(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   dist2(ixGlo1:ixGhi1,ixGlo2:ixGhi2),weight(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
      double precision  xcent1,xcent2,rounded,squared
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Coordinates of center point:'
      read(unitpar,*) xcent1,xcent2
      if(verbose)then
         write(*,*)'Center of rounded grid:',xcent1,xcent2
         write(*,*)'If the rectangle is mapped to the (-1,-1,1,1) square'
         write(*,*)'give the distances to be rounded, and to be squared:'
      endif
      read(unitpar,*)rounded,squared

* Normalized distances in the L1 and L2 norms

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         dist1(ix_1,ix_2)=max(abs((x(ix_1,ix_2,1)-xcent1)/(x(ixmax1,ixmax2,1)-
     &      xcent1)),abs((x(ix_1,ix_2,2)-xcent2)/(x(ixmax1,ixmax2,2)-xcent2)))
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         dist2(ix_1,ix_2)=sqrt(((x(ix_1,ix_2,1)-xcent1)/(x(ixmax1,ixmax2,1)-
     &      xcent1))**2+((x(ix_1,ix_2,2)-xcent2)/(x(ixmax1,ixmax2,2)-xcent2))
     &      **2)
      enddo
      enddo

* The weight increases from 0 to 1 for distance between 0 and rounded, and
* it drops back to 0 in the distance range rounded and squared. 
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         if(dist1(ix_1,ix_2).lt.rounded)then
            weight(ix_1,ix_2)=dist1(ix_1,ix_2)/rounded
         else
            weight(ix_1,ix_2)=(squared-dist1(ix_1,ix_2))/(squared-rounded)
         endif
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         weight(ix_1,ix_2)=min(one,max(zero,weight(ix_1,ix_2)))
      enddo
      enddo

* Shrink all coordinates by the factor 1 + weight*(dist1/dist2 - 1)
*      Where weight is 0 there is no distortion, where weight is 1 the grid is round

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         weight(ix_1,ix_2)=one + weight(ix_1,ix_2)*(dist1(ix_1,ix_2)/
     &      dist2(ix_1,ix_2) - one)
      enddo
      enddo

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         x(ix_1,ix_2,1)=weight(ix_1,ix_2)*(x(ix_1,ix_2,1)-xcent1)+xcent1
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         x(ix_1,ix_2,2)=weight(ix_1,ix_2)*(x(ix_1,ix_2,2)-xcent2)+xcent2
         
      enddo
      enddo

      return
      end

*=============================================================================

      subroutine polarvar(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Given r=x(:,1) and phi=2*Pi*x(:,2) 
* rotate the v_r,v_phi vector components to v_x,v_y

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,iw,jw
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wi(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),wj(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
      double precision  cosphi(ixGlo1:ixGhi1,ixGlo2:ixGhi2),
     &   sinphi(ixGlo1:ixGhi1,ixGlo2:ixGhi2),pi2
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Indices of var_r and var_phi:'
      read(unitpar,*)iw,jw

      pi2=8*atan(one)
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         cosphi(ix_1,ix_2)=cos(x(ix_1,ix_2,2)*pi2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         sinphi(ix_1,ix_2)=sin(x(ix_1,ix_2,2)*pi2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wi(ix_1,ix_2)=w(ix_1,ix_2,iw)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wj(ix_1,ix_2)=w(ix_1,ix_2,jw)
      enddo
      enddo

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,iw)=cosphi(ix_1,ix_2)*wi(ix_1,ix_2)-sinphi(ix_1,ix_2)*
     &      wj(ix_1,ix_2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,jw)=cosphi(ix_1,ix_2)*wj(ix_1,ix_2)+sinphi(ix_1,ix_2)*
     &      wi(ix_1,ix_2)
      enddo
      enddo

      return
      end

*=============================================================================

*=============================================================================

      subroutine cartesian(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Given r=x(:,1) and phi=x(:,2) on a polar grid
* rotate the v_x,v_y vector components to v_r,v_phi

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,iw,jw
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wi(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),wj(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Indices of var_x and var_y:'
      read(unitpar,*)iw,jw

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wi(ix_1,ix_2)=w(ix_1,ix_2,iw)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         wj(ix_1,ix_2)=w(ix_1,ix_2,jw)
      enddo
      enddo

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,iw)=cos(x(ix_1,ix_2,2))*wi(ix_1,ix_2)+
     &      sin(x(ix_1,ix_2,2))*wj(ix_1,ix_2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,jw)=cos(x(ix_1,ix_2,2))*wj(ix_1,ix_2)-
     &      sin(x(ix_1,ix_2,2))*wi(ix_1,ix_2)
      enddo
      enddo

      return
      end

*=============================================================================
      subroutine shiftgrid(ixmin1,ixmin2,ixmax1,ixmax2,wnew)

* Shifts grid relative to the variables. Padding is done by the border values.

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,ix,ixinside,dix,idim,iw
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw)
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Give direction and size of shift: idim, dix'
      read(unitpar,*)idim,dix

      do iw_3=1,nw
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,iw_3)=wnew(ix_1,ix_2,iw_3)
      enddo
      enddo
      enddo
        if(idim.eq.1)then
           do iw=1,nw
              do ix= ixmin1,ixmax1
                 ixinside=min(ixmax1,max(ixmin1,ix-dix))
                 do ix_1=ixmin2,ixmax2
                    wnew(ix,ix_1,iw)=w(ixinside,ix_1,iw)
                 enddo
              end do
           end do 
        else if(idim.eq.2)then
           do iw=1,nw
              do ix= ixmin2,ixmax2
                 ixinside=min(ixmax2,max(ixmin2,ix-dix))
                 do ix_1=ixmin1,ixmax1
                    wnew(ix_1,ix,iw)=w(ix_1,ixinside,iw)
                 enddo
              end do
           end do 
        else
          call die('Error in ShiftGrid: Unknown direction')
      end if

      return
      end

*=============================================================================
      subroutine sheargrid(ixmin1,ixmin2,ixmax1,ixmax2,wnew)

* Shears grid relative to the variables. Padding is done by the border values.

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,ix,ix1,ix2,dix,idim1,idim2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wnew(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2,nw),angle
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)
     &   'Give directions and angle of shear: idim1,idim2,angle'
      read(unitpar,*)idim1,idim2,angle
      if(idim1.eq.idim2)call die('Error in ShearGrid: idim1==idim2')
      angle=angle*atan(one)/45

      do iw_3=1,nw
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,iw_3)=wnew(ix_1,ix_2,iw_3)
      enddo
      enddo
      enddo
      do ix1= ixmin1,ixmax1 
      do ix2= ixmin2,ixmax2 
        if(idim2.eq.1)then
              dix=nint(-tan(angle)*(ix1-ixmin1)) 
        else if(idim2.eq.2)then
              dix=nint(-tan(angle)*(ix2-ixmin2)) 
        else
              call die('Error in ShearGrid: Unknown 2nd direction')
         end if
         if(idim1.eq.1)then
              ix=min(ixmax1,max(ixmin1,ix1-dix))
              do iw_1=1,nw
                 wnew(ix1,ix2,iw_1)=w(ix,ix2,iw_1) 
              enddo
         else if(idim1.eq.2)then
              ix=min(ixmax2,max(ixmin2,ix2-dix))
              do iw_1=1,nw
                 wnew(ix1,ix2,iw_1)=w(ix1,ix,iw_1) 
              enddo
         else
               call die('Error in ShearGrid: Unknown 1st direction')
         end if
      enddo
      enddo

      return
      end

*=============================================================================
      subroutine multiply(exponent,ixmin1,ixmin2,ixmax1,ixmax2,w)

* Multiply all variables by some function of the first coordinate

      include 'vacdef.f'

      integer  exponent,ixmin1,ixmin2,ixmax1,ixmax2,ix
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),r(ixGlo1:ixGhi1)
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Define axial symmetry: cylinder/sphere/nozzle'
      read(unitpar,'(a)')typeaxial

      do ix_1=ixmin1,ixmax1
         r(ix_1)=x(ix_1,ixmin2,1)
      enddo

      if(typeaxial.eq.'cylinder')then
            do ix_1=ixmin1,ixmax1
               area(ix_1)=r(ix_1)
            enddo
      else if(typeaxial.eq.'sphere')then
            do ix_1=ixmin1,ixmax1
               area(ix_1)=r(ix_1)**2
            enddo
      else if(typeaxial.eq.'nozzle')then
            if(verbose)write(*,*)
     &         "Warning in Multiply: This is Yee's nozzle problem!"
            do ix_1=ixmin1,ixmax1
               area(ix_1) =1.398+0.347*(1-2/(exp(1.6*r(ix_1) -8)+1))
            enddo
      else
            call die('Error in Multiply: Unknown axial symmetry type!')
      end if

      do ix= ixmin1,ixmax1
         do iw_2=1,nw
         do ix_1=ixmin2,ixmax2
            w(ix,ix_1,iw_2)=w(ix,ix_1,iw_2)*area(ix)**exponent
         enddo
         enddo
      enddo

      return
      end

*=============================================================================
      subroutine rotatevar(ixmin1,ixmin2,ixmax1,ixmax2,iwmin,iwmax,w)

* Rotate a pair of vector variables around some axis

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,iwmin,iwmax,iw1,iw2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iwmin:iwmax),
     &   w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2),w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2),angle
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Give indices of 2 variables and angle in degrees'
      read(unitpar,*)iw1,iw2,angle
      angle=angle*atan(one)/45
      if(iw1.eq.iw2.or.iw1.lt.1.or.iw1.gt.nw.or.iw2.lt.1.or.
     &   iw2.gt.nw)call die( 'Error in RotateVar: Incorrect iw1 and iw2')

      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w1(ix_1,ix_2)=w(ix_1,ix_2,iw1)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w2(ix_1,ix_2)=w(ix_1,ix_2,iw2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,iw1)=cos(angle)*w1(ix_1,ix_2)-sin(angle)*w2(ix_1,ix_2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         w(ix_1,ix_2,iw2)=sin(angle)*w1(ix_1,ix_2)+cos(angle)*w2(ix_1,ix_2)
      enddo
      enddo

      return
      end


*=============================================================================
      subroutine transposexy(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Transpose the first two coordinates of the grid (x), the variables (w), 
* and exchange the vector components as required

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,ixoldmin1,ixoldmin2,ixoldmax1,
     &   ixoldmax2,ix1,ix2,iw,ivect,idim,qnvector
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
*-----------------------------------------------------------------------------
      ixoldmin1=ixmin1
      ixoldmin2=ixmin2
      ixoldmax1=ixmax1
      ixoldmax2=ixmax2
      ixmin1=ixoldmin2
      ixmax1=ixoldmax2
       
      ixmin2=ixoldmin1
      ixmax2=ixoldmax1
* Transpose x
      do idim=1,ndim
*      !!For sake of f90tof77 x(ixmin1:ixmax1,ixmin2:ixmax2,idim)=transpose(x(ixoldmin1:ixoldmax1,ixoldmin2:ixoldmax2,idim)) is replaced:
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=x(ix_1,ix_2,idim)
         enddo
         enddo
         do ix1=ixmin1,ixmax1
      do ix2=ixmin2,ixmax2
            x(ix1,ix2,idim)=tmp(ix2,ix1)
         enddo
      enddo
      enddo

* Swap the X and Y coordinates
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         tmp(ix_1,ix_2)=x(ix_1,ix_2,1)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         x(ix_1,ix_2,1)=x(ix_1,ix_2,2)
      enddo
      enddo
      do ix_2=ixmin2,ixmax2
      do ix_1=ixmin1,ixmax1
         x(ix_1,ix_2,2)=tmp(ix_1,ix_2)
      enddo
      enddo

* Transpose w
      do iw=1,nw
*      !!For sake of f90tof77 w(ixmin1:ixmax1,ixmin2:ixmax2,iw)=transpose(w(ixoldmin1:ixoldmax1,ixoldmin2:ixoldmax2,iw)) is replaced by
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=w(ix_1,ix_2,iw)
         enddo
         enddo
         do ix1=ixmin1,ixmax1
      do ix2=ixmin2,ixmax2
            w(ix1,ix2,iw)=tmp(ix2,ix1)
         enddo
      enddo
      enddo

* Swap the vector variables
* qnvector is only used to avoid compiler warning when nvector=0
      qnvector=nvector
      do ivect=1,qnvector
         if(verbose)write(*,"(a,i1,a)")
     &      'Index of first component of vector variable #',ivect,':'
         read(unitpar,*)iw
         if(iw.ge.nw.or.iw.lt.1)call die(
     &      'Error in TransposeXY: Incorrect iw.')
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            tmp(ix_1,ix_2)=w(ix_1,ix_2,iw)
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            w(ix_1,ix_2,iw)=w(ix_1,ix_2,iw+1)
         enddo
         enddo
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            w(ix_1,ix_2,iw+1)=tmp(ix_1,ix_2)
         enddo
         enddo
      end do

      return
      end

*=============================================================================
      subroutine regrid(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Change the number of grid points and extrapolate and interpolate the
* original cell positions x and averaged values w.

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,nix1,nix2,ixnewmin1,ixnewmin2,
     &   ixnewmax1,ixnewmax2,iw,idim
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),q(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2)
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Define number of gridpoints in each direction:'
      read(unitpar,*) nix1,nix2

      ixnewmin1=ixmin1
      ixnewmin2=ixmin2
      ixnewmax1=ixmin1+nix1-1
      ixnewmax2=ixmin2+nix2-1

      do idim=1,ndim
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            q(ix_1,ix_2)=x(ix_1,ix_2,idim)
         enddo
         enddo
         call regrid1(ixmin1,ixmin2,ixmax1,ixmax2,ixnewmin1,ixnewmin2,
     &      ixnewmax1,ixnewmax2,q)
         do ix_2=ixnewmin2,ixnewmax2
         do ix_1=ixnewmin1,ixnewmax1
            x(ix_1,ix_2,idim)=q(ix_1,ix_2)
         enddo
         enddo
      enddo
      do iw=1,nw
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            q(ix_1,ix_2)=w(ix_1,ix_2,iw)
         enddo
         enddo
         call regrid1(ixmin1,ixmin2,ixmax1,ixmax2,ixnewmin1,ixnewmin2,
     &      ixnewmax1,ixnewmax2,q)
         do ix_2=ixnewmin2,ixnewmax2
         do ix_1=ixnewmin1,ixnewmax1
            w(ix_1,ix_2,iw)=q(ix_1,ix_2)
         enddo
         enddo
      enddo

      ixmax1=ixnewmax1
      ixmax2=ixnewmax2

      return
      end

*=============================================================================
      subroutine regrid1(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,
     &   ixOmax1,ixOmax2,q)

* Interpolate q from grid determined by ixI to ixO. Use the distances measured
* between the Cartesian gridpoints, i.e. distances in generalized coordinates.
* The DOMAIN ixImin-0.5..ixImax+0.5 is rediscretized by ixOmax-ixOmin+1 points.

      include 'vacdef.f'

      integer  ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2,ixI1,ixI2,ixO1,ixO2,dixI1,dixI2
      double precision  q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),qnew(ixGlo1:ixGhi1,
     &   ixGlo2:ixGhi2),dxO1,dxO2,xO1,xO2,coeff1(0:1),coeff2(0:1)
*-----------------------------------------------------------------------------

* Grid spacing of the output grid stretched onto the integer input grid
      dxO1=(ixImax1-ixImin1+one)/(ixOmax1-ixOmin1+one)
      dxO2=(ixImax2-ixImin2+one)/(ixOmax2-ixOmin2+one)

      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         qnew(ix_1,ix_2)=zero
      enddo
      enddo
      do ixO1=ixOmin1,ixOmax1
      do ixO2=ixOmin2,ixOmax2

*        Location of the output grid point
         xO1=ixImin1-half+(ixO1-ixOmin1+half)*dxO1
         xO2=ixImin2-half+(ixO2-ixOmin2+half)*dxO2

*        Index of the input grid point to the left of xO within ixImin..ixImax-1
         ixI1=min(ixImax1-1,max(ixImin1,int(xO1)))
         ixI2=min(ixImax2-1,max(ixImin2,int(xO2)))

*        Calculate bilinear interpolation/extrapolation coefficients
         coeff1(1)=xO1-ixI1
         coeff2(1)=xO2-ixI2
         coeff1(0)=1-coeff1(1)
         coeff2(0)=1-coeff2(1)

*        Interpolate q into qnew
         do dixI1=0,1
         do dixI2=0,1
            qnew(ixO1,ixO2)=qnew(ixO1,ixO2)+(coeff1(dixI1)*coeff2(dixI2))*
     &         q(ixI1+dixI1,ixI2+dixI2)
         enddo
         enddo

      enddo
      enddo

      do ix_2=ixOmin2,ixOmax2
      do ix_1=ixOmin1,ixOmax1
         q(ix_1,ix_2)=qnew(ix_1,ix_2)
         
      enddo
      enddo

      return
      end

*=============================================================================
      subroutine stretchgrid(qdomain,ixmin1,ixmin2,ixmax1,ixmax2)

* Stretch the grid logarithmically in direction idim segment by segment
* The original computational domain size is preserved if qdomain is true,
* and the first and last grid center locations are preserved if it is false.

      include 'vacdef.f'

      logical  qdomain
      integer  ixmin1,ixmin2,ixmax1,ixmax2,ix,ixL,ixR,idim,iseg,nseg
      integer  qixhi
      PARAMETER(  qixhi=10000)
      double precision  qxL,qxR,qdxL,qdxR,qdxsum,qdx(qixhi)
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Direction of stretch and number of segments'
      read(unitpar,*)idim,nseg

      if(idim.gt.ndim.or.idim.lt.1)call die(
     &   'Error in StretchGrid: Invalid direction')
      if(nseg.lt.1)call die(
     &   'Error in StretchGrid: Invalid number of segments')

      if(qdomain)then
         qxL=1.5*x(ixmin1,ixmin2,idim)-0.5*x(ixmin1+1,ixmin2+1,idim)
         qxR=1.5*x(ixmax1,ixmax2,idim)-0.5*x(ixmax1-1,ixmax2-1,idim)
         if(verbose)write(*,*)'Old domain from',qxL,' to',qxR
      else
         qxL=x(ixmin1,ixmin2,idim)
         qxR=x(ixmax1,ixmax2,idim)
         if(verbose)write(*,*)'Old centers from',qxL,' to',qxR
      endif

      if(qxL.ge.qxR)call die(
     &   'Error in StretchGrid: qxR<qxL, use "domain" or "grid" action')

      if(verbose)write(*,*)'xsecond-xfirst is scaled to 1 for segment 1'
      if(idim.eq.1)then
            if(ixmax1.gt.qixhi)call die( 
     &         'Error in StretchGrid: Too big grid, change qixhi')
            ixL=ixmin1+1
            qdxL=one
            qdx(ixL)=qdxL
            if(qdomain)then
                qdxsum=1.5D0
            else
                qdxsum=1.0D0
            endif
            do iseg=1,nseg
               if(verbose)write(*,*)'xlast-xprev for segment:',iseg
               read(unitpar,*)qdxR
               if(iseg.lt.nseg)then
                  if(verbose)write(*,*)'Cell number of end of segment:',iseg
                  read(unitpar,*)ixR
                  if(ixR.le.ixL.or.ixR.ge.ixmax1) call die(
     &               'Error in StretchGrid: bad segment position')
               else
                  if(verbose)write(*,*)'Last segment ends at right boundary'
                  ixR=ixmax1
               end if
               do ix=ixL+1,ixR
                  qdx(ix)=qdxL*(qdxR/qdxL)**(dble(ix-ixL)/dble(ixR-ixL))
                  qdxsum=qdxsum+qdx(ix)
               enddo
               qdxL=qdxR
               ixL=ixR
            end do
            if(qdomain)qdxsum=qdxsum+half*qdx(ixR)
            do ix_1=ixmin1+1,ixmax1
               qdx(ix_1)=qdx(ix_1)*(qxR-qxL)/qdxsum
            enddo
            if(qdomain)then
               do ix_1=ixmin2,ixmax2
                  x(ixmin1,ix_1,idim)=qxL+half*qdx(ixmin1+1)
               enddo
            endif
            do ix=ixmin1+1,ixmax1
               do ix_1=ixmin2,ixmax2
                  x(ix,ix_1,idim)=x(ix-1,ix_1,idim)+qdx(ix)
               enddo
            enddo
         
      else if(idim.eq.2)then
            if(ixmax2.gt.qixhi)call die( 
     &         'Error in StretchGrid: Too big grid, change qixhi')
            ixL=ixmin2+1
            qdxL=one
            qdx(ixL)=qdxL
            if(qdomain)then
                qdxsum=1.5D0
            else
                qdxsum=1.0D0
            endif
            do iseg=1,nseg
               if(verbose)write(*,*)'xlast-xprev for segment:',iseg
               read(unitpar,*)qdxR
               if(iseg.lt.nseg)then
                  if(verbose)write(*,*)'Cell number of end of segment:',iseg
                  read(unitpar,*)ixR
                  if(ixR.le.ixL.or.ixR.ge.ixmax2) call die(
     &               'Error in StretchGrid: bad segment position')
               else
                  if(verbose)write(*,*)'Last segment ends at right boundary'
                  ixR=ixmax2
               end if
               do ix=ixL+1,ixR
                  qdx(ix)=qdxL*(qdxR/qdxL)**(dble(ix-ixL)/dble(ixR-ixL))
                  qdxsum=qdxsum+qdx(ix)
               enddo
               qdxL=qdxR
               ixL=ixR
            end do
            if(qdomain)qdxsum=qdxsum+half*qdx(ixR)
            do ix_1=ixmin2+1,ixmax2
               qdx(ix_1)=qdx(ix_1)*(qxR-qxL)/qdxsum
            enddo
            if(qdomain)then
               do ix_1=ixmin1,ixmax1
                  x(ix_1,ixmin2,idim)=qxL+half*qdx(ixmin2+1)
               enddo
            endif
            do ix=ixmin2+1,ixmax2
               do ix_1=ixmin1,ixmax1
                  x(ix_1,ix,idim)=x(ix_1,ix-1,idim)+qdx(ix)
               enddo
            enddo
         
      end if

      if(qdomain)then
         qxL=1.5*x(ixmin1,ixmin2,idim)-0.5*x(ixmin1+1,ixmin2+1,idim)
         qxR=1.5*x(ixmax1,ixmax2,idim)-0.5*x(ixmax1-1,ixmax2-1,idim)
         if(verbose)write(*,*)'New domain from',qxL,' to',qxR
      else
         qxL=x(ixmin1,ixmin2,idim)
         qxR=x(ixmax1,ixmax2,idim)
         if(verbose)write(*,*)'New centers from',qxL,' to',qxR
      endif

      return
      end

*=============================================================================
      subroutine perturbvar(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Perturb a variable within limits ixP by adding the product of sine waves
* in each direction. The phases of the waves are relative to ixPmin.

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,ixPmin1,ixPmin2,ixPmax1,ixPmax2,
     &   idim,iw
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),dw,wavenum(ndim),
     &   phase(ndim)
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Variable index, amplitude and',' ixP^DL '
      read(unitpar,*)iw,dw,ixPmin1,ixPmax1,ixPmin2,ixPmax2

      if(verbose)write(*,*)'wavenumber and phase for each idim:'
      read(unitpar,*)(wavenum(idim),phase(idim),idim=1,ndim)

      do ix_2=ixPmin2,ixPmax2
      do ix_1=ixPmin1,ixPmax1
         w(ix_1,ix_2,iw)=w(ix_1,ix_2,iw)+dw* sin(phase(1)+(x(ix_1,ix_2,1)-
     &      x(ixPmin1,ixPmin2,1))*wavenum(1))*sin(phase(2)+
     &      (x(ix_1,ix_2,2)-x(ixPmin1,ixPmin2,2))*wavenum(2))
      enddo
      enddo

      return
      end

*=============================================================================
      subroutine ini_shocktube(ixmin1,ixmin2,ixmax1,ixmax2,w)

* The shocktube is divided into nseg segments in the chosen idim direction.
* Linear interpolation in segments with given left and right states.

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
      integer  ix,ixL,ixR,idim,iw,iseg,nseg,ieqpar
      double precision  wL(nw),wR(nw)

*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Normal of slab symmetry'
      read(unitpar,*)idim
      if(verbose)write(*,*)'Number of segments'
      read(unitpar,*)nseg
      if(verbose)write(*,*)'All variables at minimal position:'
      read(unitpar,*)(wL(iw_1),iw_1=1,nw)
      if(idim.eq.1)then
            ixL=ixmin1
            
            do iseg=1,nseg
               if(verbose)write(*,*)'All variables at end of segment:',iseg
               read(unitpar,*)(wR(iw_1),iw_1=1,nw)
               if(iseg.lt.nseg)then
                  if(verbose)write(*,*)'Cell number of end of segment:',iseg
                  read(unitpar,*)ixR
                  
                  if(ixR.le.ixL) call die(
     &               'Error in ini_shocktube: bad segment position')
               else
                  if(verbose)write(*,*)'Last segment ends at right boundary'
                  ixR=ixmax1      
                  
               end if
               do iw=1,nw
                  do ix=ixL,ixR
                    if(ix.ge.ixmin1.and.ix.le.ixmax1)then
                       do ix_1=ixmin2,ixmax2
                          w(ix,ix_1,iw)=((ixR-ix)*wL(iw)+(ix-ixL)*
     &                       wR(iw))/(ixR-ixL)
                       enddo
                    endif

*      w(ix,ixmin2:ixmax2,iw)=((x(ixR,ixmin2:ixmax2,1)-x(ix,ixmin2:ixmax2,1))*wL(iw)+ &
*      (x(ix,ixmin2:ixmax2,1)-x(ixL,ixmin2:ixmax2,1))*wR(iw))&
*      /(x(ixR,ixmin2:ixmax2,1)-x(ixL,ixmin2:ixmax2,1)) 
                  end do
               end do
               ixL=ixR
               do iw_1=1,nw
                  wL(iw_1)=wR(iw_1)
               enddo
            end do 
      else if(idim.eq.2)then
            ixL=ixmin2
            
            do iseg=1,nseg
               if(verbose)write(*,*)'All variables at end of segment:',iseg
               read(unitpar,*)(wR(iw_1),iw_1=1,nw)
               if(iseg.lt.nseg)then
                  if(verbose)write(*,*)'Cell number of end of segment:',iseg
                  read(unitpar,*)ixR
                  
                  if(ixR.le.ixL) call die(
     &               'Error in ini_shocktube: bad segment position')
               else
                  if(verbose)write(*,*)'Last segment ends at right boundary'
                  ixR=ixmax2      
                  
               end if
               do iw=1,nw
                  do ix=ixL,ixR
                    if(ix.ge.ixmin2.and.ix.le.ixmax2)then
                       do ix_1=ixmin1,ixmax1
                          w(ix_1,ix,iw)=((ixR-ix)*wL(iw)+(ix-ixL)*
     &                       wR(iw))/(ixR-ixL)
                       enddo
                    endif

*      w(ixmin1:ixmax1,ix,iw)=((x(ixmin1:ixmax1,ixR,2)-x(ixmin1:ixmax1,ix,2))*wL(iw)+ &
*      (x(ixmin1:ixmax1,ix,2)-x(ixmin1:ixmax1,ixL,2))*wR(iw))&
*      /(x(ixmin1:ixmax1,ixR,2)-x(ixmin1:ixmax1,ixL,2)) 
                  end do
               end do
               ixL=ixR
               do iw_1=1,nw
                  wL(iw_1)=wR(iw_1)
               enddo
            end do 
      else
            call die('Error in Ini_ShockTube: Unknown dimension')
      end if

      if(verbose)write(*,*)'Eqpar:'
      read(unitpar,*)(eqpar(ieqpar),ieqpar=1,neqpar+nspecialpar)

      return
      end

*=============================================================================
      subroutine ini_wave(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Sum of sine waves in each direction and each variable

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,ix1,ix2,idim,iw,ieqpar
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wmean(nw),
     &   wavenum(ndim),ampl(ndim),phase(ndim)
      logical  nextiw
*-----------------------------------------------------------------------------

      if(verbose)write(*,*)'Mean values:'
      read(unitpar,*)(wmean(iw_1),iw_1=1,nw)
      do iw=1,nw
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            w(ix_1,ix_2,iw)=wmean(iw)
         enddo
         enddo
      end do
      do iw=1,nw
         if(verbose)write(*,*)'Variable iw:',iw
7705     continue
            if(verbose)write(*,*)
     &         'Amplitude,wavenum,phase for each dim,nextiw=T/F:'
            read(unitpar,*)(ampl(idim),wavenum(idim),phase(idim),idim=
     &         1,ndim),nextiw
            do idim=1,ndim
                do ix_2=ixmin2,ixmax2
                do ix_1=ixmin1,ixmax1
                   w(ix_1,ix_2,iw)=w(ix_1,ix_2,iw)+ampl(idim)*
     &                sin(phase(idim)+x(ix_1,ix_2,idim)*wavenum(idim))
                enddo
                enddo
            enddo
            if(nextiw)goto 7706  
         goto 7705  
7706     continue
      enddo

      if(verbose)write(*,*)'Eqpar:'
      read(unitpar,*)(eqpar(ieqpar),ieqpar=1,neqpar+nspecialpar)

      return
      end

*=============================================================================
      subroutine wave1(ixmin1,ixmin2,ixmax1,ixmax2,w)

* Sine waves with arbitrary wave vectors and shifts using rationalized angle
* units (1.0 = full circle)

      include 'vacdef.f'

      integer  ixmin1,ixmin2,ixmax1,ixmax2,idim,iw
      double precision  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),wmean,ampl,
     &   wavenum(ndim),phase,pi2
      logical  lastiw
*-----------------------------------------------------------------------------

      pi2=8*atan(one)
      if(verbose)write(*,*)
     &   'w(iw)=wmean+dw*sin(2*pi*[x*kx+y*ky+phase]) (Note the 2*pi!)'
7707  continue
         if(verbose)write(*,*)
     &      'Give iw,wmean,dw,k(idim),phase,lastiw (quit with T):'
         read(unitpar,*)iw,wmean,ampl,(wavenum(idim),idim=1,ndim),phase,lastiw
         do ix_2=ixmin2,ixmax2
         do ix_1=ixmin1,ixmax1
            w(ix_1,ix_2,iw)=wmean+ampl*sin(pi2*(wavenum(1)*
     &         x(ix_1,ix_2,1)+wavenum(2)*x(ix_1,ix_2,2)+phase))
         enddo
         enddo
         if(lastiw)goto 7708  
      goto 7707  
7708  continue

      return
      end
*=============================================================================
* Some interface routines for subroutines often used in the VACUSR module
* to keep the compiler happy 
*=============================================================================
      subroutine gradient(realgrad,q,ixmin1,ixmin2,ixmax1,ixmax2,idir,gradq)
      logical  realgrad
      integer  ixmin1,ixmin2,ixmax1,ixmax2,idir
      double precision  q(*),gradq(*)
      call die('Error: VACINI cannot call gradient !')
      end

*=============================================================================
      subroutine laplace4(q,ixmin1,ixmin2,ixmax1,ixmax2,laplaceq)
      integer  ixmin1,ixmin2,ixmax1,ixmax2
      double precision  q(*),laplaceq(*)
      call die('Error: VACINI cannot call laplace4 !')
      end

*=============================================================================
      subroutine ensurebound(dix,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,
     &   ixOmin2,ixOmax1,ixOmax2,qt,w)
      integer  dix,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,
     &   ixOmax2
      double precision  qt,w(*)
      call die('Error: VACINI cannot call ensurebound !')
      end
*=============================================================================
* end module vacini
*##############################################################################


