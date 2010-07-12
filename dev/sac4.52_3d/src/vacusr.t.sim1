!##############################################################################
! module vacusr - sim1 ! setvac -d=22 -g=204,204 -p=hdadiab -u=sim1


INCLUDE:vacusr.gravity.t
INCLUDE:vacusr.viscosity.t

!=============================================================================
subroutine specialini(ix^L,w)

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,1:nw)

double precision:: rhoin,xcent^D,radius
double precision:: inirho,iniene
double precision:: onemor,inix,ddx
double precision:: p_1,p_2

integer:: iii_,iix_1,info,i,j
double precision:: pi,comi,eneu,sum,mode,bmax,l
character*79 atmfilename

integer:: ix_1,ix_2,ix_3

double precision:: ixc_1,ixc_2,ixc_3
double precision:: rfc,a1,a2,a3


double precision:: xc2,xc3,x_max,r_0,r_tr,b_Amp,delta,r2,r3,r_tec

!-----------------------------------------------------------------------------

{^IFMPI if (ipe .ne. 0) read(unitpar,'(a)') atmfilename}

{^IFMPI if (ipe .eq. 0) then}
write(*,*)'param load'
read(unitpar,'(a)') atmfilename
write(*,*) 'from file ',atmfilename
open(33,file=atmfilename,status='old')


!iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)

!iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)

! 1.6Mm

iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)

!iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)

!iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)

{^IFMPI endif}

{^IFMPI call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)}
{^IFMPI call MPI_BCAST(iniene,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)}

do ix_1=ixGhi1,ixGlo1,-1
{^IFMPI if (ipe .eq. 0)} read(33,*) inix,inirho

{^IFMPI if (ipe .eq. 0)} print*,'inix,inirho=',inix,inirho


{^IFMPI call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)}
{^IFMPI call MPI_BCAST(inix,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)}
{^IFMPI call MPI_BCAST(inirho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)}

 do ix_2=ixGlo2,ixGhi2
  do ix_3=ixGlo3,ixGhi3 

   x(ix_1,ix_2,ix_3,1)=inix !*1000.d0
   w(ix_1,ix_2,ix_3,rho_)=inirho
   w(ix_1,ix_2,ix_3,e_)=iniene
   w(ix_1,ix_2,ix_3,m1_)=0.0
   w(ix_1,ix_2,ix_3,m2_)=0.0
   w(ix_1,ix_2,ix_3,m3_)=0.0   

  enddo
 enddo

enddo
{^IFMPI if (ipe .eq. 0)} close(33)

{^IFMPI if (ipe .eq. 0) } print*,'grav=',eqpar(grav0_),eqpar(grav1_),eqpar(grav2_),eqpar(grav3_)



call primitive(ix^L,w)

do ix_3=ixGlo3,ixGhi3
 do ix_2=ixGlo2,ixGhi2
  do ix_1=ixGhi1-1,ixGlo1,-1 

comi=-abs(x(ix_1+1,ix_2,ix_3,1)-x(ix_1,ix_2,ix_3,1))

w(ix_1,ix_2,ix_3,p_)=w(ix_1+1,ix_2,ix_3,p_)+w(ix_1,ix_2,ix_3,rho_)*comi*1.d0*eqpar(grav1_)



  enddo
 enddo
enddo

!goto 200
do ix_3=ixGlo3,ixGhi3
 do ix_2=ixGlo2,ixGhi2
  do ix_1=ixGlo1+2,ixGhi1-2
       
       w(ix_1,ix_2,ix_3,rho_)=-(1.D0/eqpar(grav1_))*(1.D0/(12.D0*(x(ix_1+1,ix_2,ix_3,1)-x(ix_1,ix_2,ix_3,1))))&
                               *(w(ix_1+2,ix_2,ix_3,p_)-8.D0*w(ix_1+1,ix_2,ix_3,p_)+8.D0*w(ix_1-1,ix_2,ix_3,p_)&
			       -w(ix_1-2,ix_2,ix_3,p_))
               


     enddo
   enddo
 enddo  



!lower boundary
do ix_1=ixmin1+4,ixmin1+2,-1
  do ix_2=ixmin2,ixmax2
    do ix_3=ixmin3,ixmax3
        p_1=w(ix_1+2,ix_2,ix_3,p_)-8.d0*w(ix_1+1,ix_2,ix_3,p_)+8.d0*w(ix_1-1,ix_2,ix_3,p_)
        p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
        w(ix_1-2,ix_2,ix_3,p_) = 12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1-1,ix_2,ix_3,1))*p_2+p_1

!         p_1=w(ix_1+2,ix_2,p_)-8.d0*w(ix_1+1,ix_2,p_)+8.d0*w(ix_1-1,ix_2,p_)
!         p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
!         w(ix_1-2,ix_2,p_) = 12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2+p_1

       enddo
    enddo
 enddo


!upper boundary
do ix_1=ixmax1-4,ixmax1-2
   do ix_2=ixmin2,ixmax2
      do ix_3=ixmin3,ixmax3
         
          p_1=w(ix_1-2,ix_2,ix_3,p_)-8.d0*w(ix_1-1,ix_2,ix_3,p_)+8.d0*w(ix_1+1,ix_2,ix_3,p_)
          p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
          w(ix_1+2,ix_2,ix_3,p_) = -12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1-1,ix_2,ix_3,1))*p_2+p_1

!           p_1=w(ix_1-2,ix_2,p_)-8.d0*w(ix_1-1,ix_2,p_)+8.d0*w(ix_1+1,ix_2,p_)
!           p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
!           w(ix_1+2,ix_2,p_) = -12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2+p_1

      enddo
   enddo
enddo


200 continue



call conserve(ix^L,w)

w(ix^S,eb_)=w(ix^S,e_)
w(ix^S,e_)=0.d0

w(ix^S,rhob_)=w(ix^S,rho_)
w(ix^S,rho_)=0.d0


w(ix^S,b1_)=0.d0  !0.4d0



goto 90
! ********************** smooth  tube Bz magnetic field***************
xc2=4.d6  !m
xc3=4.d6  !m

x_max=8.d6 !m

r_0=1.d6
r_tr=0.2d6

b_Amp=0.08d0
delta=400.d0
 
do ix_1=ixGlo1,ixGhi1
 do ix_2=ixGlo2,ixGhi2
  do ix_3=ixGlo3,ixGhi3

    
     
     r2=(x(ix_1,ix_2,ix_3,2)-xc2)
     r3=(x(ix_1,ix_2,ix_3,3)-xc3)     
     
     
     r_tec=sqrt(r2**2.d0+r3**2.d0)         
     
!     if (r_tec .le. r_0+2.0d0*r_tr) w(ix_1,ix_2,ix_3,b1_)=b_Amp
!     if (r_tec .gt. r_0-2.0d0*r_tr) 

     w(ix_1,ix_2,ix_3,b1_)=b_Amp*(one-((atan((r_tec-(r_0+r_tr))/x_max*delta))+Pi/2.d0)/Pi)

  enddo
 enddo
enddo
    
   
! **************************************************************
90 continue


w(ix^S,eb_)=w(ix^S,eb_)-(w(ix^S,b1_)**2.d0+w(ix^S,b2_)**2.d0+w(ix^S,b3_)**2.d0)*&
            (1.d0-eqpar(gamma_)/2.d0)/(eqpar(gamma_)-1.d0)

w(ix^S,bg1_)=w(ix^S,b1_)
w(ix^S,bg2_)=w(ix^S,b2_)
w(ix^S,bg3_)=w(ix^S,b3_)
w(ix^S,b1_)=0.d0
w(ix^S,b2_)=0.d0
w(ix^S,b3_)=0.d0



return
end


!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)


include 'vacdef.f'

integer:: ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
double precision:: fdt,fdthalf2

double precision:: pre(ixG^T),tem(ixG^T),kapr(ixG^T),so(ixG^T),flux(ixG^T)
double precision:: tau(ixG^T),ine(ixG^T)

double precision:: preg(ixG^T),pret(ixG^T)

integer:: rix_1,i,j
double precision:: mol_0, rrr_

double precision:: fsokr,avgflux

integer:: iw,iiw,iix_1

integer:: ix_1,ix_2,ix_3


double precision:: s_period,xz,xc2,xc3,s_z,delta_r,s_r0,vvv
double precision:: xc1Mm,xc2Mm,xc3Mm
double precision:: xx, yy
double precision:: r(ixG^T), vvx(ixG^T), vvy(ixG^T), vvz(ixG^T)
double precision:: bbx(ixG^T), bby(ixG^T)
double precision:: Vphi, bphi, Vr, tt0, ddt, AA, max_vx, max_vy

!*****************
double precision:: t01,t02,a1,a2,s1,s2,sf,rad,rfc,sdep,tdeps,tdep, tdepc, sigma2
double precision:: s_period120, s_period90, s_period60, s_period30
double precision:: norm, am
!-----------------------------------------------------------------------------

eqpar(eta_)=0.d0

!call addsource_diff(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

eqpar(nu_)=1.0d0
!eqpar(nu_)=0.d0

call addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)



if(abs(eqpar(nu_))>smalldouble)&
   call addsource_visc(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
 
xc1Mm=0.1   !Mm        z axis
xc2Mm=1.0 !0.99d0  !Mm        x axis
!xc2Mm=0.9d0 !0.99d0  !Mm        x axis

xc3Mm=1.0 !0.99d0  !Mm        y axis
!xc3Mm=0.9 !0.99d0  !Mm        y axis

!xc2Mm=0.98333333d0  !Mm        x axis
!xc3Mm=0.98333333d0  !Mm        y axis


xz=xc1Mm*1.0d6  !m        z axis
xc2=xc2Mm*1.0d6  !m        x axis
xc3=xc3Mm*1.0d6  !m        y axis



s_z=16000.d0 !m z axis
delta_r=100000.d0   !m
s_r0=50000.d0       !m


AA=200.d0 ! for s_period 120 s
!AA=150.d0 ! for s_period 200 s
!AA=500.d0

!tt0=10.d0
!ddt=5.d0

s_period=120.d0

s_period120=120.d0
s_period90=90.d0
s_period60=60.d0
s_period30=30.d0





!s_period=200.d0

!tdep=sin(qt*2.d0*pi/s_period) 

!tdeps=sin(qt*2.d0*pi/s_period) 
!tdepc=cos(qt*2.d0*pi/s_period) 

am=2.89d0 
norm=1.d0/am ! normirovka

tdeps=norm*sin(qt*2.d0*pi/s_period120)+sin(qt*2.d0*pi/s_period90)+sin(qt*2.d0*pi/s_period60)+sin(qt*2.d0*pi/s_period30) 
tdepc=norm*cos(qt*2.d0*pi/s_period120)+cos(qt*2.d0*pi/s_period90)+cos(qt*2.d0*pi/s_period60)+cos(qt*2.d0*pi/s_period30) 

Vr=0.0d0

do ix_1=ixImin1,ixImax1
 do ix_2=ixImin2,ixImax2
  do ix_3=ixImin3,ixImax3

     
      xx=x(ix_1,ix_2,ix_3,2)-xc2
      yy=x(ix_1,ix_2,ix_3,3)-xc3
     
      r(ix_1,ix_2,ix_3)=sqrt(xx**2.d0+yy**2.d0)
      
      
!      Vphi=AA*tdep*exp(-(r(ix_1,ix_2,ix_3)-s_r0)**2.d0/delta_r**2.d0)     
      
        vvx(ix_1,ix_2,ix_3)=AA*tdeps*exp(-(r(ix_1,ix_2,ix_3)-s_r0)**2.d0/(delta_r**2.d0))     
        vvy(ix_1,ix_2,ix_3)=-AA*tdepc*exp(-(r(ix_1,ix_2,ix_3)-s_r0)**2.d0/(delta_r**2.d0))     

      
      
!      if (r(ix_1,ix_2,ix_3) .ne. 0.d0) then
	 
!        vvx(ix_1,ix_2,ix_3)=(Vr*xx-Vphi*yy)/r(ix_1,ix_2,ix_3)
!        vvy(ix_1,ix_2,ix_3)=(Vr*yy+Vphi*xx)/r(ix_1,ix_2,ix_3)
	

!      else

!        vvx(ix_1,ix_2,ix_3)=0.0d0
!        vvy(ix_1,ix_2,ix_3)=0.0d0

 
!     endif

   enddo
 enddo
enddo



!tdep=exp(-(qt-tt0)**2.d0/ddt**2.d0)



do ix_1=ixImin1,ixImax1
 do ix_2=ixImin2,ixImax2
  do ix_3=ixImin3,ixImax3


 vvv=exp(-(x(ix_1,ix_2,ix_3,1)-xz)**2.d0/s_z**2.d0)

  
 w(ix_1,ix_2,ix_3,m2_)=w(ix_1,ix_2,ix_3,m2_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*&
		      vvx(ix_1,ix_2,ix_3)*vvv*qdt

 w(ix_1,ix_2,ix_3,m3_)= w(ix_1,ix_2,ix_3,m3_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*&
		      vvy(ix_1,ix_2,ix_3)*vvv*qdt



! w(ix_1,ix_2,ix_3,b2_)=w(ix_1,ix_2,ix_3,b2_)+bbx(ix_1,ix_2,ix_3)*vvv*qdt

! w(ix_1,ix_2,ix_3,b3_)=w(ix_1,ix_2,ix_3,b3_)+bby(ix_1,ix_2,ix_3)*vvv*qdt


	      
 w(ix_1,ix_2,ix_3,e_)=w(ix_1,ix_2,ix_3,e_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*&
                      ((vvx(ix_1,ix_2,ix_3)*vvv)**2.d0+(vvy(ix_1,ix_2,ix_3)*vvv)**2.d0)*qdt/2.d0
		      
 
 !Vx[i,j]=Vx[i,j]*G/Vx_m
 !Vy[i,j]=Vy[i,j]*G/Vy_m

   enddo
 enddo
enddo


  
;max_vx=maxval(vvx)
;max_vy=maxval(vvy)
     
     
     
!     tdep=sin(qt*2.d0*pi/s_period)
!
!
!!     vvv=AA*tdep*dexp(-(r1/s_rad1)**2.d0)*dexp(-(r2/s_rad2)**2.d0)*dexp(-(r3/s_rad3)**2.d0)     
!
!     w(ix_1,ix_2,ix_3,e_)=w(ix_1,ix_2,ix_3,e_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*(vvv**2.d0)*qdt/2.d0
!     w(ix_1,ix_2,ix_3,m1_)=w(ix_1,ix_2,ix_3,m1_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*vvv*qdt
!
!  enddo
! enddo
!enddo

{^IFMPI if (ipe.eq.0)} write(*,*) '***time=',qt

end



!=============================================================================
subroutine specialbound(qt,ix^L,iw,iB,w)
include 'vacdef.f'

integer:: ix_1,ix_2

integer:: iw^LIM,idim^LIM
double precision:: qt,w(ixG^T,1:nw)
integer:: ix,ix^D,ixe,ixf,ix^L,ixpair^L,idim,iw,iB
integer:: iwv,jdim

integer:: Ns,i,j
double precision:: ki


double precision:: tmpp1,tmpp2

Ns=1

select case(iB)

case(1)

case(2)


case default
stop 'error iB!'
end select

return
end

!=============================================================================
subroutine getdt_special(w,ix^L)

! If the Coriolis force is made very strong it may require time step limiting,
! but this is not implemented here.

include 'vacdef.f'
double precision:: w(ixG^T,nw)
integer:: ix^L
!-----------------------------------------------------------------------------

!call getdt_diff(w,ix^L)

if(abs(eqpar(nu_))>smalldouble)&
   call getdt_visc(w,ix^L)


call getdt_grav(w,ix^L)

return
end


subroutine specialeta(w,ix^L,idirmin)
 
include 'vacdef.f'
 
double precision:: w(ixG^T,nw)
integer:: ix^L,idirmin
!-----------------------------------------------------------------------------
 
stop 'specialeta is not defined'
end

