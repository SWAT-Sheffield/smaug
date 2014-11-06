;programma dlia postroenija magnitnich silovich liniy
;
;02.10.2006 Fedun Viktor Appl Math Sheff univ
;stroka 24 - chtenie fila coordinat sozdannogo  vacini (VAC)
;stroka 89 - chtanie fila  so znachenijami magnitnogo polia Bx,By,Bz 
;                                 openr,1,'B_t_3000_1000_32_252.dat',/f77_unf
;stroka 169 z0_min=120000.0 uroven' na kotorom vichisliajetsia 
;raspredilenie magnitnogi polia na ploskosti
;stroka 226  Kolichestvo urovney na kotorich vichisliaetsia magnitnoe pole nn=10.0
;stroka 253 iik_max maximalnoe kolichestvo tochek na urovne gde 
;funktsija imeet maximum iik_max=10

function B_x,x_i,y_i,z_i
 RETURN,y_i
end  
function B_y,x_i,y_i,z_i
 RETURN,x_i
end  






;**********  begin data from file
close,1

openr,1,'/data/ap1vf/gorz300gmf.out',/f77_unf
for kkk=1,10 do begin


tarr=dblarr(1)
headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)
nn=0


readu,1,headline
print,'head = ', headline
readu,1,it,time,ndim,neqpar,nw
print,'it, time', it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
print, '*',ndim 
nx=lonarr(ndim)
readu,1,nx
print, 'nx = ', nx
eqpar=dblarr(neqpar)
readu,1,eqpar
print,eqpar
readu,1,varname
print,varname

n1=nx(0) 
n2=nx(1)

xx=dblarr(n1,n2,ndim)

if (nn eq 0) then w=dblarr(n2,n1,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2)
;e2=dblarr(n1,n2)
readu,1,xx
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=rotate(wi,1)
endfor

kkk=kkk+1


endfor
close,1

bbx=dblarr(n1,n2)
bbz=dblarr(n1,n2)

bbx(*,*)=reform(w(*,*,5))
bbz(*,*)=reform(w(*,*,4))





DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW


wx=400
wy=wx

window,1, xsize=wx,ysize=wy
window,2, XSIZE=600, YSIZE=600
window,3, xsize=wx,ysize=wy


wset,1

loadct, 5
!p.multi=[0,1,0,0,1]


x=dblarr(n1)
z=dblarr(n2)


;****************end data from file 




for i=0,n1-1 do begin
    x(i)=xx(i,0,0)
endfor

minimum_z=max(xx(0,0,1))
maximum_z=6.5e6

for k=0,n2-1 do begin
    z(k)=xx(0,k,1)
endfor
    
; ********** rassmatrivaem ne vsu trubtu a do urovnia maximum_z ********


  if aa eq 1 then begin
     k=0
     while (maximum_z-xx(0,k,1))*(maximum_z-xx(0,k+1,1)) gt 0.0 do begin         
     k=k+1
     endwhile
    n2=fix(k+(maximum_z-xx(0,k,1))/(xx(0,k+1,1)-xx(0,k,1)))+1
  endif
; ********** rassmatrivaem ne vsu trubtu a do urovnia maximum_z ********



bx=dblarr(n1,n2)
bz=dblarr(n1,n2)

  for i=0,n1-1 do begin
   for j=0,n2-1 do begin
      bx(i,j)=bbx(i,j)
      bz(i,j)=bbz(i,j)
   endfor
  endfor

xyz_color=dblarr(1)

mi_x=min(x)
ma_x=max(x)


mi_z=minimum_z
ma_z=maximum_z

x0= 0.0 ;mi_x

;**************************** NACHALNIY UROVEN' DLIA VICHISLENIJA POLIA
;z0_min=10.0 ;1000000.0
z0_min=35000.0
z0=z0_min ;mi_z

     k=0
     while (z0-z(k))*(z0-z(k+1)) gt 0.0 do begin         
     k=k+1
;    print,'k',k , z0,z(k),z(k+1),kk,(z0-z(k))*(z0-z(k+1))     
     endwhile
    kk=k+(z0-z(k))/(z(k+1)-z(k))


b_total=sqrt(bx(*,*)^2.0+bz(*,*)^2.0)
k=0


max_bx=max(bx)
min_bx=min(bx)
max_bz=max(bz)
min_bz=min(bz)

print,'bx bz =', min_bx,max_bx,min_bz,max_bz

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL 

xm=dblarr(1)
zm=dblarr(1)


print,'mi_x, ma_x, mi_y,ma_y,mi_z,ma_z = ', mi_x, ma_x,mi_z,ma_z




xm(0)=x0
zm(0)=z0

xyz_color(0)=1.0

color_min=60.0
color_max=256.0


dLx=4000.0; 0.6
dLz=40001.0; 0.7


df=1.0
f=1.0

;************* Kolichestvo urovney na kotorich vichisliaetsia magnitnoe pole *********************
nn=10.0

max_b=max(b_total)
min_b=min(b_total)

b_total_all=sqrt(bx(*,*)^2.0+bz(*,*)^2.0)

max_b_all=max(b_total_all)
min_b_all=min(b_total_all)

db=max_b-min_b
print,'db=',db
print, 'max_b_all, min_b_all', max_b_all, min_b_all

delta_b=db/nn


sq_b_real=dblarr(nn)

;** massivu x i y so znachenijami pixelov of windows
xr=dblarr(1)
;**
;** massivu x i y so znachenijami indexov sootvetstvuyuschich tochek
xxr=dblarr(1)
;**

;** iik - obschee kolichestvo tochek po vsey ploschadke
iik=0

;******************************************************************************
;** iik_max maximalnoe kolichestvo tochek na urovne gde funktsija imeet maximum
iik_max=20

bool_f=1

lv=dblarr(nn)

for j=0,nn-1 do begin
 lv(j)=max_b-(j+1)*delta_b

seconds=systime(1)

contour, b_total, LEVELS = lv(j),PATH_INFO=PInfo,PATH_XY=Pxy,/PATH_DATA_COORDS

!Y.STYLE = 5
!X.STYLE = 5
contour, b_total, LEVELS = lv(j), POS=[0., 0., 1.0, 1.0],/fill

print,n_elements(T)

nxy=(N_elements(Pxy))/2.0

xy=dblarr(2,nxy)

for i=0,N_elements(Pxy)/2.0-1 do begin
    xy(0,i)= Pxy(0,i)
    xy(1,i)= Pxy(1,i)    
endfor 


sq_b = POLY_AREA(xy(0,*), xy(1,*), /double)  
PRINT, 'area = ', sq_b, round(sq_b), '  square pixels' , j

; raznost' ploschadey esli ne perviy raz
if j gt 1 then begin
  koef_sq=(sq_b-sq_b_real(j-1))/sq_b_real(1)  
  sq_b_real(j)=sq_b
  koef_f=lv(j)/lv(1)
  T=tvrd()-T_pred
  T_pred=tvrd()
endif else begin
  sq_b_real(j)=sq_b
  koef_sq=1.0
  koef_f=1.0
  T=tvrd()
  T_pred=T
endelse
  

print,'koef sq f =',koef_sq,koef_f,sq_b_real(j),j
;print,T

n_sq=10000

rm=dblarr(2,n_sq)
rm=round(randomu(seconds,2,n_sq)*(wx-1))
i=0
iik_tec=0
end_koef=koef_f*koef_sq*(iik_max-1)
print,'end_koef= ',end_koef
print, 'n elements, xxr,yyr', n_elements(xxr),n_elements(xxr)


;read,aa
endfor


wset,3
lv_r=dblarr(nn)
for j=0,nn-1 do begin
  lv_r(j)=lv(nn-1-j)
endfor
!Y.STYLE = 5
!X.STYLE = 5

plot, b_total
read,aa
contour, b_total, LEVELS = lv_r,C_ANNOTATION=string(lv_r,format="(f6.2)"), POS=[0., 0., 1.0, 1.0]

 
;************** massivu dlia tochek v upper ******************
x_u=dblarr(n_elements(xxr))
z_u=dblarr(n_elements(xxr))
;*****************************************************************

;************** massivu dlia tochek v middle footpoint ******************
x_m=dblarr(n_elements(xxr))
z_m=dblarr(n_elements(xxr))
;*****************************************************************

;************** massivu dlia tochek v footpoint ******************
x_f=dblarr(n_elements(xxr))
z_f=dblarr(n_elements(xxr))
;*****************************************************************


wset,2
    

surface, DIST(5), /NODATA, /SAVE, XRANGE=[mi_x, ma_x], $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z], XSTYLE=1, $
   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.5, xtitle='x',ytitle='y',title='Field Lines', $
   ztitle='z', POSITION=[0.15, 0.2, 0.85, 0.85], az=30.0, ax=30.0



	 
print,'!X.Crange[1], !Y.Crange[1], !Z.Crange[0]',!X.Crange, !Y.Crange, !Z.Crange	 
;stop

ijk=0

while ijk le n_elements(xxr)-1 do begin 

ii=xxr(ijk)
jj=yyr(ijk)

z0=z0_min

    i=0
    j=0
   while i lt ii do begin
     i=i+1
   endwhile
     x0=x(i-1)+(ii-(i-1))*(x(i)-x(i-1))   
     x_m(ijk)=x0
     z_m(ijk)=z0
print, '****ijk=',ijk,x0,y0,z0

ijk_start=1
     
while (x0 ge mi_x and x0 le ma_x and z0 ge mi_z and z0 le ma_z) do begin 


    i=0
    j=0
    k=0
 if ijk_start ne 1 then begin        
     while (x0-x(i))*(x0-x(i+1)) gt 0.0 do begin     
      i=i+1
     endwhile    
    ii=i+(x0-x(i))/(x(i+1)-x(i))
;    print,'i',i, x0,x(i),x(i+1),ii, (x0-x(i))*(x0-x(i+1))
 
     while (y0-y(j))*(y0-y(j+1)) gt 0.0 do begin         
     j=j+1
     endwhile
    jj=j+(y0-y(j))/(y(j+1)-y(j))
   ; print,'j',j , y0,y(j),y(j+1),jj,(y0-y(j))*(y0-y(j+1))     
 endif
 
 ijk_start=2
 
     while (z0-z(k))*(z0-z(k+1)) gt 0.0 do begin         
     k=k+1
     endwhile
    kk=k+(z0-z(k))/(z(k+1)-z(k))
 ;   print,'k',k , z0,z(k),z(k+1),kk, ma_z, z(n3),n3



    Bx_inter=INTERPOLATE(bx,ii,jj,kk)
    By_inter=INTERPOLATE(by,ii,jj,kk)
    Bz_inter=INTERPOLATE(bz,ii,jj,kk)
;print,'Bx,By,Bz', Bx_inter,Bz_inter,Bz_inter
   sq=sqrt(Bx_inter^2.0+By_inter^2.0+Bz_inter^2.0)

cl=color_max-(max_b_all-sq)/(max_b_all-min_b_all)*(color_max-color_min)
   
   xyz_color=[xyz_color,cl]
   
   dx=-df*dLx*Bx_inter/sq
   dy=-df*dLy*By_inter/sq
   dz=-df*dLz*Bz_inter/sq   

;print, dx,dy,dz
   x0=x0+dx
   y0=y0+dy
   z0=z0+dz

;print,'f= ',f
 if f eq 1.0 then begin
   x0_pre=x0
   y0_pre=y0
   z0_pre=z0   
   f=-1.0
 endif 
   
   xm=[xm,x0]
   ym=[ym,y0]
   zm=[zm,z0]
;  print,xm,ym,zm
   
;PLOTS, xm,ym,zm, color=xyz_color, psym=3, /T3D
endwhile
 if df eq -1.0 then begin
   f=1.0
 endif
  
;   print,'######*****',x0_s,y0_s
     
 if f eq -1.0 then begin
  df=-1.0
;    print,'*** ijk', ijk,xm(n_elements(xm)-1),ym(n_elements(ym)-1),zm(n_elements(zm)-1) 
;    print,xm(n_elements(xm)),ym(n_elements(ym)),zm(n_elements(zm))     
    x_f(ijk)=xm(n_elements(xm)-1)
    y_f(ijk)=ym(n_elements(ym)-1)
    z_f(ijk)=zm(n_elements(zm)-1)
;    read,aa  
    x0=x0_pre
    y0=y0_pre
    z0=z0_pre  
    
 endif else begin 
    x_u(ijk)=xm(n_elements(xm)-1)
    y_u(ijk)=ym(n_elements(ym)-1)
    z_u(ijk)=zm(n_elements(zm)-1)
  df=1.0 
  ijk=ijk+1
 endelse
endwhile



PLOTS, xm,ym,zm, color=xyz_color, psym=3, /T3D 
PLOTS, x_u,y_u,z_u, color=150, psym=8, /T3D 

;PLOTS, x_m,y_m,z_m, color=120, psym=8, /T3D 
;PLOTS, x_f,y_f,z_f, color=200, psym=8, /T3D 

;********** begin additional axis

;   Axis, XAXIS=-1, !X.Crange[0], !Y.Crange[0], !Z.Crange[1],  $ 
;         XTICKNAME = name, COLOR=COLOR, XTICKLEN=0

;   Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[0], ZAXIS=-1, /T3D, $
;         ZTICKNAME = name, COLOR=COLOR, ZTICKLEN=0, /nodata      
;

;   Axis, !X.Crange[1], !Y.Crange[1], !Z.Crange[0], ZAXIS=-1, /T3D,  $ 
;         ZTICKNAME = name, COLOR=COLOR, ZTICKLEN=0, /nodata     

;   Axis, !X.Crange[0], !Y.Crange[0], !Z.Crange[0], ZAXIS=-1, /T3D,  $ 
;         ZTICKNAME = name, COLOR=COLOR, ZTICKLEN=0, /nodata     

;********** end additional axis

 xyouts,50,5, file_name, /device, color=150

;************ begin bar window
   px = !x.window * !d.x_vsize     ;Position of frame in device units
   py = !y.window * !d.y_vsize
   
   sx = px(1)-px(0)                ;Size of frame in device units
   sy = py(1)-py(0)
   sx = sx/1.25

   bx    = fix(px(0)+sx*1.3)

   by    = fix(1.6*py(0))
   bsx   = fix(sx*0.04)
   bsy   = fix(sy/2.0)
   barpos= [bx,by,bx+bsx,by+bsy]


   barim=findgen(1,bsy)/(bsy-1)*(max_b_all-min_b_all)+min_b_all


   
   barim=rebin(barim,bsx,bsy,/sample)
   
    bs=bytscl(barim, top=195)
    bs=bs+color_min
   
   print, barim
   
   tv, bs,bx,by,/device

   BRANGE=float([min_b_all,max_b_all])

   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1, ytitle=BTITLE, $
   TICKLEN=-.15,CHARSIZE=CHARSIZE,YTICKS=BTICKS,YMINOR=BMINOR
   
;************ end bar window

 ss_bx_min='bx min ='+strTrim(string(min_bx,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(4.6*py(0)), ss_bx_min, /device, color=150

 ss_bx_max='bx max ='+strTrim(string(max_bx,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(4.5*py(0)), ss_bx_max, /device, color=150

 ss_by_min='by min ='+strTrim(string(min_by,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(4.4*py(0)), ss_by_min, /device, color=150

 ss_by_max='by max ='+strTrim(string(max_by,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(4.3*py(0)), ss_by_max, /device, color=150

 ss_bz_min='bz min ='+strTrim(string(min_bz,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(4.2*py(0)), ss_bz_min, /device, color=150

 ss_bz_max='bz max ='+strTrim(string(max_bz,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(4.1*py(0)), ss_bz_max, /device, color=150
 
 ss_b_min='b min ='+strTrim(string(min_b_all,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(3.9*py(0)), ss_b_min, /device, color=120

 ss_b_max='b max ='+strTrim(string(max_b_all,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(3.8*py(0)), ss_b_max, /device, color=120

 ss_z0_b_calculation='z0 b ='+strTrim(string(z0_min,format="(f9.1)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(3.6*py(0)), ss_z0_b_calculation, /device, color=255
 
 ss_z_max='z max ='+strTrim(string(maximum_z,format="(f9.1)"),2)
 xyouts,fix(px(0)+sx*1.2),fix(3.5*py(0)), ss_z_max, /device, color=255
 
;************ vivod na PS ****************
a_sm=100	
x_shift=2.2
 
 
set_plot, 'ps'

device, filename=file_name+'.ps',/Color
surface, DIST(5), /NODATA, /SAVE, XRANGE=[mi_x, ma_x], $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z], XSTYLE=1, $
   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.5, xtitle='x',ytitle='y',title='Field Lines', $
   ztitle='z', POSITION=[0.15, 0.2, 0.85, 0.85], az=30.0, ax=30.0  

PLOTS, xm,ym,zm, color=xyz_color, psym=3, /T3D 
PLOTS, x_u,y_u,z_u, color=150, psym=8, /T3D 

;PLOTS, x_m,y_m,z_m, color=120, psym=8, /T3D 
;PLOTS, x_f,y_f,z_f, color=200, psym=8, /T3D 


;************ begin bar PS
   px = !x.window * !d.x_vsize     ;Position of frame in device units
   py = !y.window * !d.y_vsize
   print, 'px,py',px,py
   sx = px(1)-px(0)                ;Size of frame in device units
   sy = py(1)-py(0)
   sx = sx/1.25

   bx    = fix(px(0)+sx*1.3)

   by    = fix(1.0*py(0))
   bsx   = fix(sx*0.04)
   bsy   = fix(sy/2.0)
   barpos= [bx,by,bx+bsx,by+bsy]



    barim=findgen(1,bsy)/(bsy-1)*(max_b_all-min_b_all)+min_b_all

    bs=bytscl(barim, top=195)
    bs=bs+color_min
   
     n_b=n_elements(bs)   
      aaa=dblarr(2,n_b)
       for j=0, 1 do begin
        for i=0, n_b-1 do begin
         aaa(j,i)=bs(i)
        endfor
       endfor   
       
   tv, aaa,bx,by, ysize=bsy, xsize=bsx,/device

   BRANGE=float([min_b_all,max_b_all])

   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device ;,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1, ytitle=BTITLE, $
   TICKLEN=-.15,CHARSIZE=0.7,YTICKS=BTICKS,YMINOR=BMINOR
   
;************ end bar PS

 ss_bx_min='bx min ='+strTrim(string(min_bx,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(4.6*py(0)), ss_bx_min, /device, charsize=0.7

 ss_bx_max='bx max ='+strTrim(string(max_bx,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(4.45*py(0)), ss_bx_max, /device, charsize=0.7

 ss_by_min='by min ='+strTrim(string(min_by,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(4.3*py(0)), ss_by_min, /device, charsize=0.7

 ss_by_max='by max ='+strTrim(string(max_by,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(4.15*py(0)), ss_by_max, /device, charsize=0.7

 ss_bz_min='bz min ='+strTrim(string(min_bz,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(4.0*py(0)), ss_bz_min, /device,  charsize=0.7

 ss_bz_max='bz max ='+strTrim(string(max_bz,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(3.85*py(0)), ss_bz_max, /device, charsize=0.7
 
 ss_b_min='b min ='+strTrim(string(min_b_all,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(3.60*py(0)), ss_b_min, /device, color=120, charsize=0.7

 ss_b_max='b max ='+strTrim(string(max_b_all,format="(f7.2)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(3.45*py(0)), ss_b_max, /device, color=120, charsize=0.7

 ss_z0_b_calculation='z0 b ='+strTrim(string(z0_min,format="(f9.1)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(3.15*py(0)), ss_z0_b_calculation, /device, charsize=0.7

 ss_z_max='z max ='+strTrim(string(maximum_z,format="(f9.1)"),2)
 xyouts,fix(px(0)+sx*1.25),fix(2.95*py(0)), ss_z_max, /device, charsize=0.7


; xyouts,a_sm*x_Shift,a_sm, file_name, /device,  charsize=0.7

device, /close
set_plot, 'x'

 
end
