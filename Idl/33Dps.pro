

;openr,1,'/data/ap1vf/3D_cont3.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf
cname='4b0_b4'
directory='/fastdata/cs1mkg/smaug/spic'+cname+'_3d/'
;while not(eof(1)) do begin
pic=100
;for ipic=pic,pic do begin
for ipic=1,4 do begin
headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)
name='zerospic1__'
nn=0

close,1








;picid=((90*ipic)*1000L)
picid=((89+15*ipic)*1000L)
;picid=ipic*1000L
;picid=ipic*5+4
;picid=ipic;
outfile=directory+name+strtrim(string(picid),2)+'.out'
print,'ipic=',ipic
;openr,1,outfile,/f77_unf
openr,1,outfile





readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


n1=nx(0)
n2=nx(1)
n3=nx(2)

x=dblarr(n1,n2,n3,ndim)
wi=dblarr(n1,n2,n3)
w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

indexss=strtrim(picid,2)

;a = strlen(indexs)                                                  
;case a of                                                           
; 1:indexss='c300000000'+indexs                                             
; 2:indexss='c30000000'+indexs                                              
; 3:indexss='c3000000'+indexs                                               
; 4:indexss='c300000'+indexs
; 5:indexss='c30000'+indexs                                             
; 6:indexss='c3000'+indexs                                              
; 7:indexss='c300'+indexs                                               
; 8:indexss='c30'+indexs                                                  
;endcase 

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667
time=picid*0.001

T=dblarr(n1,n2,n3)

T[*,*,*]=(w[*,*,*,4]+w[*,*,*,8])
T[*,*,*]=T[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/(w[*,*,*,0]+w[*,*,*,9])/2.0
T[*,*,*]=T[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0+(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.d0

T[*,*,*]=(gamma-1.d0)*T[*,*,*]


C=dblarr(n1,n2,n3)
C[*,*,*]=sqrt(gamma*T[*,*,*]/(reform(w[*,*,*,0]+w[*,*,*,9])))
	 
;data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9))/C(*,*,*))	 
data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9)))

 hData = PTR_NEW(data, /NO_COPY) 
 ;a='Mx/'+indexss  
 a=cname+'_'+indexss                                              
 Slicer_ps, hdata, time, a


;for i=2,3 do begin
; data=reform(w(*,*,*,i)/(w(*,*,*,0)+w(*,*,*,9)))
; hData = PTR_NEW(data, /NO_COPY) 

;case i of                                                           
; 1:a='Vz/'+indexss                                          
; 2:a='Vx/'+indexss                                              
; 3:a='Vy/'+indexss                                               ;

;endcase

; Slicer_ps, hdata, time, a
 
;endfor

nn=nn+1
print, 'time - - - ', time
;end loop over each input file
;endwhile
endfor 
end
