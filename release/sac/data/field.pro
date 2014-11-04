!P.multi=[0,1,1]
nptx=100
nptz=100

bz=dblarr(nptx,nptz)
bx=dblarr(nptx,nptz)
bt=dblarr(2,nptx,nptz)

close,10
openr,10,'bbb_gen.dat',/f77_unformatted
   readu,10,bz,bx
close,10

  bt(0,*,*)=bx(*,*)
  bt(1,*,*)=bz(*,*) 
;tvframe, reform(bt(0,*,*)), /bar
;tvframe, reform(bt(1,*,*)), /bar
vector_field,bt, aa
;plot, [aa(0,*),aa(1,*)], psym=10
velovect, bx,bz, length=1
end
 
