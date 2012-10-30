it=long(1)
time=double(1)

xxs=double(1)
xxe=double(1)
zzs=double(1)
zze=double(1)
ww=dblarr(100,186)


close,1

openr,1,'vz_y50.10000',/f77_unf

readu,1,it,time,xxs,xxe,zzs,zze
readu,1,ww
close, 10

end
