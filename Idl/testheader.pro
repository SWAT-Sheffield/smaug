filename='../data/shearalfven2d_0.out'
nfile=0
npictinfile=1
str2arr,filename,filenames,nfile
gettype,filenames,filetypes,npictinfiles
openfile,10,filename,filetypes(0)

getpict,10,filetypes(0),npict,x,w,headline,phys,it,time,$
          gencoord,ndim,neqpar,nw,nx,eqpar,variables,rBody,error


close,10

end
