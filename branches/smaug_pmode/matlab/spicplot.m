   slice=48;

   val1=reshape(wd(2,:,slice,:),128,128);
   val2=reshape(wd(1,:,slice,:)+wd(10,:,1,:),128,128);
   surf(val1./val2,'LineStyle','none');
   colorbar;
   
   xsec=1:4:128;
   ysec=1:4:128;
   nx=32;
   ny=32;
   
   val1=reshape(wd(3,xsec,slice,ysec),nx,ny)/(reshape(wd(1,xsec,slice,ysec),nx,ny)+reshape(wd(10,xsec,slice,ysec),nx,ny));
   val2=reshape(wd(4,xsec,slice,ysec),nx,ny)/(reshape(wd(1,xsec,slice,ysec),nx,ny)+reshape(wd(10,xsec,slice,ysec),nx,ny));
   
   %mv1=max(max(val1));
   %mv2=max(max(val2));
   
   vmod=sqrt(val1.*val1+val2.*val2);
   mv=max(max(vmod));
   
   nval1=val1./mv;
   nval2=val2./mv;
   [x,y]=meshgrid(1:nx,1:ny);
   
   quiver(x,y,nval1,nval2);