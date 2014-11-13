%plot sections through 3d array
   %slice=48;
   x=linspace(0,4,128);
   y=linspace(0,4,128);
   z=linspace(0,6,128);
   
   
   
   
   nrange=3:126;
   
   ax=x(nrange);
   ay=y(nrange);
   az=z(nrange);
   [x1,x2,x3] = meshgrid(ax,ay,az);
   val1=reshape(wd(2,nrange,nrange,nrange),124,124,124);
   val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);
   myval=shiftdim(val1./val2,1);
   
   %P = [2 1 3];
   %x1 = permute(x1, P);
   %x2 = permute(x2, P);
   %x3 = permute(x3, P);
   %myval = permute(myval, P);
   
   
  %h= slice(myval,64, 64, 4);
  %h=figure();
  %hold on;
  h=slice(myval,80, 64,8);
  hold on;
  set(h,'EdgeColor','none','FaceColor','interp');
  %set(h,'XData',ax,'YData',ay,'ZData',az);
  hax=get(h,'Children');
  set(gca,'CameraPosition',[-606.298 -914.02 280.537]);
  set(gca,'Xlim',[0 124],'Ylim',[0 124],'Zlim',[0 124]);
  
  set(gca,'XTickLabel',{'0';'1.6';'3.2'})
  set(gca,'YTickLabel',{'0';'1.6';'3.2'})
  set(gca,'ZTickLabel',{'0.09';'0.99';'1.94';'2.88';'3.83';'4.77';'5.72'})
  cmap=colormap('Jet');
  caxis([-4000 5000]);
  colormap(cmap);
  hc=colorbar();
  set(hc,'Ylim',[-4000 5000]);
  
  
  