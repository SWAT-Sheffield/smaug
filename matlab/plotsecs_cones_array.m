   
directory='/storage2/mikeg/results/amp10k/spic4b0_1_3d/';
extension='.out';

ndirectory='/storage2/mikeg/results/amp10k/spic4b0_1_3d/imagesscones/';
nextension='.jpg';


for i=0:642
%for i=0:5:30    

id=int2str(1000*i);
filename=[directory,'zerospic1__',id,extension];
imfile=[ndirectory,'im1_',id,nextension];
disp([id filename]);
   fid=fopen(trim(filename));
   %fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   headline=trim(setstr(fread(fid,79,'char')'));
   it=fread(fid,1,'integer*4'); time=fread(fid,1,'float64');
 
   ndim=fread(fid,1,'integer*4');
   neqpar=fread(fid,1,'integer*4'); 
   nw=fread(fid,1,'integer*4');
   nx=fread(fid,3,'integer*4');
   
   nxs=nx(1)*nx(2)*nx(3);
   varbuf=fread(fid,7,'float64');
   
   gamma=varbuf(1);
   eta=varbuf(2);
   g(1)=varbuf(3);
   g(2)=varbuf(4);
   g(3)=varbuf(5);
   
   
   varnames=trim(setstr(fread(fid,79,'char')'));
   
   for idim=1:ndim
      X(:,idim)=fread(fid,nxs,'float64');
   end
   
   for iw=1:nw
      %fread(fid,4);
      w(:,iw)=fread(fid,nxs,'float64');
      %fread(fid,4);
   end
   
   nx1=nx(1);
   nx2=nx(2);
   nx3=nx(3);
   
   xx=reshape(X(:,1),nx1,nx2,nx3);
   yy=reshape(X(:,2),nx1,nx2,nx3);
   zz=reshape(X(:,3),nx1,nx2,nx3);
   
   
 
  % extract variables from w into variables named after the strings in wnames
wd=zeros(nw,nx1,nx2,nx3);
for iw=1:nw
  
     tmp=reshape(w(:,iw),nx1,nx2,nx3);
     wd(iw,:,:,:)=tmp;
end


%w=tmp(iw);
  

clear tmp; 
   
   
   fclose(fid);




%plot sections through 3d array
   %slice=48;
   x=linspace(0,4,128);
   y=linspace(0,4,128);
   z=linspace(0,6,128);
   
   xmin=3;
   xmax=126;
   ymin=3;
   ymax=126;  
    
   
   nrange=3:126;
   
   ax=x(nrange);
   ay=y(nrange);
   az=z(nrange);

      xsec=1:4:128;
   ysec=1:4:128;
   nx=32;
   ny=32;
   
 
   
   
   
    % Decide where in data space you want to plot cones. This example selects the 
    % full range of x and y in eight steps and the range 3 to 15 in four steps in z 
    % using linspace and meshgrid.
    nxrange = linspace(xmin,xmax,8);
    nyrange = linspace(ymin,ymax,8);
    nzrange = 3:4:126;
    [cx cy cz] = meshgrid(nxrange,nyrange,nzrange);
    
    

    









   %plot sections through 3d array
   %slice=48;
   x=linspace(0,4,128);
   y=linspace(0,4,128);
   z=linspace(0,6,128);
   
   
   
   
   nrange=3:126;
   
   ax=x(nrange);
   ay=y(nrange);
   az=z(nrange);
   [x1,x2,x3] = meshgrid(nrange,nrange,nrange);
   
   val1=reshape(wd(4,nrange,nrange,nrange),124,124,124);
   val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);
   uu=shiftdim(val1./val2,1);

   
   
   val1=reshape(wd(3,nrange,nrange,nrange),124,124,124);
   val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);
   uv=shiftdim(val1./val2,1);

   
   
   val1=reshape(wd(2,nrange,nrange,nrange),124,124,124);
   val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);
   uw=shiftdim(val1./val2,1);
   
   
   % Draw the cones, setting the scale factor to 5 to make the cones larger than 
    % the default size:
    hcones = coneplot(x1,x2,x3,uu,uv,uw,cx,cy,cz,5);
    % Set the coloring of each cone using FaceColor and EdgeColor:
    set(hcones,'FaceColor','red','EdgeColor','none')

    % Calculate the magnitude of the vector field (which represents wind speed) to 
    % generate scalar data for the slice command:
    hold on
   
   
   %P = [2 1 3];
   %x1 = permute(x1, P);
   %x2 = permute(x2, P);
   %x3 = permute(x3, P);
   %myval = permute(myval, P);
   
   
  %h= slice(myval,64, 64, 4);
  %h=figure();
  %hold on;
      h=slice(uw,80, 64,8);
      hold on;
      set(h,'EdgeColor','none','FaceColor','interp');
      %set(h,'XData',ax,'YData',ay,'ZData',az);
      hax=get(h,'Children');
      set(gca,'CameraPosition',[491.298 -661.02 631.537]);
      set(gca,'Xlim',[0 124],'Ylim',[0 124],'Zlim',[0 124]);

      set(gca,'XTickLabel',{'0';'1.6';'3.2'})
      set(gca,'YTickLabel',{'0';'1.6';'3.2'})
      set(gca,'ZTickLabel',{'0.09';'0.99';'1.94';'2.88';'3.83';'4.77';'5.72'})
      cmap=colormap('Jet');
      caxis([-4000 5000]);
      colormap(cmap);
      hc=colorbar();
      set(hc,'Ylim',[-4000 5000]);
      
       
  print('-djpeg', imfile);
  
  hold off
      
  
end
  