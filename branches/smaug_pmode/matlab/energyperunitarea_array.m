
%directory='/storage2/mikeg/results/spic6b0_1_3d/';
%directory='/storage2/mikeg/results/spic5b0_b1G_3d/';
directory='/storage2/mikeg/results/spic6b0_2_3d/';
%directory='/storage2/mikeg/results/spic4b0_2_3d/';
%directory='/storage2/mikeg/results/spic3p0a_0_2_3d/';
%directory='/storage2/mikeg/results/spic2p3a_0_3_3d/';
%directory='/storage2/mikeg/results/spic6p7a_0_0_3d/';
%directory='/storage2/mikeg/results/spic4p3a_0_1_3d/';
%directory='/storage2/mikeg/results/spic4p3_0_1_3d/';
extension='.out';

%ndirectory='/storage2/mikeg/results/spic5b0_b1G_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic6b0_1_3d/images_3d_vsecs/';
ndirectory='/storage2/mikeg/results/spic6b0_2_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic4b0_2_3d/images/';
%ndirectory='/storage2/mikeg/results/spic2p3a_0_3_3d/images/';
%ndirectory='/storage2/mikeg/results/spic6p7a_0_0_3d/images/';
%ndirectory='/storage2/mikeg/results/spic3p0a_0_2_3d/images/';
%ndirectory='/storage2/mikeg/results/spic4p3a_0_1_3d/images/';
%ndirectory='/storage2/mikeg/results/spic4p3_0_1_3d/images/';
nextension='.jpg';

wspacename='6b0_2_3dmatlab_perturb.mat';
%esumcorona=0;
%esumtran=0;
%esumchrom=0;

%edifcorona=0;
%ediftran=0;
%edifchrom=0;

ebsumcorona=0;
ebsumtran=0;
ebsumchrom=0; 



%esumarray=zeros(1089,124);
%esum=zeros(1,124);
%for i=1:1089
%for i=1:1385
%period=673.4;
%period=673.4;
%nt=1611;%nt=1182;
%period=231;
%nt=1182;

period=180.0;
nt=884;
for i=334:nt %1182
%for i=334:634    

id=int2str(1000*i);
filename=[directory,'zerospic1__',id,extension];
timetext=['time=',int2str(i),'s'];
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
   
   
   
   
   nrange=3:126;
   
   ax=x(nrange);
   ay=y(nrange);
   az=z(nrange);
   [x1,x2,x3] = meshgrid(ax,ay,az);
   %val1=reshape(wd(2,nrange,nrange,nrange),124,124,124);
   %val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);


   %myval=shiftdim(val1./val2,1);
   val1=reshape(wd(5,nrange,nrange,nrange),124,124,124)+reshape(wd(9,nrange,nrange,nrange),124,124,124);
   val2=reshape(wd(5,nrange,nrange,nrange),124,124,124);
   val3=reshape(wd(9,nrange,nrange,nrange),124,124,124);
   myval=shiftdim(val1,1);
   myval2=shiftdim(val2,1);
   myval3=shiftdim(val3,1);
   for ih=1:124
    %total energy
    temp=sum(myval(:,:,ih));
    esum(ih)=sum(temp);
    %perturbed energy
    temp=sum(myval2(:,:,ih));
    esumdif(ih)=sum(temp);
    %background energy
    temp=sum(myval3(:,:,ih));
    ebsum(ih)=sum(temp);
   end

   
   
   
	esumarray(i,:)=(esum);
    edifarray(i,:)=(esumdif./esum);
       
esumcorona=esumcorona+sum(esum(48:124));
esumtran=esumtran+sum(esum(39:47));
esumchrom=esumchrom+sum(esum(1:38)); 

edifcorona=edifcorona+sum(esumdif(48:124));
ediftran=ediftran+sum(esumdif(39:47));
edifchrom=edifchrom+sum(esumdif(1:38));   
       
       
end 

escor=esumcorona/nt;
estran=esumtran/nt;
eschrom=esumchrom/nt;

edcor=edifcorona/nt;
edtran=ediftran/nt;
edchrom=edifchrom/nt;

escorper=period*esumcorona/nt;
estranper=period*esumtran/nt;
eschromper=period*esumchrom/nt;




sumper=escorper+eschromper+estranper;
sume=escor+eschrom+estran;
sumdif=(edcor+edtran+edchrom);


r1=100*escor/sume;
r2=100*eschrom/sume;
r3=100*estran/sume;

r4=100*edcor/sumdif;
r5=100*edchrom/sumdif;
r6=100*edtran/sumdif;

ebsumcorona=ebsumcorona+sum(ebsum(48:124));
ebsumtran=ebsumtran+sum(ebsum(39:47));
ebsumchrom=ebsumchrom+sum(ebsum(1:38)); 

sumeb=ebsumcorona+ebsumtran+ebsumchrom;

r7=100*ebsumcorona/sumeb;
r8=100*ebsumtran/sumeb;
r9=100*ebsumchrom/sumeb;




cmap=colormap(jet(256));
colormap(cmap);

imfile=[ndirectory,'energyvtime_',id,nextension];
surf(esumarray','LineStyle','none');
hold on


set(gca,'CameraPosition',[400 45 17320.508]);
set(gca,'YTickLabel',{'0.09';'0.99';'1.94';'2.88';'3.83';'4.77';'5.72';'6.67'})
colorbar;
xlabel(gca,'Time (seconds)');
ylabel(gca,'Height (Mm)');
title(gca,'Energy Dependence for the 0,0 Mode with a 673.4s Driver'); 
print('-djpeg', imfile); 


hold off


save(wspacename);


