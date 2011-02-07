% Default settings for all the macros

%clear;
%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007
%resource = findResource('scheduler', 'configuration', 'sge');

%The set function can be used to set different propoerties of the resource
%set(resource, 'configuration', 'sge');
defaultParallelConfig('sge');
matlabpool('open','sge',16);
% Confirmation of parameters already set
global Doask;
Doask=0;

% File parameters
filename='';
physics='';
phys='';
ndir=[];
logfilename='';

% Transformation parameters
Transform='';
nxreg=[];
nxregold=[];
polar_r='';
polar_phi='';

% Function parameters
func='';
autorange=[];
fmin=[];
fmax=[];

% Plotting parameters
cut='';
plotmode='';
plottitle='default';
View=[-37.5 30];
Colorbar=0;
Shading='flat';
Contourlevel=30;
Contourstyle='g-';
Quiverscale=1;

% multiplot=[] gives the default number of subplots depending on nfile,nfunc
% multiplot=[3,2] defines 3 by 2 subplots
multiplot=[];

% Number of info items on bottom and in the header
Bottomline=2;
Headerline=2;

% Animation parameters
npict=[];
firstpict=1;
dpict=1;
npictmax=100;
doanimate=1;

% Printing parameters
dohardplot=0;
Printfile='Movie/matlab';
Device='ps';
Orient='landscape';

%mypars


% Read the npict-th picture from 1 or more files
filename='../../testcont2.out';

%filename=askstr('filename(s)  ',filename);
[filenames,nfile]=str2arr(filename);
[asciifile,filesize,pictsize,fileerror]=studyfile(filenames);
if fileerror; return; end;
dispnum('asciifile(s)  =',asciifile);
disp(' ');
npictinfile=floor(filesize./pictsize);
fprintnum('npictinfile(s)=',npictinfile);
if ~isempty(npict); fprintf('\n'); end;



parfor npict=1:48

%npict=415;
outname=sprintf('../../ascdata/testascmat_%d.out',npict);
imoutname=sprintf('../../plots/testascmat_%d.jpg',npict);    
    

%if max(npictinfile)>1
%   npict=asknum('npict(s) (eg. 2 or [3 4 5])',npict,nfile);
%else
%   npict=npictinfile;
%   dispnum('npict         =',npict);
%end

for ifile=1:nfile
   if npict(ifile)>npictinfile(ifile)
      disp(' ');
      disp(['Reducing npict=' num2str(npict(ifile)) ' for file ' ...
         num2str(ifile) ' to npictinfile=' num2str(npictinfile(ifile)) ' !']);
      npict(ifile)=npictinfile(ifile);
   end;
   fid=fopen(trim(filenames(ifile,:)));
   fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   %get_head;
   
   % Read header of a VAC data file

fileerror=0;
if asciifile(ifile)
   headline=trim(fgetl(fid));
   if ~isstr(headline);fileerror=1;return;end;
   tmp=sscanf(fgetl(fid),'%d %f %d %d %d',5)';
   it=tmp(1); time=tmp(2); ndim=tmp(3); neqpar=tmp(4); nw=tmp(5); clear tmp;
   gencoord=ndim<0; ndim=abs(ndim);
   nx=sscanf(fgetl(fid),'%d',ndim)';
   eqpar=sscanf(fgetl(fid),'%f',neqpar)';
   Variables=trim(fgetl(fid));
else
   [tmp,ntmp]=fread(fid,4);
   if ntmp<4;fileerror=1;return;end;
   headline=trim(setstr(fread(fid,79)'));
   fread(fid,4);

   fread(fid,4);
   it=fread(fid,1,'int32'); time=fread(fid,1,'float64');
   ndim=fread(fid,1,'int32');
   gencoord=ndim<0; ndim=abs(ndim);
   neqpar=fread(fid,1,'int32'); nw=fread(fid,1,'int32');
   fread(fid,4);

   fread(fid,4);
   nx=fread(fid,ndim,'int32');
   fread(fid,4);

   fread(fid,4);
   eqpar=fread(fid,neqpar,'float64');
   fread(fid,4);

   fread(fid,4);
   Variables=trim(setstr(fread(fid,79)'));
   fread(fid,4);
end

% Extract physics from headline if it is defined
i=length(headline);
while i>1 & headline(i)~='_'; i=i-1; end;
if headline(i)=='_'
   physics=headline(i+1:length(headline));
   headline=headline(1:i-1);
end
% Extraxt number of vector components NDIR from last character of physics
% and extract phys without the number of dimesions and components
if ~isempty(physics)
   ndir=str2num(physics(length(physics)));
   phys=physics(1:length(physics)-2);
end

% Extract info from names of Variables
[variables,ntmp]=str2arr(Variables);
xnames=variables(1:ndim,:);
wnames=variables(ndim+1:ndim+nw,:);

% It can optionally contain the names of the equation parameters
if ntmp==ndim+nw+neqpar
   eqparnames=variables(ndim+nw+1:ntmp,:);
   for ieqpar=1:neqpar
      eval([eqparnames(ieqpar,:) '= eqpar(ieqpar);']);
   end
end
clear ntmp;

nxs=prod(nx);
if ndim==1
   nx1=nx;
   nx2=1;
   nx3=1;
end
if ndim==2
   nx1=nx(1);
   nx2=nx(2);
   nx3=1;
end
if ndim==3
   nx1=nx(1);
   nx2=nx(2);
   nx3=nx(3);
end

   
   
   disp(' ');
   disp(['headline=''',headline,''' ']);
   if ~isempty(physics); disp(['physics =''',physics,''' ']); end;
   if fileerror
      fprintf('Error: Could not read the %d',dpict);
      disp(['. picture from data file ' filenames(ifile,:)]);
      return
   end
   fprintf('it=%d time=%f ndim=%d neqpar=%d nw=%d\n',it,time,ndim,neqpar,nw);
   dispnum('eqpar=',eqpar);
   dispnum('nx   =',nx);
   if gencoord; read_transform_par; end;
   if nfile>1
      disp(['Reading:',Variables,' (with _',num2str(ifile),')']);
   else
      disp(['Reading:',Variables]);
   end;
   %get_body;
   
   
   % Read a snapshot from a VAC data file (binary or ascii).
% First X and w are read, (X is capitalized to avoid possible conflict with 
% the name of its first component). Transform X and w according to Transform
% if generalized coordinates are found.
% Extract the variables given in the xwname string array read by get_head.

X=zeros(nxs,ndim);
w=zeros(nxs,nw);
if asciifile(ifile)
   xw=fscanf(fid,'%f',[ndim+nw,nxs]);
   X=xw(1:ndim,:)';
   w=xw(ndim+1:nw+ndim,:)';
   clear xw;
   fgetl(fid);
else
   fread(fid,4);
   for idim=1:ndim
      X(:,idim)=fread(fid,nxs,'float64');
   end
   fread(fid,4);
   for iw=1:nw
      fread(fid,4);
      w(:,iw)=fread(fid,nxs,'float64');
      fread(fid,4);
   end
end

% extract variables from X into variables named after the strings in xnames
for idim=1:ndim
  if ndim==2
     tmp=reshape(X(:,idim),nx1,nx2);
  else
     tmp=X(:,idim);
  end
  eval([xnames(idim,:),'=tmp;']);
end
if ndim==1
   xx=X;
elseif ndim==2
   xx=reshape(X(:,1),nx1,nx2);
   yy=reshape(X(:,2),nx1,nx2);
end

% extract variables from w into variables named after the strings in wnames
for iw=1:nw
  if ndim==2
     tmp=reshape(w(:,iw),nx1,nx2);
  else
     tmp=w(:,iw);
  end
  eval([wnames(iw,:),'=tmp;']);
end
clear tmp;

if gencoord & ndim==2 
   if strcmp(Transform,'polar')
      polargrid;
   elseif strcmp(Transform,'regular')
      regulargrid;
   end
end

   
   
   
   
   fclose(fid);
   if nfile>1
      % Rename the variables to rho_1,rho_2 etc for more than one file
      for i=1:size(variables,1)
        eval([trim(variables(i,:)) '_' num2str(ifile) '=' variables(i,:) ';']);
      end
   end
end

[nr,nc]=size(x);
ntot=nr*nc;

ascdat=zeros(nr,nc,12);
ascdat(:,:,1)=x(:,:);
ascdat(:,:,2)=y(:,:);
ascdat(:,:,3)=h(:,:);
ascdat(:,:,4)=m1(:,:);
ascdat(:,:,5)=m2(:,:);
ascdat(:,:,6)=e(:,:);
ascdat(:,:,7)=b1(:,:);
ascdat(:,:,8)=b2(:,:);
ascdat(:,:,9)=eb(:,:);
ascdat(:,:,10)=rhob(:,:);
ascdat(:,:,11)=bg1(:,:);
ascdat(:,:,12)=bg2(:,:);

nascdata=reshape(ascdat,ntot,12);
save( outname, 'nascdata', '-ascii')



h1=figure('Visible','off','IntegerHandle','Off');

hold on;

%Need to alter the camer position so
%hget the child property of the figure which are the axes
hax=get(h1,'Children');
set(hax,'CameraPosition',[-203.3 -1635.05 866.03]);

%mesh(xx,yy,zz);
%plot3(x,y,z,'o');
subplot(1,2,1);
surf(m1(1:8:1976,1:8:400)./(h(1:8:1976,1:8:400)+rhob(1:8:1976,1:8:400)),'EdgeColor','none');
subplot(1,2,2);
surf(m2(1:8:1976,1:8:400)./(h(1:8:1976,1:8:400)+rhob(1:8:1976,1:8:400)),'EdgeColor','none');




%print the most recent figure to a file
print( '-djpeg', imoutname);
hold off;


end
