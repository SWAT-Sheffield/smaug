

directory1='/storage2/mikeg/results/spic4p3a_0_1_3d/vapor_4p3a_0_1_3d/';
directory2='/storage2/mikeg/results/spic4p3a_0_1_3d/images_3d_vsecs/';
extension='.out';

ndirectory='/storage2/mikeg/results/spic4p3a_0_1_3d/addimages/';

nextension='.jpg';

for i=1:1:882
%for i=334:634

%i1=i+51;
i1=i-1;
id1=int2str(i1);
if i1<10
   id1=['000',int2str(i1)]; 
elseif i1<100
   id1=['00',int2str(i1)];                
elseif i1 <1000
    id1=['0',int2str(i1)];       
elseif i1<1000
    id1=[int2str(i1)];   
end


id2=int2str(1000*i);


imfile1=[directory1,'im',id1,nextension];
imfile2=[directory2,'im1_',id2,nextension];

outimfile=[ndirectory,'im_',id2,nextension];

im1=imread(imfile1);
im2=imread(imfile2);

[nr,nc,ncc]=size(im1);

im2n = imresize(im2, [nr nc]);

imout=[im1 im2n];
%hold on
%subplot(1,2,1);
%imshow(im1);
%subplot(1,2,2);
%imshow(im2);

imshow(imout);





%title(gca,'The 0,0 Mode with a 670.0s Driver'); 
print('-djpeg', outimfile); 


%hold off






end


