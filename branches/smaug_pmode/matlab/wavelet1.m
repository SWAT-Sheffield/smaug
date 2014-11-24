%wavelet analysis tests
%input file generated using
%energyperunitarea_array.m
load('6b0_3dmatlab_perturb.mat');
plot(esumarray(:,62));
hold on
plot(esumarray(:,31));
plot(esumarray(:,123));
surf(esumarray)
surf(esumarray(:,20:123))
scales=2:2:123;
wname='coif3'
cwt(esumarray(:,31),scales,wname,'abslvl');
cwt(esumarray(:,63),scales,wname,'abslvl');
cwt(esumarray(:,111),scales,wname,'abslvl');
wcoher(esumarray(:,31),esumarray(:,63),scales,wname);
wcoher(esumarray(:,31),esumarray(:,63),scales,wname,'plot');
wcoher(esumarray(:,31),esumarray(:,63),scales,wname,'plot','wcs');
cwt(esumarray(:,31),scales,'sym4');
CWTcoeffs=cwt(esumarray(:,31),scales,'sym4');
imagesc(t,1:163,abs(CWTcoeffs));
imagesc(1:1:900,1:163,abs(CWTcoeffs));