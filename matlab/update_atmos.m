ntold=nt;
nt=331;

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

nt=ntold;