%sumatmos
for ih=1:124
    temp=sum(myval3(:,:,ih));
    ebsum(ih)=sum(temp);
end

ebsumcorona=esumcorona+sum(ebsum(48:124));
ebsumtran=esumtran+sum(ebsum(39:47));
ebsumchrom=esumchrom+sum(ebsum(1:38)); 

sumeb=ebsumcorona+ebsumtran+ebsumchrom;

r8=100*ebsumcorona/sumeb;
r9=100*ebsumtran/sumeb;
r10=100*ebsumchrom/sumeb;