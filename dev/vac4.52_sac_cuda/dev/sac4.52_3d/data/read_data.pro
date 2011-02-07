
close,11

btot=dblarr(112,101)
timeout=dblarr(2400,101)
Bxx=dblarr(112,101)
Bzz=dblarr(112,101)

openr, 11, 's_data.dat'
 for i=0,111 do begin
  for j=0,100 do begin
		readf, 11, bb
		btot(i,j)=bb
  endfor		
 endfor
 for i=0,2399 do begin
  for j=0,100 do begin
		readf, 11, bb
		timeout(i,j)=bb
  endfor		
 endfor

 for i=0,111 do begin
  for j=0,100 do begin
		readf, 11, bb
		Bxx(i,j)=bb
  endfor		
 endfor

 for i=0,111 do begin
  for j=0,100 do begin
		readf, 11, bb
		Bzz(i,j)=bb
  endfor		
 endfor
 
 
close,11


end
