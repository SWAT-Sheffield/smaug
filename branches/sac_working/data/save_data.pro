
openw, 11, 's_data.dat'
 for i=0,111 do begin
  for j=0,100 do begin
		printf, 11, btot(i,j)
  endfor		
 endfor
 for i=0,2399 do begin
  for j=0,100 do begin
		printf, 11, timeout(i,j)
  endfor		
 endfor

 for i=0,111 do begin
  for j=0,100 do begin
		printf, 11, Bxx(i,j)
  endfor		
 endfor

 for i=0,111 do begin
  for j=0,100 do begin
		printf, 11, Bzz(i,j)
  endfor		
 endfor
 
 
close,11


end
