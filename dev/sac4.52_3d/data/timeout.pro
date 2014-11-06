close,1
openr,1,'../time.out'

as='      '
num=0.d0
time=dblarr(1)

while not(eof(1)) do begin
readf,1,num

time=[time,num]
endwhile

time=time[1:n_elements(time)-1]

plot,time

end

