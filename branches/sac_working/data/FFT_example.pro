; Define the number of points and the interval:  
N = 100  
T = 0.1  
  
; Midpoint+1 is the most negative frequency subscript:  
N21 = N/2 + 1  
; The array of subscripts:  
F = INDGEN(N) 
stop 
; Insert negative frequencies in elements F(N/2 +1), ..., F(N-1):  
F[N21] = N21 -N + FINDGEN(N21-2)  
  
; Compute T0 frequency:  
F = F/(N*T)  
  
; Shift so that the most negative frequency is plotted first:  
PLOT, /YLOG, SHIFT(F, -N21), SHIFT(ABS(FFT(F, -1)), -N21)  
end
