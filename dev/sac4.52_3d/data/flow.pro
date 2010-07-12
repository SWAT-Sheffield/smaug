; the field:  
vx = RANDOMU(seed, 5, 5, 5)  
vy = RANDOMU(seed, 5, 5, 5)  
vz = RANDOMU(seed, 5, 5, 5)  
  
; Set up the 3D scaling system:  
SCALE3, xr=[0,4], yr=[0,4], zr = [0,4]   
  
; Plot the vector field:  
FLOW3, vx, vy, vz 
end
