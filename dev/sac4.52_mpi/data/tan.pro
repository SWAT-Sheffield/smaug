n1=100
xmax=16000000.0
xc=8000000.0
x=findgen(n1)/n1*xmax
print,((x-xc)/xmax)


y=0.4*((atan((x-xc)/xmax*10000.d0))+!Pi/2.d0)/!Pi


print,'**',atan((x-xc)/xmax*100.d0)
plot,x-xc,y, psym=2

end
