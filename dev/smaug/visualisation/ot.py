import matplotlib.pyplot as plt
import numpy as np

# Define the smaug output variables
x=[]
y=[]
rho=[]
mx=[]
my=[]
e=[]
bx=[]
by=[]

# Read the data
data_input = open("../tmpout/zero1_0_np0202_000.out","r")

for line in data_input:
	value = line.split()
	if(len(value)== 12 and value[0]!='x'):
		x.append(float(value[0]))
		y.append(float(value[1]))
		rho.append(float(value[2]))
		mx.append(float(value[3]))
		my.append(float(value[4]))
		e.append(float(value[5]))
		bx.append(float(value[6]))
		by.append(float(value[7]))

data_input.close()

print y
#X, Y = np.meshgrid(x, y)
Z = np.sqrt(np.power(bx,2) + np.power(by,2))

plt.figure(1)
plt.contourf(X, Y, Z)
plt.show()