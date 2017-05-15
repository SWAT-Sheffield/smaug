from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os

# SETUP __________________________________________________________________

# Setup the output picture resulution
dim_x = 252	 
dim_y = 252

# The number of colors 
number_of_colors = 512

# Maximum number of frames
number_of_frames = 5000

# DPI quality
quality = 100

# True of False, black background for video of white for publication
black_background = False

# Font setup
font = {'weight' : 'light', 'size'   : 8}

# Frame title, e.g OT / Density / T = 0.003 s
title = 'OT / Density / T = '

# Create a video (True or False)
video = False

# Colormap
cm = cm.afmhot
# ________________________________________________________________________

def grid(x, y, z, resX=dim_x, resY=dim_y):
	'''
	This Function merges the different regions into a 2 dimensional grid.

	Parameters
	----------
		x - The x coordinates from the SMAUG+ result files - spatial
		y - The y coordinates from the SMAUG+ result files - spatial
		z - The z coordinates from the SMAUG+ result files - rho, e, bx, bz

	Returns
	-------
		X - meshgrid x coordinates
		Y - meshgrid y coordinates
		Z - 2D iamge of rho, e, bx, bz...
	'''

	# Create the x and y domain
	xi = np.linspace(min(x), max(x), resX)
	yi = np.linspace(min(y), max(y), resY)

	# Create the datagrid
	Z = griddata(x, y, z, xi, yi,  interp='linear')

	# Create the meshgrid
	X, Y = np.meshgrid(xi, yi)
	return X, Y, Z

def black_bg(ax):
	'''
	This fuction changes the axis and background colors.

	Parameters:
		ax - by matplotlib

	Returns:
	'''

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_linewidth(0.1)
	ax.spines['left'].set_linewidth(0.1)
	ax.spines['bottom'].set_color('white')
	ax.spines['left'].set_color('white')

	ax.title.set_color('white')
	ax.yaxis.label.set_color('white')
	ax.xaxis.label.set_color('white')
	ax.tick_params(axis='x', colors='white')
	ax.tick_params(axis='y', colors='white')

	ax.tick_params(axis='both', direction='out')
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	return 0

for i in range(0, number_of_frames):

	# Array Initialzation
	x , y= [],[]
	mx, my=[],[]
	bx, by=[],[]
	rho, e=[],[]

	frame = "../tmpout/flux1_" + str(i)
	frame = frame + "_np0101_00"

	# Merge the regions

	for j in range(0,1):
		name = frame
		name = name + str(j) + ".out"
		data_input = open(name,"r")
		print name
		for line in data_input:

			# Read line by line
			value = line.split()

			# Read the header information
			if(len(value)== 5):
				t = float(value[1])

			# Read the data
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
	plt.scatter(x, y, c=rho, alpha=0.5, edgecolors='none')
	plt.show()
	kaka
	# Create a frame

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)

	# Apply the font setup
	plt.rc('font', **font)

	# Setup black background
	if black_background: black_bg(ax)

	# Axis and image title 
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_title(title + str(t) + ' s')

	X, Y, Z = grid(x, y, rho)
	sc = ax.contourf(X, Y, Z, number_of_colors, cmap=cm)

	color_bar = plt.colorbar(sc)
	color_bar.ax.tick_params(axis='y', direction='out')

	if black_background:

		# Colorbar axis colors and properties
		color_bar.ax.yaxis.set_tick_params(color='white')
		plt.setp(plt.getp(color_bar.ax.axes, 'yticklabels'), color='white')

		# Save the figure
		plt.savefig(str(i) + '.png', format='png', dpi=quality, facecolor='black')
	else:

		# Colorbar axis colors and properties
		color_bar.ax.yaxis.set_tick_params(color='black')
		plt.setp(plt.getp(color_bar.ax.axes, 'yticklabels'), color='black')

		# Save the figure
		plt.savefig(str(i) + '.png', format='png', dpi=quality, facecolor='white')

	# Clear axis, figure and close
	plt.cla()
	plt.clf()
	plt.close() 

if video:
	os.system("ffmpeg -f image2 -i %d.png -vcodec mpeg4 -b 800k video2.avi")

