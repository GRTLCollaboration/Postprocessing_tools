# GRChombo
# Copyright 2012 The GRChombo collaboration.
# Please refer to LICENSE in GRChombo's root directory.
#

# script to convert the output of AHFinder into .png images for the AH and for the area and spin

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import glob, os
from os import path
import sys #argv
import math #exp

#######################################################################################
# INPUTS

# expected files "coords_AH1_%06d.dat" and "stats_AH1.dat" and same for "AH2" and "AH3" in mergers

location = "./"
center = [512., 512., 0.]

scale_mult = 0.0 # add ratio margin to plots
scale_add  = 0.0 # add margin to plots
jump  = 1 	# increase to jump through files and skip some
update_scale = False
scale_damping = 0.00

#######################################################################################
# CODE

if len(sys.argv)>1:
    location = sys.argv[1]
    if location[-1:]!='/':
        location += '/'
    print('Location set to: ',location)

num_AH = len(glob.glob("stats_*.dat"))

files  = [[]]*num_AH 
times  = files[:]
areas  = files[:]
spins  = files[:]
masses = files[:]
o_x	   = files[:]
o_y	   = files[:]
o_z	   = files[:]
c_x	   = files[:]
c_y	   = files[:]
c_z	   = files[:]

def readStats():

	for h in range(0, num_AH):

		stats = np.loadtxt(location + ("stats_AH%d.dat" % (h+1)))
		ncols = np.size(stats,1)

		times[h] = stats[:,0]
		files[h] = stats[:,1]
		areas[h] = stats[:,2]
		spins[h] = stats[:,3]
		masses[h] = stats[:,4]

		o_x[h] = stats[:,5]
		o_y[h] = stats[:,6]
		o_z[h] = stats[:,7]

		if ncols>8:
			c_x[h] = stats[:,8]
			c_y[h] = stats[:,9]
			c_z[h] = stats[:,10]
		else
			c_x[h] = o_x[h][:]
			c_y[h] = o_y[h][:]
			c_z[h] = o_z[h][:]

def plotAH():

	steps = [0]*num_AH
	active = [True]*num_AH # active = moving forward in 'stats' file
	step = 0

	mini = c-x[0][0] - center[0]
	maxi = mini

	while(True):

		if update_scale:
			mini = np.average([c_x[h][0] for h in range(num_AH)]) - center[0]
			maxi = mini

		fig = plt.figure(figsize=(10,10))
		ax = fig.add_subplot(111, projection='3d')
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')

		min_time = np.min([(1.e20 if (steps[h] >= len(times[h]) or not active[h]) else times[h][steps[h]]) for h in range(num_AH)])
		max_time = np.max([(0. if (steps[h] >= len(times[h]) or not active[h]) else times[h][steps[h]]) for h in range(num_AH)])
		fig.suptitle('Time = %09.4f' % min_time)

		for h in range(num_AH):

			# reached last step
			if steps[h] >= len(times[h]):
				active[h] = False
				continue

			# relevant for 'jump' > 1
			# when step catches up with merger, align time with 'jump' interval
			if not active[h]:
				if times[h][steps[h]] < max_time:
					active[h] = True
					while(steps[h] < len(times[h]) and times[h][steps[h]] < max_time):
						steps[h] += 1

			# yet at earlier step (e.g. merger first step is at a much later time)
			if times[h][steps[h]] > min_time:
				active[h] = False
				continue

			step_h = steps[h]
			time_h = times[h][step_h]
			file_h = files[h][step_h]

			steps[h] += jump

			# file inexistent (e.g. AHFinder diverged and added a 'nan' to the 'stats' file)
			if not path.exists(location + 'coords_AH%d_%06d.dat' % (h+1, file_h)):
				continue

			print("Plotting Horizon %d, step %d, file %d, time %f" % (h+1, step_h, file_h, time_h))

			out = np.loadtxt(location + 'coords_AH%d_%06d.dat' % (h+1, file_h), unpack=True)
			u = out[0]
			v = out[1]
			F = out[2]

			x_c = o_x[h][step_h] - center[0]
			y_c = o_y[h][step_h] - center[1]
			z_c = o_z[h][step_h] - center[2]

			x = F * np.sin(u) * np.cos(v) + x_c;
			y = F * np.sin(u) * np.sin(v) + y_c;
			z = F * np.cos(u)			  + z_c;

			if step == 0 or update_scale:
				mini = min(min(x), min(y), min(z), mini)
				maxi = max(max(x), max(y), max(z), maxi)
				
				#add scale to each side
				dx = (maxi - mini) * scale_mult + scale_add
				mini -= dx
				maxi += dx

			# Plot the surface
			ax.scatter(x, y, z, s=(0.1 if (num_AH > 1) else 5), color = 'black')
			# Plot the real center (not the origin)
			ax.scatter([c_x[h][step_h]],[c_y[h][step_h]],[c_z[h][step_h]], s=(5 if (num_AH > 1) else 30), color = 'red')

		if not any(active):
			break

		damp = math.exp(-step * scale_damping * jump)

		if step == 0 or update_scale or scale_damping != 0.:
			print("Using (min,max) scale = (%f, %f)" % (mini * damp, maxi * damp))

		ax.set_xlim([mini * damp, maxi * damp])
		ax.set_ylim([mini * damp, maxi * damp])
		ax.set_zlim([mini * damp, maxi * damp])

		plt.savefig("AHs_%06d.png" % step, bbox_inches = 'tight', dpi=150)
		plt.close()

		step += 1

	make_movie = input('Make video (y/n)? ')
	if make_movie=='y' or make_movie=='Y':
		os.system('ffmpeg -r 5 -s 1920x1080 -i AHs_%06d.png AHs.mp4')

def plotAreas():

	if(len(areas[0])>0):
		print("Plotting areas")

		plt.figure(figsize=(12,8))

		colors = cm.rainbow(np.linspace(0, 1, num_AH))
		for h in range(0, num_AH):
			plt.scatter(times[h],areas[h], label = "Area AH%d" % (h+1), facecolor=colors[h])
		plt.ylabel("Area ")
		plt.xlabel('t $[m]$')
		plt.legend(loc= "best")
		# #plt.ylim([34,30])
		plt.savefig("areasAHs.png", bbox_inches = 'tight')
		plt.close()

def plotSpins():

	if(len(spins[0])>0):
		print("Plotting spins")

		plt.figure(figsize=(12,8))

		colors = cm.rainbow(np.linspace(0, 1, num_AH))
		for h in range(0, num_AH):
			plt.scatter(times[h],spins[h], label = "Spin AH%d" % (h+1), facecolor=colors[h])
		plt.ylabel("Spin ")
		plt.xlabel('t $[m]$')
		plt.legend(loc= "best")
		# #plt.ylim([34,30])
		plt.savefig("spinsAHs.png", bbox_inches = 'tight')
		plt.close()

readStats()
plotAreas()
plotSpins()
# plotAH()
