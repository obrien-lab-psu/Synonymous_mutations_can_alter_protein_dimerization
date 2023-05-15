#!/bin/usr/python

# Written by Dan Nissley, February 2016

# Modified to be lean and mean by Dan Nissley, March 2017

# Calculate the mean first passage time of a protein Go model from the unfolded
# to the folded state in order to determine its mean rate of folding in bulk solution.

# April 4th 2017 - updated to print centered time bins not the time of the frame

# August 18th 2017 - updated to require that the protein remain folded for 10 frames
#                    after it initially folds before it is considered to have folded
#		     at the initial time. This version also assumes you are using
#		     fraction of native contacts, not RMSD, as your order parameter

# load modules
import sys
import os
import numpy as np

# check to make sure the correct number of commandline arguments have been entered
if len(sys.argv) != 6:

	print 'Usage is: python fpt_analysisV3.2.py [1] [2] [3] [4] [5]'
	print '[1] = file name for file containing list of props files to analyze'
	print '[2] = number of frames in each props file'
	print '[3] = the Q threshold separating the folded and unfolded states'
	print '[4] = the number of steps per frame'
	print '[5] = the time step in units of femtoseconds'
	print 'ALL COMMAND LINE ARGUMENTS ARE REQUIRED'

	# die 
	sys.exit()

else:

	pass

# list of files containing Q time series
f_list = open(sys.argv[1])           

# the total number of frames in each props_ts file; equal to (nsteps/nsavc)
nframes = int(sys.argv[2])           

# the Q threshold to differentiate F molecules from U molecule 
q_cut = float(sys.argv[3])           

# number of integration steps simulated between printing coordinates in quench phase
steps_per_frame = int(sys.argv[4]) 

# time step in units of picoseconds
tstep = float(sys.argv[5])

# make empty matrix to hold folding status of protein
num_folded = np.zeros((nframes, 1))

# counter for total number of trajectories
num_traj = 0

# counter for the final frame to print;
# this is so we don't print a bunch of unneeded
# zeroes once all trajectories have folded
max_frame = 0

# loop through the set of props input files
for file in f_list:

	# split out the new line character
	file_name = file.split('\n')[0]

	# load the contents of the file and save the
	f_q = np.loadtxt(file_name)
	
	# the starting state of the protein is unfolded
	folded = 0

	# loop over frames in the trajectory
	for t1 in range (0, len(f_q)-10):
		
		# if the threshold is reached
		if (f_q[t1] > q_cut) and (folded == 0):

			# define a new Boolean parameter;
			# this parameter checks to see
			# if the F state is kinetically stable
			still_folded = True

			# loop through a set of ten frames, starting
			# with the frame after the current frame 
			for t2 in range (t1+1, t1+11):

				# check to see if this frame is also F
				if f_q[t2] > q_cut:

					# if it is, then do nothing
					pass

				# if the frame isn't folded, change still_folded to False,
			        # indicating the folded state was not kinetically stable
				else:

					still_folded = False

				#print t1, t2, f_q[t1], f_q[t2], still_folded

			# if the protein remained folded for 10 additional frames
			if still_folded == True:

				# then consider it to have passed to F at t1
				folded = 1
				
				#print '%.9f' %float((t1*steps_per_frame*tstep)+(0.5*steps_per_frame*tstep))
				
				# if a protein just folded at this frame
				if t1 > max_frame:

					# modify the value of max_frame
					max_frame = t1

		# if the threshold hasn't been reached
		else:	
			#do nothing
			pass

		# if this traj has folded in this step or any other
		# then it contributes +1 to the number of folded trajectories
		num_folded[t1] += folded		

	# increment the total number of trajectories analyzed
	num_traj += 1

#sys.exit()

#print 'Frame   ', '\t', 'Time (ps)', '\t', 'Survival Probability of Unfolded State'
#print '---------------------------------------------'
print '%.5f' %(0), '\t', '%.9f' %float(0), '\t',  '%.9f' %(1)

for value in range (0, len(f_q)-10):

	#Assume two-state folding behavior - if the protein isn't unfolded it is folded
	# i.e., P_U = 1 - P_F

	print '%.5f' %(value+1), '\t', '%.9f' %float((value*steps_per_frame*tstep)+(0.5*steps_per_frame*tstep)), '\t',  '%.9f' %(1-(num_folded[value]/float(num_traj)))
