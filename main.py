#!/usr/bin/python3
import numpy as np
import math
import random


def FreeSpacePathLossModel(d, f):
	'''
	Free Space Path Loss Model

	Args:
		d: distance in m
		f: frequency in MHz

	Return:
		path loss in dB at d
	'''
	return 20 * math.log10(d) + 20 * math.log10(f) - 27.55

def LogDistancePathLossModel(d):
	'''
	Log Distance Path Loss Model

	Args:
		d: distance in m

	Return:
		path loss in dB at d
	'''
	return LogDistancePathLossModel.PL0 + 10 * LogDistancePathLossModel.delta * \
		math.log10(d / LogDistancePathLossModel.d0)

# Settings in the ICIOT paper
LogDistancePathLossModel.d0 = 1000 # Reference distance in m
LogDistancePathLossModel.PL0 = 130 # Reference path loss in dB
LogDistancePathLossModel.delta = 2.1 # Path loss exponent
# Settings in the IPSN paper (indoor building)
#LogDistancePathLossModel.d0 = 40 # Reference distance in m
#LogDistancePathLossModel.PL0 = 127.41 # Reference path loss in dB
#LogDistancePathLossModel.delta = 3.57 # Path loss exponent


# Important parameters
#N_x = 1000
#N_y = 1000
#N_cnt = N_x * N_y
#Unit_sr = 50
L = 50000			# Edge of analysis area
#G_x = math.floor(N_x * Unit_sr / Unit_gw)
#G_y = math.floor(N_y * Unit_sr / Unit_gw)
G_x = 6				# Number of gw potential locations on x coordinate
G_y = 6				# Number of gw potential locations on y coordinate
G_cnt = G_x * G_y   # Number of gw potential locations
Unit_gw = math.floor(L / G_x) # Unit length between gw grid points
print(Unit_gw)

# RSSI threshold of SF7, 8, 9, 10
# Copied from ICIOT paper and ICDCS paper
RSSI_k = [-123, -126, -129, -132]

# Generate the grid candidate set N and G with their x, y coordinates
# N for sensor placement and G for gateway placement
empty_position = np.zeros((1, 2))
#N = np.tile(empty_position, (N_cnt, 1))
G = np.tile(empty_position, (G_cnt, 1))
#for i in range(N_x):
#	for j in range(N_y):
#		N[i][j][:] = [i * Unit_sr, j * Unit_sr]
for p in range(G_x):
	for q in range(G_y):
		G[p * G_x + q][:] = [p * Unit_gw, q * Unit_gw]

# Randomly generate sensor positions
sr_cnt = 50000		# Number of sensors
sr_loc = []
for i in range(sr_cnt):	
	new_loc = [random.random() * L, random.random() * L]
	sr_loc.append(new_loc)
sr_loc = np.array(sr_loc)
print(sr_loc)

# Generate path loss matrix PL between sensor i and gateway j at (p,q)
PL = np.zeros((sr_cnt, G_cnt))
for i in range(sr_cnt):
	for j in range(G_cnt):
		dist = np.sqrt(np.sum((sr_loc[i][:] - G[j][:])**2))
		PL[i][j] = LogDistancePathLossModel(dist)

print(PL)
