#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

import ICIOT
import propagation


########################################
# Important parameters
########################################
class params:
	#N_x = 1000
	#N_y = 1000
	#N_cnt = N_x * N_y
	#Unit_sr = 50
	L = 50000			# Edge of analysis area
	#G_x = math.floor(N_x * Unit_sr / Unit_gw)
	#G_y = math.floor(N_y * Unit_sr / Unit_gw)
	G_x = 6     		# Number of gw potential locations on x coordinate
	G_y = 6			# Number of gw potential locations on y coordinate
	Unit_gw = math.floor(L / G_x) # Unit length between gw grid points
	desired_gw_cnt = 10 # Desired gateways to place

	T = 1200			# Sampling period in s
	Ptx_max = 23		# Maximum allowed transmission power in dBm

	# Spreading factors
	SF = ['SF7', 'SF8', 'SF9', 'SF10']
	# RSSI threshold of SF7, 8, 9, 10 in dBm
	# Copied from ICIOT paper and ICDCS paper
	RSSI_k = [-123, -126, -129, -132]
	# Air time of SF7, 8, 9, 10 in s
	# Copied from ICIOT paper under 50 Bytes payload
	AirTime_k = [0.098, 0.175, 0.329, 0.616]


def plot(sr_info, G):
	# Visualize the placement and device configuration
	# sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	plt.figure()
	colorList = ['b', 'g', 'y', 'm', 'w', 'r', 'k', 'c']
	color = [colorList[int(i)] for i in list(sr_info[:, 2])]
	plt.scatter(sr_info[:, 0], sr_info[:, 1], c=color , s=5)
	color = ['r' for i in range(gw_cnt)]
	plt.scatter(G[:, 0], G[:, 1], s=G[:, 2]*50, c=color, marker='^')
	plt.xlabel('X (m)'); plt.ylabel('Y (m)');
	# plt.legend()
	plt.savefig('vis.png')
	plt.show()

########################################
# Main Process
########################################
def main():
	# Preparation
	# Generate the grid candidate set N and G with their x, y coordinates
	# N for sensor placement and G for gateway placement
	# N = []
	G = []				# [x, y, placed or not]
	#for i in range(N_x):
	#	for j in range(N_y):
	#		new_loc = [i * Unit_sr, j * Unit_sr]
	#		N.append(new_loc)
	for p in range(params.G_x):
		for q in range(params.G_y):
			new_loc = [p * params.Unit_gw, q * params.Unit_gw, 0]
			G.append(new_loc)
	G = np.array(G)
	gw_cnt = params.G_x * params.G_y  # Number of gw potential locations

	# Randomly generate sensor positions
	sr_cnt = 1000 #50000		# Number of sensors
	sr_info = []		# [x, y, SF, Ptx]
	for i in range(sr_cnt):	
		k = random.randint(0, len(params.SF)-1)
		new_loc = [random.random() * params.L, random.random() * params.L, \
			k, params.Ptx_max]
		sr_info.append(new_loc)
	sr_info = np.array(sr_info)
	# print(sr_info)

	# Generate path loss matrix PL between sensor i and gateway j at (p,q)
	PL = np.zeros((sr_cnt, gw_cnt))
	for i in range(sr_cnt):
		for j in range(gw_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			dist = np.sqrt(np.sum((loc1 - loc2)**2))
			PL[i][j] = propagation.LogDistancePathLossModel(dist)
	# print(PL)

	gw_place, sr_info = ICIOT.ICIOTAlg(sr_info, G, PL, params)
	print(gw_place)

	# Write sensor and gateway information to file
	with open ('sr_loc.txt', 'w') as out:
		for i in range(sr_cnt):
			out.write(str(round(sr_info[i, 0], 2)) + ' ' + \
				str(round(sr_info[i, 1], 2)) + ' ' + \
				str(int(sr_info[i, 2])) + ' ' + \
				str(round(sr_info[i, 3])) + '\n')
	with open ('gw_loc.txt', 'w') as out:
		gw_cnt = G.shape[0]
		for i in range(gw_cnt):
			if G[i, 2]:
				out.write(str(round(G[i, 0], 2)) + ' ' + \
					str(round(G[i, 1], 2)) + '\n')

	plot(sr_info, G)


if __name__ == '__main__':
	main()
