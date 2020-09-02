#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import logging
import os
logging.basicConfig(level=logging.INFO)

import ICIOT
import propagation


########################################
# Important parameters
########################################
class params:
	L = 50000			# Edge of analysis area
	N_x = 1000			# Number of sensor potential locations on x coordinate
	N_y = 1000			# Number of sensor potential locations on y coordinate
	Unit_sr = math.floor(L / (N_x-1)) # Unit length between gw grid points
	#G_x = math.floor(N_x * Unit_sr / Unit_gw)
	#G_y = math.floor(N_y * Unit_sr / Unit_gw)
	G_x = 6     	 	# Number of gateway potential locations on x coordinate
	G_y = 6				# Number of gateway potential locations on y coordinate
	Unit_gw = math.floor(L / (G_x-1)) # Unit length between gw grid points
	desired_gw_cnt = 10 # Desired gateways to place

	T = 1200			# Sampling period in s
	Ptx_max = 20		# Maximum allowed transmission power in dBm

	# Spreading factors
	SF = ['SF7', 'SF8', 'SF9', 'SF10']
	# RSSI threshold of SF7, 8, 9, 10 in dBm
	# Copied from ICIOT paper and ICDCS paper
	RSSI_k = [-123, -126, -129, -132]
	# Air time of SF7, 8, 9, 10 in s
	# Copied from ICIOT paper under 50 Bytes payload
	AirTime_k = [0.098, 0.175, 0.329, 0.616]

	# Redundancy level
	K = 2   			# Coverage level
	M = 2				# Connectivity level

# which algorithm to run
class run:
	ICIOT = False

def GetDist(propFunc, params):
	'''
	Get the transmission distance for each SF under maximum transmission power
	using binary search

	Args:
		propFunc: selected propagation function model
		params: important parameters

	Return:
		maxDist: the max reachable distance for each SF, in m
	'''
	maxDist = []
	for i in range(len(params.SF)):
		minRSSI = params.RSSI_k[i]
		d_min = GetDist.d_min
		d_max = GetDist.d_max

		# start the binary search loop
		while d_max - d_min > GetDist.epsilon:
			d_mid = 0.5 * (d_min + d_max)
			PL = propFunc(d=d_mid, f=868)
			RSSI = propagation.GetRSSI(params.Ptx_max, PL)
			if RSSI > minRSSI: # distance is not long enough
				d_min = d_mid
			else:			   # distance is too long
				d_max = d_mid

		maxDist.append(0.5 * (d_min + d_max))
	return maxDist

GetDist.d_min = 0.0      # lower distance bound in m in the binary search
GetDist.d_max = 20000.0  # upper distance bound in m in the binary search
GetDist.epsilon = 50.0   # stop granularity in m in the binary search

def GetCoverage(grid, grid_y, Unit_grid, target, R, L):
	'''
	Get a binary coverage matrix, where (i,j) indicates whether grid point 
	i covers target j

	Args:
		grid: list of candidate grid locations
		grid_y: number of sensor potential locations on y coordinate
		Unit_grid: unit of grid points
		target: list of target locations to cover
		R: coverage radius
		L: range of the field

	Return:
		cov: a N_grid x N_target matrix indicating coverage situation
	'''
	N_grid = grid.shape[0]
	N_target = target.shape[0]
	cov = np.zeros((N_grid, N_target))

	for j in range(N_target):
		target_loc = target[j, :]
		x_min = max(0.0, target_loc[0] - R)
		x_idx_min = math.ceil(x_min / Unit_grid)
		x_max = min(L, target_loc[0] + R)
		x_idx_max = math.floor(x_max / Unit_grid)
		y_min = max(0.0, target_loc[1] - R)
		y_idx_min = math.ceil(y_min / Unit_grid)
		y_max = min(L, target_loc[1] + R)
		y_idx_max = math.floor(y_max / Unit_grid)

		# Start searching in the small rectangle
		for x_idx in range(x_idx_min, x_idx_max+1):
			for y_idx in range(y_idx_min, y_idx_max+1):
				idx = x_idx * grid_y + y_idx
				grid_loc = grid[idx, :]
				dist = np.sqrt(np.sum((grid_loc - target_loc)**2))
				if dist <= R:
					cov[idx, j] = 1
	return cov


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
	G = []				# [x, y, placed or not]
	for p in range(params.G_x):
		for q in range(params.G_y):
			new_loc = [p * params.Unit_gw, q * params.Unit_gw, 0]
			G.append(new_loc)
	G = np.array(G)
	G_cnt = params.G_x * params.G_y  # Number of gw potential locations

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
	PL = np.zeros((sr_cnt, G_cnt))
	for i in range(sr_cnt):
		for j in range(G_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			dist = np.sqrt(np.sum((loc1 - loc2)**2))
			PL[i][j] = propagation.LogDistancePathLossModel(d=dist)
	print(PL)
	

	maxDist = GetDist(propagation.LogDistancePathLossModel, params)
	print(maxDist)


	cov_gw_sr = GetCoverage(G[:, :2], params.G_y, params.Unit_gw, sr_info[:, :2], \
				 maxDist[len(params.SF)-1], params.L)
	print(cov_gw_sr)

	# Start greedily place gateway
	#while True:
	#	sr_info = []
	#	for idx in range(gw_cnt):
			# Try to place gateway at gateway location idx and configure the sensor
	#		new_sr = []




	if run.ICIOT:
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

		# Plot result
		plot(sr_info, G)


if __name__ == '__main__':
	main()
