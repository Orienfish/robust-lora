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

	# Version of log propagation model
	LogPropVer = 'Dongare'

	T = 1200			# Sampling period in s
	
	# Spreading Factors options
	SF = ['SF7', 'SF8', 'SF9', 'SF10']
	# RSSI threshold of SF7, 8, 9, 10 in dBm
	# Copied from ICIOT paper and ICDCS paper
	RSSI_k = [-123, -126, -129, -132]
	# Air time of SF7, 8, 9, 10 in s
	# Copied from ICIOT paper under 50 Bytes payload
	AirTime_k = [0.098, 0.175, 0.329, 0.616]
	# SNR threshold of SF7, 8, 9, 10 in dB
	# Copied from ICDCS paper
	SNR_k = [-6, -9, -12, -15]

	# Transmission power options in dBm
	Ptx = [2, 5, 8, 11, 14, 17, 20]
	Ptx_max = 20		# Maximum allowed transmission power in dBm
	# Power consumption using each transmission power options in W
	# Copied from the TOSN paper (Liando 2019) for SX1276 chipset
	PowerTx = [0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.4]

	# Channel options
	CH = [0, 1, 2, 3, 4, 5, 6, 7]

	# Redundancy level
	K = 2   			# Coverage level
	M = 2				# Connectivity level

	# Battery
	bat_cap = 2 		# Battery capacity in Ah
	bat_volt = 3.3		# Battery supply voltage in V

	# Power of MCU
	P_MCU_off = 174.65e-6 # Power of MCU (Arduino Uno) in deep sleep in mW
	P_MCU_on = 23.48e-3 # Power of MCU (Arduino Uno) on in mW
	P_R_off = 1e-4		# Power of radio in deep sleep in W

	# Power of additive white Gaussian noise with zero-mean
	N0 = 10

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
			PL = propFunc(d=d_mid, f=868, ver=params.LogPropVer)
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

		# only search in the surrounding grid space to save time
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

def GetLifetime(SF, Ptx, params):
	'''
	Calculate the lifetime of node using SF and Ptx using the formula in 
	the TOSN paper (Liando 2019) for SX1276 chip

	Args:
		SF: index of Spreading Factor used by the node in 'SF7' - 'SF10'
		Ptx: transmission power used by the node in dBm
		params: important parameters

	Return:
		Lifetime: estimated lifetime in years
	'''

	Ttx = params.AirTime_k[params.SF.index(SF)]
	PowerTx = params.PowerTx[params.Ptx.index(Ptx)]
	P_R_TX = (params.T - Ttx) * (params.P_MCU_off + params.P_R_off) + \
			 Ttx * (params.P_MCU_on + PowerTx)
	P_R_TX = P_R_TX / params.T               # Average power in W
	E_bat = params.bat_cap * params.bat_volt # Total energy of battery in W*h
	Lifetime = E_bat / P_R_TX                # Lifetime in h
	return Lifetime / 24 / 365				 # Lifetime in year

def GetPDR(sr_info, G, PL, noise, params):
	'''
	Get packet delivery ratio based on sensitivity and SNR

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensor i and potential gateways
		noise: total noise level at gateway j using SFk and channel q
		params: important parameters

	Return:
		PDR: a vector shows the PDR at each sensor i
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)

	# Use a dictionary to record the list of nodes using the SFk and channel q
	N_kq = {}
	for SFk in params.SF:
		for CHq in params.CH:
			N_kq[SFk + '_' + str(CHq)] = []
	for i in range(sr_cnt):
		SFk = params.SF[int(sr_info[i, 2])]
		CHq = params.CH[int(sr_info[i, 4])]
		N_kq[SFk + '_' + str(CHq)].append(i)

	for j in range(gw_cnt):
		if not G[j, 2]: # No gateway is placed at the current location
			continue

		# Evaluate the PDR between sensor i and gateway j
		for i in range(sr_cnt):
			# Get the probability that the sensitivity requirement is satisfied
			Prx = propagation.GetRSSI(sr_info[i, 3], PL[i, j])
			RSSI_th = params.RSSI_k[int(sr_info[i, 2])]     # Sensitivity threshold
			Prob_ss = 0.5 * (1 + math.erf((Prx - RSSI_th) / \
				(propagation.LogDistancePathLossModel.sigma * math.sqrt(2))))
			print('Prob_ss', sr_info[i, 3], PL[i, j], Prob_ss, RSSI_th - Prx)

			# Get the probability that the SNR requirement is satisfied
			SFk = params.SF[int(sr_info[i, 2])]
			CHq = params.CH[int(sr_info[i, 4])]
			# Count the number of end nodes using the same SFk and channel q
			Cnt_kq = len(N_kq[SFk + '_' + str(CHq)])
			# Sum of the expected reception power at gateway j of all nodes using
			# the same SFk and channel q
			Prx_list = list(map(lambda idx: \
				propagation.GetRSSI(sr_info[idx, 3], PL[idx, j]), \
				N_kq[SFk + '_' + str(CHq)]))
			Prx_sum = sum(Prx_list)
			# Left hand side
			temp = params.SNR_k[k] * (1 - math.exp(-2 * Cnt_kq * params.AirTime_k[k] / params.T))
			LHS = Prx - params.SNR_k[k] * params.N0 - temp * (Prx_sum - Prx)
			# Converted sigma
			Sigma_new = propagation.LogDistancePathLossModel.sigma * \
				math.sqrt(temp ** 2 * Cnt_kq + 1)
			Prob_snr = 0.5 * (1 + math.erf(LHS / (Sigma_new * math.sqrt(2))))
			print('Prob_snr', Prob_snr, LHS, Sigma_new)


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
	gw_cnt = params.G_x * params.G_y  # Number of gw potential locations

	# Sum of expected RSSI at gateway j using SFk and channel q, as noise
	noise = np.zeros((gw_cnt, len(params.SF), len(params.CH)))

	# Randomly generate sensor positions
	sr_cnt = 1000 #50000		# Number of sensors
	sr_info = []		# [x, y, SF, Ptx, CH]
	for i in range(sr_cnt):	
		k = random.randint(0, len(params.SF)-1) # SFk
		q = random.randint(0, len(params.CH)-1) # Channel q
		new_loc = [random.random() * params.L, random.random() * params.L, \
			k, params.Ptx_max, q]
		sr_info.append(new_loc)
	sr_info = np.array(sr_info)
	# print(sr_info)

	# Generate path loss and distance matrix between sensor i and gateway j
	PL = np.zeros((sr_cnt, gw_cnt))
	dist = np.zeros((sr_cnt, gw_cnt))
	for i in range(sr_cnt):
		for j in range(gw_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			dist[i, j] = np.sqrt(np.sum((loc1 - loc2)**2))
			PL[i, j] = propagation.LogDistancePathLossModel(d=dist[i, j], \
				ver=params.LogPropVer)

			# Update the expected noise level at each gateway
			Prx = propagation.GetRSSI(sr_info[i, 3], PL[i, j]) # Expected reception power
			noise[j, SF.index(sr_info[i, 2]), CH.index(sr_info[i, 3])] += Prx
	# print(PL)
	

	maxDist = GetDist(propagation.LogDistancePathLossModel, params)
	print(maxDist)


	cov_gw_sr = GetCoverage(G[:, :2], params.G_y, params.Unit_gw, sr_info[:, :2], \
				 maxDist[len(params.SF)-1], params.L)
	print(cov_gw_sr)

	# Start greedily place gateway
	#while True:
	#	obj_old = -np.inf
	#	next_idx = -1
	#	next_sr_info = None
	#	for gw_idx in range(gw_cnt):
	#		if G[gw_idx, 2]: # A gateway has been placed at this location
	#			continue
			# Try to place gateway at gateway location idx and configure the sensor
	#		dist_idx = dist[:, gw_idx] # distance to gateway idx
	#		sort_idx = np.argsort(min_dist, axis=0) # index of distance sorting
	#		for i in range(sr_cnt):
	#			sr_idx = int(sort_idx[i])
	#			if dist[sr_idx, gw_idx] > maxDist[len(params.SF)-1]:
					# the following sensors exceed the maximum communication range
					# therefore we do not need to consider them
	#				break

	print(GetLifetime('SF7', 11, params))
	GetPDR(sr_info[0, :], G, PL[0, :], 0, params)
	plt.figure()
	plt.scatter(sr_info[0, 0], sr_info[0, 1])
	plt.scatter(G[:, 0], G[:, 1], marker='^')
	plt.xlabel('X (m)'); plt.ylabel('Y (m)');
	plt.show()


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
