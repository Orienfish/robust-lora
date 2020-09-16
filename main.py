#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import logging
import os
import copy
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
	Ptx = [20, 17, 14, 11, 8, 5]
	Ptx_max = 20		# Maximum allowed transmission power in dBm
	# Power consumption using each transmission power options in W
	# Copied from the TOSN paper (Liando 2019) for SX1276 chipset
	PowerTx = [0.4, 0.4, 0.3, 0.25, 0.2, 0.15]

	# Channel options
	CH = [0, 1, 2, 3, 4, 5, 6, 7]

	# Redundancy level
	# K = 2   			# Coverage level
	M = 2				# Connectivity level

	# Battery
	bat_cap = 2 		# Battery capacity in Ah
	bat_volt = 3.3		# Battery supply voltage in V

	# Power of MCU
	P_MCU_off = 174.65e-6 # Power of MCU (Arduino Uno) in deep sleep in mW
	P_MCU_on = 23.48e-3 # Power of MCU (Arduino Uno) on in mW
	P_R_off = 1e-4		# Power of radio in deep sleep in W

	# Power of additive white Gaussian noise with zero-mean
	N0 = 1e-15 # in W

	PDR_th = 0.7		# PDR threshold at each end node
	Lifetime_th = 1		# Lifetime threshold at each end node in years

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

def GetLifetime(SFk, Ptx, PDR, params):
	'''
	Calculate the lifetime of node using SF and Ptx using the formula in 
	the TOSN paper (Liando 2019) for SX1276 chip

	Args:
		SFk: index of Spreading Factor used by the node
		Ptx: transmission power used by the node in dBm
		PDR: packet delivery ratio at the node
		params: important parameters

	Return:
		Lifetime: estimated lifetime in years
	'''

	Ttx = params.AirTime_k[SFk] / PDR
	PowerTx = params.PowerTx[params.Ptx.index(Ptx)]
	P_R_TX = (params.T - Ttx) * (params.P_MCU_off + params.P_R_off) + \
			 Ttx * (params.P_MCU_on + PowerTx)
	P_R_TX = P_R_TX / params.T               # Average power in W
	E_bat = params.bat_cap * params.bat_volt # Total energy of battery in W*h
	Lifetime = E_bat / P_R_TX                # Lifetime in h
	return Lifetime / 24 / 365				 # Lifetime in year

def GetPDR(sr_info, G, PL, N_kq, params, idx):
	'''
	Get packet delivery ratio (based on sensitivity and SNR) at sensor idx
	with specific settings in sr_info

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix
		N_kq: a dictionary recording the list of nodes using SFk and channel q
		params: important parameters
		idx: index of sensor being considered

	Return:
		PDR: packet delivery ratio at sensor idx
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

	P_idxj = [] # Packet delivery ratio between sensor idx and gateway j
	for j in range(gw_cnt):
		if not G[j, 2]: # No gateway is placed at the current location
			continue

		# Evaluate the PDR between sensor idx and gateway j
		k = int(sr_info[idx, 2])	# SFk
		Ptx = sr_info[idx, 3]		# Transmission power
		q = int(sr_info[idx, 4])	# Channel q

		# Get the probability that the sensitivity requirement is satisfied
		Prx = propagation.GetRSSI(Ptx, PL[idx, j])
		RSSI_th = params.RSSI_k[k]     # Sensitivity threshold
		Prob_ss = 0.5 * (1 + math.erf((Prx - RSSI_th) / \
			(propagation.LogDistancePathLossModel.sigma * math.sqrt(2))))
		logging.debug('Prob_ss: {} Ptx: {} PL[i,j]: {} RSSI_th - Prx: {}'.format( \
			Prob_ss, Ptx, PL[idx, j], RSSI_th - Prx))

		# Count the number of end nodes using the same SFk and channel q and
		# calculate the collision probability under pure Aloha
		Cnt_kq = len(N_kq[str(k) + '_' + str(q)])
		Coll = 1 - math.exp(-2 * Cnt_kq * params.AirTime_k[k] / params.T)
		logging.debug('Cnt of nodes using same SF and CH: {} Coll Prob: {}'.format(Cnt_kq, Coll))
		# Sum of the expected reception power at gateway j of all nodes using
		# the same SFk and channel q
		Noise_node_idx = N_kq[str(k) + '_' + str(q)].copy()
		Noise_node_idx.remove(idx)
		Prx_dB_list = list(map(lambda i: \
			propagation.GetRSSI(sr_info[i, 3], PL[i, j]), Noise_node_idx))
		Prx_W_list = list(map(lambda Prx_dB: 10 ** (Prx_dB / 10), Prx_dB_list))
		Prx_W_sum = sum(Prx_W_list)
		Prx_W_noise = Coll * Prx_W_sum + params.N0
		logging.debug('Noise from other signals in W: {} AWGN noise in W: {}'.format( \
			Coll * Prx_W_sum, params.N0))
		Prx_dB_noise = 10 * math.log10(Prx_W_noise)

		# Get the probability that the SNR requirement is satisfied	
		Val = Prx - Prx_dB_noise - params.SNR_k[k]
		logging.debug('Val: {} Prx: {} Prx_dB_noise: {}'.format(Val, Prx, Prx_dB_noise))
		# Converted sigma
		Sigma_new = propagation.LogDistancePathLossModel.sigma * \
			math.sqrt(Coll ** 2 * (Cnt_kq - 1) + 1)
		Prob_snr = 0.5 * (1 + math.erf(Val / (Sigma_new * math.sqrt(2))))
		logging.debug('Prob_snr: {} sigma_new: {}'.format(Prob_snr, Sigma_new))

		P_idxj.append(Prob_ss * Prob_snr)
		logging.debug('Prob_ss: {} Prob_snr: {}'.format(Prob_ss, Prob_snr))

	P_idxj = np.array(P_idxj)
	PDR = 1 - np.prod(1 - P_idxj)
	return PDR

def DeviceConfiguration(sr_info_cur, sr_info, G, PL, N_kq_cur, m_gateway_cur, \
	params, sr_idx, gw_idx, SF_min):
	'''
	Greedily select the SF, Ptx and CH configuration for sensor idx

	Args:
		sr_info_cur: current complete sensor/end device configuration
		sr_info: original sensor/end device configuration
		G: gateway placement
		PL: path loss matrix
		N_kq_cur: current dictionary recording the list of nodes using SFk and channel q
		m_gateway_cur: current m-gateway connectivity at each end device
		params: important parameters
		sr_idx: index of sensor being considered
		gw_idx: index of gateway being considered
		SF_min: the min SF that can used by this sensor/end devices

	Return:
		succ: whether the configuration attempt is success
		sr_info_cur: updated configuration
	'''
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)

	succ = False # whether a new assignment is success
	for k in range(SF_min, SF_cnt):

		sr_info_cur[sr_idx, 2] = k

		for q in range(CH_cnt):

			sr_info_cur[sr_idx, 4] = q
			N_kq_cur[str(k) + '_' + str(q)].append(sr_idx)

			for Ptx in params.Ptx:
				
				sr_info_cur[sr_idx, 3] = Ptx

				PDR = GetPDR(sr_info_cur, G, PL, N_kq_cur, params, sr_idx)
				Lifetime = GetLifetime(k, Ptx, PDR, params)

				logging.debug('k: {} Ptx: {} q: {} PDR: {} Lifetime: {}'.format( \
					k, Ptx, q, PDR, Lifetime))

				if PDR >= params.PDR_th and Lifetime >= params.Lifetime_th:
					# stop searching once PDR and lifetime requirements are met
					succ = True
					break

			if succ:
				break

			# Not selected to use this assignment, remove it
			N_kq_cur[str(k) + '_' + str(q)].remove(sr_idx)

		if succ:
			break

	if succ:
		m_gateway_cur[sr_idx] += 1
		logging.debug('succ: {} k: {} Ptx: {} q: {}'.format(succ, \
			sr_info_cur[sr_idx, 2], sr_info_cur[sr_idx, 3], sr_info_cur[sr_idx, 4]))
		# Only keep the assignment if it serves as the primary connection
		# In other words, it uses the smallest SF
		if sr_info[sr_idx, 2] >= 0 and sr_info_cur[sr_idx, 2] > sr_info[sr_idx, 2]: 
			# There is existing assignment and existing SF is smaller
			# Reset to the previous device settings
			sr_info_cur[sr_idx, :] = np.copy(sr_info[sr_idx, :])

	return succ, sr_info_cur, N_kq_cur, m_gateway_cur

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
	sr_cnt = 100 #50000		# Number of sensors
	sr_info = []		# [x, y, SF, Ptx, CH]
	for i in range(sr_cnt):	
		k = -1 #random.randint(0, len(params.SF)-1) # SFk
		q = -1 # random.randint(0, len(params.CH)-1) # Channel q
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
	# print(PL)

	# Use a dictionary to record the list of nodes using the SFk and channel q
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)
	N_kq = dict()
	for k in range(SF_cnt):
		for q in range(CH_cnt):
			N_kq[str(k) + '_' + str(q)] = []

	maxDist = GetDist(propagation.LogDistancePathLossModel, params)
	print(maxDist)

	#cov_gw_sr = GetCoverage(G[:, :2], params.G_y, params.Unit_gw, sr_info[:, :2], \
	#			 maxDist[len(params.SF)-1], params.L)
	#print(cov_gw_sr)

	m_gateway = np.zeros((sr_cnt, 1)) # Current gateway connectivity at each end node
	uncover_old = sr_cnt * params.M
	rounds = 0

	# Start greedily place gateway
	while True:

		# Variables to record the best gateway location in the current round
		bnft_old = -np.inf
		
		for gw_idx in range(gw_cnt):
			if G[gw_idx, 2]: # A gateway has been placed at this location
				continue

			# Try to place gateway at gateway location idx and configure the sensor
			G[gw_idx, 2] = 1
			sr_info_cur = np.copy(sr_info)
			N_kq_cur = copy.deepcopy(N_kq)
			m_gateway_cur = np.copy(m_gateway)

			dist_gw_idx = dist[:, gw_idx] # distance to gateway idx
			sort_idx = np.argsort(dist_gw_idx, axis=0) # index of distance sorting
			for i in range(sr_cnt):
				sr_idx = int(sort_idx[i])
				if dist[sr_idx, gw_idx] > maxDist[len(params.SF)-1]:
					# the following sensors exceed the maximum communication range
					# therefore we do not need to consider them
					break

				# Calculate the min SF that can used by this end device according to distance
				SF_min = maxDist.index( \
					list(filter(lambda maxD: dist[sr_idx, gw_idx] <= maxD, maxDist))[0])
				logging.debug('SF_min: {}'.format(SF_min))
				# Try to assign the minimum SF and channel based on the current sr_info
				# and collision-possible nodes N_kq
				succ, sr_info_cur, N_kq_cur, m_gateway_cur = DeviceConfiguration(sr_info_cur, \
					sr_info, G, PL, N_kq_cur, m_gateway_cur, params, sr_idx, gw_idx, SF_min)
				
				if not succ:
					# All assignment attempts are failed, no more assignment can be made
					logging.debug('All assignment attempts are failed!')
					break
						
			# Calculate benefit of placing gateway at this location
			uncover_new = np.ones((sr_cnt, 1)) * params.M - m_gateway_cur
			uncover_new = np.sum(uncover_new[uncover_new > 0])
			bnft = uncover_old - uncover_new

			# Update global best benefit value if necessary
			if bnft > bnft_old:
				bnft_old = bnft
				next_idx = gw_idx
				next_sr_info = sr_info_cur
				next_N_kq = N_kq_cur
				next_m_gateway = m_gateway_cur
				next_uncover = uncover_new

			# Logging
			logging.info('gw_idx: {} Benefit: {} Max bnft: {} Max idx: {}'.format(\
				gw_idx, bnft, bnft_old, next_idx))

			# Reset
			G[gw_idx, 2] = 0

		# Check if there is no benefit to gain, end the searching while loop
		if bnft == 0:
			logging.info('No more gateway placement can provide m-gateway connectivity benefit!')
			break

		# Place a gateway at next_idx with the max benefit and update info
		G[next_idx, 2] = 1
		sr_info = np.copy(next_sr_info)
		N_kq = next_N_kq.copy()
		m_gateway = np.copy(next_m_gateway)
		uncover_old = next_uncover

		logging.info('Placed gateway #{} at grid {} [{},{}]'.format( \
			rounds, next_idx, G[next_idx, 0], G[next_idx, 1]))
		logging.info('Uncover: {}'.format(uncover_old))

		# Check if m-gateway connectivity has been met at all end nodes
		# If so, terminate the placement process
		if uncover_old <= 0:
			break

		rounds += 1


	plot(sr_info, G)

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
