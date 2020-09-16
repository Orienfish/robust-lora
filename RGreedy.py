#!/usr/bin/python3
import numpy as np
import math
# import matplotlib.pyplot as plt
import logging
import copy

import propagation

########################################
# RGreedy Alg
########################################
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


def RGreedyAlg(sr_info_ori, G_ori, PL, dist, params):
	'''
	Call the robust gateway placement algorithm

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		dist: distance matrix between sensors and potential gateways
		params: important parameters

	Returns:
		sr_info: sensor configuration
		G: resulted gateway placement
	'''
	# Make a deep copy of the original numpy array to avoid changes
	sr_info = np.copy(sr_info_ori)
	G = np.copy(G_ori)
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

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
		bnft_best = -np.inf
		
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
			if bnft > bnft_best:
				bnft_best = bnft
				next_idx = gw_idx
				next_sr_info = sr_info_cur
				next_N_kq = N_kq_cur
				next_m_gateway = m_gateway_cur
				next_uncover = uncover_new

			# Logging
			logging.debug('gw_idx: {} Benefit: {} Max bnft: {} Max idx: {}'.format(\
				gw_idx, bnft, bnft_best, next_idx))

			# Reset
			G[gw_idx, 2] = 0

		# Check if there is no benefit to gain, end the searching while loop
		if bnft_best == 0:
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

	return sr_info, G, m_gateway