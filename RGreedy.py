#!/usr/bin/python3
import numpy as np
import math
# import matplotlib.pyplot as plt
import logging
import copy

import propagation

########################################
# Robust-Driven Greedy Alg
########################################

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
		P_idxj: packet delivery ratio between sensor idx and gatway j
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
		Noise_node_info = N_kq[str(k) + '_' + str(q)].copy()
		#print(Noise_node_info)
		#print(sr_info[idx, :])
		Noise_node_info = [x for x in Noise_node_info if not (x[0] == sr_info[idx, 0] and \
			x[1] == sr_info[idx, 1])] # Remove the same node
		#print(Noise_node_info)

		Prx_W_list = []
		for noise_node in Noise_node_info:
			loc1 = noise_node[:2]
			loc2 = G[j, :2]
			Noise_PL = propagation.LogDistancePathLossModel(d=np.sqrt(np.sum((loc1 - loc2)**2)), \
				ver=params.LogPropVer)
			Prx_dB = propagation.GetRSSI(noise_node[3], Noise_PL)
			Prx_W_list.append(10 ** (Prx_dB / 10))
		#print(Prx_W_list)

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
	return PDR, P_idxj

def DeviceConfiguration(sr_info_cur, sr_info, G, PL, N_kq_cur, m_gateway_cur, \
	params, sr_idx, gw_idx):
	'''
	Greedily select the SF, Ptx and CH configuration for sensor sr_idx considering
	a new gateway is placed at gw_idx

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

	Return:
		succ: whether the configuration attempt is success
		sr_info_cur: updated device configuration
		N_kq_cur: updated collision dictionary
		m_gateway_cur: updated m-gateway connectivity status
		PDR: packet delivery ratio due to the configuraton
		Lifetime: lifetime of the device due to the configuration
	'''
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)

	succ = False # whether a new assignment is success
	for k in range(SF_cnt):

		sr_info_cur[sr_idx, 2] = k
		for q in range(CH_cnt):

			sr_info_cur[sr_idx, 4] = q
			N_kq_cur[str(k) + '_' + str(q)].append(sr_info_cur[sr_idx, :])
			for Ptx in params.Ptx:
				
				sr_info_cur[sr_idx, 3] = Ptx

				PDR, _ = GetPDR(sr_info_cur, G, PL, N_kq_cur, params, sr_idx)
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
			N_kq_cur[str(k) + '_' + str(q)] = [x for x in N_kq_cur[str(k) + '_' + str(q)] \
				if not (x[0] == sr_info_cur[sr_idx, 0] and x[1] == sr_info_cur[sr_idx, 1])]

		if succ:
			break

	if succ:
		m_gateway_cur[sr_idx] += 1
		logging.debug('succ: {} k: {} Ptx: {} q: {}'.format(succ, \
			sr_info_cur[sr_idx, 2], sr_info_cur[sr_idx, 3], sr_info_cur[sr_idx, 4]))

		if sr_info[sr_idx, 2] >= 0:
			# If there are multple connectivities
			new_k = int(sr_info_cur[sr_idx, 2])
			old_k = int(sr_info[sr_idx, 2])
			new_q = int(sr_info_cur[sr_idx, 4])
			old_q = int(sr_info[sr_idx, 4])

			if new_k == old_k and new_q == old_q:
				# If the new configuration is the same as the old configuration
				# remove duplicates from collision dictionary

				# First remove both of them
				N_kq_cur[str(new_k) + '_' + str(new_q)] = \
					[x for x in N_kq_cur[str(new_k) + '_' + str(new_q)] \
					if not (x[0] == sr_info[sr_idx, 0] and x[1] == sr_info[sr_idx, 1])]
				# Then add back one
				N_kq_cur[str(new_k) + '_' + str(new_q)].append(sr_info_cur[sr_idx, :])

			elif new_k > old_k or \
				(new_k == old_k and sr_info_cur[sr_idx, 3] > sr_info[sr_idx, 3]): 
				# This new connectivity should be backup connectivity
				# Remove the new node from collision dictionary
				# Reset to the previous device settings 
				N_kq_cur[str(new_k) + '_' + str(new_q)] = \
					[x for x in N_kq_cur[str(new_k) + '_' + str(new_q)] \
					if not (x[0] == sr_info_cur[sr_idx, 0] and x[1] == sr_info_cur[sr_idx, 1])]
				sr_info_cur[sr_idx, :] = np.copy(sr_info[sr_idx, :])
			else:
				# The old connectivity should be backup connectivity
				# Remove the old node from collision dictionary
				# Keep the current device settings 
				N_kq_cur[str(old_k) + '_' + str(old_q)] = \
					[x for x in N_kq_cur[str(old_k) + '_' + str(old_q)] \
					if not (x[0] == sr_info[sr_idx, 0] and x[1] == sr_info[sr_idx, 1])]

	else:
		# All device configuration attempts are failed, reset current sensor info
		logging.debug('All device configuration attempts are failed!')
		sr_info_cur[sr_idx, :] = np.copy(sr_info[sr_idx, :])


	return succ, PDR, Lifetime

def UpdateConn(sr_info, G, PL, N_kq, params):
	'''
	Update current gateway connectivity at each end node, since there might be
	gateway placement already in the G

	Args:
		sr_info: original sensor/end device configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		N_kq: a dictionary recording traffic allocation
		params: important parameters

	Return:
		sr_info_cur: updated sensor/end device configuration
		m_gateway_cur: updated gateway connectivity
		N_kq_cur: updated global traffic
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

	m_gateway = np.zeros((sr_cnt, 1))
	gw_placed = G[:, 2]
	if np.sum(gw_placed) == 0:
		# No placed gateway
		return m_gateway

	# If there is placed gateway, update device configuration
	for gw_idx in range(gw_cnt):
		if not G[gw_idx, 2]: # No gateway is placed at this location
			continue

		# A gateway is placed at j	
		for sr_idx in range(sr_cnt):
			if propagation.GetRSSI(params.Ptx_max, PL[sr_idx, gw_idx]) >= params.RSSI_k[-1]:
				# If the RSSI under max tx power exceeds the minimum RSSI threshold,
				# we reckon this gateway has the probability of covering end device i 

				# Try to assign the minimum SF and channel based on the current sr_info
				# and collision-possible nodes N_kq
				succ, PDR_i, Lifetime_i = DeviceConfiguration(sr_info, \
					sr_info, G, PL, N_kq, m_gateway, params, sr_idx, gw_idx)
				logging.info('succ: {} PDR_i: {} Lifetime_i: {}'.format(succ, PDR_i, Lifetime_i))

	return m_gateway


def RGreedyAlg(sr_info_ogn, G_ogn, PL, dist, N_kq, params, GreedyParams):
	'''
	Call the robust gateway placement algorithm

	Args:
		sr_info: original read-only sensor placement and configuration
		G: original read-only gateway placement
		PL: path loss matrix between sensors and potential gateways
		dist: distance matrix between sensors and potential gateways
		N_kq: a dictionary recording traffic allocation
		params: important parameters
		GreedyParams: parameters for this greedy algorithms

	Returns:
		sr_info: sensor configuration
		G: resulted gateway placement
		m_gateway: current m-gateway connectivity
		N_kq: current traffic allocation
	'''
	# Make a deep copy of the original numpy array to avoid changes
	sr_info = np.copy(sr_info_ogn)
	G = np.copy(G_ogn)
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

	# Update current gateway connectivity at each end node, since there might be
	# gateway placement already in the G
	m_gateway = UpdateConn(sr_info, G, PL, N_kq, params)
	conn_cur = np.ones((sr_cnt, 1)) * params.M - m_gateway
	uncover_old = np.sum(conn_cur[conn_cur > 0])
	PDR_old = np.zeros((sr_cnt, 1))
	Lifetime_old = np.zeros((sr_cnt, 1))
	rounds = 0

	# Start greedily place gateway
	while True:

		# Variables to record the best gateway location in the current round
		bnft_best = -np.inf

		# If uncover count is lower than the threshold, fire end-of-exploration acceleration
		G_remain = np.copy(G)
		gw_mask = np.ones_like(G[:, 2], dtype=bool)
		PL_remain = np.copy(PL)
		if GreedyParams.end and uncover_old < params.M * sr_cnt * GreedyParams.end_thres:
			gw_mask = np.zeros_like(G[:, 2], dtype=bool) # Clear gw mask
			conn_cur = np.ones((sr_cnt, 1)) * params.M - m_gateway
			for j in range(gw_cnt):
				for i in range(sr_cnt):
					# Only keep the gateway that can cover uncovered end devices
					if conn_cur[i] > 0 and \
						propagation.GetRSSI(params.Ptx_max, PL[i, j]) >= params.RSSI_k[-1]:
						# If the RSSI under max tx power exceeds the minimum RSSI threshold,
						# we reckon this gateway has the probability of covering uncovered end devices 
						# print(propagation.GetRSSI(params.Ptx_max, PL[i, j]))
						gw_mask[j] = True
						break
			G_remain = G_remain[gw_mask, :]
			PL_remain = PL_remain[:, gw_mask]
		
		logging.info('Remained candidate sr: {} gw: {}'.format(sr_cnt, G_remain.shape[0]))
		
		# The following operation is run on sr_info, G_remain and PL_remain
		remained_cnt = G_remain.shape[0]
		for gw_idx in range(remained_cnt):
			if G_remain[gw_idx, 2]: # A gateway has been placed at this location
				continue

			# Try to place gateway at gateway location idx and configure the sensor
			# Create a deep copy of the original info for this gateway placement attempt
			G_remain[gw_idx, 2] = 1
			sr_info_cur = np.copy(sr_info)
			N_kq_cur = copy.deepcopy(N_kq)
			m_gateway_cur = np.copy(m_gateway)
			PDR_cur = np.copy(PDR_old)
			Lifetime_cur = np.copy(Lifetime_old)

			PL_gw_idx = PL_remain[:, gw_idx] # path loss to gateway idx
			sort_idx = np.argsort(PL_gw_idx, axis=0) # sort the PL from min to max

			for i in range(sr_cnt):
				sr_idx = int(sort_idx[i])

				if propagation.GetRSSI(params.Ptx_max, PL_remain[sr_idx, gw_idx]) < params.RSSI_k[-1]:
					# If the RSSI under max tx power will still below the minimum RSSI threshold
					# then the following sensors exceed the maximum communication range
					# therefore we do not need to consider them
					logging.debug('Exceeds maximum communication range!')
					break

				# Try to assign the minimum SF and channel based on the current sr_info
				# and collision-possible nodes N_kq
				succ, PDR_i, Lifetime_i = DeviceConfiguration(sr_info_cur, \
					sr_info, G_remain, PL_remain, N_kq_cur, m_gateway_cur, params, sr_idx, gw_idx)
				logging.debug('succ: {} PDR_i: {} Lifetime_i: {}'.format(succ, PDR_i, Lifetime_i))
				
				if not succ:
					# All assignment attempts are failed, no more assignment can be made
					break

				# Update the new PDR and lifetime after assignment
				PDR_cur[sr_idx] = PDR_i
				Lifetime_cur[sr_idx] = Lifetime_i
						
			# Calculate benefit of placing gateway at this location
			conn_cur = np.ones((sr_cnt, 1)) * params.M - m_gateway_cur
			uncover_new = np.sum(conn_cur[conn_cur > 0])
			# Benefit function: promote new connectivity and PDR, lifetime constraints
			bnft_cov = uncover_old - uncover_new
			bnft_pdr = GreedyParams.w_pdr * uncover_old * \
				np.sum(PDR_cur - np.ones((sr_cnt, 1)) * params.PDR_th)
			bnft_lifetime = GreedyParams.w_lifetime * uncover_old * \
				np.sum(Lifetime_cur - np.ones((sr_cnt, 1)) * params.Lifetime_th)
			bnft = bnft_cov + bnft_pdr + bnft_lifetime
			logging.info('cov: {} pdr: {} lifetime: {}'.format(bnft_cov, bnft_pdr, bnft_lifetime))

			# Update global best benefit value if necessary
			if bnft > bnft_best:
				bnft_best = bnft
				best_idx = gw_idx
				best_sr_info = sr_info_cur
				best_N_kq = N_kq_cur
				best_m_gateway = m_gateway_cur
				best_uncover = uncover_new
				best_PDR = PDR_cur
				best_Lifetime = Lifetime_cur

			# Logging
			logging.info('gw_loc: #{} {} Benefit: {} Max bnft: {} Max idx: {}'.format(\
				gw_idx, G[gw_idx, :2], bnft, bnft_best, best_idx))

			# Reset
			G_remain[gw_idx, 2] = 0

			# Early ending if already satisfy all uncovers
			#if uncover_new == 0:
			#	break

		# Check if there is no benefit to gain, end the searching while loop
		if bnft_best == 0:
			logging.info('No more gateway placement can provide m-gateway connectivity benefit!')
			break

		# Place a gateway at next_idx with the max benefit and update info
		G_remain[best_idx, 2] = 1
		G[gw_mask, :] = G_remain # Copy G_remain to G
		sr_info = best_sr_info
		N_kq = best_N_kq
		m_gateway = best_m_gateway
		uncover_old = best_uncover
		PDR_old = best_PDR
		Lifetime_old = best_Lifetime

		logging.info('Placed gateway #{} at grid {} [{},{}]'.format( \
			rounds, best_idx, G[best_idx, 0], G[best_idx, 1]))
		logging.info('Uncover: {}'.format(uncover_old))

		# Check if m-gateway connectivity has been met at all end nodes
		# If so, terminate the placement process
		if uncover_old <= 0:
			break

		rounds += 1

	return sr_info, G, m_gateway, N_kq