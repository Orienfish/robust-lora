#!/usr/bin/python3
import numpy as np
import math
# import matplotlib.pyplot as plt
import logging
import copy
from geneticalgorithm import geneticalgorithm as ga

import propagation

########################################
# Robust-Driven Genetic Alg
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
	return PDR, P_idxj

def fitness(X):
	'''
	Fitness function given the variables

	Attributes:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between end devices and candidate gateway locations
		params: important parameters
	'''

	sr_info = fitness.sr_info
	G = fitness.G
	PL = fitness.PL
	params = fitness.params
	GeneticParams = fitness.GeneticParams
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

	# Extract the corresponding variable from X and fill them to correct places
	G[:, 2] = X[:gw_cnt]
	sr_info[:, 2] = X[gw_cnt:gw_cnt+sr_cnt]
	Ptx_idx = X[gw_cnt+sr_cnt:gw_cnt+sr_cnt*2].astype(int)
	sr_info[:, 3] = list(map(lambda i: params.Ptx[i], Ptx_idx))
	sr_info[:, 4] = X[gw_cnt+sr_cnt*2:gw_cnt+sr_cnt*3]

	# Use a dictionary to record the list of nodes using the SFk and channel q
	# Note that the dictionary only records the primary connection
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)
	N_kq = dict()
	for k in range(SF_cnt):
		for q in range(CH_cnt):
			N_kq[str(k) + '_' + str(q)] = []

	# Fill in N_kq
	for idx in range(sr_cnt):
		k = int(sr_info[idx, 2])	# SFk
		q = int(sr_info[idx, 4])	# Channel q
		label = str(k) + '_' + str(q)
		N_kq[label].append(idx)

	PDR = np.zeros((sr_cnt, 1))
	Lifetime = np.zeros((sr_cnt, 1))
	m_gateway = np.zeros((sr_cnt, 1)) # Current gateway connectivity at each end node

	for idx in range(sr_cnt):
		PDR[idx], _ = GetPDR(sr_info, G, PL, N_kq, params, idx)
		k = int(sr_info[idx, 2])	# SFk
		# Find the index of closet Transmission power
		Ptx_idx = np.abs(params.Ptx - sr_info[idx, 3]).argmin()
		Ptx = params.Ptx[Ptx_idx]
		Lifetime[idx] = GetLifetime(k, Ptx, PDR[idx], params)

		# Update m gateway connectivity
		for j in range(gw_cnt):
			if not G[j, 2]: # No gateway has been placed at this location
				continue

			# If the expected RSSI meets the threshold, then we say one more gateway
			# connectivity is established
			if propagation.GetRSSI(Ptx, PL[idx, j]) > params.RSSI_k[k]:
				m_gateway[idx] += 1

	gw_cnt = np.sum(G[:, 2])
	cost_pdr = GeneticParams.w_pdr * \
		np.sum(np.ones((sr_cnt, 1)) * params.PDR_th - PDR)
	cost_lifetime = GeneticParams.w_lifetime * \
		np.sum(np.ones((sr_cnt, 1)) * params.Lifetime_th - Lifetime)
	cost_gateway = GeneticParams.w_gateway_conn * \
		np.sum(np.ones((sr_cnt, 1)) * params.M - m_gateway)
	cost = gw_cnt + cost_pdr + cost_lifetime + cost_gateway
	logging.info('gw_cnt: {} pdr: {} lifetime: {} m-gateway conn: {}'.format( \
		gw_cnt, cost_pdr, cost_lifetime, cost_gateway))

	return cost

def RGeneticAlg(sr_info, G, PL, dist, params, GeneticParams):
	'''
	Call the robust gateway placement algorithm

	Args:
		sr_info: original read-only sensor placement and configuration
		G: original read-only gateway placement
		PL: path loss matrix between sensors and potential gateways
		dist: distance matrix between sensors and potential gateways
		params: important parameters

	Returns:
		sr_info: sensor configuration
		G: resulted gateway placement
	'''
	# Set the attributes of the fitness funtion
	fitness.sr_info = sr_info
	fitness.G = G
	fitness.PL = PL
	fitness.dist = dist
	fitness.params = params
	fitness.GeneticParams = GeneticParams

	# Initialize variable constraints to input to GA
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

	di = gw_cnt + sr_cnt * 3
	#vartype_gw = np.array([['bool']] * gw_cnt)
	#vartype_dc = np.array([['int']] * sr_cnt * 3)
	#vartype = np.concatenate((vartype_gw, vartype_dc), axis=0)
	varbound_gw = np.array([[0, 1]] * gw_cnt)
	varbound_sf = np.array([[0, 3]] * sr_cnt)
	varbound_ptx = np.array([[0, 5]] * sr_cnt)
	varbound_ch = np.array([[0, 7]] * sr_cnt)
	varbound = np.concatenate((varbound_gw, varbound_sf, varbound_ptx, varbound_ch), axis=0)

	algorithm_param = {'max_num_iteration': GeneticParams.it,\
                   'population_size': GeneticParams.pop,\
                   'mutation_probability':0.1,\
                   'elit_ratio': 0.01,\
                   'crossover_probability': 0.5,\
                   'parents_portion': 0.3,\
                   'crossover_type':'uniform',\
                   'max_iteration_without_improv': None}
	X = np.ones((di))
	fitness(X)

	model=ga(function=fitness,\
	         dimension=di,\
	         variable_type='int',\
	         variable_boundaries=varbound,\
	         algorithm_parameters=algorithm_param)

	model.run()

	convergence = model.report
	solution = model.output_dict['variable']

	# Wrap up the solution
	G[:, 2] = solution[:gw_cnt]
	sr_info[:, 2] = solution[gw_cnt:gw_cnt+sr_cnt]
	Ptx_idx = solution[gw_cnt+sr_cnt:gw_cnt+sr_cnt*2].astype(int)
	sr_info[:, 3] = list(map(lambda i: params.Ptx[i], Ptx_idx))
	sr_info[:, 4] = solution[gw_cnt+sr_cnt*2:gw_cnt+sr_cnt*3]

	m_gateway = np.zeros((sr_cnt, 1)) # Current gateway connectivity at each end node

	# Update m gateway connectivity
	for idx in range(sr_cnt):
		k = int(sr_info[idx, 2])	# SFk
		Ptx_idx = np.abs(params.Ptx - sr_info[idx, 3]).argmin()
		Ptx = params.Ptx[Ptx_idx]	# Ptx

		for j in range(gw_cnt):
			if not G[j, 2]: # No gateway has been placed at this location
				continue

			# If the expected RSSI meets the threshold, then we say one more gateway
			# connectivity is established
			if propagation.GetRSSI(Ptx, PL[idx, j]) > params.RSSI_k[k]:
				m_gateway[idx] += 1

	return sr_info, G, m_gateway
