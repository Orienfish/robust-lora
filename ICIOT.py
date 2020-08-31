#!/usr/bin/python3
import numpy as np
import math
import matplotlib.pyplot as plt
import logging

import propagation

########################################
# ICIOT Alg
########################################
def GetCij(sr_info, G, PL, params):
	'''
	Get important notation Cij in the ICIOT paper showing the connection

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		params: important parameters

	Return:
		Cij: notation of whether the connection between sensor i and 
			gateway j is connectable
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	Cij = np.zeros((sr_cnt, gw_cnt))
	for j in range(gw_cnt):
		if not G[j, 2]: # no gateway is placed at j
			continue
		# a gateway is placed at j
		for i in range(sr_cnt):
			Prx = propagation.GetRSSI(sr_info[i, 3], PL[i, j])
			if Prx >= params.RSSI_k[int(sr_info[i, 2])]: # can be received
				Cij[i, j] = 1

	return Cij

def GetPDR(sr_info, G, Cij, params):
	'''
	Get packet delivery ratio
	
	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		Cij: connection indicator
		params: important parameters

	Return:
		PDR: a vector shows the PDR at each sensor i
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	SF_cnt = len(params.SF)
	
	# Calculate traffic load between sensor i and gw j, lambda_ij
	lambda_ij = np.zeros((sr_cnt, gw_cnt))
	N_jk = np.zeros((gw_cnt, SF_cnt)) # Number of nodes connecting to gw j
									  # and using the SFk
	for j in range(gw_cnt):
		# Update number of nodes that might collide
		for i in range(sr_cnt):
			if not Cij[i, j]: # Save computation since Cij is sparse
				continue
			k = int(sr_info[i, 2]) # Get SFk at sensor i
			N_jk[j, k] += Cij[i, j]
		# Calculate traffic load 
		for i in range(sr_cnt):
			if not Cij[i, j]:
				continue
			k = int(sr_info[i, 2]) # Get SFk at sensor i
			lambda_ij[i, j] = N_jk[j, k] * params.AirTime_k[k] / params.T
		# print(j, N_jk[j, :])

	# Calculate packet delivery ratio between sensor i and gw j, pi_ij
	pi_ij = np.multiply(Cij, np.exp(-2 * lambda_ij))

	# Calculate total PDR at sensor i
	PDR = 1 - np.prod(1 - pi_ij, axis=1)
	
	return PDR

def GetEnergyPerPacket(sr_info, params):
	'''
	Calculate the energy consuming in sending one packet

	Args:
		sr_info: sensor placement and configuration
		params: important parameters

	Return:
		ei: energy consumption per packet at sensor i in mJ = mW * s
	'''
	sr_cnt = len(sr_info)
	pi_mW = 10 ** (sr_info[:, 3] / 10)
	AirTime_i = [params.AirTime_k[int(sr_info[i, 2])] for i in range(sr_cnt)]
	AirTime_i = np.array(AirTime_i)
	ei = np.multiply(pi_mW, AirTime_i)
	return np.array(ei)

def DeviceConfiguration(sr_info, G, PL, params):
	'''
	The hybrid algorithm in the ICIOT paper to configure SF and transmission
	power of each end device

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		params: important parameters

	Return:
		sr_info_new: new device configuration
	'''
	# Calculate the number of end devices using each SF
	sr_cnt = sr_info.shape[0]
	N_k = sr_cnt / np.sum(1 / np.array(params.AirTime_k)) / params.AirTime_k
	N_k = np.rint(N_k)
	# print(N_k)
	# Note the sum of N_k might not equal to sr_cnt, but this is fine

	min_dist = np.ones((sr_cnt, 1)) * np.inf # the min distance of each 
											 # end device to any deployed gateway
	gw_link = np.empty((sr_cnt, 1)) # associated gateway with each end device
	# only keep the deployed gateways' row/column in G and PL
	PL = PL[:, G[:, 2] == 1]
	G = G[G[:, 2] == 1, :] 
	gw_cnt = G.shape[0]	# number of deployed gateway
	for i in range(sr_cnt):
		for j in range(gw_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			cur_dist = np.sqrt(np.sum((loc1 - loc2)**2))
			if cur_dist < min_dist[i]:
				min_dist[i] = cur_dist
				gw_link[i] = j

	# Get the index of distance sorting
	min_index = np.argsort(min_dist, axis=0)
	# print(min_index)

	# Assign SFs using EquiP strategy from the closest gateway
	k = 0 # SFk
	for i in range(sr_cnt):
		if N_k[k] <= 0:
			k += 1
		N_k[k] -= 1
		sr_index = int(min_index[i])
		gw_index = int(gw_link[sr_index])
		sr_info[sr_index, 2] = k
		pi = params.RSSI_k[k] + PL[sr_index][gw_index]
		if pi <= params.Ptx_max: # if transmission power does not violate
			sr_info[sr_index, 3] = max(pi, 0.0)
			continue
		# if transmission power violates the upperbound
		# print('tx violate')
		RSSI_max = params.Ptx_max - PL[sr_index][gw_index]
		try: # try to assign the minimum possible SF
			newRSSI = list(filter(lambda k: k <= RSSI_max, params.RSSI_k))[0]
			newSFk = RSSI_k.index(newRSSI)
			newpi = RSSI_k[newSFk] + PL[sr_index][gw_index]
			logging.debug("Reassign for tx pow violation: \
				old SF:{} pow:{} new SF:{} pow:{}".format( \
				k, pi, newSFk, newpi))
			sr_info[sr_index, 2] = newSFk
			sr_info[sr_index, 3] = newpi
		except: # if no SF could work, assign the largest SF
			sr_info[sr_index, 2] = len(params.SF) - 1
			sr_info[sr_index, 3] = params.Ptx_max
			logging.debug("Reassign failed, RSSI_max: {}".format(RSSI_max))

	return sr_info

def ICIOTAlg(sr_info, G, PL, params):
	'''
	Call the ICIOT gateway placement algorithm

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		params: important parameters

	Returns:
		gw_place: a binary vector of gateway placement decision
		sr_info: sensor configuration
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	# Start the greedy gateway placement algorithm
	for rounds in range(params.desired_gw_cnt):
		obj_old = -np.inf
		next_idx = -1
		next_sr_info = None
		for idx in range(gw_cnt):
			if G[idx, 2]: # a gateway has been placed at this location
				continue
			# try to place gateway at this location
			logging.debug("Try to place a gateway at {}".format(idx))
			G[idx, 2] = 1

			# assign powers and SFs
			sr_info = DeviceConfiguration(sr_info, G, PL, params)

			# Calculate Cij from each sensor i to gw j
			Cij = GetCij(sr_info, G, PL, params)
			logging.debug("Cij: {}".format(Cij))

			# Calculate PDR at each sensor i
			PDR = GetPDR(sr_info, G, Cij, params)
			logging.debug("PDR: {}".format(PDR))

			# Calculate energy per packet at each sensor i
			ei = GetEnergyPerPacket(sr_info, params)
			logging.debug("ei: {}".format(ei))

			# Calculate energy efficiency
			EE = np.divide(PDR, ei)
			logging.debug("EE: {}".format(EE))

			# Calculate objective value
			gw_place = G[:, 2] # a binary vector indicating gw placement
			alpha = 0.05       # weight parameter
			obj = np.sum(EE) / sr_cnt - alpha * np.sum(gw_place) / params.desired_gw_cnt
			logging.debug("obj1: {}".format(np.sum(EE) / sr_cnt))
			logging.debug("obj2: {}".format(alpha * np.sum(gw_place) / gw_cnt))
			logging.debug("obj: {}".format(obj))

			# plot(sr_info, G)

			# update objective value if necessary
			if obj > obj_old:
				obj_old = obj
				next_idx = idx
				next_sr_info = np.copy(sr_info)

			# reset
			G[idx, 2] = 0

		# place a gateway at next_idx with the max new objective value
		G[next_idx, 2] = 1
		sr_info = np.copy(next_sr_info)
		logging.info("Placed gateway #{} at grid {} [{},{}]".format( \
			rounds, next_idx, G[next_idx, 0], G[next_idx, 1]))

	gw_place = G[:, 2]
	return gw_place, sr_info