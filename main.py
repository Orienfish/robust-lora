#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import logging
import os
import time
logging.basicConfig(level=logging.INFO)

import ICIOT
import RGreedy
import propagation


########################################
# Important parameters
########################################
class params:
	L = 60000			# Edge of analysis area in m
	#N_x = 1000			# Number of sensor potential locations on x coordinate
	#N_y = 1000			# Number of sensor potential locations on y coordinate
	#Unit_sr = math.floor(L / (N_x-1)) # Unit length between gw grid points
	#G_x = math.floor(N_x * Unit_sr / Unit_gw)
	#G_y = math.floor(N_y * Unit_sr / Unit_gw)
	G_x = 12     	 	# Number of gateway potential locations on x coordinate
	G_y = 12			# Number of gateway potential locations on y coordinate
	Unit_gw = math.floor(L / G_x) # Unit length between gw grid points
	desired_gw_cnt = 10 # Desired gateways to place by ICIOT alg

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
	# Can be verified with toa.py
	AirTime_k = [0.098, 0.175, 0.329, 0.616]
	# SNR threshold of SF7, 8, 9, 10 in dB
	# Copied from ICDCS paper
	SNR_k = [-6, -9, -12, -15]

	# Transmission power options in dBm
	Ptx = [20, 17, 14, 11, 8, 5]
	Ptx_max = 20		# Maximum allowed transmission power in dBm
	Ptx_min = 5			# Minimum allowed transmission power in dBm
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

	PDR_th = 0.8		# PDR threshold at each end node
	Lifetime_th = 1		# Lifetime threshold at each end node in years
	w_pdr = 10
	w_lifetime = 1

# which algorithm to run
class run:
	iter = 1
	RGreedy = False
	RGenetic = False
	ICIOT = True

def init(sr_cnt, G, params):
	'''
	Initialize end device and path loss matrix

	Args:
		sr_cnt: number of end devices
		G: gateway locations
		params: important parameters

	Return:
		sr_info: end device settings
		PL: path loss matrix
		dist: distance matrix
	'''
	sr_info = []				# [x, y, SF, Ptx, CH]
	for i in range(sr_cnt):	
		k = -1 #random.randint(0, len(params.SF)-1) # SFk
		q = -1 # random.randint(0, len(params.CH)-1) # Channel q
		new_loc = [random.random() * params.L, random.random() * params.L, \
			k, params.Ptx_max, q]
		sr_info.append(new_loc)
	sr_info = np.array(sr_info)
	# print(sr_info)

	# Generate path loss and distance matrix between sensor i and gateway j
	gw_cnt = G.shape[0]
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

	return sr_info, PL, dist

def eval(sr_info_res, G_res, PL, params):
	'''
	Evaluate the PDR and lifetime of the solution

	Args:
		sr_info_res: generated sensor/end device configuration
		G_res: generated gateway placement
		PL: path loss matrix
		params: important parameters

	Return:
		PDR: an array of PDR at each end device
		PDR_gw: a matrix of PDR between end device-gateway pair
		lifetime: an array of lifetime at each end device 
	'''
	sr_cnt = sr_info_res.shape[0]
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)

	# First gather the number of nodes using the same SF and channel
	N_kq = dict()
	for k in range(SF_cnt): # Init N_kq
		for q in range(CH_cnt):
			N_kq[str(k) + '_' + str(q)] = []

	for idx in range(sr_cnt):		# Fill in N_kq
		k = int(sr_info_res[idx, 2])	# SFk
		q = int(sr_info_res[idx, 4])	# Channel q
		label = str(k) + '_' + str(q)
		N_kq[label].append(idx)

	PDR, PDR_gw, lifetime = [], [], []
	for idx in range(sr_cnt):
	
		newPDR, newPDR_gw = RGreedy.GetPDR(sr_info_res, G_res, PL, N_kq, params, idx)
		PDR.append(newPDR)
		PDR_gw.append(newPDR_gw)

		k = int(sr_info_res[idx, 2])	# SFk
		# Find the index of closet Transmission power
		Ptx_idx = np.abs(params.Ptx - sr_info_res[idx, 3]).argmin()
		lifetime.append(RGreedy.GetLifetime(k, params.Ptx[Ptx_idx], newPDR, params))

	print(np.array(PDR_gw))
	print(PDR)
	print(lifetime)
	return PDR, PDR_gw, lifetime

def plot(sr_info, G, version):
	# Visualize the placement and device configuration
	# sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	plt.figure()
	colorList = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
	color = [colorList[int(i)] for i in list(sr_info[:, 2])]
	plt.scatter(sr_info[:, 0], sr_info[:, 1], c=color , s=5)
	color = ['r' for i in range(gw_cnt)]
	plt.scatter(G[:, 0], G[:, 1], s=G[:, 2]*50, c=color, marker='^')
	plt.xlabel('X (m)'); plt.ylabel('Y (m)');
	# plt.legend()
	filename = 'vis_{}.png'.format(version)
	plt.savefig(filename)
	# plt.show()

def SaveInfo(sr_info, G, method):
	'''
	Save the generated solution to text file

	Args:
		sr_info: generated sensor/end device configuration solution
		G: generated gateway placement solution
		method: a string showing the algorithm, to be added to file name
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]

	# Write sensor and gateway information to file
	filename = 'sr_{}.txt'.format(method)
	with open (filename, 'w') as out:
		for i in range(sr_cnt):
			# Note: we need to convert SF to data rate (DR)
			# SF7 - DR3, SF8 - DR2, SF9 - DR1, SF10 - DR0
			out.write(str(round(sr_info[i, 0], 2)) + ' ' + \
				str(round(sr_info[i, 1], 2)) + ' ' + \
				str(int(3 - sr_info[i, 2])) + ' ' + \
				str(round(sr_info[i, 3])) + '\n')
	filename = 'gw_{}.txt'.format(method)
	with open (filename, 'w') as out:
		for i in range(gw_cnt):
			if G[i, 2]:
				out.write(str(round(G[i, 0], 2)) + ' ' + \
					str(round(G[i, 1], 2)) + '\n')

def SaveRes(sr_cnt, M, gw_cnt, time):
	# Log results
	with open('res.txt', 'a') as out:
		out.write(str(sr_cnt) + ' ' + str(M) + ' ' + \
			str(gw_cnt) + ' ' + str(time) + '\n')


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
			new_loc = [(p + 0.5) * params.Unit_gw, (q + 0.5) * params.Unit_gw, 0]
			G.append(new_loc)
	G = np.array(G)
	gw_cnt = params.G_x * params.G_y  # Number of gw potential locations

	# Randomly generate sensor positions
	#for sr_cnt in [100, 500, 1000, 5000]: # Number of sensors
	for sr_cnt in [100]:
		for it in range(run.iter):
			# Experiment iterations to evaluate diff random init
			logging.info('sr_cnt: {} iter: {}'.format(sr_cnt, it))

			sr_info, PL, dist = init(sr_cnt, G, params)

			if run.RGreedy:
				for M in [1]: #[3, 2, 1]:
					params.M = M
					logging.info('Running M = {}'.format(params.M))
					st_time = time.time()
					sr_info_res, G_res, m_gateway_res = \
						RGreedy.RGreedyAlg(sr_info, G, PL, dist, params)
					run_time = time.time() - st_time

					# Show m-gateway connectivity at each end device
					print(np.reshape(m_gateway_res, (1, -1)))

					# Print out PDR and lifetime at each end device
					eval(sr_info_res, G_res, PL, params)

					# Plot and log result
					plot(sr_info_res, G_res, 'R{}_{}'.format(params.M, sr_cnt))
					SaveRes(sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

					# Write sensor and gateway information to file
					method = 'RGreedy_{}_{}_{}'.format(sr_cnt, it, M)
					SaveInfo(sr_info_res, G_res, method)
				

			if run.ICIOT:
				st_time = time.time()
				sr_info_res, G_res = ICIOT.ICIOTAlg(sr_info, G, PL, params)
				run_time = time.time() - st_time

				# This paper assumes all end devices share one channel
				sr_info_res[:, 4] = np.ones((1, sr_cnt))
				print(sr_info_res)

				eval(sr_info_res, G_res, PL, params)

				# Plot result
				plot(sr_info_res, G_res, 'ICIOT')

				# Write sensor and gateway information to file
				method = 'ICIOT_{}_{}'.format(sr_cnt, it)
				SaveInfo(sr_info_res, G_res, method)


if __name__ == '__main__':
	main()
