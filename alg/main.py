#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import logging
import os
import time
import sys
import argparse
logging.basicConfig(level=logging.INFO)

import ICIOT
import RGreedy
import RGenetic
import propagation
import ReadData
import clustering
import optInterface

dataFile = None # './data/dataLA.csv'
origin = [33.5466, -118.7025]
PLFile = None #'./data/path_loss_mat.npy'


########################################
# Important parameters
########################################
class params:
	L = 30000			# Edge of analysis area in m
	sr_cnt = 100        # Number of end devices
	gw_dist = 6000      # Distance between two gateways in m

	# Version of log propagation model
	LogPropVer = 'ICIOT'

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
	bat_cap = 3 		# Battery capacity in Ah
	bat_volt = 3.3		# Battery supply voltage in V

	# Power of MCU
	P_MCU_off = 174.65e-6 # Power of MCU (Arduino Uno) in deep sleep in mW
	P_MCU_on = 23.48e-3 # Power of MCU (Arduino Uno) on in mW
	P_R_off = 1e-4		# Power of radio in deep sleep in W

	# Power of additive white Gaussian noise with zero-mean
	N0 = 1e-15 # in W

	PDR_th = 0.8		# PDR threshold at each end node
	Lifetime_th = 2		# Lifetime threshold at each end node in years

# Parameters for the DBSCAN clustering algorithm
class ClusterParams:
	eps = 5000			# Distance of neighborhood
	min_samples = 10	# How tolerant the algorithm is towards noise

# Parameters for the greedy algorithm
class GreedyParams:
	w_pdr = 1e-4	    # Weight for PDR
	w_lifetime = 0.5e-4	# Weight for lifetime
	cluster = False		# Whether use the clustering-based acceleration
	end = False			# Whether use the end-of-exploration acceleration
	end_thres = 0.3		# The threshold to kick off end-of-exploration acceleration

# Parameters for the ICIOT algorithm
class ICIOTParams:
	desired_gw_cnt = 12 # Number of desired gateways
	alpha = 1       	# Weight parameter in the objective function

# Parameters for the genetic algorithm
class GeneticParams:
	pop = 100			# Popupation count
	it = 1000			# Iteration count
	w_pdr = 0.5			# Weight for PDR
	w_lifetime = 0.05	# Weight for lifetime
	w_gateway_conn = 0.025 # Weight for m-gateway connectivity

# which algorithm to run
class run:
	iteration = 1
	M = [1] #[1, 2, 3] #[3, 2, 1]
	RGreedy = False  	# Pure greedy algorithm
	RGreedy_c = False	# With cluster-based acceleration
	RGreedy_e = False	# With end-of-exploration acceleration
	RGreedy_ce = False	# With both accleration techniques
	RGenetic = False
	ICIOT = True

def init(params):
	'''
	Initialize end device and path loss matrix

	Args:
		params: important parameters

	Return:
		sr_info: initial end device settings
		G: set of candidate gateway locations
		PL: path loss matrix
		dist: distance matrix
	'''
	sr_info = []				# [x, y, SF, Ptx, CH]
	# If dataFile is not provided, randomely generate
	# Else, read from the data file
	if dataFile == None:
		coor = np.random.rand(params.sr_cnt, 2) * params.L
		x_max = params.L
		y_max = params.L
	else:
		coor = ReadData.ReadFile(dataFile, origin)
		x_max = np.max(coor[:, 0])
		y_max = np.max(coor[:, 1])

	# Fill in the coordinate date to the initial sensor info
	for i in range(coor.shape[0]):	
		k = -1 #random.randint(0, len(params.SF)-1) # SFk
		q = -1 # random.randint(0, len(params.CH)-1) # Channel q
		new_loc = [coor[i, 0], coor[i, 1], k, params.Ptx_max, q]
		sr_info.append(new_loc)

	sr_info = np.array(sr_info)
	# print(sr_info)
	sr_cnt = sr_info.shape[0]
	

	# Generate the grid candidate set N and G with their x, y coordinates
	# N for sensor placement and G for gateway placement
	G = []				# [x, y, placed or not]
	x_gw = params.gw_dist / 2
	y_gw = params.gw_dist / 2
	while x_gw < x_max:
		while y_gw < y_max:
			new_loc = [x_gw, y_gw, 0]
			G.append(new_loc)
			y_gw += params.gw_dist
		x_gw += params.gw_dist
		y_gw = params.gw_dist / 2

	G = np.array(G)
	# print(G)
	gw_cnt = G.shape[0]

	# Generate path loss and distance matrix between sensor i and 
	# candidate gateway j
	PL = np.zeros((sr_cnt, gw_cnt))	
	dist = np.zeros((sr_cnt, gw_cnt))
	for i in range(sr_cnt):
		for j in range(gw_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			dist[i, j] = np.sqrt(np.sum((loc1 - loc2)**2))
			PL[i, j] = propagation.LogDistancePathLossModel(d=dist[i, j], \
				ver=params.LogPropVer)
	#np.savetxt('./data/PL_gen.csv', PL, delimiter=',')
	# If PL is provided, overwrite the isomophic one
	if PLFile != None:
		PL = np.load(PLFile).T
		logging.info('Load PL mat: {}'.format(PL.shape))
		#np.savetxt('./data/PL_import.csv', PL, delimiter=',')
		PL = PL + 20.0

	# Use a dictionary to record the list of nodes using the SFk and channel q
	# Note that the dictionary only records the primary connection
	# Initialize empty dictionary
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)
	N_kq = dict()
	for k in range(SF_cnt):
		for q in range(CH_cnt):
			N_kq[str(k) + '_' + str(q)] = []

	# Save the end device locations and candidate gateway locations
	np.savetxt('./relaxOpt/sr_loc.csv', sr_info[:, :2], delimiter=',')
	np.savetxt('./relaxOpt/gw_loc.csv', G[:, :2], delimiter=',')

	# optInterface.TestLifetime(params)
	c_ijks = optInterface.GenerateCijks(sr_info, G, PL, params)
	Ptx_cnt = len(params.Ptx)
	for k in range(SF_cnt):
		for q in range(Ptx_cnt):
			np.savetxt('./relaxOpt/cijk_{}_{}.csv'.format(k, q), \
				c_ijks[:, :, k, q], fmt='%d', delimiter=',')

	return sr_info, G, PL, dist, N_kq

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

	# Init N_kq: the number of nodes using the same SF and channel
	N_kq = dict()
	for k in range(SF_cnt): 
		for q in range(CH_cnt):
			N_kq[str(k) + '_' + str(q)] = []

	# Fill in N_kq
	for idx in range(sr_cnt):
		k = int(sr_info_res[idx, 2])	# SFk
		q = int(sr_info_res[idx, 4])	# Channel q
		if k != -1 and q != -1:
			label = str(k) + '_' + str(q)
			N_kq[label].append(sr_info_res[idx, :])

	PDR, PDR_gw, lifetime = [], [], []
	for idx in range(sr_cnt):
		k = int(sr_info_res[idx, 2])	# SFk
		q = int(sr_info_res[idx, 4])	# Channel q
		if k != -1 and q != -1:
			newPDR, newPDR_gw = RGreedy.GetPDR(sr_info_res, G_res, PL, N_kq, params, idx)
			PDR.append(newPDR)
			PDR_gw.append(newPDR_gw)

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
	colorList = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:gray']
	color = [colorList[int(i)] for i in list(sr_info[:, 2])]
	plt.scatter(sr_info[:, 0], sr_info[:, 1], c=color , s=5)
	color = ['r' for i in range(gw_cnt)]
	plt.scatter(G[:, 0], G[:, 1], s=G[:, 2]*50, c=color, marker='^')
	plt.xlabel('X (m)'); plt.ylabel('Y (m)');
	# plt.legend()
	filename = './vis/vis_{}.png'.format(version)
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
	filename = './res/sr_{}.txt'.format(method)
	with open (filename, 'w') as out:
		for i in range(sr_cnt):
			# Note: we need to convert SF to data rate (DR)
			# SF7 - DR3, SF8 - DR2, SF9 - DR1, SF10 - DR0
			out.write(str(round(sr_info[i, 0], 2)) + ' ' + \
				str(round(sr_info[i, 1], 2)) + ' ' + \
				str(int(3 - sr_info[i, 2])) + ' ' + \
				str(round(sr_info[i, 3])) + '\n')
	filename = './res/gw_{}.txt'.format(method)
	with open (filename, 'w') as out:
		for i in range(gw_cnt):
			if G[i, 2]:
				out.write(str(round(G[i, 0], 2)) + ' ' + \
					str(round(G[i, 1], 2)) + '\n')

def SaveRes(method, sr_cnt, M, gw_cnt, time):
	# Log results
	with open('./res/res.txt', 'a') as out:
		out.write(method + ' ' + str(sr_cnt) + ' ' + str(M) + ' ' + \
			str(gw_cnt) + ' ' + str(time) + '\n')


########################################
# Main Process
########################################
def main():
	# Process command line arguments
	parser = argparse.ArgumentParser(description='Run LoRa sensor deployment algs.')
	parser.add_argument("--dataFile", dest='dataFile', \
		help="path to data file of end devices locations")
	parser.add_argument("--PLFile", dest='PLFile', \
		help="path to PL file")
	parser.add_argument("--sr_cnt", dest='sr_cnt', type=int, \
		help="number of end devices")
	parser.add_argument("--desired_gw_cnt", dest='desired_gw_cnt', type=int, \
		help="desired number of gateways in the ICIOT algorithm")
	args = parser.parse_args()
	if args.sr_cnt:
		params.sr_cnt = args.sr_cnt
	elif args.desired_gw_cnt:
		ICIOTParams.desired_gw_cnt = args.desired_gw_cnt

	for it in range(run.iteration):
		# Initialization
		sr_info, G, PL, dist, N_kq = init(params)
		sr_cnt = sr_info.shape[0]
		gw_cnt = G.shape[0]
		logging.info('sr_cnt: {} gw_cnt: {}'.format(sr_cnt, gw_cnt))

		for M in run.M:
			params.M = M

			# Pure greedy algorithm
			if run.RGreedy:
				GreedyParams.end = False

				logging.info('Running Pure RGreedy M = {}'.format(params.M))
				st_time = time.time()
				sr_info_res, G_res, m_gateway_res, N_kq_res = \
					RGreedy.RGreedyAlg(sr_info, G, PL, dist, N_kq, params, GreedyParams)
				run_time = time.time() - st_time

				# Show m-gateway connectivity at each end device
				print(np.reshape(m_gateway_res, (1, -1)))

				# Print out PDR and lifetime at each end device
				eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				plot(sr_info_res, G_res, 'RGreedy_{}_{}_{}'.format(M, sr_cnt, it))
				SaveRes('G', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				method = 'RGreedy_{}_{}_{}'.format(M, sr_cnt, it)
				SaveInfo(sr_info_res, G_res, method)

			# Greedy algorithm with cluster-based acceleration
			if run.RGreedy_c:
				GreedyParams.cluster = True
				GreedyParams.end = False

				logging.info('Running RGreedy Cluster M = {}'.format(params.M))
				st_time = time.time()
				sr_info_res, G_res, m_gateway_res, N_kq_res = \
					clustering.RClusterAlg(sr_info, G, PL, dist, N_kq, params, ClusterParams, GreedyParams)
				run_time = time.time() - st_time

				# Show m-gateway connectivity at each end device
				print(np.reshape(m_gateway_res, (1, -1)))

				# Print out PDR and lifetime at each end device
				eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				plot(sr_info_res, G_res, 'RGreedyc_{}_{}_{}'.format(M, sr_cnt, it))
				SaveRes('Gc', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				method = 'RGreedyc_{}_{}_{}'.format(M, sr_cnt, it)
				SaveInfo(sr_info_res, G_res, method)

			# Greedy algorithm with end-of-exploration acceleration
			if run.RGreedy_e:
				GreedyParams.end = True

				logging.info('Running RGreedy End M = {}'.format(params.M))
				st_time = time.time()
				sr_info_res, G_res, m_gateway_res, N_kq_res = \
					RGreedy.RGreedyAlg(sr_info, G, PL, dist, N_kq, params, GreedyParams)
				run_time = time.time() - st_time

				# Show m-gateway connectivity at each end device
				print(np.reshape(m_gateway_res, (1, -1)))

				# Print out PDR and lifetime at each end device
				eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				plot(sr_info_res, G_res, 'RGreedye_{}_{}_{}'.format(M, sr_cnt, it))
				SaveRes('Ge', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				method = 'RGreedye_{}_{}_{}'.format(M, sr_cnt, it)
				SaveInfo(sr_info_res, G_res, method)

			# Greedy algorithm with both acceleration techniques
			if run.RGreedy_ce:
				GreedyParams.cluster = True
				GreedyParams.end = True

				logging.info('Running RGreedy Cluster+End M = {}'.format(params.M))
				st_time = time.time()
				sr_info_res, G_res, m_gateway_res, N_kq_res = \
					clustering.RClusterAlg(sr_info, G, PL, dist, N_kq, params, ClusterParams, GreedyParams)
				run_time = time.time() - st_time

				# Show m-gateway connectivity at each end device
				print(np.reshape(m_gateway_res, (1, -1)))

				# Print out PDR and lifetime at each end device
				eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				plot(sr_info_res, G_res, 'RGreedyce_{}_{}_{}'.format(M, sr_cnt, it))
				SaveRes('Gce', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				method = 'RGreedyce_{}_{}_{}'.format(M, sr_cnt, it)
				SaveInfo(sr_info_res, G_res, method)


			if run.RGenetic:
				logging.info('Running RGenetic M = {}'.format(params.M))
				st_time = time.time()
				sr_info_res, G_res, m_gateway_res = \
					RGenetic.RGeneticAlg(sr_info, G, PL, dist, params, GeneticParams)
				run_time = time.time() - st_time

				# Show m-gateway connectivity at each end device
				print(np.reshape(m_gateway_res, (1, -1)))

				# Print out PDR and lifetime at each end device
				eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				plot(sr_info_res, G_res, 'RGenetic{}_{}_{}'.format(M, sr_cnt, it))
				SaveRes('Genetic', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				method = 'RGenetic_{}_{}_{}'.format(M, sr_cnt, it)
				SaveInfo(sr_info_res, G_res, method)
					

		if run.ICIOT:
			st_time = time.time()
			sr_info_res, G_res = ICIOT.ICIOTAlg(sr_info, G, PL, params, ICIOTParams)
			run_time = time.time() - st_time

			# This paper assumes all end devices share one channel
			# We randomly allocate channels to make it the same as our assumption
			for i in range(sr_cnt):
				sr_info_res[i, 4] = random.randint(0, len(params.CH)-1)

			# Print out PDR and lifetime at each end device
			eval(sr_info_res, G_res, PL, params)

			# Plot result
			plot(sr_info_res, G_res, 'ICIOT_{}_{}'.format(sr_cnt, it))
			SaveRes('ICIOT', sr_cnt, 1, np.sum(G_res[:, 2]), run_time)

			# Write sensor and gateway information to file
			method = 'ICIOT_{}_{}'.format(sr_cnt, it)
			SaveInfo(sr_info_res, G_res, method)


if __name__ == '__main__':
	main()
