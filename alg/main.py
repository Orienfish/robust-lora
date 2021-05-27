#!/usr/bin/python3
# main file
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
import propagation
import optInterface
import utils

dir_path = os.path.dirname(os.path.realpath(__file__))


########################################
# Important parameters
########################################
class params:
	# indicate whch dataset to use, default random field
	dataset = 'random'

	if dataset == 'random':
		datapath = ''
		L = 20000  # Edge of analysis area in m, use if dataFile not provided
		sr_cnt = 5
		gw_dist = 5000
	elif dataset == 'HPWREN':
		datapath = '/HPWREN-dataset/'  # directory name with dataset
		sr_cnt = 1300  # Number of end devices
		gw_dist = 7250  # Distance between two gateways in m
	elif dataset == 'LA':
		datapath = '/LA-dataset/'
		sr_cnt = 100
		gw_dist = 6000
	else:
		print('Invalid dataset!')

	T = 1200  # Sampling period in s

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
	Ptx_max = 20  # Maximum allowed transmission power in dBm
	Ptx_min = 5  # Minimum allowed transmission power in dBm
	# Power consumption using each transmission power options in W
	# Copied from the TOSN paper (Liando 2019) for SX1276 chipset
	PowerTx = [0.4, 0.4, 0.3, 0.25, 0.2, 0.15]

	# Channel options
	CH = [0, 1, 2, 3, 4, 5, 6, 7]

	# Redundancy level
	# K = 2   			# Coverage level
	M = 2  # Connectivity level

	# Battery
	bat_cap = 3  # Battery capacity in Ah
	bat_volt = 3.3  # Battery supply voltage in V

	# Power of MCU
	P_MCU_off = 174.65e-6  # Power of MCU (Arduino Uno) in deep sleep in W
	P_MCU_on = 23.48e-3  # Power of MCU (Arduino Uno) on in W
	P_R_off = 1e-4  # Power of radio in deep sleep in W

	# Power of additive white Gaussian noise with zero-mean
	N0 = 1e-15  # in W

	PDR_th = 0.8  # PDR threshold at each end node
	Lifetime_th = 2  # Lifetime threshold at each end node in years

	LogPropVer = 'Dongare'  # Version of log propagation model


class DataParams:
	dataLoc = False  # Whether to use the predetermined locations
	PL = False  # Whether to use the PLFile of path loss matrix
	LogPropVer = 'Dongare'  # Version of log propagation model
	datasetPath = dir_path + '/../data' + params.datapath  # Path to the dataset
	srFile = datasetPath + 'sr_loc.csv'  # End device locations
	gwFile = datasetPath + 'gw_loc.csv'  # Candidate gateway locations
	PLFile = datasetPath + 'path_loss_mat.npy'  # Path loss between each device-gw pair
	GwAbleFile = datasetPath + 'gw_able.npy'  # Whether placing gateway at a location is allowed


# Parameters for the DBSCAN clustering algorithm
class ClusterParams:
	eps = 5000  # Distance of neighborhood
	min_samples = 10  # How tolerant the algorithm is towards noise


# Parameters for the greedy algorithm
class GreedyParams:
	w_pdr = 1e-4  # Weight for PDR
	w_lifetime = 0.5e-4  # Weight for lifetime
	cluster = False  # Whether use the clustering-based acceleration
	end = False  # Whether use the end-of-exploration acceleration
	end_thres = 0.3  # The threshold to kick off end-of-exploration acceleration


# Parameters for the ICIOT algorithm
class ICIOTParams:
	desired_gw_cnt = 2  # Number of desired gateways
	alpha = 1  # Weight parameter in the objective function


# Parameters for the genetic algorithm
class GeneticParams:
	pop = 100  # Popupation count
	it = 1000  # Iteration count
	w_pdr = 0.5  # Weight for PDR
	w_lifetime = 0.05  # Weight for lifetime
	w_gateway_conn = 0.025  # Weight for m-gateway connectivity


# which algorithm to run
class run:
	iteration = 1
	M = [1] # [3, 2, 1]
	RGreedy = True  # Pure greedy algorithm
	RGreedy_c = False  # With cluster-based acceleration
	RGreedy_e = False  # With end-of-exploration acceleration
	RGreedy_ce = False  # With both accleration techniques
	RGenetic = False
	ICIOT = True


def init(params):
	"""
	Initialize end device and path loss matrix

	Args:
		params: important parameters

	Return:
		sr_info: initial end device settings
		G: set of candidate gateway locations
		PL: path loss matrix
		dist: distance matrix
	"""
	#####################################################################
	# Initialize sensor end device information
	# If predetermined locations are provided, read from the dataLoc file
	# Else, randomly generate the end devices locations
	#####################################################################
	sr_info = []		# [x, y, SF, Ptx, CH]
	if DataParams.dataLoc:
		coor = np.loadtxt(DataParams.srFile, delimiter=',')
	else:
		coor = np.random.rand(params.sr_cnt, 2) * params.L
		x_max = params.L
		y_max = params.L

	# Fill in the sensor coordinates to the initial sensor info
	for i in range(coor.shape[0]):
		k = -1  # random.randint(0, len(params.SF)-1) # SFk
		q = -1  # random.randint(0, len(params.CH)-1) # Channel q
		new_loc = [coor[i, 0], coor[i, 1], k, params.Ptx_max, q]
		sr_info.append(new_loc)

	sr_info = np.array(sr_info)
	sr_cnt = sr_info.shape[0]

	#####################################################################
	# Initialize candidate gateway information
	# If predetermined locations are provided, read from the dataLoc file
	# Else, randomly generate the end devices locations
	#####################################################################
	G = []  # [x, y, placed or not, can place or not]
	if DataParams.dataLoc:
		coor = np.loadtxt(DataParams.gwFile, delimiter=',')
		gw_able = np.load(DataParams.GwAbleFile)
	else:
		x_gw = params.gw_dist / 2
		y_gw = params.gw_dist / 2
		gw_idx = 0
		coor = []
		while x_gw < x_max:
			while y_gw < y_max:
				coor.append([x_gw, y_gw])
				gw_idx += 1
				y_gw += params.gw_dist
			x_gw += params.gw_dist
			y_gw = params.gw_dist / 2
		coor = np.array(coor)
		gw_able = np.ones((coor.shape[0],))

	# Fill in the gateway coordinates to the initial gateway info
	for i in range(coor.shape[0]):
		new_loc = [coor[i, 0], coor[i, 1], 0, gw_able[i]]
		G.append(new_loc)

	G = np.array(G)
	gw_cnt = G.shape[0]

	#####################################################################
	# Initialize path loss matrix
	# If predetermined path losses are provided, read from the file
	# Else, randomly generate the path loss using the specified model
	#####################################################################
	PL = np.zeros((sr_cnt, gw_cnt))
	dist = np.zeros((sr_cnt, gw_cnt))
	for i in range(sr_cnt):
		for j in range(gw_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			dist[i, j] = np.sqrt(np.sum((loc1 - loc2) ** 2))
			PL[i, j] = propagation.LogDistancePathLossModel(d=dist[i, j],
															ver=DataParams.LogPropVer)
	# If PL file is provided, overwrite the isomophic one
	if DataParams.PL:
		PL = np.load(DataParams.PLFile).T
		logging.info('Load PL mat: {}'.format(PL.shape))
		PL = PL + 10.0

	# Use a dictionary to record the list of nodes using the SFk and channel q
	# Note that the dictionary only records the primary connection
	# Initialize empty dictionary
	SF_cnt = len(params.SF)
	CH_cnt = len(params.CH)
	N_kq = dict()
	for k in range(SF_cnt):
		for q in range(CH_cnt):
			N_kq[str(k) + '_' + str(q)] = []

	#####################################################################
	# Save information for relaxed problem optimization
	#####################################################################
	# Save the information as input to the relaxed problem optimization
	np.savetxt(dir_path + '/../opt/sr_loc.csv', sr_info[:, :2], delimiter=',')
	np.savetxt(dir_path + '/../opt/gw_loc.csv', G[:, :2], delimiter=',')
	np.savetxt(dir_path + '/../opt/pl.csv', PL, delimiter=',')

	# Generate the binary indicator for relaxed optimization
	sigma = propagation.LogDistancePathLossModel.sigma
	c_ijks = optInterface.GenerateCijks(sr_info, G, PL, sigma, params)
	Ptx_cnt = len(params.Ptx)
	for k in range(SF_cnt):
		for q in range(Ptx_cnt):
			np.savetxt(dir_path + '/../opt/cijk_{}_{}.csv'.format(k, q),
					   c_ijks[:, :, k, q], fmt='%d', delimiter=',')

	return sr_info, G, PL, dist, N_kq


########################################
# Main Process
########################################
def main():
	# Process command line arguments
	parser = argparse.ArgumentParser(description='Run LoRa sensor deployment algs.')
	parser.add_argument("--data", dest='data', nargs='?', const=True, default=False, type=bool,
						help="whether to use given dataFile of end devices locations")
	parser.add_argument("--PL", dest='PL', nargs='?', const=True, default=False, type=bool,
						help="whether to use given PLFile of path loss matrix")
	parser.add_argument("--sr_cnt", dest='sr_cnt', type=int,
						help="number of end devices")
	parser.add_argument("--RGreedy", dest='RGreedy', nargs='?', const=True, default=False, type=bool,
						help="whether to run the RGreedy alg")
	parser.add_argument("--ICIOT", dest='ICIOT', nargs='?', const=True, default=False, type=bool,
						help="whether to run the ICIOT alg")
	parser.add_argument("--desired_gw_cnt", dest='desired_gw_cnt', type=int,
						help="desired number of gateways in the ICIOT algorithm")
	args = parser.parse_args()
	if args.data:
		DataParams.dataLoc = args.data
	if args.PL:
		DataParams.PL = args.PL
	if args.sr_cnt:
		params.sr_cnt = args.sr_cnt
	if args.RGreedy:
		run.RGreedy = args.RGreedy
	if args.ICIOT:
		run.ICIOT = args.ICIOT
	if args.desired_gw_cnt:
		ICIOTParams.desired_gw_cnt = args.desired_gw_cnt

	# Create flag w.r.t. file usage status for labeling results
	flagData = 'd' if DataParams.dataLoc else ''
	flagPL = 'p' if DataParams.PL else ''

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
				utils.eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				method = 'RGreedy_{}_{}_{}{}{}'.format(M, sr_cnt, it, flagData, flagPL)
				utils.plot(sr_info_res, G_res, method)
				utils.SaveRes('G', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				utils.SaveInfo(sr_info_res, G_res, PL, method, params, DataParams)

				# If running ICIOT after, set the gateway number to the same
				if run.ICIOT:
					ICIOTParams.desired_gw_cnt = int(np.sum(G_res[:, 2]))

			# Greedy algorithm with cluster-based acceleration
			if run.RGreedy_c:
				import clustering
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
				utils.eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				method = 'RGreedyc_{}_{}_{}{}{}'.format(M, sr_cnt, it, flagData, flagPL)
				utils.plot(sr_info_res, G_res, method)
				utils.SaveRes('Gc', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				utils.SaveInfo(sr_info_res, G_res, PL, method, params, DataParams)

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
				utils.eval(sr_info_res, G_res, PL, params, DataParams)

				# Plot and log result
				method = 'RGreedye_{}_{}_{}{}{}'.format(M, sr_cnt, it, flagData, flagPL)
				utils.plot(sr_info_res, G_res, method)
				utils.SaveRes('Ge', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				utils.SaveInfo(sr_info_res, G_res, PL, method, params, DataParams)

			# Greedy algorithm with both acceleration techniques
			if run.RGreedy_ce:
				import clustering
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
				utils.eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				method = 'RGreedyce_{}_{}_{}{}{}'.format(M, sr_cnt, it, flagData, flagPL)
				utils.plot(sr_info_res, G_res, method)
				utils.SaveRes('Gce', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				utils.SaveInfo(sr_info_res, G_res, PL, method, params, DataParams)

			if run.RGenetic:
				import RGenetic as RGenetic
				logging.info('Running RGenetic M = {}'.format(params.M))
				st_time = time.time()
				sr_info_res, G_res, m_gateway_res = \
					RGenetic.RGeneticAlg(sr_info, G, PL, dist, params, GeneticParams)
				run_time = time.time() - st_time

				# Show m-gateway connectivity at each end device
				print(np.reshape(m_gateway_res, (1, -1)))

				# Print out PDR and lifetime at each end device
				utils.eval(sr_info_res, G_res, PL, params)

				# Plot and log result
				method = 'RGenetic_{}_{}_{}{}{}'.format(M, sr_cnt, it, flagData, flagPL)
				utils.plot(sr_info_res, G_res, method)
				utils.SaveRes('Genetic', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file
				utils.SaveInfo(sr_info_res, G_res, PL, method, params, DataParams)

			if run.ICIOT:
				st_time = time.time()
				sr_info_res, G_res = ICIOT.ICIOTAlg(sr_info, G, PL, params, ICIOTParams)
				run_time = time.time() - st_time

				# This paper assumes all end devices share one channel
				# We randomly allocate channels to make it the same as our assumption
				for i in range(sr_cnt):
					sr_info_res[i, 4] = random.randint(0, len(params.CH) - 1)

				# Print out PDR and lifetime at each end device
				utils.eval(sr_info_res, G_res, PL, params)

				# Plot result
				method = 'ICIOT_{}_{}_{}{}{}'.format(M, sr_cnt, it, flagData, flagPL)
				utils.plot(sr_info_res, G_res, method)
				utils.SaveRes('ICIOT', sr_cnt, params.M, np.sum(G_res[:, 2]), run_time)

				# Write sensor and gateway information to file

				utils.SaveInfo(sr_info_res, G_res, PL, method, params, DataParams)


if __name__ == '__main__':
	main()
