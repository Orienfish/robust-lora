#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

########################################
# Important parameters
########################################
#N_x = 1000
#N_y = 1000
#N_cnt = N_x * N_y
#Unit_sr = 50
L = 50000			# Edge of analysis area
#G_x = math.floor(N_x * Unit_sr / Unit_gw)
#G_y = math.floor(N_y * Unit_sr / Unit_gw)
G_x = 6				# Number of gw potential locations on x coordinate
G_y = 6				# Number of gw potential locations on y coordinate
Unit_gw = math.floor(L / G_x) # Unit length between gw grid points
# print(Unit_gw)

T = 1200			# Sampling period in s
Ptx_max = 23		# Maximum allowed transmission power in dBm
alpha = 0.05		# weight in objective value calculation

# Spreading factors
SF = ['SF7', 'SF8', 'SF9', 'SF10']
# RSSI threshold of SF7, 8, 9, 10 in dBm
# Copied from ICIOT paper and ICDCS paper
RSSI_k = [-123, -126, -129, -132]
# Air time of SF7, 8, 9, 10 in s
# Copied from ICIOT paper under 50 Bytes payload
AirTime_k = [0.098, 0.175, 0.329, 0.616]


def FreeSpacePathLossModel(d, f):
	'''
	Free Space Path Loss Model

	Args:
		d: distance in m
		f: frequency in MHz

	Return:
		path loss in dB at d
	'''
	PL = 20 * math.log10(d) + 20 * math.log10(f) - 27.55
	return PL

def LogDistancePathLossModel(d):
	'''
	Log Distance Path Loss Model

	Args:
		d: distance in m

	Return:
		path loss in dB at d
	'''
	PL = LogDistancePathLossModel.PL0 + 10 * LogDistancePathLossModel.delta * \
		math.log10(d / LogDistancePathLossModel.d0)
	return PL

# Settings in the ICIOT paper
LogDistancePathLossModel.d0 = 1000 # Reference distance in m
LogDistancePathLossModel.PL0 = 130 # Reference path loss in dB
LogDistancePathLossModel.delta = 2.1 # Path loss exponent
# Settings in the IPSN paper (indoor building)
#LogDistancePathLossModel.d0 = 40 # Reference distance in m
#LogDistancePathLossModel.PL0 = 127.41 # Reference path loss in dB
#LogDistancePathLossModel.delta = 3.57 # Path loss exponent

def GetRSSI(Ptx, PL):
	'''
	Calculate the RSSI at the receiver

	Args:
		Ptx: transmission power in dBm
		PL: path loss in dB

	Return:
		Prx: receiving power in dBm
	'''
	Prx = Ptx + GetRSSI.Gtx + GetRSSI.Grx - PL
	return Prx

GetRSSI.Gtx = 0 # Transmission antenna gain
GetRSSI.Grx = 0 # Reception antenna gain

# Test and tune path loss model's parameter
#d = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
#loss_log = [LogDistancePathLossModel(di) for di in d]
#Prx_log = [GetRSSI(20, lossi) for lossi in loss_log]
#loss_friis = [FreeSpacePathLossModel(di, 868) for di in d]
#Prx_friis = [GetRSSI(20, lossi) for lossi in loss_friis]
#plt.figure()
#plt.plot(d, Prx_log, d, Prx_friis)
#plt.show()

########################################
# ICIOT Alg
########################################
def GetCij(sr_info, G, PL):
	'''
	Get important notation Cij in the ICIOT paper showing the connection

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways

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
			Prx = GetRSSI(sr_info[i, 3], PL[i, j])
			if Prx > RSSI_k[int(sr_info[i, 2])]: # can be received
				Cij[i, j] = 1

	return Cij

def GetPDR(sr_info, G, Cij):
	'''
	Get packet delivery ratio
	
	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		Cij: connection indicator

	Return:
		PDR: a vector shows the PDR at each sensor i
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	SF_cnt = len(SF)
	
	# Calculate traffic load between sensor i and gw j, lambda_ij
	lambda_ij = np.zeros((sr_cnt, gw_cnt))
	N_jk = np.zeros((gw_cnt, SF_cnt)) # Number of nodes connecting to gw j
									  # and using the SFk
	for j in range(gw_cnt):
		# Update number of nodes that might collide
		for i in range(sr_cnt):
			k = int(sr_info[i, 2]) # get SFk at sensor i
			N_jk[j, k] += Cij[i, j]
		# Calculate traffic load 
		for i in range(sr_cnt):
			k = int(sr_info[i, 2]) # get SFk at sensor i
			lambda_ij[i, j] = N_jk[j, k] * AirTime_k[k] / T
		# print(j, N_jk[j, :])

	# Calculate packet delivery ratio between sensor i and gw j, pi_ij
	pi_ij = np.multiply(Cij, np.exp(-2 * lambda_ij))
	# print(pi_ij.shape)

	# Calculate total PDR at sensor i
	PDR = 1 - np.prod(1 - pi_ij, axis=1)
	
	return PDR

def GetEnergyPerPacket(sr_info):
	'''
	Calculate the energy consuming in sending one packet

	Args:
		sr_info: sensor placement and configuration

	Return:
		ei: energy consumption per packet at sensor i in mJ = mW * s
	'''
	sr_cnt = len(sr_info)
	pi_mW = 10 ** (sr_info[:, 3] / 10)
	AirTime_i = [AirTime_k[int(sr_info[i, 2])] for i in range(sr_cnt)]
	AirTime_i = np.array(AirTime_i)
	ei = np.multiply(pi_mW, AirTime_i)
	return np.array(ei)

def DeviceConfiguration(sr_info, G, PL, Ptx_max):
	'''
	The hybrid algorithm in the ICIOT paper to configure SF and transmission
	power of each end device

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		Ptx_max: maximum allowed transmission power setting

	Return:
		sr_info_new: new device configuration
	'''
	# Calculate the number of end devices using each SF
	sr_cnt = sr_info.shape[0]
	N_k = sr_cnt / np.sum(1 / np.array(AirTime_k)) / AirTime_k
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
		pi = RSSI_k[k] + PL[sr_index][gw_index]
		if pi <= Ptx_max: # if transmission power does not violate
			sr_info[sr_index, 3] = max(pi, 0.0)
			continue
		# if transmission power violates the upperbound
		# print('tx violate')
		RSSI_max = Ptx_max - PL[sr_index][gw_index]
		try: # try to assign the minimum possible SF
			newRSSI = list(filter(lambda k: k <= RSSI_max, RSSI_k))[0]
			newSFk = RSSI_k.index(newRSSI)
			newpi = RSSI_k[newSFk] + PL[sr_index][gw_index]
			logging.debug("Reassign for tx pow violation: \
				old SF:{} pow:{} new SF:{} pow:{}".format( \
				k, pi, newSFk, newpi))
			sr_info[sr_index, 2] = newSFk
			sr_info[sr_index, 3] = newpi
		except: # if no SF could work, assign the largest SF
			sr_info[sr_index, 2] = len(SF) - 1
			sr_info[sr_index, 3] = Ptx_max
			logging.debug("Reassign failed, RSSI_max: {}".format(RSSI_max))

	return sr_info

def ICIOTAlg(sr_info, G, PL, desired_gw_cnt):
	'''
	Call the ICIOT gateway placement algorithm

	Args:
		sr_info: sensor placement and configuration
		G: gateway placement
		PL: path loss matrix between sensors and potential gateways
		desired_gw_cnt: desired number of gateways to place

	Returns:
		gw_place: a binary vector of gateway placement decision
		sr_info: sensor configuration
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	# Start the greedy gateway placement algorithm
	for rounds in range(desired_gw_cnt):
		obj_old = -np.inf
		next_idx = -1
		next_sr_info = None
		for idx in range(gw_cnt):
			if G[idx, 2]: # a gateway has been placed at this location
				continue
			# try to place gateway at this location
			G[idx, 2] = 1

			# assign powers and SFs
			sr_info = DeviceConfiguration(sr_info, G, PL, Ptx_max)

			# Calculate Cij from each sensor i to gw j
			Cij = GetCij(sr_info, G, PL)
			# print(Cij)

			# Calculate PDR at each sensor i
			PDR = GetPDR(sr_info, G, Cij)
			#print('PDR:', PDR)

			# Calculate energy per packet at each sensor i
			ei = GetEnergyPerPacket(sr_info)
			#print('ei:', ei)

			# Calculate energy efficiency
			EE = np.divide(PDR, ei)
			#print('EE:', EE)

			# Calculate objective value
			gw_place = G[:, 2] # a binary vector indicating gw placement
			obj = np.sum(EE) / sr_cnt - alpha * np.sum(gw_place) / desired_gw_cnt
			#print('obj1:', np.sum(EE) / sr_cnt)
			#print('obj2:', alpha * np.sum(gw_place) / gw_cnt)
			#print('obj:', obj)

			# update objective value if necessary
			if obj > obj_old:
				obj_old = obj
				next_idx = idx
				next_sr_info = sr_info

			# reset
			G[idx, 2] = 0

		# place a gateway at next_idx with the max new objective value
		G[next_idx, 2] = 1
		sr_info = next_sr_info
		print("Placed gateway #{} at grid {} [{},{}]".format( \
			rounds, next_idx, G[next_idx, 0], G[next_idx, 1]))

	gw_place = G[:, 2]
	return gw_place, sr_info

########################################
# Main Process
########################################
def main():
	# Preparation
	# Generate the grid candidate set N and G with their x, y coordinates
	# N for sensor placement and G for gateway placement
	# N = []
	G = []				# [x, y, placed or not]
	#for i in range(N_x):
	#	for j in range(N_y):
	#		new_loc = [i * Unit_sr, j * Unit_sr]
	#		N.append(new_loc)
	for p in range(G_x):
		for q in range(G_y):
			new_loc = [p * Unit_gw, q * Unit_gw, 0]
			G.append(new_loc)
	G = np.array(G)
	gw_cnt = G_x * G_y  # Number of gw potential locations

	# Randomly generate sensor positions
	sr_cnt = 1000 #50000		# Number of sensors
	sr_info = []		# [x, y, SF, Ptx]
	for i in range(sr_cnt):	
		k = random.randint(0, len(SF)-1)
		new_loc = [random.random() * L, random.random() * L, k, Ptx_max]
		sr_info.append(new_loc)
	sr_info = np.array(sr_info)
	# print(sr_info)

	# Generate path loss matrix PL between sensor i and gateway j at (p,q)
	PL = np.zeros((sr_cnt, gw_cnt))
	for i in range(sr_cnt):
		for j in range(gw_cnt):
			loc1 = sr_info[i, :2]
			loc2 = G[j, :2]
			dist = np.sqrt(np.sum((loc1 - loc2)**2))
			PL[i][j] = LogDistancePathLossModel(dist)
	# print(PL)

	desired_gw_cnt = 10 # Desired gateways to place
	gw_place, sr_info = ICIOTAlg(sr_info, G, PL, desired_gw_cnt)
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

	# Visualize the placement and device configuration
	plt.figure()
	colorList = ['b', 'g', 'y', 'm', 'w', 'r', 'k', 'c']
	color = [colorList[int(i)] for i in list(sr_info[:, 2])]
	plt.scatter(sr_info[:, 0], sr_info[:, 1], c=color , s=5)
	color = ['r' for i in range(gw_cnt)]
	plt.scatter(G[:, 0], G[:, 1], s=G[:, 2]*50, c=color, marker='^')
	plt.xlabel('X (m)'); plt.ylabel('Y (m)');
	# plt.legend()
	plt.show()


if __name__ == '__main__':
	main()
