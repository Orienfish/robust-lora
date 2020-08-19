#!/usr/bin/python3
import numpy as np
import math
import random
import matplotlib.pyplot as plt


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
gw_cnt = G_x * G_y   # Number of gw potential locations
Unit_gw = math.floor(L / G_x) # Unit length between gw grid points
print(Unit_gw)

T = 1200			# Sampling period in s
Ptx_max = 23		# Maximum allowed transmission power in dBm

# Spreading factors
SF = ['SF7', 'SF8', 'SF9', 'SF10']
# RSSI threshold of SF7, 8, 9, 10 in dBm
# Copied from ICIOT paper and ICDCS paper
RSSI_k = {'SF7': -123, 'SF8': -126, 'SF9': -129, 'SF10': -132}
# Air time of SF7, 8, 9, 10 in s
# Copied from ICIOT paper under 50 Bytes payload
AirTime_k = {'SF7': 0.098, 'SF8': 0.175, 'SF9': 0.329, 'SF10': 0.616}

########################################
# Preparation
########################################
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
		new_loc = [p * Unit_gw, q * Unit_gw, 1]
		G.append(new_loc)

# Randomly generate sensor positions
sr_cnt = 50000		# Number of sensors
sr_info = []		# [x, y, SF, Ptx]
for i in range(sr_cnt):	
	k = random.randint(0, len(SF)-1)
	new_loc = [random.random() * L, random.random() * L, SF[k], Ptx_max]
	sr_info.append(new_loc)
# print(sr_loc)


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
#loss = [LogDistancePathLossModel(di) for di in d]
#Prx = [GetRSSI(20, lossi) for lossi in loss]
#plt.plot(d, Prx)
#plt.show()

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
	sr_cnt = len(sr_info)
	gw_cnt = len(G)
	Cij = np.zeros((sr_cnt, gw_cnt))
	for j in range(gw_cnt):
		if not G[j][2]: # no gateway is placed at j
			continue
		# a gateway is placed at j
		for i in range(sr_cnt):
			Prx = GetRSSI(sr_info[i][3], PL[i][j])
			if Prx > RSSI_k[sr_info[i][2]]: # can be received
				Cij[i][j] = 1

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
	sr_cnt = len(sr_info)
	gw_cnt = len(G)
	SF_cnt = len(SF)
	
	# Calculate traffic load between sensor i and gw j, lambda_ij
	lambda_ij = np.zeros((sr_cnt, gw_cnt))
	N_jk = np.zeros((gw_cnt, SF_cnt)) # Number of nodes connecting to gw j
									  # and using the SFk
	for j in range(gw_cnt):
		# Update number of nodes that might collide
		for i in range(sr_cnt):
			k = SF.index(sr_info[i][2]) # get the SF used by sensor i
			N_jk[j][k] += Cij[i][j]
		# Calculate traffic load 
		for i in range(sr_cnt):
			k = SF.index(sr_info[i][2]) # get the SF used by sensor i
			lambda_ij[i][j] = N_jk[j][k] * AirTime_k[SF[k]] / T
		# print(j, N_jk[j][:])

	# Calculate packet delivery ratio between sensor i and gw j, pi_ij
	pi_ij = np.multiply(Cij, np.exp(-2 * lambda_ij))
	# print(pi_ij.shape)

	# Calculate total PDR at sensor i
	PDR = 1 - np.prod(1 - pi_ij, axis=1)
	
	return PDR


########################################
# Main Process
########################################
# Generate path loss matrix PL between sensor i and gateway j at (p,q)
PL = np.zeros((sr_cnt, gw_cnt))
for i in range(sr_cnt):
	for j in range(gw_cnt):
		loc1 = np.array(sr_info[i][:2])
		loc2 = np.array(G[j][:2])
		dist = np.sqrt(np.sum((loc1 - loc2)**2))
		PL[i][j] = LogDistancePathLossModel(dist)
# print(PL)

# Calculate Cij
Cij = GetCij(sr_info, G, PL)
# print(Cij)

# Calculate PDR
PDR = GetPDR(sr_info, G, Cij)
print(PDR)