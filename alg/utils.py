#!/usr/bin/python3
import numpy as np
import os
import matplotlib.pyplot as plt
import logging

import RGreedy

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

	#print(np.array(PDR_gw))
	#print(PDR)
	#print(lifetime)
	return PDR, PDR_gw, lifetime

def plot(sr_info, G, method):
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
	if not os.path.exists('vis'):
		os.makedirs('vis')

	filename = './vis/vis_{}.png'.format(method)
	plt.savefig(filename)
	# plt.show()

def SaveInfo(sr_info, G, PL, method, params):
	'''
	Save the generated solution to text file

	Args:
		sr_info: generated sensor/end device configuration solution
		G: generated gateway placement solution
		PL: path losses between end devices and selected gateways
		method: a string showing the algorithm, to be added to file name
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	if not os.path.exists('res'):
		os.makedirs('res')

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

	filename = './res/pl_{}.txt'.format(method)
	# if PL file is not provided, we need to extract the ground-truth PL
	# like the process in the initialization function init()
	if params.data and not params.PL:
		PL = np.load(params.PLFile).T
		PL = PL + 10.0

	with open (filename, 'w') as out:
		for i in range(sr_cnt):
			for j in range(gw_cnt):
				if G[j, 2]:
					out.write(str(round(PL[i, j], 6)) + ' ')
			out.write('\n')

def SaveRes(method, sr_cnt, M, gw_cnt, time):
	# Log results
	if not os.path.exists('res'):
		os.makedirs('res')
	with open('res/res.txt', 'a') as out:
		out.write(method + ' ' + str(sr_cnt) + ' ' + str(M) + ' ' + \
			str(gw_cnt) + ' ' + str(time) + '\n')