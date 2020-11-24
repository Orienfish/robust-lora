#!/usr/bin/python3
import numpy as np
import propagation
import RGreedy

def TestLifetime(params):
	'''
	Test the lifetime under unconfirmed traffic with all SFs and Ptxs
	to trasform to the loosen linear lifetime constraint.
	'''
	SF_cnt = len(params.SF)
	for sf in range(SF_cnt):
		for ptx in params.Ptx:
			print(sf, ptx, RGreedy.GetLifetime(sf, ptx, 1.0, params))

	return

def GenerateCijks(sr_info, G, PL, params):
	'''
	Generate the c_ijkq matrix denoting the reachability from devices to
	gateways using SF k and transmission power s

	Args:
		sr_loc: sensor information
		G: gateway information
		PL: path loss matrix
		params: important parameters

	Return:
		c_ijk: a binary sr_cnt * gw_cnt * sf_cnt matrix
	'''
	sr_cnt = sr_info.shape[0]
	gw_cnt = G.shape[0]
	SF_cnt = len(params.SF)
	TP_cnt = len(params.Ptx)
	c_ijks = np.zeros((sr_cnt, gw_cnt, SF_cnt, TP_cnt))

	for i in range(sr_cnt):
		for j in range(gw_cnt):
			for k in range(SF_cnt):
				for s in range(TP_cnt):
					# Note that Ptx array in params is in reverse order
					if propagation.GetRSSI(params.Ptx[TP_cnt-s-1], PL[i, j]) > \
						params.RSSI_k[k]:
						c_ijks[i, j, k, s] = 1

	return c_ijks
