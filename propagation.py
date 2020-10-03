#!/usr/bin/python3
import numpy as np
import math
import matplotlib.pyplot as plt
import logging
import random
logging.basicConfig(level=logging.INFO)

def FreeSpacePathLossModel(**kwargs):
	'''
	Free Space Path Loss Model

	Args:
		d: distance in m
		f: frequency in MHz

	Return:
		path loss in dB at d
	'''
	PL = 20 * math.log10(kwargs['d']) + 20 * math.log10(kwargs['f']) - 27.55
	return PL

def LogDistancePathLossModel(**kwargs):
	'''
	Log Distance Path Loss Model without the uncertainty

	Args:
		d: distance in m

	Return:
		path loss in dB at d
	'''
	if kwargs['ver'] == 'ICIOT':
		# Settings in the ICIOT paper
		LogDistancePathLossModel.d0 = 1000 # Reference distance in m
		LogDistancePathLossModel.PL0 = 130 # Reference path loss in dB
		LogDistancePathLossModel.gamma = 2.1 # Path loss exponent
		LogDistancePathLossModel.sigma = 10.0 # Standard deviation
	elif kwargs['ver'] == 'Bor':
		# Settings in the paper (Bor 2016) for indoor building
		LogDistancePathLossModel.d0 = 40 # Reference distance in m
		LogDistancePathLossModel.PL0 = 127.41 # Reference path loss in dB
		LogDistancePathLossModel.gamma = 2.08 # Path loss exponent
		LogDistancePathLossModel.sigma = 3.57 # Standard deviation
	elif kwargs['ver'] == 'Dongare':
		# Settings in Dongare's thesis based on CMU campus measurements
		LogDistancePathLossModel.d0 = 140 # Reference distance in m
		LogDistancePathLossModel.PL0 = 105.5729 # Reference path loss in dB
		LogDistancePathLossModel.gamma = 2.1495 # Path loss exponent
		LogDistancePathLossModel.sigma = 10.0 # Standard deviation
	else:
		logging.error('Non-supported log propagation model!')
		return
	if kwargs['d'] <= LogDistancePathLossModel.d0:
		return LogDistancePathLossModel.PL0
	PL = LogDistancePathLossModel.PL0 + 10 * LogDistancePathLossModel.gamma * \
		math.log10(kwargs['d'] / LogDistancePathLossModel.d0)
	return PL


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

GetRSSI.Gtx = 0.0 # Transmission antenna gain
GetRSSI.Grx = 0.0 # Reception antenna gain

def main():
	# Test and tune path loss model's parameter
	d = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
	# ICIOT
	loss_log = list(map(lambda di: LogDistancePathLossModel(d=di, ver='ICIOT'), d))
	Prx_log_ICIOT = [GetRSSI(20, lossi) for lossi in loss_log]
	loss_logr = list(map(lambda di: LogDistancePathLossModel(d=di, ver='ICIOT')+random.normalvariate(.0, LogDistancePathLossModel.sigma), d))
	Prx_log_ICIOTr = [GetRSSI(20, lossi) for lossi in loss_logr]
	# Bor
	loss_bor = list(map(lambda di: LogDistancePathLossModel(d=di, ver='Bor'), d))
	Prx_log_Bor = [GetRSSI(20, lossi) for lossi in loss_bor]
	loss_borr = list(map(lambda di: LogDistancePathLossModel(d=di, ver='Bor')+random.normalvariate(.0, LogDistancePathLossModel.sigma), d))
	Prx_log_Borr = [GetRSSI(20, lossi) for lossi in loss_borr]
	# Dongare
	loss_don = list(map(lambda di: LogDistancePathLossModel(d=di, ver='Dongare'), d))
	Prx_log_Dongare = [GetRSSI(20, lossi) for lossi in loss_don]
	loss_donr = list(map(lambda di: LogDistancePathLossModel(d=di, ver='ICIOT')+random.normalvariate(.0, LogDistancePathLossModel.sigma), d))
	Prx_log_Dongarer = [GetRSSI(20, lossi) for lossi in loss_donr]
	# Friis
	loss_friis = list(map(lambda di: FreeSpacePathLossModel(d=di, f=902.3), d))
	Prx_friis= [GetRSSI(20, lossi) for lossi in loss_friis]

	plt.figure()
	plt.plot(d, Prx_log_ICIOTr, label='log ICIOT')
	plt.plot(d, Prx_log_Borr, label='log Bor')
	plt.plot(d, Prx_log_Dongarer, label='log Dongare')
	plt.plot(d, Prx_friis, label='friis')
	plt.legend()
	plt.savefig('loss.png')
	plt.show()

if __name__ == '__main__':
	main()

