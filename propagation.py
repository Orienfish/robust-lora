#!/usr/bin/python3
import numpy as np
import math
import matplotlib.pyplot as plt
import logging
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
		d0 = 1000 # Reference distance in m
		PL0 = 130 # Reference path loss in dB
		delta = 2.1 # Path loss exponent
	elif kwargs['ver'] == 'Bor':
		# Settings in the paper (Bor 2016) for indoor building
		d0 = 40 # Reference distance in m
		PL0 = 127.41 # Reference path loss in dB
		delta = 3.57 # Path loss exponent
	elif kwargs['ver'] == 'Dongare':
		# Settings in Dongare's thesis based on CMU campus measurements
		d0 = 140 # Reference distance in m
		PL0 = 105.5729 # Reference path loss in dB
		delta = 2.1495 # Path loss exponent
	else:
		logging.error('Non-supported log propagation model!')
		return
	if kwargs['d'] <= d0:
		return PL0
	PL = PL0 + 10 * delta * math.log10(kwargs['d'] / d0)
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

GetRSSI.Gtx = 2 # Transmission antenna gain
GetRSSI.Grx = 2 # Reception antenna gain

def main():
	# Test and tune path loss model's parameter
	d = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
	loss_log = [LogDistancePathLossModel(d=di, ver='ICIOT') for di in d]
	Prx_log_ICIOT = [GetRSSI(20, lossi) for lossi in loss_log]
	loss_log = [LogDistancePathLossModel(d=di, ver='Bor') for di in d]
	Prx_log_Bor = [GetRSSI(20, lossi) for lossi in loss_log]
	loss_log = [LogDistancePathLossModel(d=di, ver='Dongare') for di in d]
	Prx_log_Dongare = [GetRSSI(20, lossi) for lossi in loss_log]
	loss_friis = [FreeSpacePathLossModel(d=di, f=868) for di in d]
	Prx_friis = [GetRSSI(20, lossi) for lossi in loss_friis]
	plt.figure()
	plt.plot(d, Prx_log_ICIOT, label='log ICIOT')
	plt.plot(d, Prx_log_Bor, label='log Bor')
	plt.plot(d, Prx_log_Dongare, label='log Dongare')
	plt.plot(d, Prx_friis, label='friis')
	plt.legend()
	plt.savefig('loss.png')
	plt.show()

if __name__ == '__main__':
	main()

