#!/usr/bin/python3
import numpy as np
import math
import matplotlib.pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

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
'''
d = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
loss_log = [LogDistancePathLossModel(di) for di in d]
Prx_log = [GetRSSI(20, lossi) for lossi in loss_log]
loss_friis = [FreeSpacePathLossModel(di, 868) for di in d]
Prx_friis = [GetRSSI(20, lossi) for lossi in loss_friis]
plt.figure()
plt.plot(d, Prx_log, d, Prx_friis)
plt.savefig('loss.png')
plt.show()
'''
