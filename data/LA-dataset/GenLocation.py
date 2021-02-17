#!/usr/bin/python3
# The goal of this script is to generate the end devices locations and candidate
# gateway location in the LA dataset
import math
import matplotlib.pyplot as plt
import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

# Important parameters for the LA dataset
class params:
	gw_dist = 6000      # Distance between two gateways in m
	# the given data files and bool variables showing whether or not to use them
	dataFile = dir_path + '/dataLA.csv'			# End device locations
	origin = [33.5466, -118.7025]
	GwAbleFile = dir_path + '/gw_able.npy'	# Whether placing gateway at a location is allowed

def ReadFile(filename, origin):
	'''
	Read the [latitude longitude] from the file and convert them to location in km
	'''
	coor = []
	with open(filename, 'r') as f:
		firstline = f.readline() # ignore first label line
		lines = f.readlines()
		for line in lines:
			data = line.strip().split(',')
			lat = float(data[1])
			lon = float(data[2])
			x = distance(origin, [origin[0], lon])
			y = distance(origin, [lat, origin[1]])
			coor.append([x, y])
	coor = np.array(coor)
	return coor

def distance(origin, destination):
	'''
	Measure the distance in m from two [latitude, longitude] pairs
	Source: https://gist.github.com/rochacbruno/2883505
	'''
	lat1, lon1 = origin
	lat2, lon2 = destination
	radius = 6371 # km

	dlat = math.radians(lat2-lat1)
	dlon = math.radians(lon2-lon1)
	a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
	    * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
	c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
	d = radius * c * 1e3 # Convert from km to m

	return d

def plot(coor):
	# Visualize the end device locations from the data file
	plt.scatter(coor[:, 0], coor[:, 1])
	plt.xlabel('X (km)'); plt.ylabel('Y (km)');
	plt.show()

def main():
	# Read end device locations from dataLA.csv
	coor = ReadFile(params.dataFile, params.origin)
	x_max = np.max(coor[:, 0])
	y_max = np.max(coor[:, 1])

	# Fill in the coordinate date to the initial sensor info
	sr_info = []
	for i in range(coor.shape[0]):
		k = -1 #random.randint(0, len(params.SF)-1) # SFk
		q = -1 # random.randint(0, len(params.CH)-1) # Channel q
		new_loc = [coor[i, 0], coor[i, 1], k, 20, q]
		sr_info.append(new_loc)

	sr_info = np.array(sr_info)
	# print(sr_info)
	sr_cnt = sr_info.shape[0]

	# Generate the grid candidate set N and G with their x, y coordinates
	# N for sensor placement and G for gateway placement
	G = []				# [x, y, placed or not, can place or not]
	x_gw = params.gw_dist / 2
	y_gw = params.gw_dist / 2
	gw_able = np.load(params.GwAbleFile)

	gw_idx = 0
	while x_gw < x_max:
		while y_gw < y_max:
			new_loc = [x_gw, y_gw, 0, gw_able[gw_idx]]
			G.append(new_loc)
			gw_idx += 1
			y_gw += params.gw_dist
		x_gw += params.gw_dist
		y_gw = params.gw_dist / 2

	G = np.array(G)
	# print(G)
	gw_cnt = G.shape[0]

	# Save the end device locations and candidate gateway locations
	np.savetxt('./sr_loc.csv', sr_info[:, :2], delimiter=',')
	np.savetxt('./gw_loc.csv', G[:, :2], delimiter=',')


if __name__ == '__main__':
	main()