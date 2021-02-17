#!/usr/bin/python3
import math
import matplotlib.pyplot as plt
import numpy as np

# Important parameters or the LA dataset
class params:
	gw_dist = 6000      # Distance between two gateways in m

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
	coor = ReadFile('./data/dataLA.csv', [33.6711, -118.5911])
	print('Number of end devices: {}'.format(coor.shape[0]))
	plot(coor)

if __name__ == '__main__':
	main()