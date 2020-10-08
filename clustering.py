#!/usr/bin/python3
import numpy as np
from sklearn.cluster import DBSCAN
# from sklearn import metrics
import matplotlib.pyplot as plt
import logging

import propagation
import RGreedy
import main

def DBSCANAlg(X, ClusterParams):
	db = DBSCAN(eps=ClusterParams.eps, min_samples=ClusterParams.min_samples).fit(X)
	
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	n_noise_ = list(labels).count(-1)

	print('Estimated number of clusters: %d' % n_clusters_)
	print('Estimated number of noise points: %d' % n_noise_)
	# print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
	# print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
	# print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
	# print("Adjusted Rand Index: %0.3f"
	#       % metrics.adjusted_rand_score(labels_true, labels))
	# print("Adjusted Mutual Information: %0.3f"
	#       % metrics.adjusted_mutual_info_score(labels_true, labels))
	# print("Silhouette Coefficient: %0.3f"
	#       % metrics.silhouette_score(X, labels))

	return db

def plotClusters(X, labels, core_sample_indices, n_clusters):
	core_samples_mask = np.zeros_like(labels, dtype=bool)
	core_samples_mask[core_sample_indices] = True
	
	# Black removed and is used for noise instead.
	unique_labels = set(labels)
	colors = [plt.cm.Spectral(each)
	          for each in np.linspace(0, 1, len(unique_labels))]
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = [0, 0, 0, 1]

	    class_member_mask = (labels == k)

	    xy = X[class_member_mask & core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=14)

	    xy = X[class_member_mask & ~core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)

	plt.title('Estimated number of clusters: %d' % n_clusters)
	plt.savefig('clusters.png')
	plt.show()

def RClusterAlg(sr_info, G, PL, dist, N_kq, params, ClusterParams, GreedyParams):
	'''
	Call the cluster-based robust gateway placement algorithm

	Args:
		sr_info: original read-only sensor placement and configuration
		G: original read-only gateway placement
		PL: path loss matrix between sensors and potential gateways
		dist: distance matrix between sensors and potential gateways
		params: important parameters
		N_kq: a dictionary recording traffic allocation
		ClusterParams: parameters for clustering
		GeneticParams: parameters for this greedy algorithms

	Returns:
		sr_info: sensor configuration
		G: resulted gateway placement
	'''
	gw_cnt = G.shape[0]

	# Clustering
	db = DBSCANAlg(sr_info[:, :2], ClusterParams)
	labels = db.labels_
	n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
	plotClusters(sr_info[:, :2], labels, db.core_sample_indices_, n_clusters)

	for ic in range(n_clusters):
		# Extract the points in this cluster
		sr_info_mask = np.zeros_like(labels, dtype=bool)
		sr_info_mask[labels == ic] = True
		sr_info_blob = np.copy(sr_info[sr_info_mask, :])
		sr_cnt_blob = sr_info_blob.shape[0]

		gw_mask = np.zeros_like(G[:, 2], dtype=bool)
		for j in range(gw_cnt):
			for i in range(sr_cnt_blob):
				if propagation.GetRSSI(params.Ptx_max, PL[sr_info_mask, :][i, j]) >= params.RSSI_k[-1]:
					# If the RSSI under max tx power exceeds the minimum RSSI threshold,
					# we reckon this gateway has the probability of covering end devices 
					# in this cluster
					print(propagation.GetRSSI(params.Ptx_max, PL[i, j]))
					gw_mask[j] = True
					break
		G_blob = np.copy(G[gw_mask, :])
		PL_blob = PL[sr_info_mask, :][:, gw_mask]
		logging.info('Size of this blob: sr: {} gw: {} PL: {}'.format(sr_cnt_blob, \
			G_blob.shape[0], PL_blob.shape))

		# Call greedy algorithm for this cluster
		sr_info_blob, G_blob, m_gateway, N_kq = \
			RGreedy.RGreedyAlg(sr_info_blob, G_blob, PL_blob, dist, N_kq, params, GreedyParams)

		# Fill the cluster result back into the original arrays
		sr_info[sr_info_mask, :] = sr_info_blob
		G[gw_mask, :] = G_blob

		main.plot(sr_info, G, 'cluster_{}'.format(ic))
