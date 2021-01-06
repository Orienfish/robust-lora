from utils import *
import numpy as np

# Load segmented map 
mat_file = "./data/entire_seg_map.npy"
entire_seg_map = np.load(mat_file)


# Load gateway and sensor location data 
gw_loc_file = "./data/gw_init.txt"
sr_loc_file = "./data/sr_init.txt"

gw_loc_mat = np.loadtxt(gw_loc_file)
sr_loc_mat = np.loadtxt(sr_loc_file)


# (Lat, Lon) of upper-left and upper-right corners
LU = (34.2331, -118.7025)
RU = (34.2331, -117.4719)


# Calculate map resolution
R = lat_lon_to_distance(origin=LU, destination=RU)/entire_seg_map.shape[1]
print("Map resolution:{}".format(R))

path_loss_mat = np.empty([gw_loc_mat.shape[0], sr_loc_mat.shape[0]])


# Calculate path loss from gateway i to sensor j 
gw_num = gw_loc_mat.shape[0]
sr_num = sr_loc_mat.shape[0]

for i in range(gw_num):
    gw_x = int(gw_loc_mat[i,0]/R+0.5)
    gw_y = entire_seg_map.shape[0]-int(gw_loc_mat[i,1]/R+0.5)
    gateway_loc = min(max(gw_x, 1), entire_seg_map.shape[1]-1), min(max(gw_y, 1), entire_seg_map.shape[0]-1)

    print("Estimating PL of {}/{} gateway".format(i+1, gw_num))

    for j in range(sr_num):
        sr_x = int(sr_loc_mat[j,0]/R+0.5)
        sr_y = entire_seg_map.shape[0]-int(sr_loc_mat[j,1]/R+0.5)
        sensor_loc = min(max(sr_x, 1), entire_seg_map.shape[1]-1), min(max(sr_y, 1), entire_seg_map.shape[0]-1)

        # print("i={}, j={},\t gateway_loc={}, sensor_loc={}".format(i, j, gateway_loc, sensor_loc))
        # print("original:gw={}, {}\tsr={}, {}".format(gw_loc_mat[i,0], gw_loc_mat[i,1], sr_loc_mat[j,0], sr_loc_mat[j,1]))

        line_distance, path_loss = path_loss_estimation(segmented_land_map=entire_seg_map, map_LU=LU, map_RU=RU, gateway_loc=gateway_loc, sensor_loc=sensor_loc)
        path_loss_mat[i, j] = path_loss
        
        # print("line_distance={}m, path_loss={}\n".format(line_distance, path_loss))


mat_file = "./data/path_loss_mat.npy"
np.save(mat_file, path_loss_mat)
print("Saving results to {}\n".format(mat_file))