import random 
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

# Predefined parameters for various land types
# 0:unknown, 1:urban, 2:agriculture, 3:rangeland, 4:forest, 5:water, 6:barren
PL_d0_table = {0: 0, 1: 23.091, 2:20.478, 3:20.392, 4:21.254, 5:20.709, 6:21.254}
n_table   = {0: 0, 1: 4.499, 2:2.056, 3:2.274, 4:3.616, 5:2.158, 6:3.616}
sigma_table = {0: 0, 1: 1.887, 2:0.569, 3:0.622, 4:0.854, 5:0.599, 6:0.854}

#############################################################################
def path_loss_estimation(segmented_land_map, map_LU, map_RU, gateway_loc, sensor_loc, verbose=False):
    """
    Path loss estimation based on the segmented map.

    Args:
        segmented_land_map: Segmented land map with Cartesian pixel coordinate system with (0,0) in the upper-left corner
        map_LU: The (lat, lon) of upper-left corner for the target region
        map_RU: The (lat, lon) of upper-right corner for the target region
        gateway_loc: Gateway location (xG, yG) in meters with (0,0) in the lower-left corner
        sensor_loc: Gateway location (x, y) in meters with (0,0) in the lower-left corner
        verbose:

    Returns:
        L: Link distance from gateway to sensor node
        PL: Estimated path loss value
    """

    seg_map = segmented_land_map

    xG, yG = gateway_loc
    x, y = sensor_loc

    assert ((0< x < seg_map.shape[1]) and (0< y < seg_map.shape[0])), "Sensor coordinate out of boundary!"
    assert ((0< xG < seg_map.shape[1]) and (0< yG < seg_map.shape[0])), "Gateway coordinate out of boundary!"

    ################ Get traversed pixels from gateway to sensor ################
    if(xG==x):
        line_x_int = x*np.ones(abs(y-yG)+1).astype(int)
        if(yG>y):
            line_y_int = np.arange(yG, y-1, -1).astype(int)
        else:
            line_y_int = np.arange(yG, y+1).astype(int)

        line_y_float = line_y_int.astype(float)
    else:
        link_eq_slope = (yG-y)/(xG-x)
        link_eq_delta = yG-link_eq_slope*xG

        if(x>xG):
            line_x_int = np.arange(xG, x+1).astype(int)
        else:
            line_x_int = np.arange(xG, x-1, -1).astype(int)

        line_y_float = link_eq_slope*line_x_int+link_eq_delta
        line_y_int = np.floor(line_y_float+0.5).astype(int)

    line_land_rgb = seg_map[line_y_int, line_x_int] # stores rgb color for land over the line
    line_loc_float = np.stack([line_x_int, line_y_float], axis=1) # coordinator of line
    line_loc_int = np.stack([line_x_int, line_y_int], axis=1)

    ################ Convert RGB land type to predefined parameters ################
    line_land_type = []
    for i in range(0, len(line_land_rgb)):
        if(np.array_equal(line_land_rgb[i], [0, 255, 255])): #if land type is urban
            line_land_type.append(1)
        elif(np.array_equal(line_land_rgb[i], [255, 255, 0])): #if land type is argiculuture
            line_land_type.append(2)
        elif(np.array_equal(line_land_rgb[i], [255, 0, 255])): #if land type is rangeland
            line_land_type.append(3)
        elif(np.array_equal(line_land_rgb[i], [0, 255, 0])): #if land type is forrest
            line_land_type.append(4)
        elif(np.array_equal(line_land_rgb[i], [0, 0, 255])): #if land type is water
            line_land_type.append(5)
        elif(np.array_equal(line_land_rgb[i], [255, 255, 255])): #if land type is barren
            line_land_type.append(6)
        else: #land type is unknown
            line_land_type.append(0)
        
    ################ Calculate path loss ################
    map_resolution = lat_lon_to_distance(origin=map_LU, destination=map_RU)/seg_map.shape[1]

    if verbose:
        print("Map resolution={}m/pixel".format(map_resolution))

    #### Calculate the initial loss ####
    L = calc_distance(line_loc_float[0], line_loc_float[1], map_resolution)
    PL = calc_init_path_loss(dest_land_type=line_land_type[1], d=L)
    xy_start = line_loc_float[1]

    #### Iterative updating ####
    for i in range(1, line_loc_float.shape[0]-2):
        if(line_land_type[i]==line_land_type[i+1]):
            continue
        else:
            xy_end = line_loc_float[i]
            
            d = calc_distance(xy_start, xy_end, map_resolution)
            L, PL = update_path_loss(PL_prev=PL, dest_land_type=line_land_type[i], L_prev=L, d=d)
            
            xy_start = line_loc_float[i]

            if verbose:
                print("{}-step,\tLand-type={},\tL={},\tPL={}".format(i, line_land_type[i], L, PL))

    #### Calculate the last segment loss ####
    i_end = line_loc_float.shape[0]-1
    xy_end = line_loc_float[i_end]

    d = calc_distance(xy_start, xy_end, map_resolution)
    L, PL = update_path_loss(PL_prev=PL, dest_land_type=line_land_type[i_end], L_prev=L, d=d)
            
    if verbose:
        print("{}-step,\tLand-type={},\tL={},\tPL={}".format(i_end, line_land_type[i_end], L, PL))

    return L, PL


def lat_lon_to_distance(origin, destination):
    """"
    Convert lat and lon coordinate pair into distance(km)

    Args:
        origin: The (lat, lon) of origin point
        destination: The (lat, lon) of destination point

    Return:
        d: The distance from origin to destination in meters
    """
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371 #km

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c *1000

    return d


def calc_distance(ij_start, ij_end, R):
    """"
    Calculate distance from start to end point

    Args:
        ij_start: The coordinate of origin point with (0,0) in the upper-left corner
        ij_end: The coordinate of destination point with (0,0) in the upper-left corner
        R: Map resolution (m/pixel)

    Return:
        The distance from origin to destination in meters
    """
    return R*math.sqrt((ij_start[0]-ij_end[0])**2 + (ij_start[1]-ij_end[1])**2)

def calc_init_path_loss(dest_land_type, d, PL_d0_table=PL_d0_table, n_table=n_table, d0=1.0):
    """"
    Calculate the initial path loss

    Args:
        dest_land_type: The land type of end point
        d: The distance from start point to end point
        PL_d0_table: The dict stores PL_d0 for all land types
        n_table: The dict stores n values for all land types
        d0: Reference distance

    Return:
        The initial path loss value
    """
    assert (0<=dest_land_type<=6), "Invalid land type!"
    return PL_d0_table[dest_land_type]+10*n_table[dest_land_type]*math.log10(d/d0) #+random.gauss(0, sigma_table[land_type])
    
def update_path_loss(PL_prev, dest_land_type, L_prev, d, PL_d0_table=PL_d0_table, n_table=n_table):
    """"
    Calculate the initial path loss

    Args:
        PL_prev: The previous PL 
        dest_land_type: The land type of end point
        d: The distance from start point to end point
        PL_d0_table: The dict stores PL_d0 for all land types
        n_table: The dict stores n values for all land types
        d0: Reference distance

    Return:
        L_new: The updated distance value
        PL_new: The updated PL value
    """
    assert (0<=dest_land_type<=6), "Invalid land type!"

    L_new = L_prev + d
    # print("\nPL_prev:{}, dest_land_type:{}, n_table:{}, d:{}".format(PL_prev, dest_land_type, n_table[dest_land_type], d))
    PL_new = PL_prev + 10*n_table[dest_land_type]*math.log10(L_new/L_prev)

    return L_new, PL_new
