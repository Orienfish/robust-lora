from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import math
import random

################## USER ENTERED DATA #################
img_name = './data/HPWREN/HPWREN.jpg'



#in degrees
lat_top = 33.825275
lat_bottom = 32.45297
log_left = -118.27046
log_right = -115.55555

lat_sen = [33.16, 33.33, 32.95, 32.73, 32.70, 33.19, 32.89, 32.84, 33.01,
            32.60, 33.35, 33.32, 33.40, 33.43, 33.50, 33.70, 33.13,
            33.71, 33.61, 33.38, 33.46, 33.52, 33.61, 33.52, 32.88, 33.27]
log_sen = [-116.81, -116.92, -116.61, -116.58, -116.76, -116.76, -116.42, -116.42, -116.97,
            -116.84, -116.98, -116.68, -117.19, -117.59, -117.60, -116.94, -116.61,
            -117.53, -117.81, -116.62, -117.17, -116.43, -117.55, -117.48, -117.24, -116.64]
loc_name = ["bbm", "boucher", "cuyamaca", "lospinos", "lyons", "mesa", "mtlag", "mtlagobs", "mtwoodson",
            "otay", "pala", "lacruz", "red", "sanclem", "sanjuan", "sanmiguel", "ysabel",
            "santiago","signal", "sky", "SMER", "toro", "upperbell", "uppertalega", "UCSD", "warner"]
            
#create gw_able file
gw_able = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
            
#np.save("./data/HPWREN/gw_able.npy", gw_able)
#print(np.load("./data/HPWREN/gw_able.npy"))
print(np.load("./data/HPWREN/entire_seg_map.npy"))
            
#in meters
mapwidth = 180263.28
mapheight = 135562.06
gw_dist = 7250
sen_rad = 14000

sensor_num = 50

#buffer for range to generate potential sensors in degrees (optimizes sensor placement)
buf = 0.5
'''
######################################################
#open image
im = plt.imread(img_name)
sat_img = Image.open(img_name)
width, height = sat_img.size
ig,ax = plt.subplots(1)
ax.set_aspect('equal')

#show image
implot = plt.imshow(im)

#conversion factors
mtp_conv = width/mapwidth
lat_pix = (lat_top - lat_bottom)/height
log_pix = (log_left - log_right)/width
pix_lat = height/ (lat_top - lat_bottom)
pix_log = width/(log_left - log_right)

################# GATEWAY PLACEMENT ##################

#calculate number of gateways
width_gw_num = int(mapwidth/gw_dist) +1
height_gw_num = int(mapheight/gw_dist) +1

#calculate space between gateways
gw_pix = gw_dist * mtp_conv

latlog_gw = []
for i in range(width_gw_num):
	for j in range(height_gw_num):
		#places gateways 
		xcoord_gw = (i*gw_pix)
		ycoord_gw = (j*gw_pix)
		plt.scatter(xcoord_gw, ycoord_gw, c = 'm', s = 1)
		latlog_gw.insert(i+j, [lat_top-(ycoord_gw * lat_pix), log_left-(xcoord_gw * log_pix)])
nplatlog_gw = np.asarray(latlog_gw)
np.savetxt("./data/HPWREN/HPWREN_gw_loc.csv", nplatlog_gw, delimiter=",", header = "gateway locations")

################## SENSOR PLACEMENT ##################

#calculate sensor range
sen_range = sen_rad * mtp_conv

latlog_sen = []
for j in range(len(lat_sen)):
    print("Placing sensors around sensor #" + str(j+1))
    y_sen = (lat_top - lat_sen[j]) * pix_lat
    x_sen = (log_left - log_sen[j]) * pix_log
    i = 0
    while i < sensor_num:
        #generate random sensor location
        logcoord = random.uniform(log_sen[j] - buf, log_sen[j] + buf)
        latcoord = random.uniform(lat_sen[j] - buf, lat_sen[j] + buf)
        y_pix = (lat_top - latcoord) * pix_lat
        x_pix = (log_left - logcoord) * pix_log
        if y_pix > 0:
            mag = math.sqrt(pow((x_sen - x_pix),2) + pow((y_sen - y_pix),2))
            #if sensor is in range, place sensor
            if mag < sen_range:
                plt.scatter(x_pix, y_pix, s=1, c='y', alpha=0.85)
                latlog_sen.append([i, latcoord, logcoord])
                i = i + 1
    plt.scatter(x_sen, y_sen, s=15, c='g', alpha=0.85)
    
nplatlog_sen = np.asarray(latlog_sen)
np.savetxt("./data/HPWREN/dataHPWREN.csv", nplatlog_sen, delimiter=",", header = "sensor")

plt.show()
'''
