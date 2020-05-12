'''

  Software for the detection of eddies in model output following Chelton et al., Progress in Oceanography, 2011.

  Adapted from and using baseline code from Eric Oliver eddyTracking.

  Including additional detection designed by Jamie Atkins involving additional metrics for subsequent modified tracking methods to include a new type of eddy (i.e. anticyclonic/cyclonic/ACME).

'''

#****************************************************************
# Load required modules

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import eddy_functionsADD as eddy
import anomFunctions
from paramsADD import *
from itertools import repeat
#****************************************************************
# Load latitude and longitude vectors
fileobj = Dataset(pathroot + filename, mode = 'r')
lon = fileobj.variables['nav_lon'][0,:] # isolate only the lon values
lat = fileobj.variables['nav_lat'][:,0] # isolate only the lat values

# create empty lists
lon_eddies_a = []
lat_eddies_a = []
damp_eddies_a = []
area_eddies_a = []
scale_eddies_a = []
time_eddies_a = []
lon_eddies_c = []
lat_eddies_c = []
amp_eddies_c = []
area_eddies_c = []
scale_eddies_c = []
time_eddies_c = []

# loop over each day
print('eddy detection started')
print("number of time steps to loop over: ", T)
for tt in range(5):
    print("timestep: ",tt+1,". out of: ", T)
    # load SSH data
    eta, eta_miss = eddy.load_eta(tt)
    eta = eddy.remove_missing2(eta)
    ## Spatially filter SSH field
    eta_filt = eddy.spatial_filter(eta, lon, lat, res, cut_lon, cut_lat)
    # Detect lon and lat coordinates of ANTICYCLONIC eddies
    lon_eddies, lat_eddies, amp, area, scale = eddy.detect_eddies(eta_filt, lon, lat, ssh_crits, res, Npix_min, Npix_max, amp_thresh, d_thresh, cyc='anticyclonic')
    lon_eddies_a.append(lon_eddies)
    lat_eddies_a.append(lat_eddies)
    amp_eddies_a.append(amp)
    area_eddies_a.append(area)
    scale_eddies_a.append(scale)
    time_eddies_a.append(tt)
    # Detect lon and lat coordinates of CYCLONIC eddies
    lon_eddies, lat_eddies, amp, area, scale = eddy.detect_eddies(eta_filt, lon, lat, ssh_crits, res, Npix_min, Npix_max, amp_thresh, d_thresh, cyc='cyclonic')
    lon_eddies_c.append(lon_eddies)
    lat_eddies_c.append(lat_eddies)
    amp_eddies_c.append(amp)
    area_eddies_c.append(area)
    scale_eddies_c.append(scale)
    time_eddies_c.append(tt)

# Combine eddy information from all days into a list   
eddies = eddy.eddies_list(lon_eddies_a, lat_eddies_a, amp_eddies_a, area_eddies_a, scale_eddies_a, time_eddies_a, lon_eddies_c, lat_eddies_c, amp_eddies_c, area_eddies_c, scale_eddies_c, time_eddies_c)

#****************************************************************
# SAVE DATA [ intermediate step ] [ N.B. save here if additional metrics are not necessary ] 
## this will save a copy of the tracked eddies without additional metrics, i.e. the same output from original Eric Oliver eddyTracking

np.savez(data_dir+'eddy_detINT_'+ run, eddies=eddies)

# load data in at this point if necessary
#data = np.load('/home/bridge/ja16048/mscr/subsurfacedata_tracking/nep/detection_output/eddy_detINT_freeglorys.npz')
#eddies = data['eddies']

#****************************************************************

# PART 2) ADDITIONAL METRICS 

# Load additional metrics [ SSS, SST, MLD etc. ] 
sss, sst, mld, o2, depth, t, t_len = eddy.load_additional(pathroot, filename, filenameSAL, pathroot_o2, filename_o2)

# SSS and SST baseline calculations  [ to calculate anomalies later ]
sss_int = anomFunctions.clim_calcs(sss, num_years)
sss_baseline = anomFunctions.clim_repeat(sss_int, t_len)
sst_int =  anomFunctions.clim_calcs(sst, num_years)
sst_baseline = anomFunctions.clim_repeat(sst_int, t_len)

o2_int = anomFunctions.clim_calcs_o2(o2, num_years) # smae but for oxygen which involves a depth component, unlike sss and sst
o2_baseline = anomFunctions.clim_repeat_o2(o2_int, t_len) # smae but for oxygen which involves a depth component, unlike sss and sst

# additional metrics (SSS anomaly, SST anomaly, and MLD) for detected eddies
eddies_tempstore = eddy.add_metrics(eddies, o2, depth, sss_baseline, sss, sst_baseline, sst, mld, res, t, lat, lon)

# re-list the eddy information from all days into a list [ an updated version of 'eddy_list' function ] 
eddies = eddy.eddies_listVER2(lon_eddies_a, lat_eddies_a, amp_eddies_a, area_eddies_a, scale_eddies_a, time_eddies_a, lon_eddies_c, lat_eddies_c, amp_eddies_c, area_eddies_c, scale_eddies_c, time_eddies_c, eddies_tempstore)

# DISTINGUISH BETWEEN ACMEs AND "NORMAL" ANTICYCLONIC EDDIES
eddies = eddy.acme_distinct(eddies)

#****************************************************************

# SAVE DATA

np.savez(data_dir+'eddy_det_'+ run, eddies=eddies)
