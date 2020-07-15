'''

    A set of functions to accompany the eddy tracking software.

    Predominantly written by Eric Oliver from his eddy_functions.py script under the terms of the GNU General Public License, but also including some additions by Jamie Atkins.


'''

#****************************************************************
# load required modules

import numpy as np
import scipy as sp
import numpy.linalg as linalg
import scipy.signal as signal
import scipy.ndimage as ndimage
import scipy.interpolate as interpolate
import glob
import matplotlib
# Turn the followin on if you are running on storm sometimes - Forces matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from netCDF4 import Dataset
from itertools import repeat
import paramsADD
from paramsADD import *
import re
import math
#****************************************************************
# FUNCTIONS

def find_nearest(array, value):
    idx=(np.abs(array-value)).argmin()
    return array[idx], idx

def nanmean(array, axis=None):
    return np.mean(np.ma.masked_array(array, np.isnan(array)), axis)

def restrict_lonlat(lon, lat, lon1, lon2, lat1, lat2):
    '''
    Restricts latitude and longitude vectors given
    input limits.
    '''

    tmp, i1 = find_nearest(lon, lon1)
    tmp, i2 = find_nearest(lon, lon2)
    tmp, j1 = find_nearest(lat, lat1)
    tmp, j2 = find_nearest(lat, lat2)

    lon = lon[i1:i2+1]
    lat = lat[j1:j2+1]

    return lon, lat, i1, i2, j1, j2


def load_eta(tt):
    
    '''
    Loads sea surface height field 
    '''
    
    pathroot = './'

    z = pathroot + "day" + str(tt) + ".nc"
    fileobj = Dataset(z)

    eta = fileobj.variables['sossheig'][:][0,:,:]
    eta_miss = fileobj.variables['sossheig']._FillValue
    fileobj.close()

    return eta, eta_miss

def load_additional(pathroot, filename, filenameSAL, pathroot_o2, filename_o2):
    '''
    Loads additional metrics [ SSS, SST, MLD ] for Part 2) of eddy_detection.py
    '''
    fileobj = Dataset(pathroot + filename, mode = 'r')
    fileobjSAL = Dataset(pathroot + filenameSAL, mode = 'r')
    fileobj_o2 = Dataset(pathroot_o2 + filename_o2, mode = 'r')
    sss = fileobjSAL.variables['vosaline'][:][:,0,:,:] # isolating the surface depth level
    sst = fileobj.variables['sosstobs'][:][:,:,:]
    mld = fileobj.variables['somxl010'][:][:,:,:]
    o2 = fileobj_o2.variables['o2'][:][:,:,:,:]
    depth = fileobj_o2.variables['deptht'][:]
    t = np.arange(0,len(sss)) # range of timesteps in analysis period
    t_len = len(t) # number of timesteps in analysis period
    
    return sss, sst, mld, o2, depth, t, t_len

def remove_missing2(field):
    '''
    New version of remove 'missing'; replace masked [ land ] values with np.nan's
    '''
    field = field.filled(np.nan)

    return field


def interp_nans(data, indices):
    '''
    Linearly interpolates over missing values (np.nan's) in data
    Data is defined at locations in vector indices.
    '''

    not_nan = np.logical_not(np.isnan(data))

    return np.interp(indices, indices[not_nan], data[not_nan])


def match_missing(data1, data2):
    '''
    Make all locations that are missing in data2 also missing in data1
    Missing values are assumed to be np.nan.
    '''

    data1[np.isnan(data2)] = np.nan
    return data1


def spatial_filter(field, lon, lat, res, cut_lon, cut_lat):
    '''
    Performs a spatial filter, removing all features with
    wavelenth scales larger than cut_lon in longitude and
    cut_lat in latitude from field (defined in grid given
    by lon and lat).  Field has spatial resolution of res
    and land identified by np.nan's
    '''
    
    field_filt = np.zeros(field.shape)
    
    # see Chelton et al, Prog. Ocean., 2011 for explanation of factor of 1/5
    sig_lon = (cut_lon/5.) / res
    sig_lat = (cut_lat/5.) / res
    
    land = np.isnan(field)
    field[land] = nanmean(field) # replaces all land pixels with average SSH value of the domain over entire analysis period
    
    field_filt = field - ndimage.gaussian_filter(field, [sig_lat, sig_lon])
    
    field_filt[land] = np.nan
    
    return field_filt


def distance_matrix(lons,lats):
    '''Calculates the distances (in km) between any two cities based on the formulas
    c = sin(lati1)*sin(lati2)+cos(longi1-longi2)*cos(lati1)*cos(lati2)
    d = EARTH_RADIUS*Arccos(c)
    where EARTH_RADIUS is in km and the angles are in radians.
    Source: http://mathforum.org/library/drmath/view/54680.html
    This function returns the matrix.'''

    EARTH_RADIUS = 6378.1
    X = len(lons)
    Y = len(lats)
    assert X == Y, 'lons and lats must have same number of elements'

    d = np.zeros((X,X))

    #Populate the matrix.
    for i2 in range(len(lons)):
        lati2 = lats[i2]
        loni2 = lons[i2]
        c = np.sin(np.radians(lats)) * np.sin(np.radians(lati2)) + \
            np.cos(np.radians(lons-loni2)) * \
            np.cos(np.radians(lats)) * np.cos(np.radians(lati2))
        d[c<1,i2] = EARTH_RADIUS * np.arccos(c[c<1])

    return d


def detect_eddies(field, lon, lat, ssh_crits, res, Npix_min, Npix_max, amp_thresh, d_thresh, cyc='anticyclonic'):
    '''
    Detect eddies present in field which satisfy the criteria
    outlined in Chelton et al., Prog. ocean., 2011, App. B.2.

    Field is a 2D array specified on grid defined by lat and lon.

    ssh_crits is an array of ssh levels over which to perform
    eddy detection loop

    res is resolutin in degrees of field

    Npix_min, Npix_max, amp_thresh, d_thresh specify the constants
    used by the eddy detection algorithm (see Chelton paper for
    more details)

    cyc = 'cyclonic' or 'acme' or 'anticyclonic' [default] specifies type of
    eddies to be detected

    Function outputs lon, lat coordinates of detected eddies
    '''

    len_deg_lat = 111.325 # length of 1 degree of latitude [km]

    llon, llat = np.meshgrid(lon, lat)

    lon_eddies = np.array([])
    lat_eddies = np.array([])
    amp_eddies = np.array([])
    area_eddies = np.array([])
    scale_eddies = np.array([])

    # ssh_crits increasing for 'cyclonic', decreasing for 'anticyclonic'
    ssh_crits.sort()
    if cyc == 'anticyclonic':
        ssh_crits = np.flipud(ssh_crits)

    # loop over ssh_crits and remove interior pixels of detected eddies from subsequent loop steps
    for ssh_crit in ssh_crits:
 
    # 1. Find all regions with eta greater (less than) than ssh_crit for anticyclonic (cyclonic) eddies (Chelton et al. 2011, App. B.2, criterion 1)
        if cyc == 'anticyclonic':
            regions, nregions = ndimage.label( (field>ssh_crit).astype(int) )
        elif cyc == 'cyclonic':
            regions, nregions = ndimage.label( (field<ssh_crit).astype(int) )

        for iregion in range(nregions):
 
    # 2. Calculate number of pixels comprising detected region, reject if not within [Npix_min, Npix_max]
            region = (regions==iregion+1).astype(int)
            region_Npix = region.sum()
            eddy_area_within_limits = (region_Npix < Npix_max) * (region_Npix > Npix_min)
 
    # 3. Detect presence of local maximum (minimum) for anticylonic (cyclonic) eddies, reject if non-existent
            interior = ndimage.binary_erosion(region)
            exterior = region.astype(bool) ^ interior
            if interior.sum() == 0:
                continue
            if cyc == 'anticyclonic':
                has_internal_ext = field[interior].max() > field[exterior].max()
            elif cyc == 'cyclonic':
                has_internal_ext = field[interior].min() < field[exterior].min()
 
    # 4. Find amplitude of region, reject if < amp_thresh
            if cyc == 'anticyclonic':
                amp = field[interior].max() - field[exterior].mean()
            elif cyc == 'cyclonic':
                amp = field[exterior].mean() - field[interior].min()
            is_tall_eddy = amp >= amp_thresh
 
    # 5. Find maximum linear dimension of region, reject if < d_thresh
            if np.logical_not( eddy_area_within_limits * has_internal_ext * is_tall_eddy):
                continue
 
            lon_ext = llon[exterior]
            lat_ext = llat[exterior]
            d = distance_matrix(lon_ext, lat_ext)
            
            # latitudinal dependence adaptive maximum linear dimension threshold
            lat_bin = np.abs(np.mean(llat[interior]))
            if lat_bin >= 25:
                is_small_eddy = d.max() < d_thresh
            elif lat_bin < 25 and lat_bin >= 20:
                is_small_eddy = d.max() < d_thresh_20
            elif lat_bin < 20 and lat_bin >= 15:
                is_small_eddy = d.max() < d_thresh_15
            elif lat_bin < 15 and lat_bin >= 10:
                is_small_eddy = d.max() < d_thresh_10
            elif lat_bin < 10 and lat_bin >= 5:
                is_small_eddy = d.max() < d_thresh_5
            elif lat_bin < 5:
                is_small_eddy = d.max() < d_thresh_0

    # Detected eddies:
            if eddy_area_within_limits * has_internal_ext * is_tall_eddy * is_small_eddy:
                # find centre of mass of eddy
                eddy_object_with_mass = field * region
                eddy_object_with_mass[np.isnan(eddy_object_with_mass)] = 0
                j_cen, i_cen = ndimage.center_of_mass(eddy_object_with_mass)
                lon_cen = np.interp(i_cen, range(0,len(lon)), lon)
                lat_cen = np.interp(j_cen, range(0,len(lat)), lat)
                lon_eddies = np.append(lon_eddies, lon_cen)
                lat_eddies = np.append(lat_eddies, lat_cen)
                # assign (and calculate) amplitude, area, and scale of eddies
                amp_eddies = np.append(amp_eddies, amp)
                area = region_Npix * res**2 * len_deg_lat * len_deg_lon(lat_cen) # [km**2]
                area_eddies = np.append(area_eddies, area)
                scale = np.sqrt(area / np.pi) # [km]
                scale_eddies = np.append(scale_eddies, scale)
                # remove its interior pixels from further eddy detection
                eddy_mask = np.ones(field.shape)
                eddy_mask[interior.astype(int)==1] = np.nan
                field = field * eddy_mask

    return lon_eddies, lat_eddies, amp_eddies, area_eddies, scale_eddies


def eddies_list(lon_eddies_a, lat_eddies_a, amp_eddies_a, area_eddies_a, scale_eddies_a, time_eddies_a, lon_eddies_c, lat_eddies_c, amp_eddies_c, area_eddies_c, scale_eddies_c, time_eddies_c):
    
    '''
    Creates list detected eddies
    '''

    eddies = []

    for ed in range(len(lon_eddies_c)):
        eddy_tmp = {}
        eddy_tmp['lon'] = np.append(lon_eddies_a[ed], lon_eddies_c[ed])
        eddy_tmp['lat'] = np.append(lat_eddies_a[ed], lat_eddies_c[ed])
        eddy_tmp['amp'] = np.append(amp_eddies_a[ed], amp_eddies_c[ed])
        eddy_tmp['area'] = np.append(area_eddies_a[ed], area_eddies_c[ed])
        eddy_tmp['scale'] = np.append(scale_eddies_a[ed], scale_eddies_c[ed])
        eddy_tmp['type'] = list(repeat('anticyclonic',len(lon_eddies_a[ed]))) + list(repeat('cyclonic',len(lon_eddies_c[ed])))
        eddy_tmp['N'] = len(eddy_tmp['lon'])
        eddy_tmp['time'] = np.append(time_eddies_a[ed], time_eddies_c[ed])
        eddies.append(eddy_tmp)

    return eddies


def eddies_init(det_eddies):
    '''
    Initializes list of eddies. The ith element of output is
    a dictionary of the ith eddy containing information about
    position and size as a function of time, as well as type.
    '''
    
    eddies = []
    
    for ed in range(det_eddies[0]['N']):
        eddy_tmp = {}
        eddy_tmp['lon'] = np.array([det_eddies[0]['lon'][ed]])
        eddy_tmp['lat'] = np.array([det_eddies[0]['lat'][ed]])
        eddy_tmp['amp'] = np.array([det_eddies[0]['amp'][ed]])
        eddy_tmp['area'] = np.array([det_eddies[0]['area'][ed]])
        eddy_tmp['scale'] = np.array([det_eddies[0]['scale'][ed]])
        eddy_tmp['type'] = det_eddies[0]['type'][ed]
        eddy_tmp['time'] = np.array([1])  # temporary distinction between time and timestep from detection
        eddy_tmp['timestep'] = det_eddies[0]['time'][1]
        eddy_tmp['sss_anom'] = np.array([det_eddies[0]['sss_anom'][ed]])
        eddy_tmp['sst_anom'] = np.array([det_eddies[0]['sst_anom'][ed]])
        eddy_tmp['mld'] = np.array([det_eddies[0]['mld'][ed]])
        eddy_tmp['exist_at_start'] = True
        eddy_tmp['terminated'] = False
        eddies.append(eddy_tmp)
        
    return eddies


def load_rossrad():
    '''
    Load first baroclinic wave speed [m/s] and Rossby radius
    of deformation [km] data from rossrad.dat (Chelton et al., 1998)
    
    Also calculated is the first baroclinic Rossby wave speed [m/s]
    according to the formula:  cR = -beta rossby_rad**2
    '''
    
    #data = np.loadtxt('data/rossrad.dat')
    
    #cb
    data = np.loadtxt('./rossrad.dat')
    
    rossrad = {}
    rossrad['lat'] = data[:,0]
    rossrad['lon'] = data[:,1]
    rossrad['c1'] = data[:,2] # m/s
    rossrad['rossby_rad'] = data[:,3] # km
    
    R = 6371.e3 # Radius of Earth [m]
    Sigma = 2 * np.pi / (24*60*60) # Rotation frequency of Earth [rad/s]
    beta = (2*Sigma/R) * np.cos(rossrad['lat']*np.pi/180) # 1 / m s
    rossrad['cR'] = -beta * (1e3*rossrad['rossby_rad'])**2
    
    return rossrad


def is_in_ellipse(x0, y0, dE, d, x, y):
    '''
    Check if point (x,y) is contained in ellipse given by the equation

      (x-x1)**2     (y-y1)**2
      ---------  +  ---------  =  1
         a**2          b**2
         
    where:
    
      a = 0.5 * (dE + dW)
      b = dE
      x1 = x0 + 0.5 * (dE - dW)
      y1 = y0
    '''
    
    dW = np.max([d, dE])
    
    b = dE
    a = 0.5 * (dE + dW)
    
    x1 = x0 + 0.5*(dE - dW)
    y1 = y0
    
    return (x-x1)**2 / a**2 + (y-y1)**2 / b**2 <= 1


def len_deg_lon(lat):
    '''
    Returns the length of one degree of longitude (at latitude
    specified) in km.
    '''
    
    R = 6371. # Radius of Earth [km]
    
    return (np.pi/180.) * R * np.cos( lat * np.pi/180. )


def calculate_d(dE, lon, lat, rossrad, dt):
    '''
    Calculates length of search area to the west of central point.
    This is equal to the length of the search area to the east of
    central point (dE) unless in the tropics ( abs(lat) < 18 deg )
    in which case the distance a Rossby wave travels in one time step
    (dt, days) is used instead.
    '''
    
    if np.abs(lat) < 18 :
        # Rossby wave speed [km/day]
        c = interpolate.griddata(np.array([rossrad['lon'], rossrad['lat']]).T, rossrad['cR'], (lon, lat), method='linear') * 86400. / 1000.
        d = np.abs(1.75 * c * dt)
    else:
        d = dE
        
    return d


def track_eddies(eddies, det_eddies, tt, dt, dt_aviso, dE_aviso, rossrad, eddy_scale_min, eddy_scale_max):
    '''
    Given a map of detected eddies as a function of time (det_eddies)
    this function will update tracks of individual eddies at time step
    tt in variable eddies
    '''
    
    # List of unassigned eddies at time tt
    
    unassigned = range(det_eddies[tt]['N'])
    
    # For each existing eddy (t<tt) loop through unassigned eddies and assign to existing eddy if appropriate
    
    for ed in range(len(eddies)):
        
        # Check if eddy has already been terminated
        
        if not eddies[ed]['terminated']:
            
            # Define search region around centroid of existing eddy ed at last known position
            
            x0 = eddies[ed]['lon'][-1] # [deg. lon]
            y0 = eddies[ed]['lat'][-1] # [deg. lat]
            dE = dE_aviso/(dt_aviso/dt) # [km]
            d = calculate_d(dE, x0, y0, rossrad, dt) # [km]
            
            # Find all eddy centroids in search region at time tt
            
            is_near = is_in_ellipse(x0, y0, dE/len_deg_lon(y0), d/len_deg_lon(y0), det_eddies[tt]['lon'][unassigned], det_eddies[tt]['lat'][unassigned])
            
            # Check if eddies' amp  and area are between 0.25 and 2.5 of original eddy
            
            amp = eddies[ed]['amp'][-1]
            area = eddies[ed]['area'][-1]
            is_similar_amp = (det_eddies[tt]['amp'][unassigned] < amp*eddy_scale_max) * (det_eddies[tt]['amp'][unassigned] > amp*eddy_scale_min)
            is_similar_area = (det_eddies[tt]['area'][unassigned] < area*eddy_scale_max) * (det_eddies[tt]['area'][unassigned] > area*eddy_scale_min)
            
            # Check if eddies' type is the same as original eddy
            
            is_same_type = np.array([det_eddies[tt]['type'][i] == eddies[ed]['type'] for i in unassigned])
            
            # Possible eddies are those which are near, of the right amplitude, and of the same type
            
            possibles = is_near * is_similar_amp * is_similar_area * is_same_type
            if possibles.sum() > 0:
                
                # Of all found eddies, accept only the nearest one
                
                dist = np.sqrt((x0-det_eddies[tt]['lon'][unassigned])**2 + (y0-det_eddies[tt]['lat'][unassigned])**2)
                nearest = dist == dist[possibles].min()
                next_eddy = unassigned[np.where(nearest * possibles)[0][0]]
                
                # Add coordinats and properties of accepted eddy to trajectory of eddy ed
                
                eddies[ed]['lon'] = np.append(eddies[ed]['lon'], det_eddies[tt]['lon'][next_eddy])
                eddies[ed]['lat'] = np.append(eddies[ed]['lat'], det_eddies[tt]['lat'][next_eddy])
                eddies[ed]['amp'] = np.append(eddies[ed]['amp'], det_eddies[tt]['amp'][next_eddy])
                eddies[ed]['area'] = np.append(eddies[ed]['area'], det_eddies[tt]['area'][next_eddy])
                eddies[ed]['scale'] = np.append(eddies[ed]['scale'], det_eddies[tt]['scale'][next_eddy])
                eddies[ed]['time'] = np.append(eddies[ed]['time'], tt+1)
                eddies[ed]['sss_anom'] = np.append(eddies[ed]['sss_anom'], det_eddies[tt]['sss_anom'][next_eddy]) 
                eddies[ed]['sst_anom'] = np.append(eddies[ed]['sst_anom'], det_eddies[tt]['sst_anom'][next_eddy])
                eddies[ed]['mld'] = np.append(eddies[ed]['mld'], det_eddies[tt]['mld'][next_eddy])
                
                # Remove detected eddy from list of eddies available for assigment to existing trajectories
                
                unassigned.remove(next_eddy)
                
            # Terminate eddy otherwise
                
            else:
                
                eddies[ed]['terminated'] = True
                
    # Create "new eddies" from list of eddies not assigned to existing trajectories
                
    if len(unassigned) > 0:
        
        for un in unassigned:
            
            eddy_tmp = {}
            eddy_tmp['lon'] = np.array([det_eddies[tt]['lon'][un]])
            eddy_tmp['lat'] = np.array([det_eddies[tt]['lat'][un]])
            eddy_tmp['amp'] = np.array([det_eddies[tt]['amp'][un]])
            eddy_tmp['area'] = np.array([det_eddies[tt]['area'][un]])
            eddy_tmp['scale'] = np.array([det_eddies[tt]['scale'][un]])
            eddy_tmp['type'] = det_eddies[tt]['type'][un]
            eddy_tmp['time'] = np.array([tt+1])
            eddy_tmp['sss_anom'] = np.array([det_eddies[tt]['sss_anom'][un]])
            eddy_tmp['sst_anom'] = np.array([det_eddies[tt]['sst_anom'][un]])
            eddy_tmp['mld'] = np.array([det_eddies[tt]['mld'][un]])
            eddy_tmp['exist_at_start'] = False
            eddy_tmp['terminated'] = False
            eddies.append(eddy_tmp)
            
    return eddies

#****************************************************************

####################

### EDDY METRICS ###

####################

def eddy_metrics(eddies_tracked):
    
    ''' 
    *** companion function to a modified version of Eric Oliver's version of Chelton et al. (2011) eddy tracking code [ modified version available from jabris16 GitHub; original version available from ecjoliver GitHub ] ***
       
    Calculate and store additional eddy-specific metrics (rotational speed (U) at effective eddy radius, translation speed (c) of eddy, and nonlinearity (U/c) of eddy) to supplement those from eddy tracking.
    
    INPUT:
    eddies_tracked = array of tracked eddies
    
    OUTPUT:
    eddies =  dataset containing all eddy types and all metrics
        METRICS:
        'time' = time index of eddy feature in the analysis period
        'type' = anticyclonic or cyclonic or ACME eddy
        'age' = total age of eddy [days]
        'sss_anom' = eddying grid cell sea surface salinity (SSS) anomaly
        'sst_anom' = eddying grid cell sea surface temperature (SST) anomaly [deg C]
        'mld' = mixed layer thickness [m]
        'amp' = eddy amplitude [cm]
        'scale' = eddy radius [km]
        'rot_velocity' = rotational velocity of eddy feature [cm/s]
        'translation_speed' = translation speed of eddy feature [cm/s] (age must exceed 1 day; no value for first index as translation speed (c) must be calculated across subsequent eddy features/time steps)
        'nonlin' = instantaneous nonlinerity value of each eddy feature within its lifetime (age must exceed 1 day, no value for first index as c must be calculated across subsequent eddy features/time steps)
        'nonlin' = mean nonlinearity value of eddy across its lifetime
        'lat' = latitude of eddy feature [deg]
        'lon' = longitude of eddy feature [deg]
            
    '''
    
    grav = 981 # gravitational constant [cm/s]
    earth_radius = 637813700 # [cm]
    earth_rot = 0.000072921 # Earth rotation rate [rad/s]
    time_48hrs = 86400 * 2 # number of seconds in a day [seconds]
    radians_convert = np.pi/180 # radians to degrees conversion factor
    #e_factor = np.exp(-0.5) # multiplication factor as part of the calculation of U [ see Chelton et al. 2011 Appendix B ] 
    
    eddies = []
    
    for ed in range(len(eddies_tracked)):
        eddies_grid = {}
        eddies_grid['lon_index'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['lat_index'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['lon'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['lat'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['time'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['sss_anom'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['sst_anom'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['mld'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['age'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['amp'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['scale'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['rot_velocity'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['translation_speed'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['nonlin'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['area'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['type'] = []
        
        for x in range(len(eddies_tracked[ed]['lon'])):
            # Metrics from eddyTracking 
            eddies_grid['lat'][x] = eddies_tracked[ed]['lat'][x]
            eddies_grid['lon'][x] = eddies_tracked[ed]['lon'][x]
            eddies_grid['time'][x] = eddies_tracked[ed]['time'][x]
            eddies_grid['age'][x] = eddies_tracked[ed]['age']
            eddies_grid['amp'][x] = eddies_tracked[ed]['amp'][x]
            eddies_grid['scale'][x] = eddies_tracked[ed]['scale'][x]
            eddies_grid['sss_anom'][x] = eddies_tracked[ed]['sss_anom'][x] 
            eddies_grid['sst_anom'][x] = eddies_tracked[ed]['sst_anom'][x]
            eddies_grid['mld'][x] = eddies_tracked[ed]['mld'][x]
            eddies_grid['area'][x] = eddies_tracked[ed]['area'][x]
            eddies_grid['type'] = eddies_tracked[ed]['type']
            # rotational speed (U) calculation
            coriolis = np.absolute(2 * earth_rot * np.sin(eddies_tracked[ed]['lat'][x] * radians_convert))
            scalar = grav / coriolis
            amp_scale = (eddies_tracked[ed]['amp'][x] * 100) / (eddies_tracked[ed]['scale'][x] * 100000) # multiply by 100/100000 to conver amp/scale to cm
            rot_velocity = scalar * amp_scale # equation from Chelton et al. (2011) Appendix B when using eddy effective radius
            eddies_grid['rot_velocity'][x] = rot_velocity   
        for x in range(1,len(eddies_tracked[ed]['time'])-1):
            # translation speed (c) calculation [ calculated at time x based on distance between centroids at time x -1 and x + 1
            if eddies_tracked[ed]['age'] >= 3: # the eddy must have at least 3 track points to calc c
                angles = np.sin(eddies_tracked[ed]['lat'][x+1]) * np.sin(eddies_tracked[ed]['lat'][x-1]) + \
                    np.cos(eddies_tracked[ed]['lon'][x+1] - eddies_tracked[ed]['lon'][x-1]) * \
                    np.cos(eddies_tracked[ed]['lat'][x+1]) * np.cos(eddies_tracked[ed]['lat'][x-1])
                distance = earth_radius * np.arccos(angles) * radians_convert # distance equation from http://mathforum.org/library/drmath/view/54680.html
                speed = distance / time_48hrs # [cm/s]
                eddies_grid['translation_speed'][x] = speed
                # calculate nonlinearity (U/c)
                nonlin = rot_velocity / speed # U/c
                eddies_grid['nonlin'][x] = nonlin
                eddies_grid['nonlin'][0] = np.nan # first and last values will not have measures of nonlin (c relates to x -1 and x + 1)
                eddies_grid['nonlin'][-1] = np.nan # first and last values will not have measures of nonlin (c relates to x -1 and x + 1)
                
            #eddies_grid['nonlin_mean'][x] = np.nanmean(eddies_grid[ed]['nonlin'])
            #eddies_grid['type'] = eddies_grid[ed]['type']
            #eddies_grid['type'] = np.append(eddies_grid['type'], eddies_tracked[ed]['type'])
        eddies.append(eddies_grid)
        
    return eddies


def add_metrics(eddies, o2, depth, sss_baseline, sss, sst_baseline, sst, mld, spat_res, t, lat, lon):
    
    '''
    
    Add O2, SSS anomaly, SST anomaly and MLD to 'eddies' output from eddy_detection.py
    
    INPUTS:
    eddies = 'eddies' list from eddy_detection.py
    o2 = array of o2 data across lat, lon, depth
    depth = array of depth values extracted from model output
    sss_basline = baseline SSS to be subtraced from the absolute SSS value in order to calculate anomaly
    sst_baseline = as above but for SST
    spat_res = spatial resolution of model output data [degrees]
    t = array of time step values of the analysis period
    lat = array of latitude values of the analysis domain
    lon = array of longitude values of the analysis domain
    
    OUTPUTS:
    eddies_tempstore = temporary list of detected eddy features data with only additional metrics of SSS anomaly, SST anomaly and MLD; to later be added to a complete list with all metrics
    
    '''
    def gridcell_round(array, value):
        ''' function to assign eddy feature centroid (lon/lat) to the corresponding lon/lat index '''
        absolute_val_array = np.abs(array -value)
        grid_index = absolute_val_array.argmin()
        return grid_index

    eddies_tempstore = []
    
    for time in range(len(eddies)):
        print('add_metrics():' + str(time) + ' of ' + str(len(eddies)-1))
        
        eddies_met = {}
        eddies_met['sss_anom'] = np.zeros(len(eddies[time]['lon']))
        eddies_met['sst_anom'] = np.zeros(len(eddies[time]['lon']))
        eddies_met['mld'] = np.zeros(len(eddies[time]['lon']))
        eddies_met['o2'] = np.zeros(len(eddies[time]['lon']))
        eddies_met['o2_min'] = np.zeros(len(eddies[time]['lon']))
        eddies_met['o2_min_depth'] = np.zeros(len(eddies[time]['lon']))
        eddies_met['o2_anom'] = np.zeros(len(eddies[time]['lon']))
        for x in range(len(eddies[time]['lon'])):   
            lat_index = gridcell_round(lat, eddies[time]['lat'][x])
            lon_index = gridcell_round(lon, eddies[time]['lon'][x])
            # a) SSS
            sss_baseline_int = sss_baseline[list(t).index(eddies[time]['time'][0]), lat_index, lon_index]
            sss_abs = sss[list(t).index(eddies[time]['time'][0]), lat_index, lon_index]
            eddies_met['sss_anom'][x] = sss_abs - sss_baseline_int
            # b) SST
            sst_basline_int = sst_basline[list(t).index(eddies[time]['time'][0]), lat_index, lon_index]
            sst_abs = sst[list(t).index(eddies[time]['time'][0]), lat_index, lon_index]
            eddies_met['sst_anom'][x] = sst_abs - sst_baseline_int
            # c) MLD
            eddies_met['mld'][x] = mld[list(t).index(eddies[time]['time'][0]), lat_index, lon_index]
            # d) o2
            o2_baseline_int = o2_baseline[list(t).index(eddies[time]['time'][0]), :, lat_index, lon_index]
            eddies_met['o2'][x] = o2[list(t).index(eddies[time]['time'][0]), :, lat_index, lon_index]
            eddies_met['o2_anom'] = eddies_met['o2'][x] - o2_baseline_int
            eddies_met['o2_min'][x] = np.amin(eddies_met['o2'][x])
            eddies_met['o2_min_depth'][x] = depth[np.where(eddies_met['o2'][x] == eddies_met['o2_min'][x])]

        eddies_tempstore.append(eddies_met)
        
    return eddies_tempstore

def eddies_listVER2(lon_eddies_a, lat_eddies_a, amp_eddies_a, area_eddies_a, scale_eddies_a, time_eddies_a, lon_eddies_c, lat_eddies_c, amp_eddies_c, area_eddies_c, scale_eddies_c, time_eddies_c, eddies_tempstore):
    
    ''' 
    
    Relist the detected eddy features from eddy_detection.py to now include SSS and SST anomalies, and MLD.
    
    '''
    eddies = []
    
    for ed in range(len(lon_eddies_c)):
        eddy_tmp = {}
        eddy_tmp['lon'] = np.append(lon_eddies_a[ed], lon_eddies_c[ed])
        eddy_tmp['lat'] = np.append(lat_eddies_a[ed], lat_eddies_c[ed])
        eddy_tmp['amp'] = np.append(amp_eddies_a[ed], amp_eddies_c[ed])
        eddy_tmp['area'] = np.append(area_eddies_a[ed], area_eddies_c[ed])
        eddy_tmp['scale'] = np.append(scale_eddies_a[ed], scale_eddies_c[ed])
        eddy_tmp['type'] = list(repeat('anticyclonic',len(lon_eddies_a[ed]))) + list(repeat('cyclonic',len(lon_eddies_c[ed])))
        eddy_tmp['N'] = len(eddy_tmp['lon'])
        eddy_tmp['time'] = np.append(time_eddies_a[ed], time_eddies_c[ed])
        
        # now with additional metrics
        eddy_tmp['sss_anom'] = eddies_tempstore[ed]['sss_anom']
        eddy_tmp['sst_anom'] = eddies_tempstore[ed]['sst_anom']
        eddy_tmp['mld'] = eddies_tempstore[ed]['mld']
        eddy_tmp['o2'] = eddies_tempstore[ed]['o2']
        eddy_tmp['o2_min'] = eddies_tempstore[ed]['o2_min']
        eddy_temp['o2_min_depth'] = eddies_tempstore[ed]['o2_min_depth']
        
        eddies.append(eddy_tmp)
        
    return eddies

###

def acme_distinct(eddies):
    
    '''    
    
    Update list of detected eddies by distinguishing and reclassifying ACMEs and "typical" anticyclones by virtue of their negative SSS and SST anomalies in an ACME vs typically positive SSS and SST anomaly "typical" anticyclones.
    
    INPUTS:
    eddies = updated version of 'eddies', output of eddies_listVER2()
    
    OUTPUT:
    Modified version of 'eddies' with reclassification
    
    '''
    
    for time in range(len(eddies)):
        for x in range(len(eddies[time]['lon'])):
            if eddies[time]['type'][x] == 'anticyclonic' and eddies[time]['sss_anom'][x] < 0 and eddies[time]['sst_anom'][x] < 0:
                eddies[time]['type'][x] = 'acme'
            else:
                pass
            
    return eddies
