#****************************************************************
# load required modules
import numpy as np

# CONSTANTS [ do not adjust ]
res_intend = 0.25 # intended spatial resolution for model output, necessary if data used for tracking differs from the intended value [degrees]

#****************************************************************

# PARAMETERS [ adjust as necessary ] 

# data/plot output directories [in this case these are in the freeglorys/physical data directory]
data_dir = '/home/bridge/ja16048/mscr/subsurfacedata_tracking/etna/detection_output/' # where to store detection/tracking output
plot_dir = '/home/bridge/ja16048/mscr/subsurfacedata_tracking/etna/detection_output/' # where to store plot output

# physical data related pathwas/names
pathroot = '/home/bridge/ja16048/mscr/subsurfacedata_tracking/etna/data_load/freeglorys_data-load/' # directory containing SSH data to be loaded in for tracking
filename = 'freeglorys_etna_1992-2018.nc' # name of netcdf file for analysis [ containing SSH, SST, MLD, Lat, Lon data ]; if these come in separate files the loading of data in part 2) of eddy_detection.py will require adjusting 
filenameSAL = 'freeglorys_etna_SSS_1992-2018.nc' # name of netcdf file containing salinity [SSS] data

# biogeochem data related pathways/names
pathroot_o2 = '/home/bridge/ja16048/mscr/subsurfacedata_tracking/etna/data_load/freebiorys_data-load/'
filename_o2 = 'freebiorys_etna_1992-2018.nc'

num_years = 27 # number of years of data in the analysis period

# domain lat, lon values
lon1 = -170
lon2 = -100
lat1 = 13
lat2 = 50

# model output info
NAME = 'freeglorys' # Which dataset/model run for which to detect eddies
run = NAME
T = 1*9862 # Number of time steps to loop over, i.e. number of timesteps in SSH data
res = 0.25 # horizontal resolution of SSH field [degrees]
dt = 1. # Sample rate of detected eddies [days]

# spatial filter cutoff values
cut_lon = 20. # cutoff wavelenth in longitudinal direction (for filtering) [degrees]
cut_lat = 10. # cutoff wavelenth in latitudinal direction (for filtering) [degrees]

ssh_crit_max = 1.
dssh_crit = 0.01
ssh_crits = np.arange(-ssh_crit_max, ssh_crit_max+dssh_crit, dssh_crit)
ssh_crits = np.flipud(ssh_crits)

area_correction = res**2 / res_intend**2 # correction for different resolutions
Npix_min = np.floor(8*area_correction) # min number of eddy pixels
Npix_max = np.floor(1000*area_correction) # max number of eddy pixels

amp_thresh = 0.01 # minimum eddy amplitude [m]

''' this needs to be addressed for the latidudes approaching the tropics 
i.e. the 1200km rule as mentioned in Chelton et al. (2011)
'''

d_thresh = 400. # max linear dimension of eddy [km] ; Only valid outside Tropics (see Chelton et al. (2011), pp. 207)
# max linear dimension of eddy increases lineary between 25 deg N and equator up to 1200km (Chelton et al., 2011); 
# these increases are included by creating 5 deg bins (below) each with their own threshold value
increment = (1200 - 400) / 5.
d_thresh_20 = d_thresh + 1*increment
d_thresh_15 = d_thresh + 2*increment
d_thresh_10 = d_thresh + 3*increment
d_thresh_5 = d_thresh + 4*increment
d_thresh_0 = d_thresh + 5*increment


dt_aviso = 1. # Sample rate used in Chelton et al. (2011) [days]
dE_aviso = 150. # Length of search ellipse to East of eddy used in Chelton et al. (2011) [km]

#This is only used in eddy_tracking so, has been called explictly there. This removes the circular dependency so we can import params into eddy_functions!
#rossrad = eddy.load_rossrad() # Atlas of Rossby radius of deformation and first baroclinic wave speed (Chelton et al. 1998)

eddy_scale_min = 0.25 # min ratio of amplitude of new and old eddies
eddy_scale_max = 2.5 # max ratio of amplidude of new and old eddies

dt_save = 100 # Step increments at which to save data while tracking eddies
