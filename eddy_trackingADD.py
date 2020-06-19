'''

Software for the tracking of eddies in model output following Chelton et al., Progress in Oceanography, 2011.

Adapted from and using baseline code from Eric Oliver eddyTracking under the terms of the GNU General Public License.

Including additional detection designed by Jamie Atkins involving additional metrics of eddy dynamics (i.e. rotational and translational velocity and nonlinearity)


'''

# Load required modules

import numpy as np
import eddy_functionsMOD as eddy

# Load parameters

from paramsADD import *

# Automated eddy tracking

data = np.load(data_dir+'eddy_det_'+run+'.npz')
det_eddies = data['eddies']

# Initialize eddies discovered at first time step

eddies = eddy.eddies_init(det_eddies)

# Stitch eddy tracks together at future time steps

print('eddy tracking started')
print("number of time steps to loop over: ",T)

rossrad = eddy.load_rossrad() # Atlas of Rossby radius of deformation and first baroclinic wave speed (Chelton et al. 1998)

for tt in range(1, T):
    
    print("timestep: " ,tt+1,". out of: ", T)
    
    # Track eddies from time step tt-1 to tt and update corresponding tracks and/or create new eddies
    
    eddies = eddy.track_eddies(eddies, det_eddies, tt, dt, dt_aviso, dE_aviso, rossrad, eddy_scale_min, eddy_scale_max)
    
    # Save data incrementally
    
    if( np.mod(tt, dt_save)==0 ):
        
        np.savez(data_dir+'eddy_track_'+run, eddies=eddies)
        
# Add keys for eddy age and flag if eddy was still in existence at end of run
        
for ed in range(len(eddies)):
    
    eddies[ed]['age'] = len(eddies[ed]['lon'])
    
#****************************************************************
# save output [ containing no additional metrics ] 
np.savez(data_dir+'eddy_track_'+run, eddies=eddies)

#****************************************************************
# PART 2) ADDITIONAL METRICS [ INCLUDING U, c, NONLINEARITY ETC. ] 

# additional metrics
eddies = eddy.eddy_metrics(eddies)

# subset into cyclonic, anticyclonic and ACME
eddies_cyc, eddies_acme, eddies_anti = eddy.eddytype_subset(eddies)

#****************************************************************
# save output

np.savez(data_dir + 'eddy_trackAddMetrics_' + run, eddies = eddies)
np.savez(data_dir + 'eddy_trackAnticyclones_' + run, eddies_anti = eddies_anti)
np.savez(data_dir + 'eddy_trackCyclones_' + run, eddies_cyc = eddies_cyc)
np.savez(data_dir + 'eddy_trackACME_' + run, eddies_acme = eddies_acme)

#****************************************************************
