''' 

Functions for calculating seasonal means of variables (SST and SSS), in order to later calculate SST/SSS anomaly

Jamie Atkins.


'''

#****************************************************************

import numpy as np

#****************************************************************

#############################################

   ### MEAN CLIMATOLOGY CALCULATIONS ###

#############################################

def clim_calcs(field, num_years):

    '''

    Function to calculate mean value of a given variable ('field') for each grid cell across each day of the year.

    '''

    days_in_year = 365
    array_shape = np.array(([days_in_year, len(field[0,:,0]), len(field[0,0,:])]))
    t = num_years * days_in_year
    clim_mean = np.zeros(array_shape)
    t_len = np.array(range(t))
    t_store = np.zeros((days_in_year, int(t / days_in_year)), dtype = object)

    # subset timesteps for upcoming calculations
    for i in range(days_in_year):
        for j in range(num_years):
            if i == 0:
                t_store[0][j] = t_len[j * days_in_year]
            else:
                t_store[i][:] = t_store[0] + int(i)

    print('starting clim_calcs():')

    # calculate mean values
    for x in range(len(field[0,:,0])):
        print(str(x) + ' of ' + str(len(field[0,:,0])-1))
        for y in range(len(field[0,0,:])):
            for i in range(days_in_year):
                var_int = []
                for j in range(num_years):
                    var_int.append(field[t_store[i][j],x,y])
                var_array = np.array(var_int)
                var_mean = np.nanmean(var_array)
                clim_mean[i,x,y] = var_mean

    return clim_mean

#def clim_calcs_o2(field, num_years):
#    
#    '''
#    
#    Function to calculate mean value of oxygen ('field') [involving a depth component, unlike sss and sst] for each grid cell across depth range and across each day of the year
#    
#    '''
#
#    days_in_year = 365
#    t = num_years * days_in_year
#    array_shape = np.array(([days_in_year, len(field[0,:,0,0]), len(field[0,0,:,0]), len(field[0,0,0,:])]))
#    clim_mean = np.zeros(array_shape)
#    t_len = np.array(range(t))
#    t_store = np.zeros((days_in_year, int(t / days_in_year)), dtype = object)
#    
#    # subset timesteps for upcoming calculations
#    for i in range(days_in_year):
#        for j in range(num_years):
#            if i == 0:
#                t_store[0][j] = t_len[j * days_in_year]
#            else:
#                t_store[i][:] = t_store[0] + int(i)
#                
#    print('starting clim_calcs_o2():')
#    
#    # calculate mean values 
#    for d in range(len(field[0,:,0,0])):
#        print(str(d) + ' of ' + str(len(field[0,:,0,0])-1))
#        for x in range(len(field[0,0,:,0])):
#            for y in range(len(field[0,0,0,:])):
#                for i in range(days_in_year):
#                    var_int = []
#                    for j in range(num_years):
#                        var_int.append(field[t_store[i][j],d,x,y])
#                    var_array = np.array(var_int)
#                    var_mean = np.nanmean(var_array)
#                    clim_mean[i,d,x,y] = var_mean
#
#    return clim_mean

def clim_repeat(clim_mean, t):
    
    '''
    Repeat clim_mean across the number of timesteps in the analysis period to form a baseline which can be subtracted from the absolute values to calculate the anomaly values.
    
    *** N.B. ASSUMES THE ANALYSIS PERIOD STARTS ON JAN 1 OF A GIVEN YEAR. WORKING ON SOLUTION; FOR NOW HAVE TO TAKE THE NEAREST AVAILABLE YEAR STARTING JAN 1 ***
    
    INPUTS:
    clim_mean = output from 'clim_cals()' function above
    t = number of timesteps in the analysis period

    '''
    
    days_in_year = 365
    t_multi = int(t / days_in_year) + 1 # how many years to repeat for (+ 1 so that remainder can be added on to account for incomplete years)
    t_intermed = int(t / days_in_year) * days_in_year # floored nearest complete year, number of timesteps
    t_rem = t - t_intermed # number of remaining days
    
    lat_indexsum = len(clim_mean[0,:,0]) # number of latitude indices
    lon_indexsum = len(clim_mean[0,0,:]) # number of longitude indices
    baseline_mean_int = [clim_mean,] * (t_multi -1) 
    baseline_mean_add = clim_mean[0:t_rem,:,:]
    baseline_mean_int = np.reshape(baseline_mean_int, ((t_multi - 1) * days_in_year, lat_indexsum, lon_indexsum))
    baseline_mean = np.concatenate((baseline_mean_int, baseline_mean_add))
    
    return baseline_mean

#def clim_repeat_o2(clim_mean, t):
#    
#    '''
#    Same as above but for o2 which has a depth component, unlike sss and sst
#    '''
#    
#    days_in_year = 365
#    t_multi = int(t / days_in_year) + 1 # how many years to repeat for (+ 1 so that remainder can be added on to account for incomplete years)
#    t_intermed = int(t / days_in_year) * days_in_year # floored nearest complete year, number of timesteps
#    t_rem = t - t_intermed # number of remaining days
#    
#    lat_indexsum = len(clim_mean[0,0,:,0]) # number of latitude indices
#    lon_indexsum = len(clim_mean[0,0,0,:]) # number of longitude indices
#    depth_indexsum = len(clim_mean[0,:,0,0]) # number of depth indices
#    baseline_mean_int = [clim_mean,] * (t_multi -1) 
#    baseline_mean_add = clim_mean[0:t_rem,:,:,:]
#    baseline_mean_int = np.reshape(baseline_mean_int, ((t_multi - 1) * days_in_year, depth_indexsum, lat_indexsum, lon_indexsum))
#    baseline_mean = np.concatenate((baseline_mean_int, baseline_mean_add))
#    
#    return baseline_mean

def o2_filter(field, window):
    
    '''
    
    Annual running mean for each cell.
    
    'field' = 4D o2 field
    'window' = window size [days] for convolution; 365 days in this case
    
    '''
    
    # make storage array for filtered o2 data
    array_shape = np.array(([len(field[0,:,0,0]), len(field[0,0,:,0]), len(field[0,0,0,:])]))
    o2_filt = np.zeros(array_shape.shape)
    
    # running mean by convolution
    days_in_year = 365
    window = days_in_year
    for d in range(len(field[0,:,0,0])):
        print(str(d))
        for x in range(len(field[0,0,:,0])):
            for y in range(len(field[0,0,0,:])): 
                o2_runmean = running_mean(field[:,d,x,y], window) # creates a running mean 1d array at grid cell [d,x,y]
                o2_filt[d,x,y] = o2_runmean # each grid cell [d,x,y] has the 1d array assigned to it where the length of the 1d array is the number of timesteps in analysis period
                                            # then can refer back to the grid cell and time at which eddy occurs to access the background value to calculate anomaly
    return o2_filt
                
 
def running_mean(field_1d, window):
    
    '''
    
    'field' = 1D o2 array; i.e. each 1D array as loop through each grid cell (lat,lon,depth) across the analysis period
    'window' = window for convolution; where window is a year, i.e. an annual running mean
    
    '''
    
    return np.convolve(field_1d, np.ones(window), 'same') / window
