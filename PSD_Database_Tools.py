#!/usr/bin/env python

# Collection of codes to query the PSD database 


from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np
import pickle 
import glob

#### PSD Statistics ##############
## Inputs:
# Path to the directory containing all the PSD file directories 
# station name
# start and end times
# the percentile you want 

# Outputs
# Periods, Power in dB rel. 1 (m/s^2)^2/Hz. 

def get_PSD_statistics(DB_path, sta, starttime, endtime, prctile):

    starttime, endtime = UTCDateTime(starttime), UTCDateTime(endtime)
    All_PSD_Paths = glob.glob(DB_path + sta + '_PSD/*.pckl') 
    
    PSD_Array = []

    for path in All_PSD_Paths:
        directs = path.split('/')
        f = directs[2]
        info = f.split('_')
        year, time = info[5], info[6] 
        day, hour = time[0:3],time[3:5]
        PSDtime = UTCDateTime(year = int(year), julday = int(day), hour = int(hour))
        if PSDtime >= starttime and PSDtime <= endtime:
            pickle_file = open(path, 'rb')
            PSD_Array.append(pickle.load(pickle_file))
            
    PSD_Array = np.asarray(PSD_Array)    
    The_PSD_You_Want = np.percentile(PSD_Array, prctile,0)
    
    F_path = glob.glob(DB_path +  sta + '_PSD/*_freqs.txt')
    freqs = np.loadtxt(str(F_path[0]))
    
    return 1./freqs, The_PSD_You_Want
    
